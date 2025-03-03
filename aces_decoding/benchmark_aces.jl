using QuantumACES, LinearAlgebra, SparseArrays, StatsBase, BenchmarkTools
enter_folder("aces_decoding")
# Set up parameters
dist = 3
r_1 = 0.05 / 100
r_2 = 0.4 / 100
r_m = 0.8 / 100
r_im = r_m
r_r = 0.2 / 100
r_1_std_log = 0.5
r_2_std_log = r_1_std_log
r_m_std_log = 0.25
r_im_std_log = r_m_std_log
r_r_std_log = r_m_std_log
seed = UInt(0)
combined = true
# Generate the circuit
dep_param_load = get_dep_param(r_1, r_2, r_m)
dep_param = get_dep_param(;
    r_1 = r_1,
    r_2 = r_2,
    r_m = r_m,
    r_im = r_im,
    r_r = r_r,
    combined = combined,
)
log_param = get_log_param(;
    r_1 = r_1,
    r_2 = r_2,
    r_m = r_m,
    r_im = r_im,
    r_r = r_r,
    r_1_std_log = r_1_std_log,
    r_2_std_log = r_2_std_log,
    r_m_std_log = r_m_std_log,
    r_im_std_log = r_im_std_log,
    r_r_std_log = r_r_std_log,
    seed = seed,
    combined = combined,
)
rotated_param = get_rotated_param(dist)
# Load the design
tuple_number = 30
repeat_numbers = [49; 99; 25; 49; 49; 99]
d = load_design(rotated_param, dep_param_load, tuple_number, repeat_numbers, true, :none)
merit = calc_merit(d)
# Load the simulated ACES data at distance 25
dist_big = 25
rotated_param_big = get_rotated_param(dist_big)
rotated_planar_big = get_circuit(rotated_param_big, dep_param)
rotated_planar_big_log = get_circuit(rotated_param_big, log_param)
shots_set = [10^6; 10^7]
times_factor = sum(d.shot_weights .* d.tuple_times)
budget_set = round.(Int, times_factor * shots_set)
aces_data_big = load_aces(
    rotated_param_big,
    log_param,
    tuple_number,
    repeat_numbers,
    false,
    :none,
    budget_set,
    seed,
)
d_big = aces_data_big.d
@assert d_big.c == rotated_planar_big_log
# Print the noise estimation errors
noise_est_set = aces_data_big.noise_est_coll[1, :]
noise_error_set = aces_data_big.noise_error_coll[1, :]
shots_num = length(shots_set)
for s in 1:shots_num
    println("ACES with S' = 10^$(round(log10(shots_set[s]), digits=3)):")
    print(noise_error_set[s])
end
# Set up the overall estimation benchmark
benchmark_seconds = 600
function benchmark_estimate(
    d_big::Design,
    est_eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}},
    est_covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}},
)
    # Intialise variables
    gate_data = d_big.c.gate_data
    # Calculate estimation quantities
    tuple_number = length(d_big.tuple_set)
    est_eigenvalues =
        vcat([mean.(est_eigenvalues_experiment_ensemble[i]) for i in 1:tuple_number]...)
    est_covariance = calc_covariance(
        d_big,
        est_eigenvalues_experiment_ensemble,
        est_covariance_experiment_ensemble;
        weight_time = false,
        warning = d_big.full_covariance,
    )
    est_covariance_log = calc_covariance_log(est_covariance, est_eigenvalues)
    # Clip the eigenvalues, covariance matrix, and design matrix if appropriate
    clipped_indices = QuantumACES.get_clipped_indices(est_eigenvalues; warning = false)
    design_matrix = convert(SparseMatrixCSC{Float64, Int}, d_big.matrix[clipped_indices, :])
    est_eigenvalues = est_eigenvalues[clipped_indices]
    est_covariance_log = est_covariance_log[clipped_indices, clipped_indices]
    est_diag_covariance_log = sparse(Diagonal(est_covariance_log))
    est_diag_covariance_log_inv_factor =
        sparse(Diagonal(diag(est_diag_covariance_log) .^ (-1 / 2)))
    est_diag_covariance_log_inv =
        est_diag_covariance_log_inv_factor' * est_diag_covariance_log_inv_factor
    # Estimate the gate eigenvalues
    wls_unproj_gate_eigenvalues = QuantumACES.estimate_gate_eigenvalues(
        design_matrix,
        est_diag_covariance_log_inv_factor,
        est_eigenvalues,
    )
    # Project the gate error probabilities into the simplex
    wls_precision_matrix = calc_precision_matrix(
        design_matrix,
        wls_unproj_gate_eigenvalues,
        est_diag_covariance_log_inv,
    )
    (wls_gate_eigenvalues, wls_unproj_gate_probabilities_vec, wls_gate_probabilities_vec) =
        QuantumACES.split_project_gate_eigenvalues(
            wls_unproj_gate_eigenvalues,
            wls_precision_matrix,
            gate_data,
        )
    wls_gate_probabilities =
        QuantumACES.get_gate_probabilities_dict(wls_gate_probabilities_vec, gate_data)
    return wls_gate_probabilities::Dict{Gate, Vector{Float64}}
end
est_eigenvalues_experiment_ensemble = noise_est_set[1].eigenvalues_experiment_ensemble
est_covariance_experiment_ensemble = noise_est_set[1].covariance_experiment_ensemble
wls_gate_probabilities = noise_est_set[1].wls_gate_probabilities
est_wls_gate_probabilities = benchmark_estimate(
    d_big,
    est_eigenvalues_experiment_ensemble,
    est_covariance_experiment_ensemble,
)
@assert all(
    est_wls_gate_probabilities[gate] ≈ wls_gate_probabilities[gate] for
    gate in keys(est_wls_gate_probabilities)
)
# Benchmark the overall estimation
b_tot = @benchmarkable benchmark_estimate(
    d_big,
    est_eigenvalues_experiment_ensemble,
    est_covariance_experiment_ensemble,
) seconds = benchmark_seconds
b_tot_run = run(b_tot)
println("Total estimation benchmark:\n")
display(b_tot_run)
println("")
# Set up the linear inversion benchmark
gate_data = d_big.c.gate_data
tuple_number = length(d_big.tuple_set)
est_eigenvalues =
    vcat([mean.(est_eigenvalues_experiment_ensemble[i]) for i in 1:tuple_number]...)
est_covariance = calc_covariance(
    d_big,
    est_eigenvalues_experiment_ensemble,
    est_covariance_experiment_ensemble;
    weight_time = false,
    warning = d_big.full_covariance,
)
est_covariance_log = calc_covariance_log(est_covariance, est_eigenvalues)
clipped_indices = QuantumACES.get_clipped_indices(est_eigenvalues; warning = false)
design_matrix = convert(SparseMatrixCSC{Float64, Int}, d_big.matrix[clipped_indices, :])
est_eigenvalues = est_eigenvalues[clipped_indices]
est_covariance_log = est_covariance_log[clipped_indices, clipped_indices]
est_diag_covariance_log = sparse(Diagonal(est_covariance_log))
est_diag_covariance_log_inv_factor =
    sparse(Diagonal(diag(est_diag_covariance_log) .^ (-1 / 2)))
est_diag_covariance_log_inv =
    est_diag_covariance_log_inv_factor' * est_diag_covariance_log_inv_factor
# Benchmark the linear inversion
b_inv = @benchmarkable QuantumACES.estimate_gate_eigenvalues(
    $(design_matrix),
    $(est_diag_covariance_log_inv_factor),
    $(est_eigenvalues),
) seconds = benchmark_seconds
b_inv_run = run(b_inv)
println("Linear inversion benchmark:\n")
display(b_inv_run)
println("")
# Set up the projection benchmark
wls_unproj_gate_eigenvalues = QuantumACES.estimate_gate_eigenvalues(
    design_matrix,
    est_diag_covariance_log_inv_factor,
    est_eigenvalues,
)
wls_precision_matrix = calc_precision_matrix(
    design_matrix,
    wls_unproj_gate_eigenvalues,
    est_diag_covariance_log_inv,
)
# Benchmark the projection
b_proj = @benchmarkable QuantumACES.split_project_gate_eigenvalues(
    $(wls_unproj_gate_eigenvalues),
    $(wls_precision_matrix),
    $(gate_data),
) seconds = benchmark_seconds
b_proj_run = run(b_proj)
println("Projection benchmark:\n")
display(b_proj_run)
println("")
println("Complete: benchmark_decoding.jl")

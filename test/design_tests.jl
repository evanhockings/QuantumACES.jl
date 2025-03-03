using QuantumACES, Random, StatsBase, Test
# Make sure to delete the data folder even if previous tests failed
enter_folder("test")
rm("data"; force = true, recursive = true)
# Set up codes
dist = 3
vertical_dist = 4
horizontal_dist = 5
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
rotated_param = get_rotated_param(dist)
big_rotated_param = get_rotated_param(vertical_dist, horizontal_dist)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
rotated_planar = get_circuit(rotated_param, dep_param)
rotated_planar_log = update_noise(rotated_planar, log_param)
big_rotated_planar = get_circuit(big_rotated_param, log_param)
# Test circuit noise updating
log_gate_probabilities = init_gate_probabilities(rotated_planar.total_gates, log_param)
@testset "Circuit noise updating" begin
    rotated_planar_update = update_noise(rotated_planar, log_gate_probabilities)
    @test rotated_planar_log.gate_probabilities == rotated_planar_update.gate_probabilities
    @test rotated_planar_log.gate_eigenvalues == rotated_planar_update.gate_eigenvalues
end
# Set up optimisation parameters
ls_type = :gls
est_type = :prod
max_steps = 2
add_circuit = false
max_depth = 60
max_cycles = 1
repeat_points = 2
max_excursions = 2
excursion_length = 2
extra_tuple_number = 2
# Set up scaling parameters
dist_max = 5
precision = 1.0
min_repetitions = 3
# Set up simulation parameters
budget_set = [10^7; 2 * 10^7]
budget_number = length(budget_set)
repetitions = 2
max_samples = 10^5
z_score_cutoff_upper = 3.5
z_score_cutoff_lower = -3.5
clip_number = 100
clip_z_cutoff = 2.0
trial_number = 5
trial_z_abs_cutoff = 0.5
trial_z_mean_cutoff = 0.1
# Set up a basic design
rot_basic = get_basic_tuple_set(rotated_planar)
d_rot_basic = generate_design(rotated_planar, rot_basic; save_data = true)
rot_basic_merit = calc_merit(d_rot_basic)
# Optimise a design
d_rot_opt = optimise_design(
    rotated_planar;
    options = OptimOptions(;
        ls_type = ls_type,
        est_type = est_type,
        save_data = true,
        max_steps = max_steps,
        add_circuit = add_circuit,
        max_depth = max_depth,
        max_cycles = max_cycles,
        repeat_points = repeat_points,
        max_excursions = max_excursions,
        excursion_length = excursion_length,
        extra_tuple_number = extra_tuple_number,
        seed = seed,
    ),
)
rot_opt_merit = calc_merit(d_rot_opt)
@testset "Design optimisation" begin
    if ls_type == :gls
        @test rot_opt_merit.gls_expectation < rot_basic_merit.gls_expectation
    elseif ls_type == :wls
        @test rot_opt_merit.wls_expectation < rot_basic_merit.wls_expectation
    elseif ls_type == :ols
        @test rot_opt_merit.ols_expectation < rot_basic_merit.ols_expectation
    else
        throw(error("Unsupported least squares type $(ls_type)."))
    end
end
# Test conversion utilities
rand_eigenvalue_indices = randperm(sum(length.(d_rot_opt.mapping_ensemble)))
@testset "Conversion utilities" begin
    # Test index conversion
    @test rand_eigenvalue_indices ==
          pair_to_eig_idx(d_rot_opt, eig_to_pair_idx(d_rot_opt, rand_eigenvalue_indices))
    for (idx, circuit_tuple) in pairs(d_rot_opt.tuple_set)
        tuple_circuit = d_rot_opt.c.circuit[circuit_tuple]
        # Test mapping calculation
        mapping_set = d_rot_opt.mapping_ensemble[idx]
        initial_paulis = [m.initial for m in mapping_set]
        final_paulis = [m.final for m in mapping_set]
        test_final_paulis =
            [calc_pauli(initial, tuple_circuit) for initial in initial_paulis]
        @test final_paulis == test_final_paulis
        # Test Pauli string conversion
        string_initial_paulis =
            [string_to_pauli(pauli_to_string(initial)) for initial in initial_paulis]
        string_final_paulis =
            [string_to_pauli(pauli_to_string(final)) for final in final_paulis]
        @test initial_paulis == string_initial_paulis
        @test final_paulis == string_final_paulis
    end
end
# Test saving and loading and deletion
@testset "Design IO" begin
    d_rot_basic_load = load_design(
        d_rot_basic.c.circuit_param,
        d_rot_basic.c.noise_param,
        length(d_rot_basic.tuple_set),
        d_rot_basic.tuple_set_data.repeat_numbers,
        d_rot_basic.full_covariance,
        d_rot_basic.ls_type,
    )
    d_rot_opt_load = load_design(
        d_rot_opt.c.circuit_param,
        d_rot_opt.c.noise_param,
        length(d_rot_opt.tuple_set),
        d_rot_opt.tuple_set_data.repeat_numbers,
        d_rot_opt.full_covariance,
        d_rot_opt.ls_type,
    )
    @test d_rot_basic_load == d_rot_basic
    @test d_rot_opt_load == d_rot_opt
    delete_design(d_rot_basic)
    delete_design(d_rot_opt)
end
# Examine the scaling of the merit for depolarising and log-normal noise
merit_scaling = calc_merit_scaling(d_rot_opt, dist_max; save_data = true)
scaling_fit = get_scaling_fit(merit_scaling)
d_rot_opt_log = update_noise(d_rot_opt, log_param)
ensemble_scaling = calc_ensemble_scaling(
    d_rot_opt_log,
    dist_max;
    precision = precision,
    min_repetitions = min_repetitions,
    print_repetitions = min_repetitions - 1,
    seed = seed,
    save_data = true,
)
ensemble_fit = get_ensemble_fit(ensemble_scaling; precision = precision)
# Test saving and loading and deletion
@testset "Scaling IO" begin
    merit_scaling_load = load_scaling(d_rot_opt)
    ensemble_scaling_load = load_scaling(d_rot_opt_log)
    @test merit_scaling_load == merit_scaling
    @test ensemble_scaling_load == ensemble_scaling
    delete_scaling(merit_scaling)
    delete_scaling(ensemble_scaling)
end
# Test design noise updating
@testset "Design noise updating" begin
    d_rot_opt_update = update_noise(d_rot_opt, log_gate_probabilities)
    @test d_rot_opt_log.c.gate_probabilities == d_rot_opt_update.c.gate_probabilities
    @test d_rot_opt_log.c.gate_eigenvalues == d_rot_opt_update.c.gate_eigenvalues
end
# Transfer the tuple set to a differently-sized code
d_rot_big = generate_design(
    big_rotated_planar,
    d_rot_opt.tuple_set_data;
    shot_weights = d_rot_opt.shot_weights,
    diagnostics = true,
    save_data = true,
)
# Test that the merit is improved by better LS estimators
rot_big_merit = calc_merit(d_rot_big)
@testset "LS improvement" begin
    @test rot_big_merit.gls_expectation <= rot_big_merit.wls_expectation
    @test rot_big_merit.gls_marginal_expectation <= rot_big_merit.wls_marginal_expectation
    @test rot_big_merit.gls_relative_expectation <= rot_big_merit.wls_relative_expectation
    @test rot_big_merit.wls_expectation <= rot_big_merit.ols_expectation
    @test rot_big_merit.wls_marginal_expectation <= rot_big_merit.ols_marginal_expectation
    @test rot_big_merit.wls_relative_expectation <= rot_big_merit.ols_relative_expectation
    @test rot_big_merit.gls_variance <= rot_big_merit.wls_variance
    @test rot_big_merit.gls_marginal_variance <= rot_big_merit.wls_marginal_variance
    @test rot_big_merit.gls_relative_variance <= rot_big_merit.wls_relative_variance
    @test rot_big_merit.wls_variance <= rot_big_merit.ols_variance
    @test rot_big_merit.wls_marginal_variance <= rot_big_merit.ols_marginal_variance
    @test rot_big_merit.wls_relative_variance <= rot_big_merit.ols_relative_variance
end
# Simulate the design
aces_data_rot_big = simulate_aces(
    d_rot_big,
    budget_set;
    seed = seed,
    max_samples = max_samples,
    detailed_diagnostics = true,
    save_data = true,
    clear_design = true,
)
aces_data_rot_big = simulate_aces(
    d_rot_big,
    budget_set;
    repetitions = repetitions,
    seed = seed,
    split = true,
    save_data = true,
    save_interval = 1,
)
pretty_print(aces_data_rot_big, rot_big_merit)
noise_score_coll_big = get_noise_score(aces_data_rot_big, rot_big_merit)
model_violation_coll_big = get_model_violation(aces_data_rot_big)
noise_error = aces_data_rot_big.noise_error_coll[end, budget_number]
noise_score = noise_score_coll_big[end, budget_number]
display(noise_error)
display(noise_score)
delete_aces(aces_data_rot_big)
# Estimate the noise with a combined design
d_rot_comb = get_combined_design(d_rot_big)
comb_noise_est_coll_big = [
    estimate_gate_noise(d_rot_comb, noise_est) for
    noise_est in aces_data_rot_big.noise_est_coll
]
comb_model_violation_coll_big = get_model_violation(d_rot_comb, comb_noise_est_coll_big)
# Simulate a diagonal design
d_rot_diag = get_diag_design(d_rot_big)
rot_diag_merit = calc_merit(d_rot_diag)
aces_data_rot_diag = simulate_aces(d_rot_diag, budget_set; seed = seed)
noise_score_coll_diag = get_noise_score(aces_data_rot_diag, rot_diag_merit)
model_violation_coll_diag = get_model_violation(aces_data_rot_diag)
pretty_print(aces_data_rot_diag, rot_diag_merit)
# Test that the simulations agree sufficiently with the predicted distributions
@testset "Simulation merit agreement" begin
    for noise_score in noise_score_coll_big
        @test is_score_expected(noise_score, z_score_cutoff_lower, z_score_cutoff_upper)
    end
    for idx in eachindex(model_violation_coll_big)
        model_violation = model_violation_coll_big[idx]
        @test z_score_cutoff_lower <= model_violation &&
              model_violation <= z_score_cutoff_upper
        # Ensure the combined design has worse model violation scores
        comb_model_violation = comb_model_violation_coll_big[idx]
        @test model_violation <= comb_model_violation
    end
    for noise_est in aces_data_rot_diag.noise_est_coll
        @test noise_est.gls_unproj_gate_eigenvalues ≈ noise_est.wls_unproj_gate_eigenvalues
        @test noise_est.gls_gate_eigenvalues ≈ noise_est.wls_gate_eigenvalues
    end
    for noise_score in noise_score_coll_diag
        @test is_score_expected(noise_score, z_score_cutoff_lower, z_score_cutoff_upper)
    end
    for model_violation in model_violation_coll_diag
        @test z_score_cutoff_lower <= model_violation &&
              model_violation <= z_score_cutoff_upper
    end
end
# Ensure the estimation is robust to many eigenvalues being clipped
tuple_number = length(d_rot_big.tuple_set)
noise_est = aces_data_rot_big.noise_est_coll[end, budget_number]
shot_budget = noise_est.shot_budget
mapping_lengths = length.(d_rot_big.mapping_ensemble)
mapping_lower = cumsum([0; mapping_lengths[1:(end - 1)]])
M_rot = sum(mapping_lengths)
Random.seed!(seed)
rand_indices = randperm(M_rot)
Random.seed!()
clip_indices = rand_indices[1:clip_number]
clipped_eigenvalues_experiment_ensemble =
    deepcopy(noise_est.eigenvalues_experiment_ensemble)
for i in 1:clip_number
    clip_idx = clip_indices[i]
    tuple_clip_idx = findlast(clip_idx .> mapping_lower)
    mapping_clip_idx = clip_idx - mapping_lower[tuple_clip_idx]
    clipped_eigenvalues_experiment_ensemble[tuple_clip_idx][mapping_clip_idx] .= 0.0
end
covariance_experiment_ensemble = noise_est.covariance_experiment_ensemble
clipped_noise_est = estimate_gate_noise(
    d_rot_big,
    clipped_eigenvalues_experiment_ensemble,
    covariance_experiment_ensemble,
    shot_budget,
)
# Clipping so many eigenvalues probably messes with the marginal and relative z-scores, so they are omitted here
clipped_noise_error = get_noise_error(d_rot_big, clipped_noise_est)
clipped_noise_score = get_noise_score(clipped_noise_error, rot_big_merit)
rot_gls_z_score = noise_score_coll_big[end, budget_number].gls_z_score
clipped_gls_z_score = clipped_noise_score.gls_z_score
# Ensure the estimation is not substantially affected by clipping single eigenvalues
trial_clip_indices = rand_indices[(M_rot - trial_number + 1):M_rot]
trial_gls_z_scores = Vector{Float64}(undef, trial_number)
trial_gls_marginal_z_scores = Vector{Float64}(undef, trial_number)
trial_gls_relative_z_scores = Vector{Float64}(undef, trial_number)
for idx in 1:trial_number
    trial_eigenvalues_experiment_ensemble =
        deepcopy(noise_est.eigenvalues_experiment_ensemble)
    clip_idx = trial_clip_indices[idx]
    tuple_clip_idx = findlast(clip_idx .> mapping_lower)
    mapping_clip_idx = clip_idx - mapping_lower[tuple_clip_idx]
    trial_eigenvalues_experiment_ensemble[tuple_clip_idx][mapping_clip_idx] .= 0.0
    trial_noise_est = estimate_gate_noise(
        d_rot_big,
        trial_eigenvalues_experiment_ensemble,
        covariance_experiment_ensemble,
        shot_budget;
        split = true,
    )
    # Check agreement with the merits
    trial_noise_error = get_noise_error(d_rot_big, trial_noise_est)
    trial_noise_score = get_noise_score(trial_noise_error, rot_big_merit)
    trial_gls_z_scores[idx] = trial_noise_score.gls_z_score
    trial_gls_marginal_z_scores[idx] = trial_noise_score.gls_marginal_z_score
    trial_gls_relative_z_scores[idx] = trial_noise_score.gls_relative_z_score
end
gls_z_score_abs_deviation = abs.(trial_gls_z_scores .- rot_gls_z_score)
rot_gls_marginal_z_score = noise_score_coll_big[end, budget_number].gls_marginal_z_score
gls_marginal_z_score_abs_deviation =
    abs.(trial_gls_marginal_z_scores .- rot_gls_marginal_z_score)
rot_gls_relative_z_score = noise_score_coll_big[end, budget_number].gls_relative_z_score
gls_relative_z_score_abs_deviation =
    abs.(trial_gls_relative_z_scores .- rot_gls_relative_z_score)
@testset "Clipping eigenvalues" begin
    # Clipping many eigenvalues robustness
    @test rot_gls_z_score < clipped_gls_z_score
    @test clipped_gls_z_score < z_score_cutoff_upper + clip_z_cutoff
    # Single eigenvalue clipping robustness
    @test all(gls_z_score_abs_deviation .< trial_z_abs_cutoff)
    @test all(gls_marginal_z_score_abs_deviation .< trial_z_abs_cutoff)
    @test all(gls_relative_z_score_abs_deviation .< trial_z_abs_cutoff)
    @test mean(gls_z_score_abs_deviation) < trial_z_mean_cutoff
    @test mean(gls_marginal_z_score_abs_deviation) < trial_z_mean_cutoff
    @test mean(gls_relative_z_score_abs_deviation) < trial_z_mean_cutoff
end
# Remove the data folder generated for these tests
try
    rm("data")
catch
    @warn "The data folder was not empty; this is not unexpected for CI tests."
    rm("data"; force = true, recursive = true)
end
exit_folder("test")

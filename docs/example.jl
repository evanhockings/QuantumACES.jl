using QuantumACES
# Set up noise models
r_1 = 0.05 / 100
r_2 = 0.4 / 100
r_m = 0.8 / 100
total_std_log = 0.5
seed = UInt(0)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
# Generate the circuit
dist = 3
rotated_param = get_rotated_param(dist)
circuit_dep = get_circuit(rotated_param, dep_param)
circuit_log = get_circuit(rotated_param, log_param)
# Generate a design for this circuit
d = generate_design(circuit_dep)
display(d)
# Optimise a design for this circuit
d = optimise_design(circuit_dep; options = OptimOptions(; seed = seed))
display(d)
# Simulate repeated noise characterisation experiments
d_log = generate_design(circuit_log, d)
d_log = update_noise(d, log_param)
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_log, budget_set; repetitions = repetitions, seed = seed)
# Compare the performance to predictions
merit_log = calc_merit(d_log)
pretty_print(aces_data, merit_log)
# Generate a design for a large circuit
dist_big = 13
rotated_param_big = get_rotated_param(dist_big)
circuit_big = get_circuit(rotated_param_big, dep_param)
circuit_big_log = get_circuit(rotated_param_big, log_param)
d_big = generate_design(circuit_big_log, d; full_covariance = false, diagnostics = true)
# Simulate a large noise characterisation experiment
aces_data_big = simulate_aces(d_big, budget_set; seed = seed, split = true)
# Examine the performance scaling
dist_max = 7
merit_scaling = calc_merit_scaling(d, dist_max)
scaling_fit = get_scaling_fit(merit_scaling)
ensemble_scaling = calc_ensemble_scaling(d_log, dist_max; seed = seed)
ensemble_fit = get_ensemble_fit(ensemble_scaling)
# Compare the performance to predictions
wls_pred_expectation = ensemble_fit.wls_expectation_model(dist_big)
wls_pred_variance = ensemble_fit.wls_variance_model(dist_big)
wls_scores_big = [
    (noise_error.wls_nrmse .- wls_pred_expectation) / sqrt(wls_pred_variance) for
    noise_error in aces_data_big.noise_error_coll[1, :]
]
display(wls_scores_big)
wls_pred_relative_expectation = ensemble_fit.wls_relative_expectation_model(dist_big)
wls_pred_relative_variance = ensemble_fit.wls_relative_variance_model(dist_big)
wls_relative_scores_big = [
    (noise_error.wls_relative_nrmse .- wls_pred_relative_expectation) /
    sqrt(wls_pred_relative_variance) for
    noise_error in aces_data_big.noise_error_coll[1, :]
]
display(wls_relative_scores_big)
# Simulate a memory experiment
big_rounds = dist_big
big_shots = 5 * 10^6
decoder_gate_probabilities = [
    circuit_big_log.gate_probabilities
    circuit_big.gate_probabilities
    [noise_est.wls_gate_probabilities for noise_est in aces_data_big.noise_est_coll[1, :]]
]
decoder_labels = [
    "True"
    "Depolarising"
    ["ACES S=$(budget)" for budget in budget_set]
]
big_memory_data = simulate_memory(
    circuit_big_log,
    big_rounds,
    big_shots;
    seed = seed,
    decoder_gate_probabilities = decoder_gate_probabilities,
    decoder_labels = decoder_labels,
    diagnostics = true,
)
big_memory_summary = get_memory_summary(big_memory_data)
# Create a randomised design
min_randomisations = 64
target_shot_budget = 5 * 10^6
experiment_shots = 64
d_rand = generate_rand_design(
    d_log,
    min_randomisations,
    target_shot_budget,
    experiment_shots;
    seed = seed,
)
# Calculate the figure of merit of the randomised design
d_shot = get_design(d_rand)
merit_shot = calc_merit(d_shot)
# Generate circuits
qiskit_qubit_num = 17
qiskit_qubit_map = collect(0:(qiskit_qubit_num - 1))
(stim_ensemble, qiskit_ensemble) =
    get_stim_qiskit_ensemble(d_rand, qiskit_qubit_num, qiskit_qubit_map)
# Simulate the randomised circuits
simulate_stim_ensemble(d_rand, stim_ensemble, experiment_shots; seed = seed)
rand_noise_est = estimate_stim_ensemble(d_rand, experiment_shots; simulation_seed = seed)
rand_noise_error = get_noise_error(d_rand, rand_noise_est)
rand_noise_score = get_noise_score(rand_noise_error, merit_shot)
# Example hypothetical Qiskit results folder name
backend = "backend"
d_rand_filename = rand_design_filename(d_rand)
@assert d_rand_filename[(end - 4):end] == ".jld2"
qiskit_results_folder = "data/$(backend)_$(d_rand_filename[1:(end - 5)])"
# Example hypothetical Qiskit job name
prefix = "job"
example_job_1_filename = "$(qiskit_results_folder)/$(prefix)_1.pickle"
# Hypothetically process and estimate the Qiskit data
#=
process_qiskit_ensemble(
    d_rand,
    qiskit_qubit_num,
    qiskit_qubit_map,
    experiment_shots;
    backend = backend,
    prefix = prefix,
)
noise_est =
    estimate_qiskit_ensemble(d_rand, qiskit_qubit_map, experiment_shots; backend = backend)
=#
# Hypothetically examine the model violation for the noise model and combined noise model
#=
model_violation = get_model_violation(d_shot, noise_est)
d_comb = get_combined_design(d_shot)
comb_noise_est = estimate_gate_noise(d_comb, noise_est)
comb_model_violation = get_model_violation(d_comb, comb_noise_est)
aic = get_aic(d_shot, noise_est)
bic = get_bic(d_shot, noise_est)
comb_aic = get_aic(d_comb, comb_noise_est)
comb_bic = get_bic(d_comb, comb_noise_est)
=#

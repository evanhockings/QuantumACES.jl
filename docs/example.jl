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
d_log = update_noise(d, log_param)
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_log, budget_set; repetitions = repetitions, seed = seed)
# Compare the performance to predictions
merit_log = calc_merit(d_log)
display(merit_log)
pretty_print(aces_data, merit_log)
noise_score_coll = get_noise_score(aces_data, merit_log)
gls_z_scores = [noise_score.gls_z_score for noise_score in noise_score_coll]
display(gls_z_scores)
# Examine the performance scaling
dist_max = 7
merit_scaling = calc_merit_scaling(d, dist_max)
ensemble_scaling = calc_ensemble_scaling(d_log, dist_max; seed = seed)
# Generate a design for a large circuit
dist_big = 13
rotated_param_big = get_rotated_param(dist_big)
circuit_big = get_circuit(rotated_param_big, log_param)
d_big = generate_design(
    circuit_big,
    d.tuple_set_data;
    shot_weights = d.shot_weights,
    full_covariance = false,
    diagnostics = true,
)
# Simulate a large noise characterisation experiment
aces_data_big = simulate_aces(d_big, budget_set; seed = seed, split = true)
# Compare the performance to predictions
scaling_fit = get_scaling_fit(merit_scaling)
ensemble_fit = get_ensemble_fit(ensemble_scaling)
wls_pred_expectation = ensemble_fit.wls_expectation_model(dist_big)
wls_pred_variance = ensemble_fit.wls_variance_model(dist_big)
wls_z_scores_big = [
    (noise_error.wls_nrmse .- wls_pred_expectation) / sqrt(wls_pred_variance) for
    noise_error in aces_data_big.noise_error_coll[1, :]
]
display(wls_z_scores_big)
wls_pred_relative_expectation = ensemble_fit.wls_relative_expectation_model(dist_big)
wls_pred_relative_variance = ensemble_fit.wls_relative_variance_model(dist_big)
wls_relative_z_scores_big = [
    (noise_error.wls_relative_nrmse .- wls_pred_relative_expectation) /
    sqrt(wls_pred_relative_variance) for
    noise_error in aces_data_big.noise_error_coll[1, :]
]
display(wls_relative_z_scores_big)
# Create a randomised design
min_randomisations = 40
target_shot_budget = 10^7
experiment_shots = 512
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
display(d_shot)
display(merit_shot)
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
# Simulate a memory experiment using the true, estimated, and depolarising noise
rounds = dist
shots = 10^6
decoder_gate_probabilities = [
    circuit_log.gate_probabilities,
    rand_noise_est.gls_gate_probabilities,
    circuit_dep.gate_probabilities,
]
memory_data = simulate_memory(
    circuit_log,
    rounds,
    shots;
    seed = seed,
    decoder_gate_probabilities = decoder_gate_probabilities,
)
display(memory_data)
# Check the code distances
(z_dist, x_dist) = calc_memory_distances(circuit_dep)

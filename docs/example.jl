using QuantumACES

# Set up noise models
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
# Generate the circuit
dist = 3
rotated_param = get_rotated_param(dist)
rotated_planar = get_circuit(rotated_param, dep_param)
# Optimise a design for this circuit
d = optimise_design(rotated_planar; options = OptimOptions(; ls_type = :wls, seed = seed))
# Simulate repeated noise characterisation experiments
d_log = update_noise(d, log_param)
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_log, budget_set; repetitions = repetitions)
# Compare the performance to predictions for the largest measurement budget
wls_merit_log = calc_wls_merit(d_log)
fgls_z_scores =
    (aces_data.fgls_gate_norm_coll[:, 3] .- wls_merit_log.expectation) /
    sqrt(wls_merit_log.variance)
# Examine the performance scaling
dist_max = 9
dep_planar_scaling = calc_depolarising_planar_scaling(d, dist_max; ls_type = :wls)
log_planar_scaling =
    calc_lognormal_planar_scaling(d_log, dist_max; ls_type = :wls, seed = seed)
# Generate a design for a large circuit
dist_big = 13
rotated_param_big = get_rotated_param(dist_big)
rotated_planar_big = get_circuit(rotated_param_big, log_param)
d_big = generate_design(rotated_planar_big, d.tuple_set_data)
# Simulate a large noise characterisation experiment
budget_set_big = [10^6; 10^7; 10^8; 10^9]
aces_data_big = simulate_aces(d_big, budget_set_big)
# Compare the performance to predictions
pred_expectation = log_planar_scaling.expectation_fit(dist_big)
pred_variance = log_planar_scaling.variance_fit(dist_big)
fgls_z_scores_big =
    (aces_data_big.fgls_gate_norm_coll[1, :] .- pred_expectation) / sqrt(pred_variance)

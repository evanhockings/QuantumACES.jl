using AveragedCircuitEigenvalueSampling, Test
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
unrotated_param = get_unrotated_param(dist)
big_rotated_param = get_rotated_param(vertical_dist, horizontal_dist)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
rotated_planar = get_circuit(rotated_param, dep_param)
unrotated_planar = get_circuit(unrotated_param, dep_param)
big_rotated_planar = get_circuit(big_rotated_param, log_param)
# Set up optimisation parameters
rot_ls_type = :wls
unrot_ls_type = :gls
rot_max_steps = 10
rot_max_cycles = 0
unrot_max_cycles = 1
excursion_number = 2
excursion_length = 2
rot_tuple_number = 2 * length(rotated_planar.unique_layer_indices) + 3
unrot_tuple_number = 2 * length(unrotated_planar.unique_layer_indices) + 2
# Set up scaling parameters
dist_max = 5
precision = 1.0
min_repetitions = 5
# Set up simulation parameters
shots_set = [10^7; 2 * 10^7]
repetitions = 2
max_samples = 10^6
z_score_cutoff = 4.0
# Optimise and simulate a design for a rotated planar code
@testset "Rotated planar design" begin
    # Set up a basic design
    rot_basic = get_basic_tuple_set(rotated_planar)
    d_rot_basic = generate_design(rotated_planar, rot_basic; save_data = true)
    # Test saving and loading and deletion
    d_rot_basic_load = load_design(
        d_rot_basic.c.circuit_param,
        d_rot_basic.c.noise_param,
        length(d_rot_basic.tuple_set),
        d_rot_basic.tuple_set_data.repeat_numbers,
        d_rot_basic.full_covariance,
        d_rot_basic.ls_type,
    )
    @test d_rot_basic_load == d_rot_basic
    delete_design(d_rot_basic)
    rot_basic_merit = calc_wls_merit(d_rot_basic)
    # Optimise a design
    d_rot_opt = optimise_design(
        rotated_planar;
        options = OptimOptions(;
            ls_type = rot_ls_type,
            save_data = true,
            max_steps = rot_max_steps,
            max_cycles = rot_max_cycles,
            excursion_number = excursion_number,
            excursion_length = excursion_length,
            max_tuple_number = rot_tuple_number,
            seed = seed,
        ),
    )
    rot_opt_merit = calc_wls_merit(d_rot_opt)
    @test rot_opt_merit.expectation < rot_basic_merit.expectation
    # Test saving and loading and deletion
    d_rot_opt_load = load_design(
        d_rot_opt.c.circuit_param,
        d_rot_opt.c.noise_param,
        length(d_rot_opt.tuple_set),
        d_rot_opt.tuple_set_data.repeat_numbers,
        d_rot_opt.full_covariance,
        d_rot_opt.ls_type,
    )
    @test d_rot_opt_load == d_rot_opt
    delete_design(d_rot_opt)
    # Examine the scaling of the merit for depolarising and log-normal noise
    dep_planar_scaling =
        calc_depolarising_planar_scaling(d_rot_opt, dist_max; save_data = true)
    d_rot_opt_log = update_noise(d_rot_opt, log_param)
    log_planar_scaling = calc_lognormal_planar_scaling(
        d_rot_opt_log,
        dist_max;
        precision = precision,
        min_repetitions = min_repetitions,
        save_data = true,
    )
    # Test saving and loading and deletion
    dep_planar_scaling_load = load_scaling(d_rot_opt, d_rot_opt.ls_type)
    @test dep_planar_scaling_load.merit_scaling == dep_planar_scaling.merit_scaling
    delete_scaling(dep_planar_scaling)
    log_planar_scaling_load = load_scaling(d_rot_opt_log, d_rot_opt_log.ls_type)
    @test log_planar_scaling_load.expectation_scaling ==
          log_planar_scaling.expectation_scaling
    @test log_planar_scaling_load.variance_scaling == log_planar_scaling.variance_scaling
    @test log_planar_scaling_load.eigenvalues_scaling ==
          log_planar_scaling.eigenvalues_scaling
    delete_scaling(log_planar_scaling)
    # Transfer the tuple set to a differently-sized code
    d_rot_big = generate_design(
        big_rotated_planar,
        d_rot_opt.tuple_set_data;
        shot_weights = d_rot_opt.shot_weights,
        diagnostics = true,
        save_data = true,
    )
    # Test that the merit is improved by better LS estimators
    rot_merit_set = calc_merit_set(d_rot_big)
    (rot_gls_merit, rot_wls_merit, rot_ols_merit) = rot_merit_set
    @test rot_gls_merit.expectation <= rot_wls_merit.expectation
    @test rot_wls_merit.expectation <= rot_ols_merit.expectation
    # Simulate the design
    aces_data_rot_big = simulate_aces(
        d_rot_big,
        shots_set;
        seed = seed,
        max_samples = max_samples,
        save_data = true,
        clear_design = true,
    )
    aces_data_rot_big = simulate_aces(
        d_rot_big,
        shots_set;
        repetitions = repetitions,
        seed = seed,
        max_samples = max_samples,
        detailed_diagnostics = true,
        save_data = false,
    )
    pretty_print(aces_data_rot_big, rot_merit_set)
    delete_aces(aces_data_rot_big)
    # Test that the simulations agree sufficiently with the predicted distributions
    rot_big_gls_z_scores =
        (aces_data_rot_big.fgls_gate_norm_coll .- rot_gls_merit.expectation) ./
        sqrt(rot_gls_merit.variance)
    rot_big_wls_z_scores =
        (aces_data_rot_big.wls_gate_norm_coll .- rot_wls_merit.expectation) ./
        sqrt(rot_wls_merit.variance)
    rot_big_ols_z_scores =
        (aces_data_rot_big.ols_gate_norm_coll .- rot_ols_merit.expectation) ./
        sqrt(rot_ols_merit.variance)
    @test all(abs.(rot_big_gls_z_scores) .< z_score_cutoff)
    @test all(abs.(rot_big_wls_z_scores) .< z_score_cutoff)
    @test all(abs.(rot_big_ols_z_scores) .< z_score_cutoff)
    # Remove the data folder generated for these tests
    # This tests the removal functions, as it will fail if the folder is not empty
    rm("data")
end
# Optimise and simulate a design for a unrotated planar code
@testset "Unrotated planar design" begin
    # Set up a basic design
    unrot_basic = get_basic_tuple_set(unrotated_planar)
    d_unrot_basic = generate_design(unrotated_planar, unrot_basic)
    unrot_basic_merit = calc_gls_merit(d_unrot_basic)
    # Optimise a design
    d_unrot_opt = optimise_design(
        unrotated_planar;
        options = OptimOptions(;
            ls_type = unrot_ls_type,
            save_data = true,
            max_cycles = unrot_max_cycles,
            excursion_number = excursion_number,
            excursion_length = excursion_length,
            max_tuple_number = unrot_tuple_number,
            seed = seed,
        ),
    )
    unrot_opt_merit = calc_gls_merit(d_unrot_opt)
    @test unrot_opt_merit.expectation < unrot_basic_merit.expectation
    # Update the noise on the design
    d_unrot_opt_log = update_noise(d_unrot_opt, log_param)
    unrot_merit_set = calc_merit_set(d_unrot_opt_log)
    (unrot_gls_merit, unrot_wls_merit, unrot_ols_merit) = unrot_merit_set
    @test unrot_gls_merit.expectation <= unrot_wls_merit.expectation
    @test unrot_wls_merit.expectation <= unrot_ols_merit.expectation
    @test unrot_gls_merit.expectation != unrot_opt_merit.expectation
    @test unrot_gls_merit.variance != unrot_opt_merit.variance
    @test unrot_gls_merit.cond_num ≈ unrot_opt_merit.cond_num
    @test unrot_gls_merit.pinv_norm ≈ unrot_opt_merit.pinv_norm
    # Examine the scaling of the merit for depolarising and log-normal noise
    dep_planar_scaling = calc_depolarising_planar_scaling(d_unrot_opt, dist_max)
    log_planar_scaling = calc_lognormal_planar_scaling(
        d_unrot_opt_log,
        dist_max;
        precision = precision,
        min_repetitions = min_repetitions,
    )
    # Simulate the design
    aces_data_unrot_log =
        simulate_aces(d_unrot_opt_log, shots_set; repetitions = repetitions, seed = seed)
    pretty_print(aces_data_unrot_log, unrot_merit_set)
    # Test that the simulations agree sufficiently with the predicted distributions
    unrot_log_gls_z_scores =
        (aces_data_unrot_log.fgls_gate_norm_coll .- unrot_gls_merit.expectation) ./
        sqrt(unrot_gls_merit.variance)
    unrot_log_wls_z_scores =
        (aces_data_unrot_log.wls_gate_norm_coll .- unrot_wls_merit.expectation) ./
        sqrt(unrot_wls_merit.variance)
    unrot_log_ols_z_scores =
        (aces_data_unrot_log.ols_gate_norm_coll .- unrot_ols_merit.expectation) ./
        sqrt(unrot_ols_merit.variance)
    @test all(abs.(unrot_log_gls_z_scores) .< z_score_cutoff)
    @test all(abs.(unrot_log_wls_z_scores) .< z_score_cutoff)
    @test all(abs.(unrot_log_ols_z_scores) .< z_score_cutoff)
end
# Make sure to delete the data folder even if previous tests failed
rm("data"; force = true, recursive = true)
exit_folder("test")

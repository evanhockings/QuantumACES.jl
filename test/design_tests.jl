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
rotated_param = RotatedPlanarParameters(dist)
unrotated_param = UnrotatedPlanarParameters(dist)
big_rotated_param = RotatedPlanarParameters(vertical_dist, horizontal_dist)
dep_param = DepolarisingParameters(r_1, r_2, r_m)
log_param = LogNormalParameters(r_1, r_2, r_m, total_std_log; seed = seed)
rotated_planar = Code(rotated_param, dep_param)
unrotated_planar = Code(unrotated_param, dep_param)
big_rotated_planar = Code(big_rotated_param, log_param)
# Set up optimisation parameters
rot_ls_type = :wls
unrot_ls_type = :gls
rot_max_cycles = 0
unrot_max_cycles = 1
rot_tuple_num = 2 * length(rotated_planar.unique_layer_indices) + 3
unrot_tuple_num = 2 * length(unrotated_planar.unique_layer_indices) + 2
excursion_len = 2
excursion_num = 2
rot_max_steps = 5
# Set up scaling parameters
dist_max = 4
precision = 1.0
min_repetitions = 5
# Set up simulation parameters
shots_set = [10^7; 2 * 10^7]
repetitions = 2
max_samples = 10^6
z_score_cutoff = 4.0
# Optimise and simulate a design for a rotated planar code
@testset "Rotated planar design" begin
    # Set up a trivial design
    rot_trivial = TrivialTupleSet(rotated_planar)
    d_rot_trivial = GenerateDesign(rotated_planar, rot_trivial; save_data = true)
    # Test saving and loading and deletion
    d_rot_trivial_load = load_design(
        d_rot_trivial.code.code_param,
        d_rot_trivial.code.noise_param,
        length(d_rot_trivial.tuple_set),
        d_rot_trivial.tuple_set_data.repeat_numbers,
        d_rot_trivial.full_covariance,
    )
    @test d_rot_trivial_load == d_rot_trivial
    delete_design(d_rot_trivial)
    rot_trivial_merit = WLSMerit(d_rot_trivial)
    # Optimise a design
    d_rot_opt = OptimiseDesign(
        rotated_planar;
        ls_type = rot_ls_type,
        max_cycles = rot_max_cycles,
        tuple_num = rot_tuple_num,
        excursion_len = excursion_len,
        excursion_num = excursion_num,
        max_steps = rot_max_steps,
        seed = seed,
        save_data = true,
    )
    rot_opt_merit = WLSMerit(d_rot_opt)
    @test rot_opt_merit.expectation < rot_trivial_merit.expectation
    # Test saving and loading and deletion
    d_rot_opt_load = load_design(
        d_rot_opt.code.code_param,
        d_rot_opt.code.noise_param,
        length(d_rot_opt.tuple_set),
        d_rot_opt.tuple_set_data.repeat_numbers,
        d_rot_opt.full_covariance,
    )
    @test d_rot_opt_load == d_rot_opt
    delete_design(d_rot_opt)
    # Examine the scaling of the merit for depolarising and log-normal noise
    dep_scaling_data = DepolarisingScaling(d_rot_opt, dist_max; save_data = true)
    d_rot_opt_log = Update(d_rot_opt, log_param)
    log_scaling_data = LogNormalScaling(
        d_rot_opt_log,
        dist_max;
        precision = precision,
        min_repetitions = min_repetitions,
        save_data = true,
    )
    # Test saving and loading and deletion
    dep_scaling_data_load = load_scaling(d_rot_opt, d_rot_opt.ls_type, :dep)
    @test dep_scaling_data_load == dep_scaling_data
    delete_scaling(dep_scaling_data)
    log_scaling_data_load = load_scaling(d_rot_opt_log, d_rot_opt_log.ls_type, :log)
    @test log_scaling_data_load == log_scaling_data
    delete_scaling(log_scaling_data)
    # Transfer the tuple set to a differently-sized code
    d_rot_big = GenerateDesign(
        big_rotated_planar,
        d_rot_opt.tuple_set_data;
        shot_weights = d_rot_opt.shot_weights,
        save_data = true,
    )
    # Test that the merit is improved by better LS estimators
    rot_merit_set = MeritSet(d_rot_big)
    (rot_gls_merit, rot_wls_merit, rot_ols_merit) = rot_merit_set
    @test rot_gls_merit.expectation <= rot_wls_merit.expectation
    @test rot_wls_merit.expectation <= rot_ols_merit.expectation
    # Simulate the design
    aces_data_rot_big = SimulateACES(
        d_rot_big,
        shots_set;
        seed = seed,
        max_samples = max_samples,
        save_data = true,
        clear_design = true,
    )
    aces_data_rot_big = SimulateACES(
        d_rot_big,
        shots_set;
        repetitions = repetitions,
        seed = seed,
        max_samples = max_samples,
        save_data = false,
    )
    PrettyPrint(aces_data_rot_big, rot_merit_set)
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
    # Set up a trivial design
    unrot_trivial = TrivialTupleSet(unrotated_planar)
    d_unrot_trivial = GenerateDesign(unrotated_planar, unrot_trivial)
    unrot_trivial_merit = GLSMerit(d_unrot_trivial)
    # Optimise a design
    d_unrot_opt = OptimiseDesign(
        unrotated_planar;
        ls_type = unrot_ls_type,
        max_cycles = unrot_max_cycles,
        tuple_num = unrot_tuple_num,
        excursion_len = excursion_len,
        excursion_num = excursion_num,
        seed = seed,
    )
    unrot_opt_merit = GLSMerit(d_unrot_opt)
    @test unrot_opt_merit.expectation < unrot_trivial_merit.expectation
    # Update the noise on the design
    d_unrot_opt_log = Update(d_unrot_opt, log_param)
    unrot_merit_set = MeritSet(d_unrot_opt_log)
    (unrot_gls_merit, unrot_wls_merit, unrot_ols_merit) = unrot_merit_set
    @test unrot_gls_merit.expectation <= unrot_wls_merit.expectation
    @test unrot_wls_merit.expectation <= unrot_ols_merit.expectation
    @test unrot_gls_merit.expectation != unrot_opt_merit.expectation
    @test unrot_gls_merit.variance != unrot_opt_merit.variance
    @test unrot_gls_merit.cond_num ≈ unrot_opt_merit.cond_num
    @test unrot_gls_merit.pinv_norm ≈ unrot_opt_merit.pinv_norm
    # Examine the scaling of the merit for depolarising and log-normal noise
    dep_scaling_data = DepolarisingScaling(d_unrot_opt, dist_max)
    log_scaling_data = LogNormalScaling(
        d_unrot_opt_log,
        dist_max;
        precision = precision,
        min_repetitions = min_repetitions,
    )
    # Simulate the design
    aces_data_unrot_log =
        SimulateACES(d_unrot_opt_log, shots_set; repetitions = repetitions, seed = seed)
    PrettyPrint(aces_data_unrot_log, unrot_merit_set)
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

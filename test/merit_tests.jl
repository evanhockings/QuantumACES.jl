using QuantumACES, LinearAlgebra, ForwardDiff, Random, Distributions, Test
# Set up codes
dist = 3
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
rotated_param = get_rotated_param(dist)
unrotated_param = get_unrotated_param(dist)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
rotated_planar = get_circuit(rotated_param, log_param)
unrotated_planar = get_circuit(unrotated_param, dep_param)
# Set up designs
rot_basic = get_basic_tuple_set(rotated_planar)
rot_tuple_set = [[rotated_planar.circuit_tuple]; rot_basic]
d_rot = generate_design(rotated_planar, rot_tuple_set)
unrot_basic = get_basic_tuple_set(unrotated_planar)
unrot_tuple_set = [[unrotated_planar.circuit_tuple]; unrot_basic]
d_unrot = generate_design(unrotated_planar, unrot_tuple_set)
# Test that the basic experiment numbers are correctly generated
@testset "Basic experiment numbers" begin
    d_rot_basic = generate_design(rotated_planar, rot_basic)
    d_unrot_basic = generate_design(unrotated_planar, unrot_basic)
    @test d_rot_basic.experiment_numbers ==
          QuantumACES.get_basic_experiment_numbers(rotated_planar)
    @test d_unrot_basic.experiment_numbers ==
          QuantumACES.get_basic_experiment_numbers(unrotated_planar)
end
# Test that we can generate codes with a range of different parameters
test_param_1 = get_rotated_param(
    dist + 2,
    dist;
    dynamically_decouple = false,
    single_qubit_time = 1.0,
    two_qubit_time = 2.0,
    dynamical_decoupling_time = 1.0,
    meas_reset_time = 20.0,
)
test_code_1 = get_circuit(test_param_1, log_param)
test_param_2 = get_rotated_param(
    dist,
    dist + 1;
    check_type = :standard,
    gate_type = :cx,
    dynamically_decouple = false,
    pad_identity = false,
)
test_code_2 = get_circuit(test_param_2, log_param)
test_param_3 = get_unrotated_param(
    dist,
    dist + 1;
    pad_identity = false,
    single_qubit_time = 1.0,
    two_qubit_time = 2.0,
    meas_reset_time = 20.0,
)
test_code_3 = get_circuit(test_param_3, dep_param)
# Set up gradient descent parameters
max_steps = 10
rot_covariance_log = calc_covariance_log(d_rot)
N_rot = rotated_planar.N
C_rot = length(d_rot.tuple_set)
rot_mapping_lengths = length.(d_rot.mapping_ensemble)
rot_gate_eigenvalues_diag = Diagonal(d_rot.c.gate_eigenvalues)
unrot_covariance_log = calc_covariance_log(d_unrot)
N_unrot = unrotated_planar.N
C_unrot = length(d_unrot.tuple_set)
unrot_mapping_lengths = length.(d_unrot.mapping_ensemble)
unrot_gate_eigenvalues_diag = Diagonal(d_unrot.c.gate_eigenvalues)
# Test the merit gradient for GLS
@testset "GLS merit gradient" begin
    # Test GLS gradient descent
    (d_rot_gls, rot_covariance_log_gls, rot_merit_descent_gls) = gls_optimise_weights(
        d_rot,
        rot_covariance_log;
        options = OptimOptions(; ls_type = :gls, max_steps = max_steps),
    )
    @test d_rot_gls.shot_weights != d_rot.shot_weights
    @test rot_covariance_log_gls != rot_covariance_log
    # Test the covariance matrix output by gradient descent is correct
    rot_gls_opt_merit = calc_gls_merit(d_rot_gls)
    rot_gls_opt_cov_eigvals =
        eigvals(QuantumACES.calc_gls_covariance(d_rot_gls, rot_covariance_log_gls))
    (rot_gls_opt_expectation, rot_gls_opt_variance) =
        QuantumACES.nrmse_moments(rot_gls_opt_cov_eigvals)
    @test rot_gls_opt_merit.eigenvalues ≈ rot_gls_opt_cov_eigvals
    @test rot_gls_opt_merit.expectation ≈ rot_gls_opt_expectation
    @test rot_gls_opt_merit.variance ≈ rot_gls_opt_variance
    # Check that gradient descent improved the figure of merit
    rot_gls_unopt_expectation = calc_gls_moments(d_rot, rot_covariance_log)[1]
    @test rot_gls_opt_expectation < rot_gls_unopt_expectation
    # Test the GLS gradient
    gls_1 = time()
    rot_gls_shot_weights = d_rot_gls.shot_weights
    rot_gls_shot_weights_factor_inv = QuantumACES.get_shot_weights_factor_inv(
        rot_gls_shot_weights,
        d_rot_gls.tuple_times,
        rot_mapping_lengths,
    )
    rot_covariance_log_gls_unweighted =
        rot_covariance_log_gls * rot_gls_shot_weights_factor_inv
    rot_covariance_log_gls_unweighted_inv = QuantumACES.sparse_covariance_inv(
        rot_covariance_log_gls_unweighted,
        rot_mapping_lengths,
    )
    (gls_expectation_grad_log, gls_expectation) = QuantumACES.calc_gls_merit_grad_log(
        d_rot_gls,
        rot_gls_shot_weights,
        rot_covariance_log_gls_unweighted_inv,
    )
    gls_2 = time()
    # Compare against ForwardDiff
    function DifferentiableGLSExpectation(log_shot_weights)
        # Scale the covariance matrix by the shot weights
        shot_weights = exp.(-log_shot_weights) / sum(exp.(-log_shot_weights))
        tuple_times_factor = sum(shot_weights .* d_rot_gls.tuple_times)
        rot_shot_weights_factor_inv =
            (1 / tuple_times_factor) * Diagonal(
                vcat(
                    [
                        shot_weights[idx] * ones(rot_mapping_lengths[idx]) for
                        idx in 1:C_rot
                    ]...,
                ),
            )
        rot_covariance_log_inv_scaled =
            rot_covariance_log_gls_unweighted_inv * rot_shot_weights_factor_inv
        # Calculate the covariance matrix of the gate eigenvalues, including the first-order Taylor approximation
        rot_gls_gate_eigenvalues_cov_scaled = Symmetric(
            rot_gate_eigenvalues_diag *
            inv(
                cholesky(
                    Symmetric(
                        Array(d_rot.matrix' * rot_covariance_log_inv_scaled * d_rot.matrix),
                    ),
                ),
            ) *
            rot_gate_eigenvalues_diag,
        )
        # Calculate the expectation
        cov_trace = tr(rot_gls_gate_eigenvalues_cov_scaled)
        cov_sq_trace = tr(rot_gls_gate_eigenvalues_cov_scaled^2)
        expectation = sqrt(cov_trace) * (1 - cov_sq_trace / (4 * cov_trace^2)) / sqrt(N_rot)
        return expectation
    end
    rot_gls_log_shot_weights = -log.(rot_gls_shot_weights)
    rot_gls_log_shot_weights .-= minimum(rot_gls_log_shot_weights)
    gls_expectation_grad_log_test =
        ForwardDiff.gradient(DifferentiableGLSExpectation, rot_gls_log_shot_weights)
    gls_3 = time()
    @test gls_3 - gls_2 > gls_2 - gls_1
    @test gls_expectation_grad_log ≈ gls_expectation_grad_log_test
    println(
        "The GLS gradient took $(round(gls_2 - gls_1, digits=3)) s to compute, whereas ForwardDiff took $(round(gls_3 - gls_2, digits=3)) s.",
    )
    println(
        "The optimised GLS figure of merit is $(round(rot_gls_opt_expectation, digits=3)), whereas the unoptimised GLS figure of merit is $(round(rot_gls_unopt_expectation, digits=3)).",
    )
end
# Test the merit gradient for WLS
@testset "WLS merit gradient" begin
    # Test WLS gradient descent
    (d_rot_wls, rot_covariance_log_wls, rot_merit_descent_wls) = wls_optimise_weights(
        d_rot,
        rot_covariance_log;
        options = OptimOptions(; ls_type = :wls, max_steps = max_steps),
    )
    @test d_rot_wls.shot_weights != d_rot.shot_weights
    @test rot_covariance_log_wls != rot_covariance_log
    # Check the covariance matrix output by gradient descent yields the correct quantities
    rot_wls_opt_merit = calc_wls_merit(d_rot_wls)
    rot_wls_opt_cov_eigvals =
        eigvals(QuantumACES.calc_wls_covariance(d_rot_wls, rot_covariance_log_wls))
    (rot_wls_opt_expectation, rot_wls_opt_variance) =
        QuantumACES.nrmse_moments(rot_wls_opt_cov_eigvals)
    @test rot_wls_opt_merit.eigenvalues ≈ rot_wls_opt_cov_eigvals
    @test rot_wls_opt_merit.expectation ≈ rot_wls_opt_expectation
    @test rot_wls_opt_merit.variance ≈ rot_wls_opt_variance
    # Check that gradient descent improved the figure of merit
    rot_wls_unopt_expectation = calc_wls_moments(d_rot, rot_covariance_log)[1]
    @test rot_wls_opt_expectation < rot_wls_unopt_expectation
    # Test the WLS gradient
    wls_1 = time()
    rot_wls_shot_weights = d_rot_wls.shot_weights
    rot_wls_shot_weights_factor_inv = QuantumACES.get_shot_weights_factor_inv(
        rot_wls_shot_weights,
        d_rot_wls.tuple_times,
        rot_mapping_lengths,
    )
    rot_covariance_log_wls_unweighted =
        rot_covariance_log_wls * rot_wls_shot_weights_factor_inv
    (wls_expectation_grad_log, wls_expectation) = QuantumACES.calc_wls_merit_grad_log(
        d_rot_wls,
        rot_wls_shot_weights,
        rot_covariance_log_wls_unweighted,
    )
    wls_2 = time()
    # Compare against ForwardDiff
    function DifferentiableWLSExpectation(log_shot_weights)
        # Scale the covariance matrix by the shot weights
        shot_weights = exp.(-log_shot_weights) / sum(exp.(-log_shot_weights))
        tuple_times_factor = sum(shot_weights .* d_rot_wls.tuple_times)
        rot_shot_weights_factor =
            tuple_times_factor * Diagonal(
                vcat(
                    [
                        (1 / shot_weights[idx]) * ones(rot_mapping_lengths[idx]) for
                        idx in 1:C_rot
                    ]...,
                ),
            )
        rot_covariance_log_scaled =
            rot_covariance_log_wls_unweighted * rot_shot_weights_factor
        rot_covariance_log_diag_inv = Diagonal(rot_covariance_log_scaled)^(-1)
        # Calculate the covariance matrix of the gate eigenvalues, including the first-order Taylor approximation
        rot_wls_estimator =
            rot_gate_eigenvalues_diag *
            inv(
                cholesky(
                    Symmetric(
                        Array(d_rot.matrix' * rot_covariance_log_diag_inv * d_rot.matrix),
                    ),
                ),
            ) *
            d_rot.matrix' *
            rot_covariance_log_diag_inv
        rot_wls_gate_eigenvalues_cov_scaled =
            Symmetric(rot_wls_estimator * rot_covariance_log_scaled * rot_wls_estimator')
        # Calculate the expectation
        cov_trace = tr(rot_wls_gate_eigenvalues_cov_scaled)
        cov_sq_trace = tr(rot_wls_gate_eigenvalues_cov_scaled^2)
        expectation = sqrt(cov_trace) * (1 - cov_sq_trace / (4 * cov_trace^2)) / sqrt(N_rot)
        return expectation
    end
    rot_wls_log_shot_weights = -log.(rot_wls_shot_weights)
    rot_wls_log_shot_weights .-= minimum(rot_wls_log_shot_weights)
    wls_expectation_grad_log_test =
        ForwardDiff.gradient(DifferentiableWLSExpectation, rot_wls_log_shot_weights)
    wls_3 = time()
    @test wls_3 - wls_2 > wls_2 - wls_1
    @test wls_expectation_grad_log ≈ wls_expectation_grad_log_test
    println(
        "The WLS gradient took $(round(wls_2 - wls_1, digits=3)) s to compute, whereas ForwardDiff took $(round(wls_3 - wls_2, digits=3)) s.",
    )
    println(
        "The optimised WLS figure of merit is $(round(rot_wls_opt_expectation, digits=3)), whereas the unoptimised WLS figure of merit is $(round(rot_wls_unopt_expectation, digits=3)).",
    )
end
# Test the merit gradient for OLS
@testset "OLS merit gradient" begin
    # Test OLS gradient descent
    (d_unrot_ols, unrot_covariance_log_ols, unrot_merit_descent_ols) = ols_optimise_weights(
        d_unrot,
        unrot_covariance_log;
        options = OptimOptions(; ls_type = :ols, max_steps = max_steps),
    )
    @test d_unrot_ols.shot_weights != d_unrot.shot_weights
    @test unrot_covariance_log_ols != unrot_covariance_log
    # Check the covariance matrix output by gradient descent yields the correct quantities
    unrot_ols_opt_merit = calc_ols_merit(d_unrot_ols)
    unrot_ols_opt_cov_eigvals =
        eigvals(QuantumACES.calc_ols_covariance(d_unrot_ols, unrot_covariance_log_ols))
    (unrot_ols_opt_expectation, unrot_ols_opt_variance) =
        QuantumACES.nrmse_moments(unrot_ols_opt_cov_eigvals)
    @test unrot_ols_opt_merit.eigenvalues ≈ unrot_ols_opt_cov_eigvals
    @test unrot_ols_opt_merit.expectation ≈ unrot_ols_opt_expectation
    @test unrot_ols_opt_merit.variance ≈ unrot_ols_opt_variance
    # Check that gradient descent improved the figure of merit
    unrot_ols_unopt_expectation = calc_ols_moments(d_unrot, unrot_covariance_log)[1]
    @test unrot_ols_opt_expectation < unrot_ols_unopt_expectation
    # Test the OLS gradient
    ols_1 = time()
    unrot_ols_shot_weights = d_unrot_ols.shot_weights
    unrot_ols_shot_weights_factor_inv = QuantumACES.get_shot_weights_factor_inv(
        unrot_ols_shot_weights,
        d_unrot_ols.tuple_times,
        unrot_mapping_lengths,
    )
    unrot_covariance_log_ols_unweighted =
        unrot_covariance_log_ols * unrot_ols_shot_weights_factor_inv
    unrot_ols_estimator =
        unrot_gate_eigenvalues_diag *
        inv(bunchkaufman(Symmetric(Array(d_unrot.matrix' * d_unrot.matrix)))) *
        d_unrot.matrix'
    unrot_ols_estimator_covariance =
        unrot_ols_estimator * unrot_covariance_log_ols_unweighted
    unrot_ols_gram_covariance = unrot_ols_estimator' * unrot_ols_estimator_covariance
    (ols_expectation_grad_log, ols_expectation) = QuantumACES.calc_ols_merit_grad_log(
        d_unrot_ols,
        unrot_ols_shot_weights,
        unrot_ols_estimator,
        unrot_ols_estimator_covariance,
        unrot_ols_gram_covariance,
    )
    ols_2 = time()
    # Compare against ForwardDiff
    function DifferentiableOLSExpectation(log_shot_weights)
        # Scale the covariance matrix by the shot weights
        shot_weights = exp.(-log_shot_weights) / sum(exp.(-log_shot_weights))
        tuple_times_factor = sum(shot_weights .* d_unrot_ols.tuple_times)
        unrot_shot_weights_factor =
            tuple_times_factor * Diagonal(
                vcat(
                    [
                        (1 / shot_weights[idx]) * ones(unrot_mapping_lengths[idx]) for
                        idx in 1:C_unrot
                    ]...,
                ),
            )
        unrot_covariance_log_scaled =
            unrot_covariance_log_ols_unweighted * unrot_shot_weights_factor
        # Calculate the covariance matrix of the gate eigenvalues, including the first-order Taylor approximation
        unrot_ols_gate_eigenvalues_cov_scaled = Symmetric(
            unrot_ols_estimator * unrot_covariance_log_scaled * unrot_ols_estimator',
        )
        # Calculate the expectation
        cov_trace = tr(unrot_ols_gate_eigenvalues_cov_scaled)
        cov_sq_trace = tr(unrot_ols_gate_eigenvalues_cov_scaled^2)
        expectation =
            sqrt(cov_trace) * (1 - cov_sq_trace / (4 * cov_trace^2)) / sqrt(N_unrot)
        return expectation
    end
    unrot_ols_log_shot_weights = -log.(unrot_ols_shot_weights)
    unrot_ols_log_shot_weights .-= minimum(unrot_ols_log_shot_weights)
    ols_expectation_grad_log_test =
        ForwardDiff.gradient(DifferentiableOLSExpectation, unrot_ols_log_shot_weights)
    ols_3 = time()
    @test ols_3 - ols_2 > ols_2 - ols_1
    @test ols_expectation_grad_log ≈ ols_expectation_grad_log_test
    println(
        "The OLS gradient took $(round(ols_2 - ols_1, digits=3)) s to compute, whereas ForwardDiff took $(round(ols_3 - ols_2, digits=3)) s.",
    )
    println(
        "The optimised OLS figure of merit is $(round(unrot_ols_opt_expectation, digits=3)), whereas the unoptimised OLS figure of merit is $(round(unrot_ols_unopt_expectation, digits=3)).",
    )
end
# Compare the merits for the three different LS estimators
@testset "LS merit comparison" begin
    # Compare the merits when optimised for the three LS estimators
    (d_rot_set, rot_covariance_log_set, rot_merit_descent_set, rot_merit_array) =
        compare_ls_optimise_weights(d_rot, rot_covariance_log)
    # Test that gradient descent converges quickly and to merits within expected bounds
    max_steps = 35
    merit_min = 2.4
    merit_max = 4.0
    @test all(length.(rot_merit_descent_set) .<= max_steps)
    @test all([rot_merit_descent_set[idx][end] .> merit_min for idx in 1:3])
    @test all([rot_merit_descent_set[idx][end] .< merit_max for idx in 1:3])
    # Test that the GLS merit is better than the WLS merit
    @test rot_merit_array[1, 1] < rot_merit_array[2, 2]
    # Test that the WLS merit is better than the OLS merit
    @test rot_merit_array[2, 2] < rot_merit_array[3, 3]
    # Test that the GLS merit is best with shot weights optimised for GLS
    @test rot_merit_array[1, 1] <= rot_merit_array[2, 1]
    @test rot_merit_array[1, 1] <= rot_merit_array[3, 1]
    # Test that the WLS merit is best with shot weights optimised for WLS
    @test rot_merit_array[2, 2] <= rot_merit_array[1, 2]
    @test rot_merit_array[2, 2] <= rot_merit_array[3, 2]
    # Test that the OLS merit is best with shot weights optimised for OLS
    @test rot_merit_array[3, 3] <= rot_merit_array[1, 3]
    @test rot_merit_array[3, 3] <= rot_merit_array[2, 3]
    # Pretty print the results
    pretty_print(rot_merit_array)
end
# Test the functions for growing and pruning designs
@testset "Growing and pruning designs" begin
    # Generate the design
    d_rot_basic = generate_design(rotated_planar, rot_basic)
    rot_covariance_log_basic = calc_covariance_log(d_rot_basic)
    # Grow the design
    (d_rot_grow, rot_covariance_log_grow) = QuantumACES.grow_design(
        d_rot_basic,
        rot_covariance_log_basic,
        rotated_planar.circuit_tuple,
    )
    # Test printing functionality
    pretty_print(d_rot_grow)
    println(get_mapping_string(d_rot_grow.mapping_ensemble[end][end], d_rot_grow.c))
    # Test that the grown design and covariance matrix are correct
    rot_grow = [rot_basic; [rotated_planar.circuit_tuple]]
    d_rot_grow_test = generate_design(rotated_planar, rot_grow)
    rot_covariance_log_test = calc_covariance_log(d_rot_grow)
    @test d_rot_grow.c == d_rot_grow_test.c
    @test d_rot_grow.full_covariance == d_rot_grow_test.full_covariance
    @test d_rot_grow.matrix == d_rot_grow_test.matrix
    @test d_rot_grow.tuple_set == d_rot_grow_test.tuple_set
    @test d_rot_grow.tuple_set_data == d_rot_grow_test.tuple_set_data
    @test d_rot_grow.mapping_ensemble == d_rot_grow_test.mapping_ensemble
    @test d_rot_grow.experiment_ensemble == d_rot_grow_test.experiment_ensemble
    @test d_rot_grow.covariance_dict_ensemble == d_rot_grow_test.covariance_dict_ensemble
    @test d_rot_grow.prep_ensemble == d_rot_grow_test.prep_ensemble
    @test d_rot_grow.meas_ensemble == d_rot_grow_test.meas_ensemble
    @test d_rot_grow.tuple_times ≈ d_rot_grow_test.tuple_times
    @test d_rot_grow.shot_weights ≈ d_rot_grow_test.shot_weights
    @test d_rot_grow.experiment_numbers == d_rot_grow_test.experiment_numbers
    @test d_rot_grow.experiment_number == d_rot_grow_test.experiment_number
    @test d_rot_grow.ls_type == d_rot_grow_test.ls_type
    @test rot_covariance_log_grow ≈ rot_covariance_log_test
    # Test completing the design
    d_rot_grow_completed =
        complete_design(generate_design(rotated_planar, rot_grow; full_covariance = false))
    @test d_rot_grow.c == d_rot_grow_completed.c
    @test d_rot_grow.full_covariance == d_rot_grow_completed.full_covariance
    @test d_rot_grow.matrix == d_rot_grow_completed.matrix
    @test d_rot_grow.tuple_set == d_rot_grow_completed.tuple_set
    @test d_rot_grow.tuple_set_data == d_rot_grow_completed.tuple_set_data
    @test d_rot_grow.mapping_ensemble == d_rot_grow_completed.mapping_ensemble
    @test d_rot_grow.experiment_ensemble == d_rot_grow_completed.experiment_ensemble
    @test d_rot_grow.covariance_dict_ensemble ==
          d_rot_grow_completed.covariance_dict_ensemble
    @test d_rot_grow.prep_ensemble == d_rot_grow_completed.prep_ensemble
    @test d_rot_grow.meas_ensemble == d_rot_grow_completed.meas_ensemble
    @test d_rot_grow.tuple_times ≈ d_rot_grow_completed.tuple_times
    @test d_rot_grow.shot_weights ≈ d_rot_grow_completed.shot_weights
    @test d_rot_grow.experiment_numbers == d_rot_grow_completed.experiment_numbers
    @test d_rot_grow.experiment_number == d_rot_grow_completed.experiment_number
    @test d_rot_grow.ls_type == d_rot_grow_completed.ls_type
    # Prune the design
    T = length(rot_grow)
    (d_rot_prune, rot_covariance_log_prune) =
        QuantumACES.prune_design(d_rot_grow, rot_covariance_log_grow, T)
    # Test that the pruned design and covariance matrix are correct
    @test d_rot_prune.c == d_rot_basic.c
    @test d_rot_prune.full_covariance == d_rot_basic.full_covariance
    @test d_rot_prune.matrix == d_rot_basic.matrix
    @test d_rot_prune.tuple_set == d_rot_basic.tuple_set
    @test d_rot_prune.tuple_set_data == d_rot_basic.tuple_set_data
    @test d_rot_prune.mapping_ensemble == d_rot_basic.mapping_ensemble
    @test d_rot_prune.experiment_ensemble == d_rot_basic.experiment_ensemble
    @test d_rot_prune.covariance_dict_ensemble == d_rot_basic.covariance_dict_ensemble
    @test d_rot_prune.prep_ensemble == d_rot_basic.prep_ensemble
    @test d_rot_prune.meas_ensemble == d_rot_basic.meas_ensemble
    @test d_rot_prune.tuple_times ≈ d_rot_basic.tuple_times
    @test d_rot_prune.shot_weights ≈ d_rot_basic.shot_weights
    @test d_rot_prune.experiment_numbers == d_rot_basic.experiment_numbers
    @test d_rot_prune.experiment_number == d_rot_basic.experiment_number
    @test d_rot_prune.ls_type == d_rot_basic.ls_type
    @test rot_covariance_log_prune ≈ rot_covariance_log_basic
    # Test that the sparse covariance inverse works correctly
    rot_covariance_log_grow_inv = QuantumACES.sparse_covariance_inv(
        rot_covariance_log_grow,
        length.(d_rot_grow.mapping_ensemble),
    )
    @test Array(rot_covariance_log_grow_inv) ≈ inv(cholesky(Array(rot_covariance_log_grow)))
end
# Test the NRMSE probability distribution
@testset "NRMSE probability distribution" begin
    # Randomly sample eigenvalues according to the WLS estimator covariance matrix
    S = 10^9
    repetitions = 1000
    N = d_rot.c.N
    gate_eigenvalues = d_rot.c.gate_eigenvalues
    (rot_eigenvalues, rot_covariance) = QuantumACES.calc_eigenvalues_covariance(d_rot)
    # Sample the eigenvalues according to the calculated WLS estimator covariance matrix
    est_eigenvalues_distribution =
        MvNormal(rot_eigenvalues, Array((1 / S) * rot_covariance))
    Random.seed!(seed)
    est_eigenvalues_matrix = rand(est_eigenvalues_distribution, repetitions)
    Random.seed!()
    # Calculate the WLS gate eigenvalues and the NRMSE
    est_eigenvalues_coll = Vector{Vector{Float64}}(undef, repetitions)
    wls_gate_eigenvalues_coll = Vector{Vector{Float64}}(undef, repetitions)
    wls_gate_norm_coll = Vector{Float64}(undef, repetitions)
    for idx in 1:repetitions
        est_eigenvalues_coll[idx] = est_eigenvalues_matrix[:, idx]
        wls_gate_eigenvalues_coll[idx] =
            QuantumACES.wls_estimate_gate_eigenvalues(d_rot, est_eigenvalues_coll[idx])
        wls_gate_norm_coll[idx] =
            sqrt(S / N) * norm(wls_gate_eigenvalues_coll[idx] - gate_eigenvalues, 2)
    end
    # Generate the NRMSE probability distribution
    nrmse_min = 2.0
    nrmse_max = 5.0
    pdf_int = 0.01
    nrmse_values = collect(nrmse_min:pdf_int:nrmse_max)
    (gls_rot_merit, wls_rot_merit, ols_rot_merit) = calc_merit_set(d_rot)
    wls_rot_nrmse_pdf = nrmse_pdf(wls_rot_merit.eigenvalues, nrmse_values)
    # Test that the predicted and simulated values are consistent enough
    @test minimum(wls_gate_norm_coll) > nrmse_values[findfirst(wls_rot_nrmse_pdf .> 0.0)]
    @test maximum(wls_gate_norm_coll) < nrmse_values[findlast(wls_rot_nrmse_pdf .> 0.0)]
    @test isapprox(
        sort(wls_gate_norm_coll)[ceil(Int, repetitions / 2)],
        nrmse_values[findmax(wls_rot_nrmse_pdf)[2]];
        atol = 2e-2,
    )
    # The start and end values of the probability distribution should be zero
    @test wls_rot_nrmse_pdf[1] ≈ 0.0
    @test wls_rot_nrmse_pdf[end] ≈ 0.0
    # This simplifies our trapezoidal rule tests of the integral, expectation, and variance
    wls_rot_nrmse_pdf_integral = sum(wls_rot_nrmse_pdf) * pdf_int
    wls_rot_nrmse_pdf_expectation = sum(wls_rot_nrmse_pdf .* nrmse_values) * pdf_int
    wls_rot_nrmse_pdf_variance =
        sum(wls_rot_nrmse_pdf .* nrmse_values .^ 2) * pdf_int -
        wls_rot_nrmse_pdf_expectation^2
    @test isapprox(wls_rot_nrmse_pdf_integral, 1.0, rtol = 1e-4)
    @test isapprox(wls_rot_nrmse_pdf_expectation, wls_rot_merit.expectation, rtol = 1e-4)
    @test isapprox(wls_rot_nrmse_pdf_variance, wls_rot_merit.variance, rtol = 1e-2)
end

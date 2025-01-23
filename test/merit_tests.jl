using QuantumACES, LinearAlgebra, SparseArrays, ForwardDiff, Random, Distributions, Test
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
rotated_planar =
    get_circuit(rotated_param, log_param; noisy_prep = true, noisy_meas = false)
unrotated_planar = get_circuit(unrotated_param, dep_param)
# Set up designs
rot_basic = get_basic_tuple_set(rotated_planar)
rot_tuple = get_circuit_tuple(rotated_planar)
rot_tuple_set = [[rot_tuple]; rot_basic]
rot_basic_tuple_set_data = get_tuple_set_data(rotated_planar, rot_basic)
rot_tuple_set_data = get_tuple_set_data(rotated_planar, rot_tuple_set)
d_rot = generate_design(rotated_planar, rot_tuple_set)
unrot_basic = get_basic_tuple_set(unrotated_planar)
unrot_tuple = get_circuit_tuple(unrotated_planar)
unrot_tuple_set = [[unrot_tuple]; unrot_basic]
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
# Set up gradient descent parameters
max_steps = 10
rot_covariance_log = calc_covariance_log(d_rot)
N_rot = rotated_planar.gate_data.N
C_rot = length(d_rot.tuple_set)
rot_mapping_lengths = length.(d_rot.mapping_ensemble)
rot_gate_eigenvalues_diag = sparse(Diagonal(d_rot.c.gate_eigenvalues))
unrot_covariance_log = calc_covariance_log(d_unrot)
N_unrot = unrotated_planar.gate_data.N
C_unrot = length(d_unrot.tuple_set)
unrot_mapping_lengths = length.(d_unrot.mapping_ensemble)
unrot_gate_eigenvalues_diag = sparse(Diagonal(d_unrot.c.gate_eigenvalues))
gls_ls_type = :gls
gls_est_type = :ordinary
rot_gls_transform = get_transform(d_rot, gls_est_type)
rot_gls_gate_transform = rot_gls_transform * rot_gate_eigenvalues_diag
wls_ls_type = :wls
wls_est_type = :relative
rot_wls_transform = get_transform(d_rot, wls_est_type)
rot_wls_gate_transform = rot_wls_transform * rot_gate_eigenvalues_diag
ols_ls_type = :ols
ols_est_type = :marginal
unrot_ols_transform = get_transform(d_unrot, ols_est_type)
unrot_ols_gate_transform = unrot_ols_transform * unrot_gate_eigenvalues_diag
# Test the merit gradient for GLS
@testset "GLS merit gradient" begin
    # Test GLS gradient descent
    (d_rot_gls, rot_covariance_log_gls, rot_merit_descent_gls) = gls_optimise_weights(
        d_rot,
        rot_covariance_log;
        options = OptimOptions(;
            ls_type = gls_ls_type,
            est_type = gls_est_type,
            max_steps = max_steps,
        ),
    )
    @test d_rot_gls.shot_weights != d_rot.shot_weights
    @test rot_covariance_log_gls != rot_covariance_log
    # Test that optimising the sum is equivalent to optimising for the ordinary figure of merit with the appropriate factor
    (d_rot_gls_sum, rot_covariance_log_gls_sum, rot_merit_descent_gls_sum) =
        gls_optimise_weights(
            d_rot,
            rot_covariance_log;
            options = OptimOptions(;
                ls_type = gls_ls_type,
                est_type = :sum,
                est_weight = 1.0,
                max_steps = max_steps,
            ),
        )
    @test d_rot_gls_sum.shot_weights ≈ d_rot_gls.shot_weights
    @test rot_covariance_log_gls_sum ≈ rot_covariance_log_gls
    @test rot_merit_descent_gls_sum ≈ rot_merit_descent_gls
    # Test the covariance matrix output by gradient descent is correct
    rot_gls_opt_merit = calc_merit(d_rot_gls)
    rot_gls_opt_cov_eigvals = eigvals(
        rot_gls_transform *
        calc_gls_covariance(d_rot_gls, rot_covariance_log_gls) *
        rot_gls_transform',
    )
    (rot_gls_opt_expectation, rot_gls_opt_variance) =
        QuantumACES.nrmse_moments(rot_gls_opt_cov_eigvals)
    if gls_est_type == :ordinary
        @test rot_gls_opt_merit.gls_cov_eigenvalues ≈ rot_gls_opt_cov_eigvals
        @test rot_gls_opt_merit.gls_expectation ≈ rot_gls_opt_expectation
        @test rot_gls_opt_merit.gls_variance ≈ rot_gls_opt_variance
    elseif gls_est_type == :marginal
        @test rot_gls_opt_merit.gls_marginal_cov_eigenvalues ≈ rot_gls_opt_cov_eigvals
        @test rot_gls_opt_merit.gls_marginal_expectation ≈ rot_gls_opt_expectation
        @test rot_gls_opt_merit.gls_marginal_variance ≈ rot_gls_opt_variance
    elseif gls_est_type == :relative
        @test rot_gls_opt_merit.gls_relative_cov_eigenvalues ≈ rot_gls_opt_cov_eigvals
        @test rot_gls_opt_merit.gls_relative_expectation ≈ rot_gls_opt_expectation
        @test rot_gls_opt_merit.gls_relative_variance ≈ rot_gls_opt_variance
    else
        throw(
            error(
                "The estimator type $(gls_est_type) must be either :ordinary, :marginal, or :relative.",
            ),
        )
    end
    # Check that gradient descent improved the figure of merit
    rot_gls_unopt_expectation =
        calc_gls_moments(d_rot, rot_covariance_log; est_type = gls_est_type)[1]
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
    rot_covariance_log_gls_unweighted_inv =
        sparse_covariance_inv(rot_covariance_log_gls_unweighted, rot_mapping_lengths)
    (gls_expectation_grad_log, gls_expectation) = QuantumACES.calc_gls_merit_grad_log(
        d_rot_gls,
        rot_gls_shot_weights,
        rot_covariance_log_gls_unweighted_inv,
        rot_gls_gate_transform,
    )
    gls_2 = time()
    # Compare against ForwardDiff
    function DifferentiableGLSExpectation(log_shot_weights)
        # Scale the covariance matrix by the shot weights
        shot_weights = exp.(-log_shot_weights) / sum(exp.(-log_shot_weights))
        times_factor = sum(shot_weights .* d_rot_gls.tuple_times)
        rot_shot_weights_factor_inv =
            (1 / times_factor) * Diagonal(
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
        # Marginalise the covariance matrix
        rot_gls_gate_eigenvalues_cov = Symmetric(
            rot_gls_transform * rot_gls_gate_eigenvalues_cov_scaled * rot_gls_transform',
        )
        # Calculate the expectation
        cov_trace = tr(rot_gls_gate_eigenvalues_cov)
        cov_sq_trace = tr(rot_gls_gate_eigenvalues_cov^2)
        N_cov = size(rot_gls_gate_eigenvalues_cov, 1)
        expectation = sqrt(cov_trace) * (1 - cov_sq_trace / (4 * cov_trace^2)) / sqrt(N_cov)
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
        "The $(gls_ls_type) $(gls_est_type) gradient took $(round(gls_2 - gls_1, digits = 3)) s to compute, whereas ForwardDiff took $(round(gls_3 - gls_2, digits = 3)) s.",
    )
    println(
        "The optimised $(gls_ls_type) $(gls_est_type) figure of merit is $(round(rot_gls_opt_expectation, digits = 3)), whereas the unoptimised $(gls_ls_type) $(gls_est_type) figure of merit is $(round(rot_gls_unopt_expectation, digits = 3)).",
    )
end
# Test the merit gradient for WLS
@testset "WLS merit gradient" begin
    # Test WLS gradient descent
    (d_rot_wls, rot_covariance_log_wls, rot_merit_descent_wls) = wls_optimise_weights(
        d_rot,
        rot_covariance_log;
        options = OptimOptions(;
            ls_type = :wls,
            est_type = wls_est_type,
            max_steps = max_steps,
        ),
    )
    @test d_rot_wls.shot_weights != d_rot.shot_weights
    @test rot_covariance_log_wls != rot_covariance_log
    # Test that optimising the product is equivalent to optimising for the relative figure of merit with the appropriate factor
    (d_rot_wls_prod, rot_covariance_log_wls_prod, rot_merit_descent_wls_prod) =
        wls_optimise_weights(
            d_rot,
            rot_covariance_log;
            options = OptimOptions(;
                ls_type = :wls,
                est_type = :prod,
                est_weight = 0.0,
                max_steps = max_steps,
            ),
        )
    @test d_rot_wls_prod.shot_weights ≈ d_rot_wls.shot_weights
    @test rot_covariance_log_wls_prod ≈ rot_covariance_log_wls
    @test rot_merit_descent_wls_prod ≈ rot_merit_descent_wls
    # Check the covariance matrix output by gradient descent yields the correct quantities
    rot_wls_opt_merit = calc_merit(d_rot_wls)
    rot_wls_opt_cov_eigvals = eigvals(
        rot_wls_transform *
        calc_wls_covariance(d_rot_wls, rot_covariance_log_wls) *
        rot_wls_transform',
    )
    (rot_wls_opt_expectation, rot_wls_opt_variance) =
        QuantumACES.nrmse_moments(rot_wls_opt_cov_eigvals)
    if wls_est_type == :ordinary
        @test rot_wls_opt_merit.wls_cov_eigenvalues ≈ rot_wls_opt_cov_eigvals
        @test rot_wls_opt_merit.wls_expectation ≈ rot_wls_opt_expectation
        @test rot_wls_opt_merit.wls_variance ≈ rot_wls_opt_variance
    elseif wls_est_type == :marginal
        @test rot_wls_opt_merit.wls_marginal_cov_eigenvalues ≈ rot_wls_opt_cov_eigvals
        @test rot_wls_opt_merit.wls_marginal_expectation ≈ rot_wls_opt_expectation
        @test rot_wls_opt_merit.wls_marginal_variance ≈ rot_wls_opt_variance
    elseif wls_est_type == :relative
        @test rot_wls_opt_merit.wls_relative_cov_eigenvalues ≈ rot_wls_opt_cov_eigvals
        @test rot_wls_opt_merit.wls_relative_expectation ≈ rot_wls_opt_expectation
        @test rot_wls_opt_merit.wls_relative_variance ≈ rot_wls_opt_variance
    else
        throw(
            error(
                "The estimator type $(wls_est_type) must be either :ordinary, :marginal, or :relative.",
            ),
        )
    end
    # Check that gradient descent improved the figure of merit
    rot_wls_unopt_expectation =
        calc_wls_moments(d_rot, rot_covariance_log; est_type = wls_est_type)[1]
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
        rot_wls_gate_transform,
    )
    wls_2 = time()
    # Compare against ForwardDiff
    function DifferentiableWLSExpectation(log_shot_weights)
        # Scale the covariance matrix by the shot weights
        shot_weights = exp.(-log_shot_weights) / sum(exp.(-log_shot_weights))
        times_factor = sum(shot_weights .* d_rot_wls.tuple_times)
        rot_shot_weights_factor =
            times_factor * Diagonal(
                vcat(
                    [
                        (1 / shot_weights[idx]) * ones(rot_mapping_lengths[idx]) for
                        idx in 1:C_rot
                    ]...,
                ),
            )
        rot_covariance_log_scaled =
            rot_covariance_log_wls_unweighted * rot_shot_weights_factor
        rot_covariance_log_diag_inv = sparse(Diagonal(rot_covariance_log_scaled)^(-1))
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
        # Marginalise the covariance matrix
        rot_wls_gate_eigenvalues_cov = Symmetric(
            rot_wls_transform * rot_wls_gate_eigenvalues_cov_scaled * rot_wls_transform',
        )
        # Calculate the expectation
        cov_trace = tr(rot_wls_gate_eigenvalues_cov)
        cov_sq_trace = tr(rot_wls_gate_eigenvalues_cov^2)
        N_cov = size(rot_wls_gate_eigenvalues_cov, 1)
        expectation = sqrt(cov_trace) * (1 - cov_sq_trace / (4 * cov_trace^2)) / sqrt(N_cov)
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
        "The $(wls_ls_type) $(wls_est_type) gradient took $(round(wls_2 - wls_1, digits = 3)) s to compute, whereas ForwardDiff took $(round(wls_3 - wls_2, digits = 3)) s.",
    )
    println(
        "The optimised $(wls_ls_type) $(wls_est_type) figure of merit is $(round(rot_wls_opt_expectation, digits = 3)), whereas the unoptimised $(wls_ls_type) $(wls_est_type) figure of merit is $(round(rot_wls_unopt_expectation, digits = 3)).",
    )
end
# Test the merit gradient for OLS
@testset "OLS merit gradient" begin
    # Test OLS gradient descent
    (d_unrot_ols, unrot_covariance_log_ols, unrot_merit_descent_ols) = ols_optimise_weights(
        d_unrot,
        unrot_covariance_log;
        options = OptimOptions(;
            ls_type = :ols,
            est_type = ols_est_type,
            max_steps = max_steps,
        ),
    )
    @test d_unrot_ols.shot_weights != d_unrot.shot_weights
    @test unrot_covariance_log_ols != unrot_covariance_log
    # Check the covariance matrix output by gradient descent yields the correct quantities
    unrot_ols_opt_merit = calc_merit(d_unrot_ols)
    unrot_ols_opt_cov_eigvals = eigvals(
        unrot_ols_transform *
        calc_ols_covariance(d_unrot_ols, unrot_covariance_log_ols) *
        unrot_ols_transform',
    )
    (unrot_ols_opt_expectation, unrot_ols_opt_variance) =
        QuantumACES.nrmse_moments(unrot_ols_opt_cov_eigvals)
    if ols_est_type == :ordinary
        @test unrot_ols_opt_merit.ols_cov_eigenvalues ≈ unrot_ols_opt_cov_eigvals
        @test unrot_ols_opt_merit.ols_expectation ≈ unrot_ols_opt_expectation
        @test unrot_ols_opt_merit.ols_variance ≈ unrot_ols_opt_variance
    elseif ols_est_type == :marginal
        @test unrot_ols_opt_merit.ols_marginal_cov_eigenvalues ≈ unrot_ols_opt_cov_eigvals
        @test unrot_ols_opt_merit.ols_marginal_expectation ≈ unrot_ols_opt_expectation
        @test unrot_ols_opt_merit.ols_marginal_variance ≈ unrot_ols_opt_variance
    elseif ols_est_type == :relative
        @test unrot_ols_opt_merit.ols_relative_cov_eigenvalues ≈ unrot_ols_opt_cov_eigvals
        @test unrot_ols_opt_merit.ols_relative_expectation ≈ unrot_ols_opt_expectation
        @test unrot_ols_opt_merit.ols_relative_variance ≈ unrot_ols_opt_variance
    else
        throw(
            error(
                "The estimator type $(ols_est_type) must be either :ordinary, :marginal, or :relative.",
            ),
        )
    end
    # Check that gradient descent improved the figure of merit
    unrot_ols_unopt_expectation =
        calc_ols_moments(d_unrot, unrot_covariance_log; est_type = ols_est_type)[1]
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
        unrot_ols_gate_transform *
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
        times_factor = sum(shot_weights .* d_unrot_ols.tuple_times)
        unrot_shot_weights_factor =
            times_factor * Diagonal(
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
        # This is already marginalised by the definition of unrot_ols_estimator
        unrot_ols_gate_eigenvalues_cov = Symmetric(
            unrot_ols_estimator * unrot_covariance_log_scaled * unrot_ols_estimator',
        )
        # Calculate the expectation
        cov_trace = tr(unrot_ols_gate_eigenvalues_cov)
        cov_sq_trace = tr(unrot_ols_gate_eigenvalues_cov^2)
        N_cov = size(unrot_ols_gate_eigenvalues_cov, 1)
        expectation = sqrt(cov_trace) * (1 - cov_sq_trace / (4 * cov_trace^2)) / sqrt(N_cov)
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
        "The $(ols_ls_type) $(ols_est_type) gradient took $(round(ols_2 - ols_1, digits = 3)) s to compute, whereas ForwardDiff took $(round(ols_3 - ols_2, digits = 3)) s.",
    )
    println(
        "The optimised $(ols_ls_type) $(ols_est_type) figure of merit is $(round(unrot_ols_opt_expectation, digits = 3)), whereas the unoptimised $(ols_ls_type) $(ols_est_type) figure of merit is $(round(unrot_ols_unopt_expectation, digits = 3)).",
    )
end
# Test the functions for growing and pruning designs
@testset "Growing and pruning designs" begin
    # Generate the design
    d_rot_basic = generate_design(rotated_planar, rot_basic_tuple_set_data)
    rot_covariance_log_basic = calc_covariance_log(d_rot_basic)
    # Grow the design
    (d_rot_grow, rot_covariance_log_grow) =
        QuantumACES.grow_design(d_rot_basic, rot_covariance_log_basic, rot_tuple)
    # Test printing functionality
    display(d_rot_grow)
    println(get_mapping_string(d_rot_grow.mapping_ensemble[end][end], d_rot_grow.c))
    # Test that the grown design and covariance matrix are correct
    d_rot_grow_test = generate_design(rotated_planar, rot_tuple_set_data)
    d_rot_grow_unfull =
        generate_design(rotated_planar, rot_tuple_set_data; full_covariance = false)
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
    # Test completing the design covariance matrix
    d_rot_grow_full = get_full_design(d_rot_grow_unfull)
    @test d_rot_grow.c == d_rot_grow_full.c
    @test d_rot_grow.full_covariance == d_rot_grow_full.full_covariance
    @test d_rot_grow.matrix == d_rot_grow_full.matrix
    @test d_rot_grow.tuple_set == d_rot_grow_full.tuple_set
    @test d_rot_grow.tuple_set_data == d_rot_grow_full.tuple_set_data
    @test d_rot_grow.mapping_ensemble == d_rot_grow_full.mapping_ensemble
    @test d_rot_grow.experiment_ensemble == d_rot_grow_full.experiment_ensemble
    @test d_rot_grow.covariance_dict_ensemble == d_rot_grow_full.covariance_dict_ensemble
    @test d_rot_grow.prep_ensemble == d_rot_grow_full.prep_ensemble
    @test d_rot_grow.meas_ensemble == d_rot_grow_full.meas_ensemble
    @test d_rot_grow.tuple_times ≈ d_rot_grow_full.tuple_times
    @test d_rot_grow.shot_weights ≈ d_rot_grow_full.shot_weights
    @test d_rot_grow.experiment_numbers == d_rot_grow_full.experiment_numbers
    @test d_rot_grow.experiment_number == d_rot_grow_full.experiment_number
    @test d_rot_grow.ls_type == d_rot_grow_full.ls_type
    # Test making the design covariance matrix diagonal
    d_rot_grow_diag = get_diag_design(d_rot_grow_test)
    @test d_rot_grow_unfull.c == d_rot_grow_diag.c
    @test d_rot_grow_unfull.full_covariance == d_rot_grow_diag.full_covariance
    @test d_rot_grow_unfull.matrix == d_rot_grow_diag.matrix
    @test d_rot_grow_unfull.tuple_set == d_rot_grow_diag.tuple_set
    @test d_rot_grow_unfull.tuple_set_data == d_rot_grow_diag.tuple_set_data
    @test d_rot_grow_unfull.mapping_ensemble == d_rot_grow_diag.mapping_ensemble
    @test d_rot_grow_unfull.experiment_ensemble == d_rot_grow_diag.experiment_ensemble
    @test d_rot_grow_unfull.covariance_dict_ensemble ==
          d_rot_grow_diag.covariance_dict_ensemble
    @test d_rot_grow_unfull.prep_ensemble == d_rot_grow_diag.prep_ensemble
    @test d_rot_grow_unfull.meas_ensemble == d_rot_grow_diag.meas_ensemble
    @test d_rot_grow_unfull.tuple_times ≈ d_rot_grow_diag.tuple_times
    @test d_rot_grow_unfull.shot_weights ≈ d_rot_grow_diag.shot_weights
    @test d_rot_grow_unfull.experiment_numbers == d_rot_grow_diag.experiment_numbers
    @test d_rot_grow_unfull.experiment_number == d_rot_grow_diag.experiment_number
    @test d_rot_grow_unfull.ls_type == d_rot_grow_diag.ls_type
    # Prune the design
    rot_tuple_idx =
        findfirst(rot_tuple == grow_tuple for grow_tuple in d_rot_grow.tuple_set)
    (d_rot_prune, rot_covariance_log_prune) =
        QuantumACES.prune_design(d_rot_grow, rot_covariance_log_grow, rot_tuple_idx)
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
    rot_grow_mapping_lengths = length.(d_rot_grow.mapping_ensemble)
    rot_covariance_log_grow_inv_factor = QuantumACES.sparse_covariance_inv_factor(
        rot_covariance_log_grow,
        rot_grow_mapping_lengths,
    )
    rot_covariance_log_grow_inv =
        sparse_covariance_inv(rot_covariance_log_grow, rot_grow_mapping_lengths)
    @test rot_covariance_log_grow_inv ≈
          rot_covariance_log_grow_inv_factor' * rot_covariance_log_grow_inv_factor
    cholesky_grow = cholesky(Array(rot_covariance_log_grow))
    cholesky_grow_inv = inv(cholesky_grow)
    cholesky_grow_inv_factor = inv(Array(cholesky_grow.L))
    @test cholesky_grow_inv ≈ cholesky_grow_inv_factor' * cholesky_grow_inv_factor
    @test rot_covariance_log_grow_inv ≈ cholesky_grow_inv
    M_grow = size(rot_covariance_log_grow, 1)
    test_data = rand(M_grow)
    test_matrix = rand(M_grow, M_grow)
    @test (rot_covariance_log_grow_inv_factor * test_matrix) \
          (rot_covariance_log_grow_inv_factor * test_data) ≈
          (cholesky_grow_inv_factor * test_matrix) \ (cholesky_grow_inv_factor * test_data)
    # Test creation of the gate eigenvalue estimator covariance matrix
    unpad_rot_gate_probabilities_cov_grow =
        calc_gate_probabilities_covariance(d_rot_grow; unpad = true)
    rot_marg_gate_probabilities_cov_grow =
        calc_gate_probabilities_covariance(d_rot_grow; est_type = :marginal)
    rot_rel_gate_probabilities_cov_grow =
        calc_gate_probabilities_covariance(d_rot_grow; est_type = :relative)
    N_pad_relative = d_rot_grow.c.gate_data.N_pad_relative
    @test rot_marg_gate_probabilities_cov_grow[1:N_pad_relative, 1:N_pad_relative] ≈
          rot_rel_gate_probabilities_cov_grow
    # Test the precision matrix
    precision_matrix = calc_precision_matrix(d_rot_grow)
    precision_matrix_diag = calc_precision_matrix(d_rot_grow; diagonal = true)
    gls_gate_eigenvalues_cov = calc_gls_covariance(d_rot_grow, rot_covariance_log_test)
    wls_gate_eigenvalues_cov_diag =
        calc_wls_covariance(d_rot_grow, sparse(Diagonal(rot_covariance_log_test)))
    @test precision_matrix ≈ inv(cholesky(gls_gate_eigenvalues_cov))
    @test precision_matrix_diag ≈ inv(cholesky(wls_gate_eigenvalues_cov_diag))
    # Test that the sparse covariance inverse works correctly when clipping certain eigenvalues
    clipped_indices = QuantumACES.get_clipped_indices(test_data)
    clipped_mapping_lengths =
        QuantumACES.get_clipped_mapping_lengths(rot_grow_mapping_lengths, clipped_indices)
    clipped_covariance_log = rot_covariance_log_grow[clipped_indices, clipped_indices]
    clipped_covariance_log_inv_factor = QuantumACES.sparse_covariance_inv_factor(
        clipped_covariance_log,
        clipped_mapping_lengths,
    )
    clipped_covariance_log_inv =
        sparse_covariance_inv(clipped_covariance_log, clipped_mapping_lengths)
    @test clipped_covariance_log_inv ≈
          clipped_covariance_log_inv_factor' * clipped_covariance_log_inv_factor
    cholesky_clipped = cholesky(Array(clipped_covariance_log))
    cholesky_clipped_inv = inv(cholesky_clipped)
    cholesky_clipped_inv_factor = inv(Array(cholesky_clipped.L))
    @test cholesky_clipped_inv ≈ cholesky_clipped_inv_factor' * cholesky_clipped_inv_factor
    @test clipped_covariance_log_inv ≈ cholesky_clipped_inv
    clipped_data = test_data[clipped_indices]
    clipped_matrix = test_matrix[clipped_indices, clipped_indices]
    @test (clipped_covariance_log_inv_factor * clipped_matrix) \
          (clipped_covariance_log_inv_factor * clipped_data) ≈
          (cholesky_clipped_inv_factor * clipped_matrix) \
          (cholesky_clipped_inv_factor * clipped_data)
end
# Test the NRMSE probability distribution
@testset "NRMSE probability distribution" begin
    # Initialise data
    S = 10^9
    repetitions = 1000
    gate_eigenvalues = get_gate_eigenvalues(d_rot)
    N = d_rot.c.gate_data.N
    marginal_gate_eigenvalues = get_marginal_gate_eigenvalues(d_rot)
    N_marginal = d_rot.c.gate_data.N_marginal
    relative_gate_eigenvalues = get_relative_gate_eigenvalues(d_rot)
    N_relative = d_rot.c.gate_data.N_relative
    rot_eigenvalues = get_eigenvalues(d_rot)
    rot_covariance = calc_covariance(d_rot)
    rot_covariance_log = calc_covariance_log(d_rot)
    rot_mapping_lengths = length.(d_rot.mapping_ensemble)
    rot_diag_covariance_log_inv_factor = QuantumACES.sparse_covariance_inv_factor(
        sparse(Diagonal(rot_covariance_log)),
        rot_mapping_lengths,
    )
    rot_merit = calc_merit(d_rot)
    # Sample the eigenvalues according to the calculated WLS estimator covariance matrix
    est_eigenvalues_distribution =
        MvNormal(rot_eigenvalues, Array((1 / S) * rot_covariance))
    Random.seed!(seed)
    est_eigenvalues_matrix = rand(est_eigenvalues_distribution, repetitions)
    Random.seed!()
    # Calculate the WLS gate eigenvalues and the NRMSE
    design_matrix = convert(SparseMatrixCSC{Float64, Int}, d_rot.matrix)
    est_eigenvalues_coll = Vector{Vector{Float64}}(undef, repetitions)
    wls_gate_eigenvalues_coll = Vector{Vector{Float64}}(undef, repetitions)
    wls_marginal_gate_eigenvalues_coll = Vector{Vector{Float64}}(undef, repetitions)
    wls_relative_gate_eigenvalues_coll = Vector{Vector{Float64}}(undef, repetitions)
    wls_gate_norm_coll = Vector{Float64}(undef, repetitions)
    wls_marginal_gate_norm_coll = Vector{Float64}(undef, repetitions)
    wls_relative_gate_norm_coll = Vector{Float64}(undef, repetitions)
    for idx in 1:repetitions
        est_eigenvalues_coll[idx] = est_eigenvalues_matrix[:, idx]
        wls_gate_eigenvalues_coll[idx] = QuantumACES.estimate_gate_eigenvalues(
            design_matrix,
            rot_diag_covariance_log_inv_factor,
            est_eigenvalues_coll[idx],
        )
        wls_gate_norm_coll[idx] =
            sqrt(S / N) * norm(wls_gate_eigenvalues_coll[idx] - gate_eigenvalues, 2)
        wls_marginal_gate_eigenvalues_coll[idx] =
            get_marginal_gate_eigenvalues(wls_gate_eigenvalues_coll[idx], d_rot.c.gate_data)
        wls_marginal_gate_norm_coll[idx] =
            sqrt(S / N_marginal) *
            norm(wls_marginal_gate_eigenvalues_coll[idx] - marginal_gate_eigenvalues, 2)
        wls_relative_gate_eigenvalues_coll[idx] =
            get_relative_gate_eigenvalues(wls_gate_eigenvalues_coll[idx], d_rot.c.gate_data)
        wls_relative_gate_norm_coll[idx] =
            sqrt(S / N_relative) *
            norm(wls_relative_gate_eigenvalues_coll[idx] - relative_gate_eigenvalues, 2)
    end
    # Generate the NRMSE probability distribution for the ordinary gate eigenvalues
    nrmse_min = 0.0
    nrmse_max = 5.0
    pdf_int = 0.01
    nrmse_values = collect(nrmse_min:pdf_int:nrmse_max)
    wls_rot_nrmse_pdf = nrmse_pdf(rot_merit.wls_cov_eigenvalues, nrmse_values)
    # Test that the predicted and simulated values are consistent enough
    @test minimum(wls_gate_norm_coll) > nrmse_values[findfirst(wls_rot_nrmse_pdf .> 0.0)]
    @test maximum(wls_gate_norm_coll) < nrmse_values[findlast(wls_rot_nrmse_pdf .> 0.0)]
    @test isapprox(
        sort(wls_gate_norm_coll)[ceil(Int, repetitions / 2)],
        nrmse_values[findmax(wls_rot_nrmse_pdf)[2]];
        atol = 0.02,
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
    @test isapprox(wls_rot_nrmse_pdf_expectation, rot_merit.wls_expectation, rtol = 1e-4)
    @test isapprox(wls_rot_nrmse_pdf_variance, rot_merit.wls_variance, rtol = 1e-2)
    # Generate the NRMSE probability distribution for the marginal gate eigenvalues
    marg_pdf_int = 0.005
    marg_nrmse_values = collect(nrmse_min:marg_pdf_int:nrmse_max)
    wls_rot_marginal_nrmse_pdf =
        nrmse_pdf(rot_merit.wls_marginal_cov_eigenvalues, marg_nrmse_values)
    # Test that the predicted and simulated values are consistent enough
    @test minimum(wls_marginal_gate_norm_coll) >
          marg_nrmse_values[findfirst(wls_rot_marginal_nrmse_pdf .> 0.0)]
    @test maximum(wls_marginal_gate_norm_coll) <
          marg_nrmse_values[findlast(wls_rot_marginal_nrmse_pdf .> 0.0)]
    @test isapprox(
        sort(wls_marginal_gate_norm_coll)[ceil(Int, repetitions / 2)],
        marg_nrmse_values[findmax(wls_rot_marginal_nrmse_pdf)[2]];
        atol = 0.02,
    )
    # The start and end values of the probability distribution should be zero
    @test wls_rot_marginal_nrmse_pdf[1] ≈ 0.0
    @test wls_rot_marginal_nrmse_pdf[end] ≈ 0.0
    # This simplifies our trapezoidal rule tests of the integral, expectation, and variance
    wls_rot_marginal_nrmse_pdf_integral = sum(wls_rot_marginal_nrmse_pdf) * marg_pdf_int
    wls_rot_marginal_nrmse_pdf_expectation =
        sum(wls_rot_marginal_nrmse_pdf .* marg_nrmse_values) * marg_pdf_int
    wls_rot_marginal_nrmse_pdf_variance =
        sum(wls_rot_marginal_nrmse_pdf .* marg_nrmse_values .^ 2) * marg_pdf_int -
        wls_rot_marginal_nrmse_pdf_expectation^2
    @test isapprox(wls_rot_marginal_nrmse_pdf_integral, 1.0, rtol = 1e-4)
    @test isapprox(
        wls_rot_marginal_nrmse_pdf_expectation,
        rot_merit.wls_marginal_expectation,
        rtol = 1e-4,
    )
    @test isapprox(
        wls_rot_marginal_nrmse_pdf_variance,
        rot_merit.wls_marginal_variance,
        rtol = 1e-2,
    )
    # Generate the NRMSE probability distribution for the relative gate eigenvalues
    wls_rot_relative_nrmse_pdf =
        nrmse_pdf(rot_merit.wls_relative_cov_eigenvalues, marg_nrmse_values)
    # Test that the predicted and simulated values are consistent enough
    @test minimum(wls_relative_gate_norm_coll) >
          marg_nrmse_values[findfirst(wls_rot_relative_nrmse_pdf .> 0.0)]
    @test maximum(wls_relative_gate_norm_coll) <
          marg_nrmse_values[findlast(wls_rot_relative_nrmse_pdf .> 0.0)]
    @test isapprox(
        sort(wls_relative_gate_norm_coll)[ceil(Int, repetitions / 2)],
        marg_nrmse_values[findmax(wls_rot_relative_nrmse_pdf)[2]];
        atol = 0.02,
    )
    # The start and end values of the probability distribution should be zero
    @test wls_rot_relative_nrmse_pdf[1] ≈ 0.0
    @test wls_rot_relative_nrmse_pdf[end] ≈ 0.0
    # This simplifies our trapezoidal rule tests of the integral, expectation, and variance
    wls_rot_relative_nrmse_pdf_integral = sum(wls_rot_relative_nrmse_pdf) * marg_pdf_int
    wls_rot_relative_nrmse_pdf_expectation =
        sum(wls_rot_relative_nrmse_pdf .* marg_nrmse_values) * marg_pdf_int
    wls_rot_relative_nrmse_pdf_variance =
        sum(wls_rot_relative_nrmse_pdf .* marg_nrmse_values .^ 2) * marg_pdf_int -
        wls_rot_relative_nrmse_pdf_expectation^2
    @test isapprox(wls_rot_relative_nrmse_pdf_integral, 1.0, rtol = 1e-4)
    @test isapprox(
        wls_rot_relative_nrmse_pdf_expectation,
        rot_merit.wls_relative_expectation,
        rtol = 1e-4,
    )
    @test isapprox(
        wls_rot_relative_nrmse_pdf_variance,
        rot_merit.wls_relative_variance,
        rtol = 1e-2,
    )
end

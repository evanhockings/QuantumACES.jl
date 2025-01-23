"""
    get_shot_weights_factor(shot_weights::Vector{Float64}, tuple_times::Vector{Float64}, mapping_lengths::Vector{Int})

Returns the shot weights factor for the sparse block diagonal circuit (log-)eigenvalue estimator covariance matrix for the shot weights `shot_weights` and tuple times `tuple_times`, where the block sizes are specified by `mapping_lengths`.
"""
function get_shot_weights_factor(
    shot_weights::Vector{Float64},
    tuple_times::Vector{Float64},
    mapping_lengths::Vector{Int},
)
    # Initialise parameters
    tuple_number = length(mapping_lengths)
    @assert length(shot_weights) == tuple_number
    @assert length(tuple_times) == tuple_number
    # Calculate the shot weights factor
    times_factor = sum(shot_weights .* tuple_times)
    shot_weights_factor =
        times_factor * Diagonal(
            vcat(
                [
                    (1 / shot_weights[i]) * ones(mapping_lengths[i]) for i in 1:tuple_number
                ]...,
            ),
        )
    return shot_weights_factor::Diagonal{Float64, Vector{Float64}}
end

"""
    get_shot_weights_factor_inv(shot_weights::Vector{Float64}, tuple_times::Vector{Float64}, mapping_lengths::Vector{Int})

Returns the shot weights inverse factor for the sparse block diagonal circuit (log-)eigenvalue estimator covariance matrix for the shot weights `shot_weights` and tuple times `tuple_times`, where the block sizes are specified by `mapping_lengths`.
"""
function get_shot_weights_factor_inv(
    shot_weights::Vector{Float64},
    tuple_times::Vector{Float64},
    mapping_lengths::Vector{Int},
)
    # Initialise parameters
    tuple_number = length(mapping_lengths)
    @assert length(shot_weights) == tuple_number
    @assert length(tuple_times) == tuple_number
    # Calculate the shot weights factor inverse
    times_factor = sum(shot_weights .* tuple_times)
    shot_weights_factor_inv =
        (1 / times_factor) * Diagonal(
            vcat([shot_weights[i] * ones(mapping_lengths[i]) for i in 1:tuple_number]...),
        )
    return shot_weights_factor_inv::Diagonal{Float64, Vector{Float64}}
end

"""
    get_shot_weights_local_grad(shot_weights::Vector{Float64}, tuple_times::Vector{Float64})

Returns the local gradient factor corresponding to each tuple's block in the covariance matrix, at shot weights `shot_weights` and tuple times `tuple_times`.
"""
function get_shot_weights_local_grad(
    shot_weights::Vector{Float64},
    tuple_times::Vector{Float64},
)
    # The local gradient applies only tuple's portion of the covariance matrix
    times_factor = sum(shot_weights .* tuple_times)
    shot_weights_local_grad = -times_factor ./ (shot_weights .^ 2)
    return shot_weights_local_grad::Vector{Float64}
end

"""
    get_merit_grad(sigma_tr::Float64, sigma_tr_grad::Vector{Float64}, sigma_sq_tr::Float64, sigma_sq_tr_grad::Vector{Float64}, N::Integer)

Returns the gradient of the figure of merit with respect to the shot weights given the trace of the gate eigenvalue estimator covariance matrix `sigma_tr` and its square `sigma_sq_tr`, the gradient with respect to the shot weights `sigma_tr_grad` and the gradient of the square `sigma_sq_tr_grad`, and the number of gate eigenvalues `N`.
"""
function get_merit_grad(
    sigma_tr::Float64,
    sigma_tr_grad::Vector{Float64},
    sigma_sq_tr::Float64,
    sigma_sq_tr_grad::Vector{Float64},
    N::Integer,
)
    # Calculate the gradient of the figure of merit, or expected NRMSE, with respect to the shot weights
    merit_grad =
        (
            sigma_tr_grad *
            (((sigma_tr^(2) / 2) + (3 * sigma_sq_tr / 8)) / (sigma_tr^(5 / 2))) -
            sigma_sq_tr_grad * (1 / (4 * sigma_tr^(3 / 2)))
        ) / sqrt(N)
    return merit_grad::Vector{Float64}
end

"""
    get_shot_weights_log_matrix(shot_weights::Vector{Float64})

Returns the matrix that transforms the gradient of the figure of merit with respect to the shot weights `shot_weights` into the gradient of the figure of merit with respect to the logarithms of the shot weights.
"""
function get_shot_weights_log_matrix(shot_weights::Vector{Float64})
    shot_weights_log_matrix = shot_weights * shot_weights' - Diagonal(shot_weights)
    return shot_weights_log_matrix::Matrix{Float64}
end

"""
    calc_gls_merit_grad_log(d::Design, shot_weights::Vector{Float64}, covariance_log_unweighted_inv::SparseMatrixCSC{Float64, Int}, gate_transform_matrix::SparseMatrixCSC{Float64, Int})

Returns the gradient of the generalised least squares (GLS) figure of merit for the design `d` with respect to the logarithms of the shot weights `shot_weights`, using the inverse of the unweighted (by the shot weights factor) covariance matrix of the circuit log-eigenvalue estimator `covariance_log_unweighted_inv`, with the estimator type implicitly specified by the gate eigenvalue transform matrix `gate_transform_matrix`.
"""
function calc_gls_merit_grad_log(
    d::Design,
    shot_weights::Vector{Float64},
    covariance_log_unweighted_inv::SparseMatrixCSC{Float64, Int},
    gate_transform_matrix::SparseMatrixCSC{Float64, Int},
)
    # Initialise data
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    mapping_indices =
        [mapping_lower[idx]:mapping_upper[idx] for idx in 1:length(d.tuple_set)]
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(shot_weights, d.tuple_times, mapping_lengths)
    shot_weights_local_grad = get_shot_weights_local_grad(shot_weights, d.tuple_times)
    shot_weights_log_matrix = get_shot_weights_log_matrix(shot_weights)
    # Precompute useful quantities
    covariance_log_inv = covariance_log_unweighted_inv * shot_weights_factor_inv
    sigma_prime = Symmetric(
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_inv * d.matrix)))),
    )
    design_weights = d.matrix' * shot_weights_factor_inv
    covariance_sigma_prime = covariance_log_inv * d.matrix * sigma_prime
    # Calculate gate transform quantities
    N = size(gate_transform_matrix, 1)
    sigma = Symmetric(gate_transform_matrix * sigma_prime * gate_transform_matrix')
    gate_sigma_prime = gate_transform_matrix' * gate_transform_matrix * sigma_prime
    gls_matrix = covariance_sigma_prime * gate_sigma_prime
    # Compute the trace of the covariance matrix and its square
    sigma_tr = tr(sigma)
    sigma_sq_tr = tr(sigma^2)
    merit = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    # Compute the gradient of the trace of the covariance matrix
    sigma_tr_diag = diag(gls_matrix * design_weights)
    sigma_tr_partial = [sum(sigma_tr_diag[indices]) for indices in mapping_indices]
    sigma_tr_grad =
        sum(sigma_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_tr_partial .* shot_weights_local_grad
    # Compute the gradient of the trace of the square of the covariance matrix
    sigma_sq_tr_diag = 2 * diag(gls_matrix * gate_sigma_prime * design_weights)
    sigma_sq_tr_partial = [sum(sigma_sq_tr_diag[indices]) for indices in mapping_indices]
    sigma_sq_tr_grad =
        sum(sigma_sq_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_sq_tr_partial .* shot_weights_local_grad
    # Calculate the gradient of the figure of merit with respect to the shot weights
    merit_grad = get_merit_grad(sigma_tr, sigma_tr_grad, sigma_sq_tr, sigma_sq_tr_grad, N)
    # Calculate the gradient of the figure of merit with respect to the log shot weights
    merit_grad_log = shot_weights_log_matrix * merit_grad
    return (merit_grad_log::Vector{Float64}, merit::Float64)
end

"""
    gls_optimise_weights(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; options::OptimOptions = OptimOptions())

Returns versions of the design `d` and circuit log-eigenvalue estimator covariance matrix `covariance_log` after optimising the shot weights with respect to the generalised least squares (GLS) figure of merit, alongside the figure of merit values at each step.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function gls_optimise_weights(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    ls_type = options.ls_type
    @assert ls_type == :gls "Inappropriate least squares optimisation type $(ls_type) supplied."
    est_type = options.est_type
    est_weight = options.est_weight
    learning_rate = options.learning_rate
    momentum = options.momentum
    learning_rate_scale_factor = options.learning_rate_scale_factor
    max_steps = options.max_steps
    convergence_threshold = options.convergence_threshold
    convergence_steps = options.convergence_steps
    diagnostics = options.grad_diagnostics
    # Initialise data
    tuple_number = length(d.tuple_set)
    shot_weights = project_simplex(d.shot_weights)
    mapping_lengths = length.(d.mapping_ensemble)
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
    covariance_log_unweighted = covariance_log * shot_weights_factor_inv
    covariance_log_unweighted_inv =
        sparse_covariance_inv(covariance_log_unweighted, mapping_lengths)
    gate_eigenvalues_diag = sparse(Diagonal(get_gate_eigenvalues(d)))
    if est_type == :sum || est_type == :prod
        ord_gate_transform_matrix = get_transform(d, :ordinary) * gate_eigenvalues_diag
        rel_gate_transform_matrix = get_transform(d, :relative) * gate_eigenvalues_diag
    else
        gate_transform_matrix = get_transform(d, est_type) * gate_eigenvalues_diag
    end
    # Perform gradient descent
    stepping = true
    step = 1
    recently_pruned = 0
    recently_zeroed = 0
    scaled_learning_rate = learning_rate
    old_shot_weights = deepcopy(shot_weights)
    velocity = zeros(tuple_number)
    merit_descent = Vector{Float64}()
    while stepping
        # Calculate the gradient of the figure of merit
        nesterov_log_shot_weights = -log.(shot_weights) + momentum * velocity
        nesterov_shot_weights =
            exp.(-nesterov_log_shot_weights) / sum(exp.(-nesterov_log_shot_weights))
        if est_type == :sum || est_type == :prod
            (ord_merit_grad_log, ord_merit) = calc_gls_merit_grad_log(
                d,
                nesterov_shot_weights,
                covariance_log_unweighted_inv,
                ord_gate_transform_matrix,
            )
            (rel_merit_grad_log, rel_merit) = calc_gls_merit_grad_log(
                d,
                nesterov_shot_weights,
                covariance_log_unweighted_inv,
                rel_gate_transform_matrix,
            )
            if est_type == :sum
                merit = est_weight * ord_merit + (1 - est_weight) * rel_merit
                merit_grad_log =
                    est_weight * ord_merit_grad_log + (1 - est_weight) * rel_merit_grad_log
            elseif est_type == :prod
                merit = ord_merit^(est_weight) * rel_merit^(1 - est_weight)
                merit_grad_log =
                    merit * (
                        est_weight * ord_merit_grad_log / ord_merit +
                        (1 - est_weight) * rel_merit_grad_log / rel_merit
                    )
            else
                throw(error("Unsupported estimator type $(est_type)."))
            end
        else
            (merit_grad_log, merit) = calc_gls_merit_grad_log(
                d,
                nesterov_shot_weights,
                covariance_log_unweighted_inv,
                gate_transform_matrix,
            )
        end
        push!(merit_descent, merit)
        # Update the shot weights, ensuring both they and the gradient are appropriately normalised
        velocity = momentum * velocity - scaled_learning_rate * merit_grad_log
        log_shot_weights_update = -log.(shot_weights) + velocity
        shot_weights_update =
            exp.(-log_shot_weights_update) / sum(exp.(-log_shot_weights_update))
        if (length(merit_descent) >= 2) && (merit_descent[end] > merit_descent[end - 1])
            shot_weights = old_shot_weights
            velocity = zeros(tuple_number)
            if recently_zeroed > 0
                scaled_learning_rate /= learning_rate_scale_factor
            end
            if diagnostics
                println(
                    "The updated shot weights worsened the figure of merit in step $(step); the velocity and shot weights have been reset$(recently_zeroed > 0 ? " and the learning rate has been reduced." : ".")",
                )
            end
            recently_pruned += 1
            recently_zeroed = convergence_steps
        else
            old_shot_weights = deepcopy(nesterov_shot_weights)
            shot_weights = shot_weights_update
        end
        # Check the step count and convergence
        if step >= convergence_steps
            merit_converged =
                all([
                    abs(1 - merit_descent[idx_2] / merit_descent[idx_1]) <
                    convergence_steps * convergence_threshold for
                    idx_1 in (step - convergence_steps + 1):step for
                    idx_2 in (step - convergence_steps + 1):step
                ]) && (
                    abs(1 - merit_descent[end - 1] / merit_descent[end]) <
                    convergence_threshold
                )
        else
            merit_converged = false
        end
        if (merit_converged || step >= max_steps) && recently_pruned == 0
            stepping = false
            if diagnostics
                if merit_converged
                    println(
                        "Converged after $(step) steps. The $(ls_type) $(est_type) figure of merit is $(round(merit_descent[step], sigdigits = 5)).",
                    )
                else
                    println(
                        "The maximum number of steps $(max_steps) has been reached without convergence. The $(ls_type) $(est_type) figure of merit is $(round(merit_descent[step], sigdigits = 5))$(max_steps > 1 ? ", which differs from the previous step by $(round(abs(merit_descent[step] - merit_descent[step - 1]), sigdigits = 5)), whereas the threshold for convergence is $(convergence_threshold)" : "").",
                    )
                end
            end
        else
            if diagnostics
                println(
                    "The $(ls_type) $(est_type) figure of merit was $(round(merit_descent[step], sigdigits = 5)) prior to step $(step).",
                )
            end
            step += 1
            if recently_pruned > 0
                recently_pruned -= 1
            end
            if recently_zeroed > 0
                recently_zeroed -= 1
            end
        end
    end
    # Update the covariance matrix
    shot_weights_factor =
        get_shot_weights_factor(shot_weights, d.tuple_times, mapping_lengths)
    covariance_log = covariance_log_unweighted * shot_weights_factor
    # Update the design
    @reset d.shot_weights = project_simplex(shot_weights)
    @reset d.ls_type = :gls
    return (
        d::Design,
        covariance_log::SparseMatrixCSC{Float64, Int},
        merit_descent::Vector{Float64},
    )
end

"""
    calc_wls_merit_grad_log(d::Design, shot_weights::Vector{Float64}, covariance_log_unweighted::SparseMatrixCSC{Float64, Int}, gate_transform_matrix::SparseMatrixCSC{Float64, Int})

Returns the gradient of the weighted least squares (WLS) figure of merit for the design `d` with respect to the logarithms of the shot weights `shot_weights`, using the unweighted (by the shot weights factor) covariance matrix of the circuit log-eigenvalue estimator `covariance_log_unweighted`, with the estimator type implicitly specified by the gate eigenvalue transform matrix `gate_transform_matrix`.
"""
function calc_wls_merit_grad_log(
    d::Design,
    shot_weights::Vector{Float64},
    covariance_log_unweighted::SparseMatrixCSC{Float64, Int},
    gate_transform_matrix::SparseMatrixCSC{Float64, Int},
)
    # Initialise data
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    mapping_indices =
        [mapping_lower[idx]:mapping_upper[idx] for idx in 1:length(d.tuple_set)]
    shot_weights_factor =
        get_shot_weights_factor(shot_weights, d.tuple_times, mapping_lengths)
    shot_weights_local_grad = get_shot_weights_local_grad(shot_weights, d.tuple_times)
    shot_weights_log_matrix = get_shot_weights_log_matrix(shot_weights)
    # Precompute useful quantities
    covariance_log_unweighted_matrix = Symmetric(Array(covariance_log_unweighted))
    covariance_log = covariance_log_unweighted * shot_weights_factor
    covariance_log_matrix = Symmetric(Array(covariance_log))
    covariance_log_diag_inv = sparse(Diagonal(covariance_log)^(-1))
    covariance_combined = covariance_log * covariance_log_diag_inv
    conjugated_inv = Symmetric(
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_diag_inv * d.matrix)))),
    )
    wls_estimator_prime = conjugated_inv * Array(d.matrix' * covariance_log_diag_inv)
    wls_estimator_prime_commutant =
        covariance_combined * d.matrix * wls_estimator_prime - covariance_combined
    # Calculate gate transform quantities
    N = size(gate_transform_matrix, 1)
    wls_estimator = gate_transform_matrix * wls_estimator_prime
    sigma = Symmetric(wls_estimator * covariance_log_matrix * wls_estimator')
    wls_gram = Symmetric(wls_estimator' * wls_estimator)
    wls_sq_gram = Symmetric(wls_gram * covariance_log_matrix * wls_gram)
    # Compute the trace of the covariance matrix and its square
    sigma_tr = tr(sigma)
    sigma_sq_tr = tr(sigma^2)
    merit = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    # Compute the gradient of the trace of the covariance matrix
    sigma_tr_diag = diag(
        Symmetric(wls_gram + 2 * Diagonal(wls_gram * wls_estimator_prime_commutant)) *
        covariance_log_unweighted_matrix,
    )
    sigma_tr_partial = [sum(sigma_tr_diag[indices]) for indices in mapping_indices]
    sigma_tr_grad =
        sum(sigma_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_tr_partial .* shot_weights_local_grad
    # Compute the gradient of the trace of the square of the covariance matrix
    sigma_sq_tr_diag =
        2 * diag(
            Symmetric(
                wls_sq_gram + 2 * Diagonal(wls_sq_gram * wls_estimator_prime_commutant),
            ) * covariance_log_unweighted_matrix,
        )
    sigma_sq_tr_partial = [sum(sigma_sq_tr_diag[indices]) for indices in mapping_indices]
    sigma_sq_tr_grad =
        sum(sigma_sq_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_sq_tr_partial .* shot_weights_local_grad
    # Calculate the gradient of the figure of merit with respect to the shot weights
    merit_grad = get_merit_grad(sigma_tr, sigma_tr_grad, sigma_sq_tr, sigma_sq_tr_grad, N)
    # Calculate the gradient of the figure of merit with respect to the log shot weights
    merit_grad_log = shot_weights_log_matrix * merit_grad
    return (merit_grad_log::Vector{Float64}, merit::Float64)
end

"""
    wls_optimise_weights(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; options::OptimOptions = OptimOptions())

Returns versions of the design `d` and circuit log-eigenvalue estimator covariance matrix `covariance_log` after optimising the shot weights with respect to the weighted least squares (WLS) figure of merit, alongside the figure of merit values at each step.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function wls_optimise_weights(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    ls_type = options.ls_type
    @assert ls_type == :wls "Inappropriate least squares optimisation type $(ls_type) supplied."
    est_type = options.est_type
    est_weight = options.est_weight
    learning_rate = options.learning_rate
    momentum = options.momentum
    learning_rate_scale_factor = options.learning_rate_scale_factor
    max_steps = options.max_steps
    convergence_threshold = options.convergence_threshold
    convergence_steps = options.convergence_steps
    diagnostics = options.grad_diagnostics
    # Initialise data
    tuple_number = length(d.tuple_set)
    shot_weights = project_simplex(d.shot_weights)
    mapping_lengths = length.(d.mapping_ensemble)
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
    covariance_log_unweighted = covariance_log * shot_weights_factor_inv
    gate_eigenvalues_diag = sparse(Diagonal(get_gate_eigenvalues(d)))
    if est_type == :sum || est_type == :prod
        ord_gate_transform_matrix = get_transform(d, :ordinary) * gate_eigenvalues_diag
        rel_gate_transform_matrix = get_transform(d, :relative) * gate_eigenvalues_diag
    else
        gate_transform_matrix = get_transform(d, est_type) * gate_eigenvalues_diag
    end
    # Perform gradient descent
    stepping = true
    step = 1
    recently_pruned = 0
    recently_zeroed = 0
    scaled_learning_rate = learning_rate
    old_shot_weights = deepcopy(shot_weights)
    velocity = zeros(tuple_number)
    merit_descent = Vector{Float64}()
    while stepping
        # Calculate the gradient of the figure of merit
        nesterov_log_shot_weights = -log.(shot_weights) + momentum * velocity
        nesterov_shot_weights =
            exp.(-nesterov_log_shot_weights) / sum(exp.(-nesterov_log_shot_weights))
        if est_type == :sum || est_type == :prod
            (ord_merit_grad_log, ord_merit) = calc_wls_merit_grad_log(
                d,
                nesterov_shot_weights,
                covariance_log_unweighted,
                ord_gate_transform_matrix,
            )
            (rel_merit_grad_log, rel_merit) = calc_wls_merit_grad_log(
                d,
                nesterov_shot_weights,
                covariance_log_unweighted,
                rel_gate_transform_matrix,
            )
            if est_type == :sum
                merit = est_weight * ord_merit + (1 - est_weight) * rel_merit
                merit_grad_log =
                    est_weight * ord_merit_grad_log + (1 - est_weight) * rel_merit_grad_log
            elseif est_type == :prod
                merit = ord_merit^(est_weight) * rel_merit^(1 - est_weight)
                merit_grad_log =
                    merit * (
                        est_weight * ord_merit_grad_log / ord_merit +
                        (1 - est_weight) * rel_merit_grad_log / rel_merit
                    )
            else
                throw(error("Unsupported estimator type $(est_type)."))
            end
        else
            (merit_grad_log, merit) = calc_wls_merit_grad_log(
                d,
                nesterov_shot_weights,
                covariance_log_unweighted,
                gate_transform_matrix,
            )
        end
        push!(merit_descent, merit)
        # Update the shot weights, ensuring both they and the gradient are appropriately normalised
        velocity = momentum * velocity - scaled_learning_rate * merit_grad_log
        log_shot_weights_update = -log.(shot_weights) + velocity
        shot_weights_update =
            exp.(-log_shot_weights_update) / sum(exp.(-log_shot_weights_update))
        if (length(merit_descent) >= 2) && (merit_descent[end] > merit_descent[end - 1])
            shot_weights = old_shot_weights
            velocity = zeros(tuple_number)
            if recently_zeroed > 0
                scaled_learning_rate /= learning_rate_scale_factor
            end
            if diagnostics
                println(
                    "The updated shot weights worsened the figure of merit in step $(step); the velocity and shot weights have been reset$(recently_zeroed > 0 ? " and the learning rate has been reduced." : ".")",
                )
            end
            recently_pruned += 1
            recently_zeroed = convergence_steps
        else
            old_shot_weights = deepcopy(nesterov_shot_weights)
            shot_weights = shot_weights_update
        end
        # Check the step count and convergence
        if step >= convergence_steps
            merit_converged =
                all([
                    abs(1 - merit_descent[idx_2] / merit_descent[idx_1]) <
                    convergence_steps * convergence_threshold for
                    idx_1 in (step - convergence_steps + 1):step for
                    idx_2 in (step - convergence_steps + 1):step
                ]) && (
                    abs(1 - merit_descent[end - 1] / merit_descent[end]) <
                    convergence_threshold
                )
        else
            merit_converged = false
        end
        if (merit_converged || step >= max_steps) && recently_pruned == 0
            stepping = false
            if diagnostics
                if merit_converged
                    println(
                        "Converged after $(step) steps. The $(ls_type) $(est_type) figure of merit is $(round(merit_descent[step], sigdigits = 5)).",
                    )
                else
                    println(
                        "The maximum number of steps $(max_steps) has been reached without convergence. The $(ls_type) $(est_type) figure of merit is $(round(merit_descent[step], sigdigits = 5))$(max_steps > 1 ? ", which differs from the previous step by $(round(abs(merit_descent[step] - merit_descent[step - 1]), sigdigits = 5)), whereas the threshold for convergence is $(convergence_threshold)" : "").",
                    )
                end
            end
        else
            if diagnostics
                println(
                    "The $(ls_type) $(est_type) figure of merit was $(round(merit_descent[step], sigdigits = 5)) prior to step $(step).",
                )
            end
            step += 1
            if recently_pruned > 0
                recently_pruned -= 1
            end
            if recently_zeroed > 0
                recently_zeroed -= 1
            end
        end
    end
    # Update the covariance matrix
    shot_weights_factor =
        get_shot_weights_factor(shot_weights, d.tuple_times, mapping_lengths)
    covariance_log = covariance_log_unweighted * shot_weights_factor
    # Update the design
    @reset d.shot_weights = project_simplex(shot_weights)
    @reset d.ls_type = :wls
    return (
        d::Design,
        covariance_log::SparseMatrixCSC{Float64, Int},
        merit_descent::Vector{Float64},
    )
end

"""
    calc_ols_merit_grad_log(d::Design, shot_weights::Vector{Float64}, ols_estimator::Matrix{Float64}, ols_estimator_covariance::Matrix{Float64}, ols_gram_covariance::Matrix{Float64})

Returns the gradient of the ordinary least squares (OLS) figure of merit for the design `d` with respect to the logarithms of the shot weights `shot_weights`, using the OLS estimator matrix `ols_estimator`, scaled by the unweighted covariance matrix in `ols_estimator_covariance`, and the OLS Gram matrix also scaled by the unweighted covariance matrix in `ols_gram_covariance`.
"""
function calc_ols_merit_grad_log(
    d::Design,
    shot_weights::Vector{Float64},
    ols_estimator::Matrix{Float64},
    ols_estimator_covariance::Matrix{Float64},
    ols_gram_covariance::Matrix{Float64},
)
    # Initialise data
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    mapping_indices =
        [mapping_lower[idx]:mapping_upper[idx] for idx in 1:length(d.tuple_set)]
    shot_weights_factor =
        get_shot_weights_factor(shot_weights, d.tuple_times, mapping_lengths)
    shot_weights_local_grad = get_shot_weights_local_grad(shot_weights, d.tuple_times)
    shot_weights_log_matrix = get_shot_weights_log_matrix(shot_weights)
    # Calculate gate transform quantities
    N = size(ols_estimator, 1)
    # Compute the trace of the covariance matrix and its square
    sigma = ols_estimator_covariance * shot_weights_factor * ols_estimator'
    sigma_tr = tr(sigma)
    sigma_sq_tr = tr(sigma^2)
    merit = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    # Compute the gradient of the trace of the covariance matrix
    sigma_tr_diag = diag(ols_gram_covariance)
    sigma_tr_partial = [sum(sigma_tr_diag[indices]) for indices in mapping_indices]
    sigma_tr_grad =
        sum(sigma_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_tr_partial .* shot_weights_local_grad
    # Compute the gradient of the trace of the square of the covariance matrix
    sigma_sq_tr_diag =
        2 * diag(ols_gram_covariance * shot_weights_factor * ols_gram_covariance)
    sigma_sq_tr_partial = [sum(sigma_sq_tr_diag[indices]) for indices in mapping_indices]
    sigma_sq_tr_grad =
        sum(sigma_sq_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_sq_tr_partial .* shot_weights_local_grad
    # Calculate the gradient of the figure of merit with respect to the shot weights
    merit_grad = get_merit_grad(sigma_tr, sigma_tr_grad, sigma_sq_tr, sigma_sq_tr_grad, N)
    # Calculate the gradient of the figure of merit with respect to the log shot weights
    merit_grad_log = shot_weights_log_matrix * merit_grad
    return (merit_grad_log::Vector{Float64}, merit::Float64)
end

"""
    ols_optimise_weights(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; options::OptimOptions = OptimOptions())

Returns versions of the design `d` and circuit log-eigenvalue estimator covariance matrix `covariance_log` after optimising the shot weights with respect to the ordinary least squares (OLS) figure of merit, alongside the figure of merit values at each step.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function ols_optimise_weights(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    ls_type = options.ls_type
    @assert ls_type == :ols "Inappropriate least squares optimisation type $(ls_type) supplied."
    est_type = options.est_type
    est_weight = options.est_weight
    learning_rate = options.learning_rate
    momentum = options.momentum
    learning_rate_scale_factor = options.learning_rate_scale_factor
    max_steps = options.max_steps
    convergence_threshold = options.convergence_threshold
    convergence_steps = options.convergence_steps
    diagnostics = options.grad_diagnostics
    # Initialise data
    tuple_number = length(d.tuple_set)
    shot_weights = project_simplex(d.shot_weights)
    mapping_lengths = length.(d.mapping_ensemble)
    gate_eigenvalues_diag = sparse(Diagonal(get_gate_eigenvalues(d)))
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
    covariance_log_unweighted = covariance_log * shot_weights_factor_inv
    if est_type == :sum || est_type == :prod
        ord_gate_transform_matrix = get_transform(d, :ordinary) * gate_eigenvalues_diag
        ord_ols_estimator =
            ord_gate_transform_matrix *
            inv(bunchkaufman(Symmetric(Array(d.matrix' * d.matrix)))) *
            d.matrix'
        ord_ols_estimator_covariance = ord_ols_estimator * covariance_log_unweighted
        ord_ols_gram_covariance = ord_ols_estimator' * ord_ols_estimator_covariance
        rel_gate_transform_matrix = get_transform(d, :relative) * gate_eigenvalues_diag
        rel_ols_estimator =
            rel_gate_transform_matrix *
            inv(bunchkaufman(Symmetric(Array(d.matrix' * d.matrix)))) *
            d.matrix'
        rel_ols_estimator_covariance = rel_ols_estimator * covariance_log_unweighted
        rel_ols_gram_covariance = rel_ols_estimator' * rel_ols_estimator_covariance
    else
        gate_transform_matrix = get_transform(d, est_type) * gate_eigenvalues_diag
        ols_estimator =
            gate_transform_matrix *
            inv(bunchkaufman(Symmetric(Array(d.matrix' * d.matrix)))) *
            d.matrix'
        ols_estimator_covariance = ols_estimator * covariance_log_unweighted
        ols_gram_covariance = ols_estimator' * ols_estimator_covariance
    end
    # Perform gradient descent
    stepping = true
    step = 1
    recently_pruned = 0
    recently_zeroed = 0
    scaled_learning_rate = learning_rate
    old_shot_weights = deepcopy(shot_weights)
    velocity = zeros(tuple_number)
    merit_descent = Vector{Float64}()
    while stepping
        # Calculate the gradient of the figure of merit
        nesterov_log_shot_weights = -log.(shot_weights) + momentum * velocity
        nesterov_shot_weights =
            exp.(-nesterov_log_shot_weights) / sum(exp.(-nesterov_log_shot_weights))
        if est_type == :sum || est_type == :prod
            (ord_merit_grad_log, ord_merit) = calc_ols_merit_grad_log(
                d,
                nesterov_shot_weights,
                ord_ols_estimator,
                ord_ols_estimator_covariance,
                ord_ols_gram_covariance,
            )
            (rel_merit_grad_log, rel_merit) = calc_ols_merit_grad_log(
                d,
                nesterov_shot_weights,
                rel_ols_estimator,
                rel_ols_estimator_covariance,
                rel_ols_gram_covariance,
            )
            if est_type == :sum
                merit = est_weight * ord_merit + (1 - est_weight) * rel_merit
                merit_grad_log =
                    est_weight * ord_merit_grad_log + (1 - est_weight) * rel_merit_grad_log
            elseif est_type == :prod
                merit = ord_merit^(est_weight) * rel_merit^(1 - est_weight)
                merit_grad_log =
                    merit * (
                        est_weight * ord_merit_grad_log / ord_merit +
                        (1 - est_weight) * rel_merit_grad_log / rel_merit
                    )
            else
                throw(error("Unsupported estimator type $(est_type)."))
            end
        else
            (merit_grad_log, merit) = calc_ols_merit_grad_log(
                d,
                nesterov_shot_weights,
                ols_estimator,
                ols_estimator_covariance,
                ols_gram_covariance,
            )
        end
        push!(merit_descent, merit)
        # Update the shot weights, ensuring both they and the gradient are appropriately normalised
        velocity = momentum * velocity - scaled_learning_rate * merit_grad_log
        log_shot_weights_update = -log.(shot_weights) + velocity
        shot_weights_update =
            exp.(-log_shot_weights_update) / sum(exp.(-log_shot_weights_update))
        if (length(merit_descent) >= 2) && (merit_descent[end] > merit_descent[end - 1])
            shot_weights = old_shot_weights
            velocity = zeros(tuple_number)
            if recently_zeroed > 0
                scaled_learning_rate /= learning_rate_scale_factor
            end
            if diagnostics
                println(
                    "The updated shot weights worsened the figure of merit in step $(step); the velocity and shot weights have been reset$(recently_zeroed > 0 ? " and the learning rate has been reduced." : ".")",
                )
            end
            recently_pruned += 1
            recently_zeroed = convergence_steps
        else
            old_shot_weights = deepcopy(nesterov_shot_weights)
            shot_weights = shot_weights_update
        end
        # Check the step count and convergence
        if step >= convergence_steps
            merit_converged =
                all([
                    abs(1 - merit_descent[idx_2] / merit_descent[idx_1]) <
                    convergence_steps * convergence_threshold for
                    idx_1 in (step - convergence_steps + 1):step for
                    idx_2 in (step - convergence_steps + 1):step
                ]) && (
                    abs(1 - merit_descent[end - 1] / merit_descent[end]) <
                    convergence_threshold
                )
        else
            merit_converged = false
        end
        if (merit_converged || step >= max_steps) && recently_pruned == 0
            stepping = false
            if diagnostics
                if merit_converged
                    println(
                        "Converged after $(step) steps. The $(ls_type) $(est_type) figure of merit is $(round(merit_descent[step], sigdigits = 5)).",
                    )
                else
                    println(
                        "The maximum number of steps $(max_steps) has been reached without convergence. The $(ls_type) $(est_type) figure of merit is $(round(merit_descent[step], sigdigits = 5))$(max_steps > 1 ? ", which differs from the previous step by $(round(abs(merit_descent[step] - merit_descent[step - 1]), sigdigits = 5)), whereas the threshold for convergence is $(convergence_threshold)" : "").",
                    )
                end
            end
        else
            if diagnostics
                println(
                    "The $(ls_type) $(est_type) figure of merit was $(round(merit_descent[step], sigdigits = 5)) prior to step $(step).",
                )
            end
            step += 1
            if recently_pruned > 0
                recently_pruned -= 1
            end
            if recently_zeroed > 0
                recently_zeroed -= 1
            end
        end
    end
    # Update the covariance matrix
    shot_weights_factor =
        get_shot_weights_factor(shot_weights, d.tuple_times, mapping_lengths)
    covariance_log = covariance_log_unweighted * shot_weights_factor
    # Update the design
    @reset d.shot_weights = project_simplex(shot_weights)
    @reset d.ls_type = :ols
    return (
        d::Design,
        covariance_log::SparseMatrixCSC{Float64, Int},
        merit_descent::Vector{Float64},
    )
end

"""
    optimise_weights(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; options::OptimOptions = OptimOptions())

Returns versions of the design `d` and circuit log-eigenvalue estimator covariance matrix `covariance_log` after optimising the shot weights with respect to the figure of merit, alongside the figure of merit values at each step.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`, which in particular specifies the least squares estimator type for which the figure of merit is calculated.
"""
function optimise_weights(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    ls_type = options.ls_type
    max_steps = options.max_steps
    # Perform gradient descent
    if max_steps == 0
        merit_descent = Float64[]
    elseif ls_type == :gls
        (d, covariance_log, merit_descent) =
            gls_optimise_weights(d, covariance_log; options = options)
    elseif ls_type == :wls
        (d, covariance_log, merit_descent) =
            wls_optimise_weights(d, covariance_log; options = options)
    elseif ls_type == :ols
        (d, covariance_log, merit_descent) =
            ols_optimise_weights(d, covariance_log; options = options)
    else
        error("Must use a valid least squares type.")
    end
    return (
        d::Design,
        covariance_log::SparseMatrixCSC{Float64, Int},
        merit_descent::Vector{Float64},
    )
end

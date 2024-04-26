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
    tuple_times_factor = sum(shot_weights .* tuple_times)
    shot_weights_factor =
        tuple_times_factor * Diagonal(
            vcat(
                [
                    (1 / shot_weights[idx]) * ones(mapping_lengths[idx]) for
                    idx in 1:tuple_number
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
    tuple_times_factor = sum(shot_weights .* tuple_times)
    shot_weights_factor_inv =
        (1 / tuple_times_factor) * Diagonal(
            vcat(
                [
                    shot_weights[idx] * ones(mapping_lengths[idx]) for idx in 1:tuple_number
                ]...,
            ),
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
    tuple_times_factor = sum(shot_weights .* tuple_times)
    shot_weights_local_grad = -tuple_times_factor ./ (shot_weights .^ 2)
    return shot_weights_local_grad::Vector{Float64}
end

"""
    get_merit_grad(sigma_tr::Float64, sigma_tr_grad::Vector{Float64}, sigma_sq_tr::Float64, sigma_sq_tr_grad::Vector{Float64}, N::Int)

Returns the gradient of the figure of merit with respect to the shot weights given the trace of the gate eigenvalue estimator covariance matrix `sigma_tr` and its square `sigma_sq_tr`, the gradient with respect to the shot weights `sigma_tr_grad` and the gradient of the square `sigma_sq_tr_grad`, and the number of gate eigenvalues `N`.
"""
function get_merit_grad(
    sigma_tr::Float64,
    sigma_tr_grad::Vector{Float64},
    sigma_sq_tr::Float64,
    sigma_sq_tr_grad::Vector{Float64},
    N::Int,
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
    calc_gls_merit_grad_log(d::Design, shot_weights::Vector{Float64}, covariance_log_unweighted_inv::SparseMatrixCSC{Float64, Int})

Returns the gradient of the generalised least squares (GLS) figure of merit for the design `d` with respect to the logarithms of the shot weights `shot_weights`, using the inverse of the unweighted (by the shot weights factor) covariance matrix of the circuit log-eigenvalue estimator `covariance_log_unweighted_inv`.
"""
function calc_gls_merit_grad_log(
    d::Design,
    shot_weights::Vector{Float64},
    covariance_log_unweighted_inv::SparseMatrixCSC{Float64, Int},
)
    # Initialise data
    N = d.c.N
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    mapping_indices =
        [mapping_lower[idx]:mapping_upper[idx] for idx in 1:length(d.tuple_set)]
    gate_eigenvalues = d.c.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(shot_weights, d.tuple_times, mapping_lengths)
    shot_weights_local_grad = get_shot_weights_local_grad(shot_weights, d.tuple_times)
    # Compute the trace of the covariance matrix and its square
    covariance_log_inv = covariance_log_unweighted_inv * shot_weights_factor_inv
    sigma_prime = Symmetric(
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_inv * d.matrix)))),
    )
    sigma = Symmetric(gate_eigenvalues_diag * sigma_prime * gate_eigenvalues_diag)
    sigma_tr = tr(sigma)
    sigma_sq_tr = tr(sigma^2)
    merit = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    # Compute the gradient of the trace of the covariance matrix
    gate_sigma_prime = gate_eigenvalues_diag^2 * sigma_prime
    gls_matrix = covariance_log_inv * d.matrix * sigma_prime * gate_sigma_prime
    sigma_tr_diag = diag(gls_matrix * d.matrix' * shot_weights_factor_inv)
    sigma_tr_partial = [sum(sigma_tr_diag[indices]) for indices in mapping_indices]
    sigma_tr_grad =
        sum(sigma_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_tr_partial .* shot_weights_local_grad
    # Compute the gradient of the trace of the square of the covariance matrix
    sigma_sq_tr_diag =
        2 * diag(gls_matrix * gate_sigma_prime * d.matrix' * shot_weights_factor_inv)
    sigma_sq_tr_partial = [sum(sigma_sq_tr_diag[indices]) for indices in mapping_indices]
    sigma_sq_tr_grad =
        sum(sigma_sq_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_sq_tr_partial .* shot_weights_local_grad
    # Calculate the gradient of the figure of merit with respect to the shot weights
    merit_grad = get_merit_grad(sigma_tr, sigma_tr_grad, sigma_sq_tr, sigma_sq_tr_grad, N)
    # Calculate the gradient of the figure of merit with respect to the log shot weights
    shot_weights_log_matrix = get_shot_weights_log_matrix(shot_weights)
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
    @assert options.ls_type == :gls "Inappropriate least squares optimisation type $(options.ls_type) supplied."
    learning_rate = options.learning_rate
    momentum = options.momentum
    learning_rate_scale_factor = options.learning_rate_scale_factor
    shot_weights_clip = options.shot_weights_clip
    max_steps = options.max_steps
    convergence_threshold = options.convergence_threshold
    convergence_steps = options.convergence_steps
    diagnostics = options.grad_diagnostics
    # Initialise data
    N = d.c.N
    tuple_number = length(d.tuple_set)
    shot_weights = project_simplex(d.shot_weights)
    mapping_lengths = length.(d.mapping_ensemble)
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
    covariance_log_unweighted = covariance_log * shot_weights_factor_inv
    covariance_log_unweighted_inv =
        sparse_covariance_inv(covariance_log_unweighted, mapping_lengths)
    # Perform gradient descent
    stepping = true
    step = 1
    recently_pruned = 0
    recently_zeroed = 0
    scaled_learning_rate = learning_rate
    unprunable = Int[]
    old_shot_weights = deepcopy(shot_weights)
    velocity = zeros(tuple_number)
    merit_descent = Vector{Float64}(undef, 0)
    while stepping
        # Calculate the gradient of the figure of merit
        nesterov_log_shot_weights = -log.(shot_weights) + momentum * velocity
        nesterov_shot_weights =
            exp.(-nesterov_log_shot_weights) / sum(exp.(-nesterov_log_shot_weights))
        (merit_grad_log, merit) =
            calc_gls_merit_grad_log(d, nesterov_shot_weights, covariance_log_unweighted_inv)
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
        # Prune a tuple from the design if one has a shot weight below the minimum
        weights_below_min = findall(shot_weights .< shot_weights_clip)
        if length(weights_below_min) > 0 && recently_pruned == 0
            min_weight_idx = findmin(shot_weights[weights_below_min])[2]
            prune_idx = weights_below_min[min_weight_idx]
            if prune_idx ∉ unprunable
                (d_prune, covariance_log_unweighted_prune) = prune_design(
                    d,
                    covariance_log_unweighted,
                    prune_idx;
                    update_weights = false,
                )
                # The pruned design must still achieve full rank
                if rank(Array(d_prune.matrix)) == N
                    prune_indices = setdiff(1:tuple_number, prune_idx)
                    prune_shot_weights_factor = get_shot_weights_factor(
                        d_prune.shot_weights,
                        d_prune.tuple_times,
                        mapping_lengths[prune_indices],
                    )
                    covariance_log_prune =
                        covariance_log_unweighted_prune * prune_shot_weights_factor
                    prune_merit = calc_gls_moments(d_prune, covariance_log_prune)[1]
                    prune_improve = all([
                        prune_merit < merit_descent[idx] for
                        idx in max(1, step - convergence_steps):step
                    ])
                    # Prune the tuple only if doing so improves the figure of merit
                    if prune_improve
                        # Prune the quantities
                        prune_indices = setdiff(1:tuple_number, prune_idx)
                        mapping_lengths = mapping_lengths[prune_indices]
                        d = d_prune
                        covariance_log_unweighted = covariance_log_unweighted_prune
                        covariance_log_unweighted_inv = sparse_covariance_inv(
                            covariance_log_unweighted,
                            mapping_lengths,
                        )
                        # Update the shot weights and reset the velocity
                        tuple_number -= 1
                        shot_weights =
                            shot_weights[prune_indices] / sum(shot_weights[prune_indices])
                        old_shot_weights = deepcopy(shot_weights)
                        velocity = velocity[prune_indices]
                        recently_pruned = convergence_steps
                        if diagnostics
                            println(
                                "Pruned tuple $(prune_idx) of $(tuple_number + 1) in step $(step); the velocity has been reset.",
                            )
                        end
                    end
                else
                    push!(unprunable, prune_idx)
                end
            end
        end
        # Check the step count and convergence
        if step > convergence_steps
            merit_converged =
                all([
                    abs(merit_descent[idx_1] - merit_descent[idx_2]) <
                    convergence_steps * convergence_threshold for
                    idx_1 in (step - convergence_steps):step for
                    idx_2 in (step - convergence_steps):step
                ]) &&
                (abs(merit_descent[end] - merit_descent[end - 1]) < convergence_threshold)
        else
            merit_converged = false
        end
        if (merit_converged || step >= max_steps) && recently_pruned == 0
            stepping = false
            if diagnostics
                if merit_converged
                    println(
                        "Converged after $(step) steps. The GLS figure of merit is $(round(merit_descent[step], sigdigits = 5)).",
                    )
                else
                    println(
                        "The maximum number of steps $(max_steps) has been reached without convergence. The GLS figure of merit is $(round(merit_descent[step], sigdigits = 5)), which differs from the previous step by $(round(abs(merit_descent[step] - merit_descent[step-1]), sigdigits=5)), whereas the threshold for convergence is $(convergence_threshold).",
                    )
                end
            end
        else
            if diagnostics
                println(
                    "The GLS figure of merit was $(round(merit_descent[step], sigdigits = 5)) prior to step $(step).",
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
    @reset d.shot_weights = shot_weights
    @reset d.ls_type = :gls
    return (
        d::Design,
        covariance_log::SparseMatrixCSC{Float64, Int},
        merit_descent::Vector{Float64},
    )
end

"""
    calc_wls_merit_grad_log(d::Design, shot_weights::Vector{Float64}, covariance_log_unweighted::SparseMatrixCSC{Float64, Int})

Returns the gradient of the weighted least squares (WLS) figure of merit for the design `d` with respect to the logarithms of the shot weights `shot_weights`, using the unweighted (by the shot weights factor) covariance matrix of the circuit log-eigenvalue estimator `covariance_log_unweighted`.
"""
function calc_wls_merit_grad_log(
    d::Design,
    shot_weights::Vector{Float64},
    covariance_log_unweighted::SparseMatrixCSC{Float64, Int},
)
    # Initialise data
    N = d.c.N
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    mapping_indices =
        [mapping_lower[idx]:mapping_upper[idx] for idx in 1:length(d.tuple_set)]
    gate_eigenvalues = d.c.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    shot_weights_factor =
        get_shot_weights_factor(shot_weights, d.tuple_times, mapping_lengths)
    shot_weights_local_grad = get_shot_weights_local_grad(shot_weights, d.tuple_times)
    # Compute the trace of the covariance matrix and its square
    covariance_log = covariance_log_unweighted * shot_weights_factor
    covariance_log_matrix = Symmetric(Array(covariance_log))
    covariance_log_diag_inv = Diagonal(covariance_log)^(-1)
    covariance_combined = Array(covariance_log * covariance_log_diag_inv)
    conjugated_inv = Symmetric(
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_diag_inv * d.matrix)))),
    )
    wls_estimator_prime = conjugated_inv * Array(d.matrix' * covariance_log_diag_inv)
    wls_estimator_prime_commutant = d.matrix * wls_estimator_prime - I
    wls_estimator = gate_eigenvalues_diag * wls_estimator_prime
    sigma = Symmetric(wls_estimator * covariance_log_matrix * wls_estimator')
    sigma_tr = tr(sigma)
    sigma_sq_tr = tr(sigma^2)
    merit = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    # Compute the gradient of the trace of the covariance matrix
    wls_gram = Symmetric(wls_estimator' * wls_estimator)
    wls_tr_matrix = Symmetric(
        wls_gram +
        2 * Diagonal(wls_gram * covariance_combined * wls_estimator_prime_commutant),
    )
    sigma_tr_diag = diag(wls_tr_matrix * Array(covariance_log_unweighted))
    sigma_tr_partial = [sum(sigma_tr_diag[indices]) for indices in mapping_indices]
    sigma_tr_grad =
        sum(sigma_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_tr_partial .* shot_weights_local_grad
    # Compute the gradient of the trace of the square of the covariance matrix
    wls_sq_gram = Symmetric(wls_gram * covariance_log_matrix * wls_gram)
    wls_sq_tr_matrix = Symmetric(
        2 * wls_sq_gram +
        4 * Diagonal(wls_sq_gram * covariance_combined * wls_estimator_prime_commutant),
    )
    sigma_sq_tr_diag = diag(wls_sq_tr_matrix * Array(covariance_log_unweighted))
    sigma_sq_tr_partial = [sum(sigma_sq_tr_diag[indices]) for indices in mapping_indices]
    sigma_sq_tr_grad =
        sum(sigma_sq_tr_partial ./ shot_weights) * d.tuple_times +
        sigma_sq_tr_partial .* shot_weights_local_grad
    # Calculate the gradient of the figure of merit with respect to the shot weights
    merit_grad = get_merit_grad(sigma_tr, sigma_tr_grad, sigma_sq_tr, sigma_sq_tr_grad, N)
    # Calculate the gradient of the figure of merit with respect to the log shot weights
    shot_weights_log_matrix = get_shot_weights_log_matrix(shot_weights)
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
    @assert options.ls_type == :wls "Inappropriate least squares optimisation type $(options.ls_type) supplied."
    learning_rate = options.learning_rate
    momentum = options.momentum
    learning_rate_scale_factor = options.learning_rate_scale_factor
    shot_weights_clip = options.shot_weights_clip
    max_steps = options.max_steps
    convergence_threshold = options.convergence_threshold
    convergence_steps = options.convergence_steps
    diagnostics = options.grad_diagnostics
    # Initialise data
    N = d.c.N
    tuple_number = length(d.tuple_set)
    shot_weights = project_simplex(d.shot_weights)
    mapping_lengths = length.(d.mapping_ensemble)
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
    covariance_log_unweighted = covariance_log * shot_weights_factor_inv
    # Perform gradient descent
    stepping = true
    step = 1
    recently_pruned = 0
    recently_zeroed = 0
    scaled_learning_rate = learning_rate
    unprunable = Int[]
    old_shot_weights = deepcopy(shot_weights)
    velocity = zeros(tuple_number)
    merit_descent = Vector{Float64}(undef, 0)
    while stepping
        # Calculate the gradient of the figure of merit
        nesterov_log_shot_weights = -log.(shot_weights) + momentum * velocity
        nesterov_shot_weights =
            exp.(-nesterov_log_shot_weights) / sum(exp.(-nesterov_log_shot_weights))
        (merit_grad_log, merit) =
            calc_wls_merit_grad_log(d, nesterov_shot_weights, covariance_log_unweighted)
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
        # Prune a tuple from the design if one has a shot weight below the minimum
        weights_below_min = findall(shot_weights .< shot_weights_clip)
        if length(weights_below_min) > 0 && recently_pruned == 0
            min_weight_idx = findmin(shot_weights[weights_below_min])[2]
            prune_idx = weights_below_min[min_weight_idx]
            if prune_idx ∉ unprunable
                (d_prune, covariance_log_unweighted_prune) = prune_design(
                    d,
                    covariance_log_unweighted,
                    prune_idx;
                    update_weights = false,
                )
                # The pruned design must still achieve full rank
                if rank(Array(d_prune.matrix)) == N
                    prune_indices = setdiff(1:tuple_number, prune_idx)
                    prune_shot_weights_factor = get_shot_weights_factor(
                        d_prune.shot_weights,
                        d_prune.tuple_times,
                        mapping_lengths[prune_indices],
                    )
                    covariance_log_prune =
                        covariance_log_unweighted_prune * prune_shot_weights_factor
                    prune_merit = calc_wls_moments(d_prune, covariance_log_prune)[1]
                    prune_improve = all([
                        prune_merit < merit_descent[idx] for
                        idx in max(1, step - convergence_steps):step
                    ])
                    # Prune the tuple only if doing so improves the figure of merit
                    if prune_improve
                        # Prune the quantities
                        mapping_lengths = mapping_lengths[prune_indices]
                        d = d_prune
                        covariance_log_unweighted = covariance_log_unweighted_prune
                        # Update the shot weights and reset the velocity
                        tuple_number -= 1
                        shot_weights =
                            shot_weights[prune_indices] / sum(shot_weights[prune_indices])
                        old_shot_weights = deepcopy(shot_weights)
                        velocity = velocity[prune_indices]
                        recently_pruned = convergence_steps
                        if diagnostics
                            println(
                                "Pruned tuple $(prune_idx) of $(tuple_number + 1) in step $(step); the velocity has been reset.",
                            )
                        end
                    end
                else
                    push!(unprunable, prune_idx)
                end
            end
        end
        # Check the step count and convergence
        if step > convergence_steps
            merit_converged =
                all([
                    abs(merit_descent[idx_1] - merit_descent[idx_2]) <
                    convergence_steps * convergence_threshold for
                    idx_1 in (step - convergence_steps):step for
                    idx_2 in (step - convergence_steps):step
                ]) &&
                (abs(merit_descent[end] - merit_descent[end - 1]) < convergence_threshold)
        else
            merit_converged = false
        end
        if (merit_converged || step >= max_steps) && recently_pruned == 0
            stepping = false
            if diagnostics
                if merit_converged
                    println(
                        "Converged after $(step) steps. The WLS figure of merit is $(round(merit_descent[step], sigdigits = 5)).",
                    )
                else
                    println(
                        "The maximum number of steps $(max_steps) has been reached without convergence. The WLS figure of merit is $(round(merit_descent[step], sigdigits = 5)), which differs from the previous step by $(round(abs(merit_descent[step] - merit_descent[step-1]), sigdigits=5)), whereas the threshold for convergence is $(convergence_threshold).",
                    )
                end
            end
        else
            if diagnostics
                println(
                    "The WLS figure of merit was $(round(merit_descent[step], sigdigits = 5)) prior to step $(step).",
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
    @reset d.shot_weights = shot_weights
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
    N = d.c.N
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    mapping_indices =
        [mapping_lower[idx]:mapping_upper[idx] for idx in 1:length(d.tuple_set)]
    shot_weights_factor =
        get_shot_weights_factor(shot_weights, d.tuple_times, mapping_lengths)
    shot_weights_local_grad = get_shot_weights_local_grad(shot_weights, d.tuple_times)
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
    shot_weights_log_matrix = get_shot_weights_log_matrix(shot_weights)
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
    @assert options.ls_type == :ols "Inappropriate least squares optimisation type $(options.ls_type) supplied."
    learning_rate = options.learning_rate
    momentum = options.momentum
    learning_rate_scale_factor = options.learning_rate_scale_factor
    shot_weights_clip = options.shot_weights_clip
    max_steps = options.max_steps
    convergence_threshold = options.convergence_threshold
    convergence_steps = options.convergence_steps
    diagnostics = options.grad_diagnostics
    # Initialise data
    N = d.c.N
    tuple_number = length(d.tuple_set)
    shot_weights = project_simplex(d.shot_weights)
    mapping_lengths = length.(d.mapping_ensemble)
    gate_eigenvalues = d.c.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
    covariance_log_unweighted = covariance_log * shot_weights_factor_inv
    ols_estimator =
        gate_eigenvalues_diag *
        inv(bunchkaufman(Symmetric(Array(d.matrix' * d.matrix)))) *
        d.matrix'
    ols_estimator_covariance = ols_estimator * covariance_log_unweighted
    ols_gram_covariance = ols_estimator' * ols_estimator_covariance
    # Perform gradient descent
    stepping = true
    step = 1
    recently_pruned = 0
    recently_zeroed = 0
    scaled_learning_rate = learning_rate
    unprunable = Int[]
    old_shot_weights = deepcopy(shot_weights)
    velocity = zeros(tuple_number)
    merit_descent = Vector{Float64}(undef, 0)
    while stepping
        # Calculate the gradient of the figure of merit
        nesterov_log_shot_weights = -log.(shot_weights) + momentum * velocity
        nesterov_shot_weights =
            exp.(-nesterov_log_shot_weights) / sum(exp.(-nesterov_log_shot_weights))
        (merit_grad_log, merit) = calc_ols_merit_grad_log(
            d,
            nesterov_shot_weights,
            ols_estimator,
            ols_estimator_covariance,
            ols_gram_covariance,
        )
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
        # Prune a tuple from the design if one has a shot weight below the minimum
        weights_below_min = findall(shot_weights .< shot_weights_clip)
        if length(weights_below_min) > 0 && recently_pruned == 0
            min_weight_idx = findmin(shot_weights[weights_below_min])[2]
            prune_idx = weights_below_min[min_weight_idx]
            if prune_idx ∉ unprunable
                (d_prune, covariance_log_unweighted_prune) = prune_design(
                    d,
                    covariance_log_unweighted,
                    prune_idx;
                    update_weights = false,
                )
                # The pruned design must still achieve full rank
                if rank(Array(d_prune.matrix)) == N
                    prune_indices = setdiff(1:tuple_number, prune_idx)
                    prune_shot_weights_factor = get_shot_weights_factor(
                        d_prune.shot_weights,
                        d_prune.tuple_times,
                        mapping_lengths[prune_indices],
                    )
                    covariance_log_prune =
                        covariance_log_unweighted_prune * prune_shot_weights_factor
                    prune_merit = calc_ols_moments(d_prune, covariance_log_prune)[1]
                    prune_improve = all([
                        prune_merit < merit_descent[idx] for
                        idx in max(1, step - convergence_steps):step
                    ])
                    # Prune the tuple only if doing so improves the figure of merit
                    if prune_improve
                        # Prune the quantities
                        mapping_lengths = mapping_lengths[prune_indices]
                        d = d_prune
                        covariance_log_unweighted = covariance_log_unweighted_prune
                        ols_estimator =
                            gate_eigenvalues_diag *
                            inv(bunchkaufman(Symmetric(Array(d.matrix' * d.matrix)))) *
                            d.matrix'
                        ols_estimator_covariance = ols_estimator * covariance_log_unweighted
                        ols_gram_covariance = ols_estimator' * ols_estimator_covariance
                        # Update the shot weights and reset the velocity
                        tuple_number -= 1
                        shot_weights =
                            shot_weights[prune_indices] / sum(shot_weights[prune_indices])
                        old_shot_weights = deepcopy(shot_weights)
                        velocity = velocity[prune_indices]
                        recently_pruned = convergence_steps
                        if diagnostics
                            println(
                                "Pruned tuple $(prune_idx) of $(tuple_number + 1) in step $(step); the velocity has been reset.",
                            )
                        end
                    end
                else
                    push!(unprunable, prune_idx)
                end
            end
        end
        # Check the step count and convergence
        if step > convergence_steps
            merit_converged =
                all([
                    abs(merit_descent[idx_1] - merit_descent[idx_2]) <
                    convergence_steps * convergence_threshold for
                    idx_1 in (step - convergence_steps):step for
                    idx_2 in (step - convergence_steps):step
                ]) &&
                (abs(merit_descent[end] - merit_descent[end - 1]) < convergence_threshold)
        else
            merit_converged = false
        end
        if (merit_converged || step >= max_steps) && recently_pruned == 0
            stepping = false
            if diagnostics
                if merit_converged
                    println(
                        "Converged after $(step) steps. The OLS figure of merit is $(round(merit_descent[step], sigdigits = 5)).",
                    )
                else
                    println(
                        "The maximum number of steps $(max_steps) has been reached without convergence. The OLS figure of merit is $(round(merit_descent[step], sigdigits = 5)), which differs from the previous step by $(round(abs(merit_descent[step] - merit_descent[step-1]), sigdigits=5)), whereas the threshold for convergence is $(convergence_threshold).",
                    )
                end
            end
        else
            if diagnostics
                println(
                    "The OLS figure of merit was $(round(merit_descent[step], sigdigits = 5)) prior to step $(step).",
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
    @reset d.shot_weights = shot_weights
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

"""
    compare_ls_optimise_weights(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; gls_options::OptimOptions = OptimOptions(; ls_type = :gls), wls_options::OptimOptions = OptimOptions(; ls_type = :wls), ols_options::OptimOptions = OptimOptions(; ls_type = :ols))

Returns versions of the design `d` and circuit log-eigenvalue estimator covariance matrix `covariance_log` after optimising the shot weights with respect to the figure of merit, alongside the figure of merit values at each step and the final merits, for all three least squares estimator types.
The optimisation for each is parameterised by the [`OptimOptions`](@ref) objects `gls_options`, `wls_options`, and `ols_options`, respectively.
"""
function compare_ls_optimise_weights(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    gls_options::OptimOptions = OptimOptions(; ls_type = :gls),
    wls_options::OptimOptions = OptimOptions(; ls_type = :wls),
    ols_options::OptimOptions = OptimOptions(; ls_type = :ols),
)
    # Perform gradient descent for all LS estimators
    (d_gls, covariance_log_gls, merit_descent_gls) =
        gls_optimise_weights(d, covariance_log; options = gls_options)
    (d_wls, covariance_log_wls, merit_descent_wls) =
        wls_optimise_weights(d, covariance_log; options = wls_options)
    (d_ols, covariance_log_ols, merit_descent_ols) =
        ols_optimise_weights(d, covariance_log; options = ols_options)
    d_set = (d_gls, d_wls, d_ols)
    covariance_log_set = (covariance_log_gls, covariance_log_wls, covariance_log_ols)
    merit_descent_set = (merit_descent_gls, merit_descent_wls, merit_descent_ols)
    # Calculate the expectation and variance for all combinations of the optimised shot weights and LS estimators
    expectations = Matrix{Float64}(undef, 3, 3)
    variances = Matrix{Float64}(undef, 3, 3)
    (expectations[1, 1], variances[1, 1]) = calc_gls_moments(d_gls, covariance_log_gls)
    (expectations[2, 1], variances[2, 1]) = calc_gls_moments(d_wls, covariance_log_wls)
    (expectations[3, 1], variances[3, 1]) = calc_gls_moments(d_ols, covariance_log_ols)
    (expectations[1, 2], variances[1, 2]) = calc_wls_moments(d_gls, covariance_log_gls)
    (expectations[2, 2], variances[2, 2]) = calc_wls_moments(d_wls, covariance_log_wls)
    (expectations[3, 2], variances[3, 2]) = calc_wls_moments(d_ols, covariance_log_ols)
    (expectations[1, 3], variances[1, 3]) = calc_ols_moments(d_gls, covariance_log_gls)
    (expectations[2, 3], variances[2, 3]) = calc_ols_moments(d_wls, covariance_log_wls)
    (expectations[3, 3], variances[3, 3]) = calc_ols_moments(d_ols, covariance_log_ols)
    merit_array = hcat(expectations, sqrt.(variances))
    return (
        d_set::Tuple{Design, Design, Design},
        covariance_log_set::Tuple{
            SparseMatrixCSC{Float64, Int},
            SparseMatrixCSC{Float64, Int},
            SparseMatrixCSC{Float64, Int},
        },
        merit_descent_set::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}},
        merit_array::Matrix{Float64},
    )
end

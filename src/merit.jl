"""
    Merit

Merit parameters for an experimental design.

# Fields

  - `circuit_param::AbstractCircuitParameters`: Circuit parameters.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.
  - `ls_type::Symbol`: Type of least squares estimator for which the merit is reported.
  - `tuple_set::Vector{Vector{Int}}`: Set of tuples which arrange the circuit layers.
  - `tuple_set_data::TupleSetData`: [`TupleSetData`](@ref) object that generates the tuple set.
  - `tuple_times::Vector{Float64}`: Time taken to implement the circuit corresponding to each tuple, normalised according to the basic tuple set.
  - `shot_weights::Vector{Float64}`: Shot weights for each tuple in the set, which add to 1.
  - `experiment_numbers::Vector{Int}`: Number of experiments for each tuple in the set.
  - `experiment_number::Int`: Total number of experiments.
  - `G::Int`: Total number of gates.
  - `N::Int`: Number of gate eigenvalues.
  - `expectation::Float64`: Expected normalised RMS error for the gate eigenvalue estimator vector.
  - `variance::Float64`: Normalised RMS error variance for the gate eigenvalue estimator vector.
  - `eigenvalues::Vector{Float64}`: Eigenvalues of the gate eigenvalue estimator covariance matrix.
  - `cond_num::Float64`: Condition number of the design matrix, the ratio of the largest and smallest singular values.
  - `pinv_norm::Float64`: Pseudoinverse norm of the design matrix, the inverse of the smallest singular value.
"""
struct Merit
    circuit_param::AbstractCircuitParameters
    noise_param::AbstractNoiseParameters
    ls_type::Symbol
    tuple_set::Vector{Vector{Int}}
    tuple_set_data::TupleSetData
    tuple_times::Vector{Float64}
    shot_weights::Vector{Float64}
    experiment_numbers::Vector{Int}
    experiment_number::Int
    G::Int
    N::Int
    expectation::Float64
    variance::Float64
    eigenvalues::Vector{Float64}
    cond_num::Float64
    pinv_norm::Float64
end

function Base.show(io::IO, merit::Merit)
    return show(io, round.((merit.expectation, sqrt(merit.variance)), digits = 5))
end

@struct_hash_equal_isequal Merit

"""
    calc_covariance(d::Design, eigenvalues::Vector{Float64}, gate_eigenvalues::Vector{Float64}; warning::Bool = true)
    calc_covariance(d::Design, eigenvalues::Vector{Float64}; warning::Bool = true)
    calc_covariance(d::Design; warning::Bool = true)

Returns the circuit eigenvalue estimator covariance matrix for the design `d` with circuit eigenvalues `eigenvalues` and gate eigenvalues `gate_eigenvalues`.
If `warning` is `true`, warns that if `d.full_covariance` is `false` this will only generate the diagonal of the covariance matrix.
"""
function calc_covariance(
    d::Design,
    eigenvalues::Vector{Float64},
    gate_eigenvalues::Vector{Float64};
    warning::Bool = true,
)
    # Warn if the design does not have the full covariance matrix
    if warning && ~d.full_covariance
        @warn "The design does not include the full covariance matrix; only the diagonal will be calculated."
    end
    # Initialise some variables
    tuple_number = length(d.tuple_set)
    # Initialise a trivial mapping to check things are working properly
    trivial_pauli = Pauli(Bool[0], 0)
    trivial_row = SparseVector{Int32, Int32}(1, Int32[], Int32[])
    trivial_track = Vector{Int16}[]
    trivial_mapping = Mapping(trivial_pauli, trivial_pauli, trivial_row, trivial_track)
    # Generate the covariance matrix for the circuit eigenvalues
    M = size(d.matrix, 1)
    tuple_times_factor = sum(d.shot_weights .* d.tuple_times)
    covariance = spzeros(Float64, Int32, M, M)
    reentrant_lock = ReentrantLock()
    for idx in 1:tuple_number
        tuple_covariance_dict = d.covariance_dict_ensemble[idx]
        tuple_experiment_number = d.experiment_numbers[idx]
        shot_weight = d.shot_weights[idx]
        if idx == 1
            tuple_offset = 0
        else
            tuple_offset = sum([length(d.mapping_ensemble[j]) for j in 1:(idx - 1)])
        end
        # Generate the relevant term in the covariance matrix for each key in the dictionary
        @threads :static for key in collect(keys(tuple_covariance_dict))
            # If the key lies on the diagonal, calculate the variance, otherwise calculate the covariance, which is a little more complicated
            upper_key = CartesianIndex(tuple_offset + key[1], tuple_offset + key[2])
            lower_key = CartesianIndex(tuple_offset + key[2], tuple_offset + key[1])
            if upper_key == lower_key
                @assert key[1] == key[2]
                @assert tuple_covariance_dict[key][1] == trivial_mapping
                # Calculate the eigenvalue
                eigenvalue = eigenvalues[upper_key[1]]
                # Calculate the variance scaling factor
                eigenvalue_experiments = tuple_covariance_dict[key][2]
                scaling_var =
                    (tuple_times_factor / shot_weight) *
                    (tuple_experiment_number / eigenvalue_experiments)
                # Calculate the variance
                if eigenvalue < 1.0
                    eigenvalue_variance = scaling_var * (1 - eigenvalue^2)
                else
                    eigenvalue_variance =
                        scaling_var * (1 - maximum(eigenvalues[eigenvalues .< 1.0])^2)
                end
                # Set the covariance matrix term
                lock(reentrant_lock) do
                    covariance[upper_key] = eigenvalue_variance
                end
            else
                prod_mapping = tuple_covariance_dict[key][1]
                # Calculate the eigenvalues
                eigenvalue_1 = eigenvalues[upper_key[1]]
                eigenvalue_2 = eigenvalues[upper_key[2]]
                eigenvalue_both =
                    exp(-(prod_mapping.design_row' * (-log.(gate_eigenvalues))))
                # Calculate the covariance scaling factor
                key_1 = CartesianIndex(key[1], key[1])
                key_2 = CartesianIndex(key[2], key[2])
                eigenvalue_1_experiments = tuple_covariance_dict[key_1][2]
                eigenvalue_2_experiments = tuple_covariance_dict[key_2][2]
                eigenvalue_both_experiments = tuple_covariance_dict[key][2]
                scaling_cov =
                    (tuple_times_factor / shot_weight) * (
                        (tuple_experiment_number * eigenvalue_both_experiments) /
                        (eigenvalue_1_experiments * eigenvalue_2_experiments)
                    )
                # Calculate the covariance
                eigenvalue_both_covariance =
                    scaling_cov * (eigenvalue_both - eigenvalue_1 * eigenvalue_2)
                if eigenvalue_both_covariance < 0.0
                    @debug "The negative off-diagonal covariance matrix entry $(round(eigenvalue_both_covariance, sigdigits = 3)) has been set to 0; consider taking more shots."
                    eigenvalue_both_covariance = 0.0
                end
                # Set the covariance matrix entries
                lock(reentrant_lock) do
                    covariance[upper_key] = eigenvalue_both_covariance
                    covariance[lower_key] = eigenvalue_both_covariance
                end
            end
        end
    end
    return covariance::SparseMatrixCSC{Float64, Int32}
end
function calc_covariance(d::Design, eigenvalues::Vector{Float64}; warning::Bool = true)
    # Calculate the covariance matrix
    gate_eigenvalues = d.c.gate_eigenvalues
    covariance = calc_covariance(d, eigenvalues, gate_eigenvalues; warning = warning)
    return covariance::SparseMatrixCSC{Float64, Int32}
end
function calc_covariance(d::Design; warning::Bool = true)
    # Generate the circuit eigenvalues
    gate_eigenvalues = d.c.gate_eigenvalues
    eigenvalues = exp.(-(d.matrix * (-log.(gate_eigenvalues))))
    # Calculate the covariance matrix
    covariance = calc_covariance(d, eigenvalues; warning = warning)
    return covariance::SparseMatrixCSC{Float64, Int32}
end

"""
    calc_eigenvalues_covariance(d::Design; warning::Bool = true)

Returns the circuit eigenvalues and circuit eigenvalue estimator covariance matrix for the design `d`.
If `warning` is `true`, warns that if `d.full_covariance` is `false` this will only generate the diagonal of the covariance matrix.
"""
function calc_eigenvalues_covariance(d::Design; warning::Bool = true)
    # Generate the circuit eigenvalues
    gate_eigenvalues = d.c.gate_eigenvalues
    eigenvalues = exp.(-(d.matrix * (-log.(gate_eigenvalues))))
    # Calculate the covariance matrix
    covariance = calc_covariance(d, eigenvalues; warning = warning)
    return (eigenvalues::Vector{Float64}, covariance::SparseMatrixCSC{Float64, Int32})
end

"""
    calc_covariance_log(d::Design; warning::Bool = true)

Returns the covariance matrix of the circuit log-eigenvalues for the design `d`.
If `warning` is `true`, warns that if `d.full_covariance` is `false` this will only generate the diagonal of the covariance matrix.
"""
function calc_covariance_log(d::Design; warning::Bool = true)
    # Generate the circuit eigenvalues and covariance matrix
    (eigenvalues, covariance) = calc_eigenvalues_covariance(d; warning = warning)
    # Use a first-order Taylor approximation to estimate the covariance matrix of the circuit log-eigenvalues
    eigenvalues_inv_diag = Diagonal(eigenvalues .^ (-1))
    covariance_log = convert(
        SparseMatrixCSC{Float64, Int},
        Symmetric(eigenvalues_inv_diag * covariance * eigenvalues_inv_diag),
    )
    return covariance_log::SparseMatrixCSC{Float64, Int}
end

"""
    sparse_covariance_inv(covariance_log::SparseMatrixCSC{Float64, Int}, mapping_lengths::Vector{Int})

Returns the inverse of the sparse block diagonal circuit log-eigenvalue estimator covariance matrix `covariance_log`, where the block sizes are specified by `mapping_lengths`.
"""
function sparse_covariance_inv(
    covariance_log::SparseMatrixCSC{Float64, Int},
    mapping_lengths::Vector{Int},
)
    # Use the block diagonal structure of the covariance matrix to speed up computation of the covariance inverse
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    M = sum(mapping_lengths)
    covariance_log_inv = spzeros(Float64, Int, M, M)
    for i in eachindex(mapping_lengths)
        mapping_range = mapping_lower[i]:mapping_upper[i]
        mapping_length = mapping_lengths[i]
        # The sparse Cholesky decomposition of A actually computes a decomposition of a permuted matrix PAP'=LL'
        # Hence A=P'LL'P, and our GLS factor is the inverse of P'L, L^(-1)P
        # Note also that A^(-1) = P'L^(-1)'L^(-1)P
        block_chol = cholesky(covariance_log[mapping_range, mapping_range])
        block_L_inv = sparse(inv(LowerTriangular(sparse(block_chol.L))))
        block_perm = sparse(1:mapping_length, block_chol.p, ones(mapping_length))
        covariance_log_inv[mapping_range, mapping_range] =
            block_perm' * block_L_inv' * block_L_inv * block_perm
    end
    return covariance_log_inv::SparseMatrixCSC{Float64, Int}
end

"""
    calc_gls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the gate eigenvalue estimator covariance matrix for the generalised least squares (GLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_gls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the inverse of the covariance matrix of the circuit log-eigenvalues
    mapping_lengths = length.(d.mapping_ensemble)
    gate_eigenvalues = d.c.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    covariance_log_inv = sparse_covariance_inv(covariance_log, mapping_lengths)
    # Calculate the covariance matrix of the gate eigenvalues using a first-order Taylor approximation
    gls_gate_eigenvalues_cov = Symmetric(
        gate_eigenvalues_diag *
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_inv * d.matrix)))) *
        gate_eigenvalues_diag,
    )
    return gls_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    calc_wls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the gate eigenvalue estimator covariance matrix for the weighted least squares (WLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_wls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the WLS estimator
    gate_eigenvalues = d.c.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    covariance_log_diag_inv = Diagonal(covariance_log)^(-1)
    wls_estimator =
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_diag_inv * d.matrix)))) *
        d.matrix' *
        covariance_log_diag_inv
    # Calculate the covariance matrix of the gate eigenvalues using a first-order Taylor approximation
    wls_gate_eigenvalues_cov = Symmetric(
        gate_eigenvalues_diag *
        wls_estimator *
        covariance_log *
        wls_estimator' *
        gate_eigenvalues_diag,
    )
    return wls_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    calc_ols_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the gate eigenvalue estimator covariance matrix for the ordinary least squares (OLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_ols_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the OLS estimator
    gate_eigenvalues = d.c.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    ols_estimator = inv(bunchkaufman(Symmetric(Array(d.matrix' * d.matrix)))) * d.matrix'
    # Calculate the covariance matrix of the gate eigenvalues using a first-order Taylor approximation
    ols_gate_eigenvalues_cov = Symmetric(
        gate_eigenvalues_diag *
        ols_estimator *
        covariance_log *
        ols_estimator' *
        gate_eigenvalues_diag,
    )
    return ols_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    calc_ls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}, ls_type::Symbol)

Returns the gate eigenvalue estimator covariance matrix for the least squares estimator specified by `ls_type` (`:gls`, `:wls`, or `:ols`) corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_ls_covariance(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int},
    ls_type::Symbol,
)
    if ls_type == :gls
        gate_eigenvalues_cov = calc_gls_covariance(d, covariance_log)
    elseif ls_type == :wls
        gate_eigenvalues_cov = calc_wls_covariance(d, covariance_log)
    elseif ls_type == :ols
        gate_eigenvalues_cov = calc_ols_covariance(d, covariance_log)
    else
        throw(error("The estimator type $(ls_type) must be either :gls, :wls, or :ols."))
    end
    return gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    nrmse_moments(cov_eigenvalues::Vector{Float64})

Returns the expectation and variance of the normalised RMS error, as determined by the eigenvalues `cov_eigenvalues` of the gate eigenvalue estimator covariance matrix.
"""
function nrmse_moments(cov_eigenvalues::Vector{Float64})
    # Calculate the trace of the gate eigenvalue estimator covariance matrix, and its square
    sigma_tr = sum(cov_eigenvalues)
    sigma_sq_tr = sum(cov_eigenvalues .^ 2)
    # Calculate the expectation and variance of the NRMSE
    N = length(cov_eigenvalues)
    expectation = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    variance = (sigma_sq_tr / (2 * sigma_tr)) * (1 - sigma_sq_tr / (8 * sigma_tr^2)) / N
    return (expectation::Float64, variance::Float64)
end

"""
    calc_gls_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the expectation and variance of the normalised RMS error for the generalised least squares (GLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_gls_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the GLS moments from the GLS gate eigenvalue estimator covariance matrix
    (gls_expectation, gls_variance) =
        nrmse_moments(eigvals(calc_gls_covariance(d, covariance_log)))
    return (gls_expectation::Float64, gls_variance::Float64)
end

"""
    calc_wls_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the expectation and variance of the normalised RMS error for the weighted least squares (WLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_wls_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the WLS moments from the WLS gate eigenvalue estimator covariance matrix
    (wls_expectation, wls_variance) =
        nrmse_moments(eigvals(calc_wls_covariance(d, covariance_log)))
    return (wls_expectation::Float64, wls_variance::Float64)
end

"""
    calc_ols_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the expectation and variance of the normalised RMS error for the ordinary least squares (OLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_ols_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the OLS moments from the OLS gate eigenvalue estimator covariance matrix
    (ols_expectation, ols_variance) =
        nrmse_moments(eigvals(calc_ols_covariance(d, covariance_log)))
    return (ols_expectation::Float64, ols_variance::Float64)
end

"""
    calc_ls_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}, ls_type::Symbol)

Returns the expectation and variance of the normalised RMS error for the least squares estimator specified by `ls_type` (`:gls`, `:wls`, or `:ols`) corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_ls_moments(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int},
    ls_type::Symbol,
)
    if ls_type == :gls
        (expectation, variance) = calc_gls_moments(d, covariance_log)
    elseif ls_type == :wls
        (expectation, variance) = calc_wls_moments(d, covariance_log)
    elseif ls_type == :ols
        (expectation, variance) = calc_ols_moments(d, covariance_log)
    else
        throw(error("The estimator type $(ls_type) must be either :gls, :wls, or :ols."))
    end
    return (expectation::Float64, variance::Float64)
end

"""
    calc_gls_merit(d::Design)

Returns the [`Merit`](@ref) object for the generalised least squares (GLS) estimator corresponding to the design `d`.
"""
function calc_gls_merit(d::Design)
    # Calculate the circuit log-eigenvalue covariance matrix and the gate eigenvalues
    covariance_log = calc_covariance_log(d)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_singular_vals = zeros(Float64, d.c.N)
    try
        design_matrix_singular_vals = svd(convert(Matrix{Float64}, Array(d.matrix))).S
    catch
        @warn "The default singular value decomposition algorithm failed; falling back to QR iteration."
        design_matrix_singular_vals =
            svd(
                convert(Matrix{Float64}, Array(d.matrix));
                alg = LinearAlgebra.QRIteration(),
            ).S
    end
    design_matrix_norm = design_matrix_singular_vals[1]
    design_matrix_pinv_norm = 1 / design_matrix_singular_vals[end]
    design_matrix_cond_num = design_matrix_norm * design_matrix_pinv_norm
    # Calculate the GLS merit from the GLS gate eigenvalue estimator covariance matrix
    gls_cov_eigenvalues = eigvals(calc_gls_covariance(d, covariance_log))
    (gls_expectation, gls_variance) = nrmse_moments(gls_cov_eigenvalues)
    gls_merit = Merit(
        d.c.circuit_param,
        d.c.noise_param,
        :gls,
        d.tuple_set,
        d.tuple_set_data,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
        length(d.c.total_gates),
        d.c.N,
        gls_expectation,
        gls_variance,
        gls_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
    )
    return gls_merit::Merit
end

"""
    calc_wls_merit(d::Design)

Returns the [`Merit`](@ref) object for the weighted least squares (WLS) estimator corresponding to the design `d`.
"""
function calc_wls_merit(d::Design)
    # Calculate the circuit log-eigenvalue covariance matrix and the gate eigenvalues
    covariance_log = calc_covariance_log(d)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_singular_vals = zeros(Float64, d.c.N)
    try
        design_matrix_singular_vals = svd(convert(Matrix{Float64}, Array(d.matrix))).S
    catch
        @warn "The default singular value decomposition algorithm failed; falling back to QR iteration."
        design_matrix_singular_vals =
            svd(
                convert(Matrix{Float64}, Array(d.matrix));
                alg = LinearAlgebra.QRIteration(),
            ).S
    end
    design_matrix_norm = design_matrix_singular_vals[1]
    design_matrix_pinv_norm = 1 / design_matrix_singular_vals[end]
    design_matrix_cond_num = design_matrix_norm * design_matrix_pinv_norm
    # Calculate the WLS merit from the WLS gate eigenvalue estimator covariance matrix
    wls_cov_eigenvalues = eigvals(calc_wls_covariance(d, covariance_log))
    (wls_expectation, wls_variance) = nrmse_moments(wls_cov_eigenvalues)
    wls_merit = Merit(
        d.c.circuit_param,
        d.c.noise_param,
        :wls,
        d.tuple_set,
        d.tuple_set_data,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
        length(d.c.total_gates),
        d.c.N,
        wls_expectation,
        wls_variance,
        wls_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
    )
    return wls_merit::Merit
end

"""
    calc_ols_merit(d::Design)

Returns the [`Merit`](@ref) object for the ordinary least squares (OLS) estimator corresponding to the design `d`.
"""
function calc_ols_merit(d::Design)
    # Calculate the circuit log-eigenvalue covariance matrix and the gate eigenvalues
    covariance_log = calc_covariance_log(d)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_singular_vals = zeros(Float64, d.c.N)
    try
        design_matrix_singular_vals = svd(convert(Matrix{Float64}, Array(d.matrix))).S
    catch
        @warn "The default singular value decomposition algorithm failed; falling back to QR iteration."
        design_matrix_singular_vals =
            svd(
                convert(Matrix{Float64}, Array(d.matrix));
                alg = LinearAlgebra.QRIteration(),
            ).S
    end
    design_matrix_norm = design_matrix_singular_vals[1]
    design_matrix_pinv_norm = 1 / design_matrix_singular_vals[end]
    design_matrix_cond_num = design_matrix_norm * design_matrix_pinv_norm
    # Calculate the OLS merit from the OLS gate eigenvalue estimator covariance matrix
    ols_cov_eigenvalues = eigvals(calc_ols_covariance(d, covariance_log))
    (ols_expectation, ols_variance) = nrmse_moments(ols_cov_eigenvalues)
    ols_merit = Merit(
        d.c.circuit_param,
        d.c.noise_param,
        :ols,
        d.tuple_set,
        d.tuple_set_data,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
        length(d.c.total_gates),
        d.c.N,
        ols_expectation,
        ols_variance,
        ols_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
    )
    return ols_merit::Merit
end

"""
    calc_ls_merit(d::Design, ls_type::Symbol)

Returns the [`Merit`](@ref) object for the least squares estimator specified by `ls_type` (`:gls`, `:wls`, or `:ols`) corresponding to the design `d`.
"""
function calc_ls_merit(d::Design, ls_type::Symbol)
    if ls_type == :gls
        merit = calc_gls_merit(d)
    elseif ls_type == :wls
        merit = calc_wls_merit(d)
    elseif ls_type == :ols
        merit = calc_ols_merit(d)
    else
        throw(error("The estimator type $(ls_type) must be either :gls, :wls, or :ols."))
    end
    return merit::Merit
end

"""
    calc_ls_merit(d::Design, ls_type::Symbol)

Returns [`Merit`](@ref) objects for all three least squares estimators (generalised, weighted, and ordinary least squares) corresponding to the design `d`.
"""
function calc_merit_set(d::Design)
    # Calculate the circuit log-eigenvalue covariance matrix and the gate eigenvalues
    covariance_log = calc_covariance_log(d)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_singular_vals = zeros(Float64, d.c.N)
    try
        design_matrix_singular_vals = svd(convert(Matrix{Float64}, Array(d.matrix))).S
    catch
        @warn "The default singular value decomposition algorithm failed; falling back to QR iteration."
        design_matrix_singular_vals =
            svd(
                convert(Matrix{Float64}, Array(d.matrix));
                alg = LinearAlgebra.QRIteration(),
            ).S
    end
    design_matrix_norm = design_matrix_singular_vals[1]
    design_matrix_pinv_norm = 1 / design_matrix_singular_vals[end]
    design_matrix_cond_num = design_matrix_norm * design_matrix_pinv_norm
    # Calculate the GLS merit from the GLS gate eigenvalue estimator covariance matrix
    gls_cov_eigenvalues = eigvals(calc_gls_covariance(d, covariance_log))
    (gls_expectation, gls_variance) = nrmse_moments(gls_cov_eigenvalues)
    gls_merit = Merit(
        d.c.circuit_param,
        d.c.noise_param,
        :gls,
        d.tuple_set,
        d.tuple_set_data,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
        length(d.c.total_gates),
        d.c.N,
        gls_expectation,
        gls_variance,
        gls_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
    )
    # Calculate the WLS merit from the WLS gate eigenvalue estimator covariance matrix
    wls_cov_eigenvalues = eigvals(calc_wls_covariance(d, covariance_log))
    (wls_expectation, wls_variance) = nrmse_moments(wls_cov_eigenvalues)
    wls_merit = Merit(
        d.c.circuit_param,
        d.c.noise_param,
        :wls,
        d.tuple_set,
        d.tuple_set_data,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
        length(d.c.total_gates),
        d.c.N,
        wls_expectation,
        wls_variance,
        wls_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
    )
    # Calculate the OLS merit from the OLS gate eigenvalue estimator covariance matrix
    ols_cov_eigenvalues = eigvals(calc_ols_covariance(d, covariance_log))
    (ols_expectation, ols_variance) = nrmse_moments(ols_cov_eigenvalues)
    ols_merit = Merit(
        d.c.circuit_param,
        d.c.noise_param,
        :ols,
        d.tuple_set,
        d.tuple_set_data,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
        length(d.c.total_gates),
        d.c.N,
        ols_expectation,
        ols_variance,
        ols_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
    )
    return (gls_merit::Merit, wls_merit::Merit, ols_merit::Merit)
end

"""
    nrmse_pdf_integrand(u::Float64, x::Float64, norm_cov_eigenvalues::Vector{Float64})

Returns the integrand of the Imhof method CDF for the distribution of the NRMSE.

Calculation follows Eq. 3.2 of `Computing the distribution of quadratic forms in normal variables` by J. P. Imhof (1961).
"""
function nrmse_pdf_integrand(u::Float64, x::Float64, norm_cov_eigenvalues::Vector{Float64})
    # Both theta and rho have been simplified as, with reference to Eq. 1.1 we see that h and delta are 1 and 0, respectively
    theta = sum(atan(eigenvalue * u) for eigenvalue in norm_cov_eigenvalues) / 2 - x * u / 2
    rho = prod((1 + (eigenvalue * u)^2) for eigenvalue in norm_cov_eigenvalues)^(1 / 4)
    integrand = sin(theta) / (u * rho)
    return integrand::Float64
end

"""
    nrmse_pdf(cov_eigenvalues::Vector{Float64}, nrmse_values::Vector{Float64}; epsilon::Float64 = 1e-5)

Returns the probability density function (PDF) for the normalised RMS error (NRMSE) of the gate eigenvalue estimator vector, which follows a generalised chi-squared distribution and whose covariance matrix has eigenvalues `cov_eigenvalues`, at the coordinates specified by `nrmse_values`.
Does not calculate values when the normal approximation to the PDF is less than a factor of `epsilon` of its maximum value.

Calculation follows Eq. 3.2 of `Computing the distribution of quadratic forms in normal variables` by J. P. Imhof (1961).
"""
function nrmse_pdf(
    cov_eigenvalues::Vector{Float64},
    nrmse_values::Vector{Float64};
    epsilon::Float64 = 1e-5,
)
    # Normalise the gate eigenvalue estimator covariance matrix eigenvalues
    N = length(cov_eigenvalues)
    norm_cov_eigenvalues = cov_eigenvalues / N
    # Calculate the moments 
    (expectation, variance) = nrmse_moments(cov_eigenvalues)
    # Calculate the normal approximation of the PDF
    nrmse_pdf_normal =
        1 / (sqrt(2 * π) * sqrt(variance)) *
        exp.(-(nrmse_values .- expectation) .^ 2 / (2 * variance))
    nrmse_pdf_normal_max = maximum(nrmse_pdf_normal)
    # Perform Imhof's method to calculate the PDF
    function imhof_cdf(x)
        return 0.5 -
               quadgk(u -> nrmse_pdf_integrand(u, x, norm_cov_eigenvalues), 0, Inf)[1] / π
    end
    nrmse_pdf = zeros(length(nrmse_values))
    for (idx, x) in enumerate(nrmse_values)
        # Only perform the calculation for sufficiently large PDF values
        if nrmse_pdf_normal[idx] / nrmse_pdf_normal_max > epsilon
            # Differentiate the CDF to obtain the PDF
            # Note that Imhof's method describes the square, whereas we seek the square root
            nrmse_pdf[idx] = 2 * x * central_fdm(5, 1)(imhof_cdf, x^2)
        end
    end
    return nrmse_pdf::Vector{Float64}
end

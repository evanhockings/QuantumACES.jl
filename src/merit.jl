"""
    Merit

Merit parameters for an experimental design.

# Fields

  - `circuit_param::AbstractCircuitParameters`: Circuit parameters.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.
  - `tuple_set_data::TupleSetData`: [`TupleSetData`](@ref) object that generates the tuple set.
  - `tuple_times::Vector{Float64}`: Time taken to implement the circuit corresponding to each tuple, normalised according to the basic tuple set.
  - `shot_weights::Vector{Float64}`: Shot weights for each tuple in the set, which add to 1.
  - `experiment_numbers::Vector{Int}`: Number of experiments for each tuple in the set.
  - `experiment_number::Int`: Total number of experiments.
  - `N::Int`: Number of gate eigenvalues.
  - `N_marginal::Int`: Number of marginal gate eigenvalues.
  - `N_relative::Int`: Number of relative precision marginal gate eigenvalues.
  - `G::Int`: Total number of gates.
  - `gls_expectation::Float64`: Expected normalised RMS error for the generalised least squares (GLS) gate eigenvalue estimator vector.
  - `gls_variance::Float64`: Normalised RMS error variance for the GLS gate eigenvalue estimator vector.
  - `gls_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the GLS gate eigenvalue estimator covariance matrix.
  - `gls_marginal_expectation::Float64`: Expected normalised RMS error for the marginal GLS gate eigenvalue estimator vector.
  - `gls_marginal_variance::Float64`: Normalised RMS error variance for the marginal GLS gate eigenvalue estimator vector.
  - `gls_marginal_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the marginal GLS gate eigenvalue estimator covariance matrix.
  - `gls_relative_expectation::Float64`: Expected normalised RMS error for the relative precision marginal GLS gate eigenvalue estimator vector.
  - `gls_relative_variance::Float64`: Normalised RMS error variance for the relative precision marginal GLS gate eigenvalue estimator vector.
  - `gls_relative_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the relative precision marginal GLS gate eigenvalue estimator covariance matrix.
  - `wls_expectation::Float64`: Expected normalised RMS error for the weighted least squares (WLS) gate eigenvalue estimator vector.
  - `wls_variance::Float64`: Normalised RMS error variance for the WLS gate eigenvalue estimator vector.
  - `wls_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the WLS gate eigenvalue estimator covariance matrix.
  - `wls_marginal_expectation::Float64`: Expected normalised RMS error for the marginal WLS gate eigenvalue estimator vector.
  - `wls_marginal_variance::Float64`: Normalised RMS error variance for the marginal WLS gate eigenvalue estimator vector.
  - `wls_marginal_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the marginal WLS gate eigenvalue estimator covariance matrix.
  - `wls_relative_expectation::Float64`: Expected normalised RMS error for the relative precision marginal WLS gate eigenvalue estimator vector.
  - `wls_relative_variance::Float64`: Normalised RMS error variance for the relative precision marginal WLS gate eigenvalue estimator vector.
  - `wls_relative_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the relative precision marginal WLS gate eigenvalue estimator covariance matrix.
  - `ols_expectation::Float64`: Expected normalised RMS error for the ordinary least squares (OLS) gate eigenvalue estimator vector.
  - `ols_variance::Float64`: Normalised RMS error variance for the OLS gate eigenvalue estimator vector.
  - `ols_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the OLS gate eigenvalue estimator covariance matrix.
  - `ols_marginal_expectation::Float64`: Expected normalised RMS error for the marginal OLS gate eigenvalue estimator vector.
  - `ols_marginal_variance::Float64`: Normalised RMS error variance for the marginal OLS gate eigenvalue estimator vector.
  - `ols_marginal_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the marginal OLS gate eigenvalue estimator covariance matrix.
  - `ols_relative_expectation::Float64`: Expected normalised RMS error for the relative precision marginal OLS gate eigenvalue estimator vector.
  - `ols_relative_variance::Float64`: Normalised RMS error variance for the relative precision marginal OLS gate eigenvalue estimator vector.
  - `ols_relative_cov_eigenvalues::Vector{Float64}`: Eigenvalues of the relative precision marginal OLS gate eigenvalue estimator covariance matrix.
  - `cond_num::Float64`: Condition number of the design matrix, the ratio of the largest and smallest singular values.
  - `pinv_norm::Float64`: Pseudoinverse norm of the design matrix, the inverse of the smallest singular value.
"""
struct Merit
    circuit_param::AbstractCircuitParameters
    noise_param::AbstractNoiseParameters
    tuple_set_data::TupleSetData
    tuple_times::Vector{Float64}
    shot_weights::Vector{Float64}
    experiment_numbers::Vector{Int}
    experiment_number::Int
    N::Int
    N_marginal::Int
    N_relative::Int
    G::Int
    gls_expectation::Float64
    gls_variance::Float64
    gls_cov_eigenvalues::Vector{Float64}
    gls_marginal_expectation::Float64
    gls_marginal_variance::Float64
    gls_marginal_cov_eigenvalues::Vector{Float64}
    gls_relative_expectation::Float64
    gls_relative_variance::Float64
    gls_relative_cov_eigenvalues::Vector{Float64}
    wls_expectation::Float64
    wls_variance::Float64
    wls_cov_eigenvalues::Vector{Float64}
    wls_marginal_expectation::Float64
    wls_marginal_variance::Float64
    wls_marginal_cov_eigenvalues::Vector{Float64}
    wls_relative_expectation::Float64
    wls_relative_variance::Float64
    wls_relative_cov_eigenvalues::Vector{Float64}
    ols_expectation::Float64
    ols_variance::Float64
    ols_cov_eigenvalues::Vector{Float64}
    ols_marginal_expectation::Float64
    ols_marginal_variance::Float64
    ols_marginal_cov_eigenvalues::Vector{Float64}
    ols_relative_expectation::Float64
    ols_relative_variance::Float64
    ols_relative_cov_eigenvalues::Vector{Float64}
    cond_num::Float64
    pinv_norm::Float64
end

function Base.show(io::IO, merit::Merit)
    # Initialise variables
    digits = 4
    # GLS merit
    gls_merit = "$(round(merit.gls_expectation, digits = digits)) ± $(round(sqrt(merit.gls_variance), digits = digits))"
    gls_marginal_merit = "$(round(merit.gls_marginal_expectation, digits = digits)) ± $(round(sqrt(merit.gls_marginal_variance), digits = digits))"
    gls_relative_merit = "$(round(merit.gls_relative_expectation, digits = digits)) ± $(round(sqrt(merit.gls_relative_variance), digits = digits))"
    # WLS merit
    wls_merit = "$(round(merit.wls_expectation, digits = digits)) ± $(round(sqrt(merit.wls_variance), digits = digits))"
    wls_marginal_merit = "$(round(merit.wls_marginal_expectation, digits = digits)) ± $(round(sqrt(merit.wls_marginal_variance), digits = digits))"
    wls_relative_merit = "$(round(merit.wls_relative_expectation, digits = digits)) ± $(round(sqrt(merit.wls_relative_variance), digits = digits))"
    # OLS merit
    ols_merit = "$(round(merit.ols_expectation, digits = digits)) ± $(round(sqrt(merit.ols_variance), digits = digits))"
    ols_marginal_merit = "$(round(merit.ols_marginal_expectation, digits = digits)) ± $(round(sqrt(merit.ols_marginal_variance), digits = digits))"
    ols_relative_merit = "$(round(merit.ols_relative_expectation, digits = digits)) ± $(round(sqrt(merit.ols_relative_variance), digits = digits))"
    # Geometric mean
    gls_geo_mean = string(
        round(
            (sqrt(merit.gls_expectation * merit.gls_relative_expectation));
            digits = digits,
        ),
    )
    wls_geo_mean = string(
        round(
            (sqrt(merit.wls_expectation * merit.wls_relative_expectation));
            digits = digits,
        ),
    )
    ols_geo_mean = string(
        round(
            (sqrt(merit.ols_expectation * merit.ols_relative_expectation));
            digits = digits,
        ),
    )
    # Display the merit
    merit_array = [
        ["GLS" gls_merit gls_marginal_merit gls_relative_merit gls_geo_mean]
        ["WLS" wls_merit wls_marginal_merit wls_relative_merit wls_geo_mean]
        ["OLS" ols_merit ols_marginal_merit ols_relative_merit ols_geo_mean]
    ]
    header = ["Type"; "Ordinary merit"; "Marginal merit"; "Relative merit"; "Geo. mean"]
    pretty_table(io, merit_array; header = header, alignment = :l)
    return nothing
end

@struct_hash_equal_isequal Merit

"""
    get_covariance_mapping_data(d::Design)

Returns covariance matrix mapping index data for the design `d`.
"""
function get_covariance_mapping_data(d::Design)
    # Initialise parameters
    tuple_number = length(d.tuple_set)
    mapping_lengths = length.(d.mapping_ensemble)
    # Covariance circuit eigenvalue measurement supports
    covariance_keys_ensemble = [
        setdiff(
            sort(collect(keys(d.covariance_dict_ensemble[i]))),
            [CartesianIndex(pauli_idx, pauli_idx) for pauli_idx in 1:mapping_lengths[i]],
        ) for i in 1:tuple_number
    ]
    covariance_key_indices_ensemble = [
        Dict(covariance_key => idx for (idx, covariance_key) in pairs(covariance_keys))
        for covariance_keys in covariance_keys_ensemble
    ]
    return (
        covariance_keys_ensemble::Vector{Vector{CartesianIndex{2}}},
        covariance_key_indices_ensemble::Vector{Dict{CartesianIndex{2}, Int}},
    )
end

"""
    calc_covariance_eigenvalues(d::Design)
    calc_covariance_eigenvalues(d::Design, gate_eigenvalues::Vector{Float64})

Returns the covariance matrix eigenvalue ensemble for the design `d` with gate eigenvalues `gate_eigenvalues`.
"""
function calc_covariance_eigenvalues(d::Design, gate_eigenvalues::Vector{Float64})
    # Initialise variables
    tuple_number = length(d.tuple_set)
    (covariance_keys_ensemble, covariance_key_indices_ensemble) =
        get_covariance_mapping_data(d)
    gate_log_eigenvalues = log.(gate_eigenvalues)
    #
    covariance_eigenvalues_ensemble = Vector{Vector{Float64}}(undef, tuple_number)
    for i in 1:tuple_number
        covariance_dict = d.covariance_dict_ensemble[i]
        covariance_keys = covariance_keys_ensemble[i]
        covariance_key_indices = covariance_key_indices_ensemble[i]
        covariance_eigenvalues_ensemble[i] = Vector{Float64}(undef, length(covariance_keys))
        for covariance_key in covariance_keys
            covariance_eigenvalues_ensemble[i][covariance_key_indices[covariance_key]] =
                exp(covariance_dict[covariance_key][1].design_row' * gate_log_eigenvalues)
        end
    end
    return covariance_eigenvalues_ensemble::Vector{Vector{Float64}}
end
function calc_covariance_eigenvalues(d::Design)
    covariance_eigenvalues_ensemble =
        calc_covariance_eigenvalues(d, get_gate_eigenvalues(d))
    return covariance_eigenvalues_ensemble::Vector{Vector{Float64}}
end

"""
    calc_covariance(d::Design; epsilon::Real = 1e-10, weight_time::Bool = true, warning::Bool = true)
    calc_covariance(d::Design, gate_eigenvalues::Vector{Float64}; epsilon::Real = 1e-10, weight_time::Bool = true, warning::Bool = true)
    calc_covariance(d::Design, eigenvalues_ensemble::Vector{Vector{Float64}}, covariance_ensemble::Vector{Vector{Float64}}; epsilon::Real = 1e-10, weight_time::Bool = true, warning::Bool = true)
    calc_covariance(d::Design, eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}}, covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}}; epsilon::Real = 1e-10, weight_time::Bool = true, warning::Bool = true)

Returns the circuit eigenvalue estimator covariance matrix for the design `d` with gate eigenvalues `gate_eigenvalues`, eigenvalue ensemble `eigenvalues_ensemble`, and covariance eigenvalue ensemble `covariance_ensemble`.
Eigenvalue variances are set to a minimum value `epsilon` and the covariance matrix is adjusted by the times factor if `weight_time` is `true`, and if `warning` is `true`, warns that if `d.full_covariance` is `false` this will only generate the diagonal of the covariance matrix.
"""
function calc_covariance(
    d::Design,
    eigenvalues_ensemble::Vector{Vector{Float64}},
    covariance_ensemble::Vector{Vector{Float64}};
    epsilon::Real = 1e-10,
    weight_time::Bool = true,
    warning::Bool = true,
)
    # Warn if the design does not have the full covariance matrix
    if warning && ~d.full_covariance
        @warn "The design does not include the full covariance matrix; only the diagonal will be calculated."
    end
    # Initialise variables
    tuple_number = length(d.tuple_set)
    mapping_lengths = length.(d.mapping_ensemble)
    trivial_pauli = Pauli(Bool[0], 0)
    trivial_row = SparseVector{Int32, Int32}(1, Int32[], Int32[])
    trivial_track = Vector{Int16}[]
    trivial_mapping = Mapping(trivial_pauli, trivial_pauli, trivial_row, trivial_track)
    # Generate the covariance matrix for the circuit eigenvalues
    covariance_key_indices_ensemble = get_covariance_mapping_data(d)[2]
    covariance_blocks = Vector{SparseMatrixCSC{Float64, Int}}(undef, tuple_number)
    @threads :static for i in 1:tuple_number
        experiment_number = d.experiment_numbers[i]
        shot_weight = d.shot_weights[i]
        covariance_dict = d.covariance_dict_ensemble[i]
        mapping_length = mapping_lengths[i]
        eigenvalues = eigenvalues_ensemble[i]
        covariance_eigenvalues = covariance_ensemble[i]
        covariance_key_indices = covariance_key_indices_ensemble[i]
        covariance_block = spzeros(Float64, mapping_length, mapping_length)
        # Generate the relevant term in the covariance matrix for each key in the dictionary
        for covariance_key in keys(covariance_dict)
            pauli_idx_1 = covariance_key[1]
            pauli_idx_2 = covariance_key[2]
            # If the key lies on the diagonal, calculate the variance, otherwise calculate the covariance, which is a little more complicated
            if pauli_idx_1 == pauli_idx_2
                @assert covariance_dict[covariance_key][1] == trivial_mapping
                # Calculate the eigenvalue
                eigenvalue = eigenvalues[pauli_idx_1]
                # Calculate the shot fraction for the variance
                eigenvalue_experiments = covariance_dict[covariance_key][2]
                variance_shot_factor =
                    (1 / shot_weight) * (experiment_number / eigenvalue_experiments)
                # Calculate the variance
                eigenvalue_variance =
                    variance_shot_factor * max((1 - eigenvalue^2), epsilon)
                # Set the covariance matrix term
                covariance_block[covariance_key] = eigenvalue_variance
            else
                # Generate covariance keys
                covariance_key_1 = CartesianIndex(pauli_idx_1, pauli_idx_1)
                covariance_key_2 = CartesianIndex(pauli_idx_2, pauli_idx_2)
                covariance_key_transpose = CartesianIndex(pauli_idx_2, pauli_idx_1)
                # Calculate the eigenvalues
                eigenvalue_1 = eigenvalues[pauli_idx_1]
                eigenvalue_2 = eigenvalues[pauli_idx_2]
                eigenvalue_12 =
                    covariance_eigenvalues[covariance_key_indices[covariance_key]]
                # Calculate the shot fraction for the variance
                eigenvalue_1_experiments = covariance_dict[covariance_key_1][2]
                eigenvalue_2_experiments = covariance_dict[covariance_key_2][2]
                eigenvalue_12_experiments = covariance_dict[covariance_key][2]
                covariance_shot_factor =
                    (1 / shot_weight) * (
                        (experiment_number * eigenvalue_12_experiments) /
                        (eigenvalue_1_experiments * eigenvalue_2_experiments)
                    )
                # Calculate the covariance
                eigenvalue_both_covariance =
                    covariance_shot_factor * (eigenvalue_12 - eigenvalue_1 * eigenvalue_2)
                # Set the covariance matrix entries
                covariance_block[covariance_key] = eigenvalue_both_covariance
                covariance_block[covariance_key_transpose] = eigenvalue_both_covariance
            end
        end
        covariance_blocks[i] = covariance_block
    end
    covariance = blockdiag(covariance_blocks...)
    # Adjust the covariance matrix by the times factor
    if weight_time
        times_factor = sum(d.shot_weights .* d.tuple_times)
        covariance = times_factor * covariance
    end
    return covariance::SparseMatrixCSC{Float64, Int}
end
function calc_covariance(
    d::Design,
    eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}},
    covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}};
    epsilon::Real = 1e-10,
    weight_time::Bool = true,
    warning::Bool = true,
)
    tuple_number = length(d.tuple_set)
    eigenvalues_ensemble =
        [mean.(eigenvalues_experiment_ensemble[i]) for i in 1:tuple_number]
    covariance_ensemble = [mean.(covariance_experiment_ensemble[i]) for i in 1:tuple_number]
    covariance = calc_covariance(
        d,
        eigenvalues_ensemble,
        covariance_ensemble;
        epsilon = epsilon,
        weight_time = weight_time,
        warning = warning,
    )
    return covariance::SparseMatrixCSC{Float64, Int}
end
function calc_covariance(
    d::Design,
    gate_eigenvalues::Vector{Float64};
    epsilon::Real = 1e-10,
    weight_time::Bool = true,
    warning::Bool = true,
)
    eigenvalues_ensemble = get_eigenvalues_ensemble(d, gate_eigenvalues)
    covariance_ensemble = calc_covariance_eigenvalues(d, gate_eigenvalues)
    covariance = calc_covariance(
        d,
        eigenvalues_ensemble,
        covariance_ensemble;
        epsilon = epsilon,
        weight_time = weight_time,
        warning = warning,
    )
    return covariance::SparseMatrixCSC{Float64, Int}
end
function calc_covariance(
    d::Design;
    epsilon::Real = 1e-10,
    weight_time::Bool = true,
    warning::Bool = true,
)
    gate_eigenvalues = get_gate_eigenvalues(d)
    covariance = calc_covariance(
        d,
        gate_eigenvalues;
        epsilon = epsilon,
        weight_time = weight_time,
        warning = warning,
    )
    return covariance::SparseMatrixCSC{Float64, Int}
end

"""
    calc_covariance_log(covariance::SparseMatrixCSC{Float64, Int}, eigenvalues::Vector{Float64})
    calc_covariance_log(d::Design; warning::Bool = true)

Returns the covariance matrix of the circuit log-eigenvalues corresponding to the circuit eigenvalues `eigenvalues` with covariance matrix `covariance`, using values generated from the design `d` and, if `warning` is `true`, warns that if `d.full_covariance` is `false` then this will only generate the diagonal of the covariance matrix.
"""
function calc_covariance_log(
    covariance::SparseMatrixCSC{Float64, Int},
    eigenvalues::Vector{Float64},
)
    # Use a first-order Taylor approximation to estimate the covariance matrix of the circuit log-eigenvalues
    eigenvalues_inv_diag = sparse(Diagonal(eigenvalues .^ (-1)))
    covariance_log =
        sparse(Symmetric(eigenvalues_inv_diag * covariance * eigenvalues_inv_diag))
    return covariance_log::SparseMatrixCSC{Float64, Int}
end
function calc_covariance_log(d::Design; warning::Bool = true)
    eigenvalues = get_eigenvalues(d)
    covariance = calc_covariance(d; warning = warning)
    covariance_log = calc_covariance_log(covariance, eigenvalues)
    return covariance_log::SparseMatrixCSC{Float64, Int}
end

"""
    sparse_covariance_inv_factor(covariance_log::SparseMatrixCSC{Float64, Int}, mapping_lengths::Vector{Int}; epsilon::Real = 1e-12, warning::Bool = true)

Returns the inverse of the lower Cholesky factor of the sparse block diagonal circuit log-eigenvalue estimator covariance matrix `covariance_log`, where the block sizes are specified by `mapping_lengths`, ensuring a smallest eigenvalue `epsilon` if the Cholesky factorisation fails, and warning if `warning` is `true`.
"""
function sparse_covariance_inv_factor(
    covariance_log::SparseMatrixCSC{Float64, Int},
    mapping_lengths::Vector{Int};
    epsilon::Real = 1e-12,
    warning::Bool = true,
)
    # Use the block diagonal structure of the covariance matrix to speed up computation of the covariance inverse Cholesky factor
    M = sum(mapping_lengths)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    covariance_log_inv_factor = spzeros(Float64, Int, M, M)
    for i in eachindex(mapping_lengths)
        mapping_range = mapping_lower[i]:mapping_upper[i]
        mapping_length = mapping_lengths[i]
        block_covariance_log = covariance_log[mapping_range, mapping_range]
        # The sparse Cholesky decomposition of A actually computes a decomposition of a permuted matrix PAP'=LL'
        # Hence A=P'LL'P, and our GLS factor is the inverse of P'L, L^(-1)P
        # Note also that A^(-1) = P'L^(-1)'L^(-1)P
        # Use the sparse Cholesky factorisation to compute the inverse
        try
            block_chol = cholesky(block_covariance_log)
            block_L_inv = sparse(inv(LowerTriangular(sparse(block_chol.L))))
            block_perm = sparse(1:mapping_length, block_chol.p, ones(mapping_length))
            covariance_log_inv_factor[mapping_range, mapping_range] =
                block_L_inv * block_perm
        catch
            block_covariance_log = sparse(Symmetric(abs.(block_covariance_log)))
            block_eig_min = Arpack.eigs(
                block_covariance_log;
                nev = 1,
                ncv = mapping_length,
                which = :SR,
            )[1][1]
            if block_eig_min < epsilon
                block_covariance_log += (epsilon - block_eig_min) * I
                if warning
                    println(
                        "WARNING: After taking the elementwise absolute value, the smallest eigenvalue $(round(block_eig_min, sigdigits = 4)) of the covariance matrix block $(i) remains small; setting to $(epsilon).",
                    )
                end
            end
            block_chol = cholesky(block_covariance_log)
            block_L_inv = sparse(inv(LowerTriangular(sparse(block_chol.L))))
            block_perm = sparse(1:mapping_length, block_chol.p, ones(mapping_length))
            covariance_log_inv_factor[mapping_range, mapping_range] =
                block_L_inv * block_perm
        end
    end
    return covariance_log_inv_factor::SparseMatrixCSC{Float64, Int}
end

"""
    sparse_covariance_inv(covariance_log::SparseMatrixCSC{Float64, Int}, mapping_lengths::Vector{Int}; epsilon::Real = 1e-12, warning::Bool = true)

Returns the inverse of the sparse block diagonal circuit log-eigenvalue estimator covariance matrix `covariance_log`, where the block sizes are specified by `mapping_lengths`, ensuring a smallest eigenvalue `epsilon` if the Cholesky factorisation fails, and warning if `warning` is `true`.
"""
function sparse_covariance_inv(
    covariance_log::SparseMatrixCSC{Float64, Int},
    mapping_lengths::Vector{Int};
    epsilon::Real = 1e-12,
    warning::Bool = true,
)
    covariance_log_inv_factor = sparse_covariance_inv_factor(
        covariance_log,
        mapping_lengths;
        epsilon = epsilon,
        warning = warning,
    )
    covariance_log_inv = covariance_log_inv_factor' * covariance_log_inv_factor
    return covariance_log_inv::SparseMatrixCSC{Float64, Int}
end

"""
    get_transform(d::Design, est_type::Symbol)

Returns a transform matrix that maps gate eigenvalues to gate eigenvalues of type depending on the estimator type `est_type` (`:ordinary`, `:marginal`, or `:relative`), calculated using the gate data of the design `d`.
"""
function get_transform(d::Design, est_type::Symbol)
    if est_type == :ordinary
        transform_matrix = get_ordinary_transform(d.c.gate_data)
    elseif est_type == :marginal
        transform_matrix = get_marginal_transform(d.c.gate_data)
    elseif est_type == :relative
        transform_matrix = get_relative_transform(d.c.gate_data)
    else
        throw(
            error(
                "The estimator type $(est_type) must be either :ordinary, :marginal, or :relative.",
            ),
        )
    end
    return transform_matrix::SparseMatrixCSC{Float64, Int}
end

"""
    get_marginal_gate_covariance(d::Design, gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}})

Returns the marginal gate eigenvalue estimator covariance matrix corresponding to the gate eigenvalue estimator covariance matrix `gate_eigenvalues_cov` for the design `d`.
"""
function get_marginal_gate_covariance(
    d::Design,
    gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}},
)
    # Get the transform matrix
    marginal_transform = get_transform(d, :marginal)
    # Transform the gate eigenvalue estimator covariance matrix
    marginal_gate_eigenvalues_cov =
        Symmetric(marginal_transform * gate_eigenvalues_cov * marginal_transform')
    return marginal_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    get_relative_gate_covariance(d::Design, gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}})

Returns the relative gate eigenvalue estimator covariance matrix corresponding to the gate eigenvalue estimator covariance matrix `gate_eigenvalues_cov` for the design `d`.
"""
function get_relative_gate_covariance(
    d::Design,
    gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}},
)
    # Get the transform matrix
    relative_transform = get_transform(d, :relative)
    # Transform the gate eigenvalue estimator covariance matrix
    relative_gate_eigenvalues_cov =
        Symmetric(relative_transform * gate_eigenvalues_cov * relative_transform')
    return relative_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    get_pad_transform(d::Design, est_type::Symbol; inverse::Bool = false, probabilities::Bool = false)

Returns a transform matrix that pads gate eigenvalues, or gate error probabilities if `probabilities` is `true`, whose type depends on the estimator type `est_type` (`:ordinary`, `:marginal`, or `:relative`), with identity eigenvaleus or error probabilities, respectively, and the transpose if `inverse` is `true`, calculated using the gate data of the design `d`.
"""
function get_pad_transform(
    d::Design,
    est_type::Symbol;
    inverse::Bool = false,
    probabilities::Bool = false,
)
    if est_type == :ordinary
        pad_transform = get_pad_transform(
            d.c.gate_data;
            inverse = inverse,
            probabilities = probabilities,
        )
    elseif est_type == :marginal
        pad_transform = get_marginal_pad_transform(
            d.c.gate_data;
            inverse = inverse,
            probabilities = probabilities,
        )
    elseif est_type == :relative
        pad_transform = get_relative_pad_transform(
            d.c.gate_data;
            inverse = inverse,
            probabilities = probabilities,
        )
    else
        throw(
            error(
                "The estimator type $(est_type) must be either :ordinary, :marginal, or :relative.",
            ),
        )
    end
end

"""
    get_wht_transform(d::Design, est_type::Symbol; inverse::Bool = false)

Returns a transform matrix that maps padded gate error probabilities to padded gate eigenvalues, or the inverse transform if `inverse` is `true`, whose type depends on the estimator type `est_type` (`:ordinary`, :marginal`, or `:relative`), calculated using the gate data of the design `d`.
"""
function get_wht_transform(d::Design, est_type::Symbol; inverse::Bool = false)
    if est_type == :ordinary
        wht_transform = get_wht_transform(d.c.gate_data; inverse = inverse)
    elseif est_type == :marginal
        wht_transform = get_marginal_wht_transform(d.c.gate_data; inverse = inverse)
    elseif est_type == :relative
        wht_transform = get_relative_wht_transform(d.c.gate_data; inverse = inverse)
    else
        throw(
            error(
                "The estimator type $(est_type) must be either :ordinary, :marginal, or :relative.",
            ),
        )
    end
    return wht_transform::SparseMatrixCSC{Float64, Int}
end

"""
    calc_gls_covariance(d::Design)
    calc_gls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the gate eigenvalue estimator covariance matrix for the generalised least squares (GLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_gls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the inverse of the covariance matrix of the circuit log-eigenvalues
    mapping_lengths = length.(d.mapping_ensemble)
    gate_eigenvalues_diag = Diagonal(get_gate_eigenvalues(d))
    covariance_log_inv = sparse_covariance_inv(covariance_log, mapping_lengths)
    # Calculate the covariance matrix of the gate eigenvalues using a first-order Taylor approximation
    gls_gate_eigenvalues_cov = Symmetric(
        gate_eigenvalues_diag *
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_inv * d.matrix)))) *
        gate_eigenvalues_diag,
    )
    return gls_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end
function calc_gls_covariance(d::Design)
    covariance_log = calc_covariance_log(d)
    gls_gate_eigenvalues_cov = calc_gls_covariance(d, covariance_log)
    return gls_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    calc_wls_covariance(d::Design)
    calc_wls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the gate eigenvalue estimator covariance matrix for the weighted least squares (WLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_wls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the WLS estimator
    gate_eigenvalues_diag = Diagonal(get_gate_eigenvalues(d))
    covariance_log_diag_inv = sparse(Diagonal(covariance_log)^(-1))
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
function calc_wls_covariance(d::Design)
    covariance_log = calc_covariance_log(d)
    wls_gate_eigenvalues_cov = calc_wls_covariance(d, covariance_log)
    return wls_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    calc_ols_covariance(d::Design)
    calc_ols_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})

Returns the gate eigenvalue estimator covariance matrix for the ordinary least squares (OLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`.
"""
function calc_ols_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the OLS estimator
    gate_eigenvalues_diag = Diagonal(get_gate_eigenvalues(d))
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
function calc_ols_covariance(d::Design)
    covariance_log = calc_covariance_log(d)
    ols_gate_eigenvalues_cov = calc_ols_covariance(d, covariance_log)
    return ols_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    calc_ls_covariance(d::Design; ls_type::Symbol = :gls, est_type::Symbol = :ordinary)
    calc_ls_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; ls_type::Symbol = :gls, est_type::Symbol = :ordinary)

Returns the gate eigenvalue estimator covariance matrix for the least squares estimator specified by `ls_type` (`:gls`, `:wls`, or `:ols`) corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`, with estimator type `est_type` (`:ordinary`, `:marginal`, or `:relative`).
"""
function calc_ls_covariance(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    ls_type::Symbol = :gls,
    est_type::Symbol = :ordinary,
)
    transform_matrix = get_transform(d, est_type)
    if ls_type == :gls
        gate_eigenvalues_cov = calc_gls_covariance(d, covariance_log)
    elseif ls_type == :wls
        gate_eigenvalues_cov = calc_wls_covariance(d, covariance_log)
    elseif ls_type == :ols
        gate_eigenvalues_cov = calc_ols_covariance(d, covariance_log)
    else
        throw(
            error("The least squares type $(ls_type) must be either :gls, :wls, or :ols."),
        )
    end
    gate_eigenvalues_cov =
        Symmetric(transform_matrix * gate_eigenvalues_cov * transform_matrix')
    return gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end
function calc_ls_covariance(d::Design; ls_type::Symbol = :gls, est_type::Symbol = :ordinary)
    covariance_log = calc_covariance_log(d)
    gate_eigenvalues_cov =
        calc_ls_covariance(d, covariance_log; ls_type = ls_type, est_type = est_type)
    return gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    nrmse_moments(gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}})    
    nrmse_moments(gate_eigenvalues_cov_eigenvalues::Vector{Float64})

Returns the expectation and variance of the normalised RMS error, as determined either by the gate eigenvalue estimator covariance matrix `gate_eigenvalues_cov` or its eigenvalues `gate_eigenvalues_cov_eigenvalues`.
"""
function nrmse_moments(gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}})
    # Calculate the trace of the gate eigenvalue estimator covariance matrix, and its square
    sigma_tr = tr(gate_eigenvalues_cov)
    sigma_sq_tr = tr(gate_eigenvalues_cov^2)
    # Calculate the expectation and variance of the NRMSE
    N = size(gate_eigenvalues_cov, 1)
    @assert N == size(gate_eigenvalues_cov, 2) "The gate eigenvalue estimator covariance matrix is not square."
    expectation = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    variance = (sigma_sq_tr / (2 * sigma_tr)) * (1 - sigma_sq_tr / (8 * sigma_tr^2)) / N
    return (expectation::Float64, variance::Float64)
end
function nrmse_moments(gate_eigenvalues_cov_eigenvalues::Vector{Float64})
    # Calculate the trace of the gate eigenvalue estimator covariance matrix, and its square
    sigma_tr = sum(gate_eigenvalues_cov_eigenvalues)
    sigma_sq_tr = sum(gate_eigenvalues_cov_eigenvalues .^ 2)
    # Calculate the expectation and variance of the NRMSE
    N = length(gate_eigenvalues_cov_eigenvalues)
    expectation = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    variance = (sigma_sq_tr / (2 * sigma_tr)) * (1 - sigma_sq_tr / (8 * sigma_tr^2)) / N
    return (expectation::Float64, variance::Float64)
end

"""
    calc_gate_probabilities_covariance(d::Design; ls_type::Symbol = :gls, est_type::Symbol = :ordinary, unpad::Bool = false)
    calc_gate_probabilities_covariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; ls_type::Symbol = :gls, est_type::Symbol = :ordinary, unpad::Bool = false)

Returns the padded gate error probabilities estimator covariance matrix for the least squares estimator specified by `ls_type` (`:gls`, `:wls`, or `:ols`) corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`, with estimator type `est_type` (`:ordinary`, `:marginal`, or `:relative`).
Instead returns the unpadded gate error probabilities estimator covariance matrix if `unpad` is `true`, which is stripped of the covariance with the identity error probabilities.
"""
function calc_gate_probabilities_covariance(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    ls_type::Symbol = :gls,
    est_type::Symbol = :ordinary,
    unpad::Bool = false,
)
    # Calculate the gate eigenvalue estimator covariance matrix
    gate_eigenvalues_cov =
        calc_ls_covariance(d, covariance_log; ls_type = ls_type, est_type = est_type)
    # Get the transform matrices
    pad_transform = get_pad_transform(d, est_type)
    wht_transform_inv = get_wht_transform(d, est_type; inverse = true)
    pad_transform_probs = get_pad_transform(d, est_type; probabilities = true)
    wht_transform = get_wht_transform(d, est_type)
    # Check the transform matrices
    @assert dropzeros!(pad_transform' * pad_transform) ≈ I
    @assert dropzeros!(pad_transform' * pad_transform_probs) ≈ I
    @assert dropzeros!(pad_transform_probs' * pad_transform) ≈ I
    probs_transform = pad_transform_probs' * wht_transform * pad_transform_probs
    @assert probs_transform ≈ pad_transform_probs' * wht_transform * pad_transform
    @assert probs_transform ≈ pad_transform' * wht_transform * pad_transform_probs
    probs_transform_inv = pad_transform' * wht_transform_inv * pad_transform
    @assert dropzeros!(probs_transform_inv * probs_transform) ≈ I
    @assert dropzeros!(probs_transform * probs_transform_inv) ≈ I
    # Calculate the gate probabilities estimator covariance matrix
    gate_probabilities_cov = Symmetric(
        wht_transform_inv *
        pad_transform *
        gate_eigenvalues_cov *
        pad_transform' *
        wht_transform_inv',
    )
    # Unpad the gate probability distributions if appropriate
    if unpad
        gate_probabilities_cov =
            Symmetric(pad_transform' * gate_probabilities_cov * pad_transform)
    end
    return gate_probabilities_cov::Symmetric{Float64, Matrix{Float64}}
end
function calc_gate_probabilities_covariance(
    d::Design;
    ls_type::Symbol = :gls,
    est_type::Symbol = :ordinary,
    unpad::Bool = false,
)
    covariance_log = calc_covariance_log(d)
    gate_probabilities_cov = calc_gate_probabilities_covariance(
        d,
        covariance_log;
        ls_type = ls_type,
        est_type = est_type,
        unpad = unpad,
    )
    return gate_probabilities_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    calc_precision_matrix(design_matrix::SparseMatrixCSC{Float64, Int}, gate_eigenvalues::Vector{Float64}, covariance_log_inv::SparseMatrixCSC{Float64, Int})
    calc_precision_matrix(d::Design, gate_eigenvalues::Vector{Float64}; diagonal::Bool = false)
    calc_precision_matrix(d::Design; diagonal::Bool = false)

Returns the precision matrix, namely the inverse of the generalised least squares gate eigenvalue estimator covariance matrix, corresponding to the design `d` with design matrix `design_matrix`, gate eigenvalues `gate_eigenvalues` and circuit log-eigenvalue estimator covariance matrix inverse `covariance_log_inv`.
Ordinarily this yields the precision matrix for generalise least squares, but if `diagonal` is `true` it calculates the inverse of the diagonal of the covariance matrix, rather than the inverse of the full matrix, yielding the precision matrix for both generalised and weighted least squares if the circuit eigenvalue estimators were uncorrelated.
"""
function calc_precision_matrix(
    design_matrix::SparseMatrixCSC{Float64, Int},
    gate_eigenvalues::Vector{Float64},
    covariance_log_inv::SparseMatrixCSC{Float64, Int},
)
    gate_eigenvalues_diag_inv = sparse(Diagonal(gate_eigenvalues .^ (-1)))
    precision_matrix = sparse(
        Symmetric(
            gate_eigenvalues_diag_inv *
            design_matrix' *
            covariance_log_inv *
            design_matrix *
            gate_eigenvalues_diag_inv,
        ),
    )
    return precision_matrix::SparseMatrixCSC{Float64, Int64}
end
function calc_precision_matrix(
    d::Design,
    gate_eigenvalues::Vector{Float64};
    diagonal::Bool = false,
)
    # Calculate the covariance matrix of the circuit log-eigenvalues
    design_matrix = convert(SparseMatrixCSC{Float64, Int}, d.matrix)
    eigenvalues = get_eigenvalues(d, gate_eigenvalues)
    covariance = calc_covariance(d, gate_eigenvalues; warning = ~diagonal)
    covariance_log = calc_covariance_log(covariance, eigenvalues)
    mapping_lengths = length.(d.mapping_ensemble)
    if diagonal
        covariance_log_inv = sparse(Diagonal(diag(covariance_log) .^ (-1)))
    else
        covariance_log_inv = sparse_covariance_inv(covariance_log, mapping_lengths)
    end
    # Calculate the precision matrix
    precision_matrix =
        calc_precision_matrix(design_matrix, gate_eigenvalues, covariance_log_inv)
    return precision_matrix::SparseMatrixCSC{Float64, Int64}
end
function calc_precision_matrix(d::Design; diagonal::Bool = false)
    precision_matrix =
        calc_precision_matrix(d, get_gate_eigenvalues(d); diagonal = diagonal)
    return precision_matrix::SparseMatrixCSC{Float64, Int64}
end

"""
    calc_gls_moments(d::Design; est_type::Symbol = :ordinary)
    calc_gls_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; est_type::Symbol = :ordinary)

Returns the expectation and variance of the normalised RMS error for the generalised least squares (GLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`, with estimator type `est_type` (`:ordinary`, `:marginal`, or `:relative`).
"""
function calc_gls_moments(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    est_type::Symbol = :ordinary,
)
    # Calculate the GLS moments from the GLS gate eigenvalue estimator covariance matrix
    gls_gate_eigenvalues_cov =
        calc_ls_covariance(d, covariance_log; ls_type = :gls, est_type = est_type)
    (gls_expectation, gls_variance) = nrmse_moments(gls_gate_eigenvalues_cov)
    return (gls_expectation::Float64, gls_variance::Float64)
end
function calc_gls_moments(d::Design; est_type::Symbol = :ordinary)
    covariance_log = calc_covariance_log(d)
    (gls_expectation, gls_variance) =
        calc_gls_moments(d, covariance_log; est_type = est_type)
    return (gls_expectation::Float64, gls_variance::Float64)
end

"""
    calc_wls_moments(d::Design; est_type::Symbol = :ordinary)
    calc_wls_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; est_type::Symbol = :ordinary)

Returns the expectation and variance of the normalised RMS error for the weighted least squares (WLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`, with estimator type `est_type` (`:ordinary`, `:marginal`, or `:relative`).
"""
function calc_wls_moments(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    est_type::Symbol = :ordinary,
)
    # Calculate the WLS moments from the WLS gate eigenvalue estimator covariance matrix
    wls_gate_eigenvalues_cov =
        calc_ls_covariance(d, covariance_log; ls_type = :wls, est_type = est_type)
    (wls_expectation, wls_variance) = nrmse_moments(wls_gate_eigenvalues_cov)
    return (wls_expectation::Float64, wls_variance::Float64)
end
function calc_wls_moments(d::Design; est_type::Symbol = :ordinary)
    covariance_log = calc_covariance_log(d)
    (wls_expectation, wls_variance) =
        calc_wls_moments(d, covariance_log; est_type = est_type)
    return (wls_expectation::Float64, wls_variance::Float64)
end

"""
    calc_ols_moments(d::Design; est_type::Symbol = :ordinary)
    calc_ols_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; est_type::Symbol = :ordinary)

Returns the expectation and variance of the normalised RMS error for the ordinary least squares (OLS) estimator corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`, with estimator type `est_type` (`:ordinary`, `:marginal`, or `:relative`).
"""
function calc_ols_moments(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    est_type::Symbol = :ordinary,
)
    # Calculate the OLS moments from the OLS gate eigenvalue estimator covariance matrix
    ols_gate_eigenvalues_cov =
        calc_ls_covariance(d, covariance_log; ls_type = :ols, est_type = est_type)
    (ols_expectation, ols_variance) = nrmse_moments(ols_gate_eigenvalues_cov)
    return (ols_expectation::Float64, ols_variance::Float64)
end
function calc_ols_moments(d::Design; est_type::Symbol = :ordinary)
    covariance_log = calc_covariance_log(d)
    (ols_expectation, ols_variance) =
        calc_ols_moments(d, covariance_log; est_type = est_type)
    return (ols_expectation::Float64, ols_variance::Float64)
end

"""
    calc_ls_moments(d::Design; ls_type::Symbol = :gls, est_type::Symbol = :ordinary, est_weight::Float64 = 0.5)
    calc_ls_moments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; ls_type::Symbol = :gls, est_type::Symbol = :ordinary, est_weight::Float64 = 0.5)

Returns the expectation and variance of the normalised RMS error for the least squares estimator specified by `ls_type` (`:gls`, `:wls`, or `:ols`) corresponding to the design `d` with circuit log-eigenvalue estimator covariance matrix `covariance_log`, with estimator type `est_type` (`:ordinary`, `:marginal`, or `:relative`).
If `est_type` is instead set to `:sum` or `:prod`, computes the arithmetic or geometric mean of the ordinary and relative precision moments, weighted by `est_weight`, the weighting allocated to the ordinary precision moments.
"""
function calc_ls_moments(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    ls_type::Symbol = :gls,
    est_type::Symbol = :ordinary,
    est_weight::Float64 = 0.5,
)
    @assert 0 <= est_weight <= 1 "The estimator weight must be between 0 and 1."
    if est_type == :sum || est_type == :prod
        ord_gate_eigenvalues_cov =
            calc_ls_covariance(d, covariance_log; ls_type = ls_type, est_type = :ordinary)
        rel_gate_eigenvalues_cov = get_relative_gate_covariance(d, ord_gate_eigenvalues_cov)
        (ord_expectation, ord_variance) = nrmse_moments(ord_gate_eigenvalues_cov)
        (rel_expectation, rel_variance) = nrmse_moments(rel_gate_eigenvalues_cov)
        if est_type == :sum
            expectation = est_weight * ord_expectation + (1 - est_weight) * rel_expectation
            variance = est_weight * ord_variance + (1 - est_weight) * rel_variance
        elseif est_type == :prod
            expectation = ord_expectation^(est_weight) * rel_expectation^(1 - est_weight)
            variance = ord_variance^(est_weight) * rel_variance^(1 - est_weight)
        else
            throw(error("Unsupported estimator type $(est_type)."))
        end
    else
        gate_eigenvalues_cov =
            calc_ls_covariance(d, covariance_log; ls_type = ls_type, est_type = est_type)
        (expectation, variance) = nrmse_moments(gate_eigenvalues_cov)
    end
    return (expectation::Float64, variance::Float64)
end
function calc_ls_moments(
    d::Design;
    ls_type::Symbol = :gls,
    est_type::Symbol = :ordinary,
    est_weight::Float64 = 0.5,
)
    covariance_log = calc_covariance_log(d)
    (expectation, variance) = calc_ls_moments(
        d,
        covariance_log;
        ls_type = ls_type,
        est_type = est_type,
        est_weight = est_weight,
    )
    return (expectation::Float64, variance::Float64)
end

"""
    calc_merit(d::Design; warning::Bool = true)
    calc_merit(d_rand::RandDesign; warning::Bool = true)

Returns a [`Merit`](@ref) object corresponding to the design `d` or randomised design `d_rand`, displaying a warning if the design does not have full covariance matrix data and `warning` is `true`.
"""
function calc_merit(d::Design; warning::Bool = true)
    if ~d.full_covariance
        if warning
            println(
                "The design does not have full covariance matrix data. Generating this data.\nBE CAREFUL: this may consume extreme amounts of memory and be extremely slow.",
            )
        end
        d_full = get_full_design(d)
        covariance_log = calc_covariance_log(d_full)
        wls_gate_eigenvalues_cov = calc_wls_covariance(d_full, covariance_log)
        gls_gate_eigenvalues_cov = deepcopy(wls_gate_eigenvalues_cov)
        ols_gate_eigenvalues_cov = calc_ols_covariance(d_full, covariance_log)
    else
        covariance_log = calc_covariance_log(d)
        gls_gate_eigenvalues_cov = calc_gls_covariance(d, covariance_log)
        wls_gate_eigenvalues_cov = calc_wls_covariance(d, covariance_log)
        ols_gate_eigenvalues_cov = calc_ols_covariance(d, covariance_log)
    end
    # Calculate the GLS merit quantities
    gls_marginal_gate_eigenvalues_cov =
        get_marginal_gate_covariance(d, gls_gate_eigenvalues_cov)
    gls_relative_gate_eigenvalues_cov =
        get_relative_gate_covariance(d, gls_gate_eigenvalues_cov)
    # Calculate the eigenvalues of the covariance matrices
    gls_cov_eigenvalues = eigvals(gls_gate_eigenvalues_cov)
    gls_marginal_cov_eigenvalues = eigvals(gls_marginal_gate_eigenvalues_cov)
    gls_relative_cov_eigenvalues = eigvals(gls_relative_gate_eigenvalues_cov)
    # Calculate the expected NRMSE and NRMSE variances
    (gls_expectation, gls_variance) = nrmse_moments(gls_cov_eigenvalues)
    (gls_marginal_expectation, gls_marginal_variance) =
        nrmse_moments(gls_marginal_cov_eigenvalues)
    (gls_relative_expectation, gls_relative_variance) =
        nrmse_moments(gls_relative_cov_eigenvalues)
    # Calculate the WLS merit quantities
    wls_marginal_gate_eigenvalues_cov =
        get_marginal_gate_covariance(d, wls_gate_eigenvalues_cov)
    wls_relative_gate_eigenvalues_cov =
        get_relative_gate_covariance(d, wls_gate_eigenvalues_cov)
    # Calculate the eigenvalues of the covariance matrices
    wls_cov_eigenvalues = eigvals(wls_gate_eigenvalues_cov)
    wls_marginal_cov_eigenvalues = eigvals(wls_marginal_gate_eigenvalues_cov)
    wls_relative_cov_eigenvalues = eigvals(wls_relative_gate_eigenvalues_cov)
    # Calculate the expected NRMSE and NRMSE variances
    (wls_expectation, wls_variance) = nrmse_moments(wls_cov_eigenvalues)
    (wls_marginal_expectation, wls_marginal_variance) =
        nrmse_moments(wls_marginal_cov_eigenvalues)
    (wls_relative_expectation, wls_relative_variance) =
        nrmse_moments(wls_relative_cov_eigenvalues)
    # Calculate the OLS merit quantities
    ols_marginal_gate_eigenvalues_cov =
        get_marginal_gate_covariance(d, ols_gate_eigenvalues_cov)
    ols_relative_gate_eigenvalues_cov =
        get_relative_gate_covariance(d, ols_gate_eigenvalues_cov)
    # Calculate the eigenvalues of the covariance matrices
    ols_cov_eigenvalues = eigvals(ols_gate_eigenvalues_cov)
    ols_marginal_cov_eigenvalues = eigvals(ols_marginal_gate_eigenvalues_cov)
    ols_relative_cov_eigenvalues = eigvals(ols_relative_gate_eigenvalues_cov)
    # Calculate the expected NRMSE and NRMSE variances
    (ols_expectation, ols_variance) = nrmse_moments(ols_cov_eigenvalues)
    (ols_marginal_expectation, ols_marginal_variance) =
        nrmse_moments(ols_marginal_cov_eigenvalues)
    (ols_relative_expectation, ols_relative_variance) =
        nrmse_moments(ols_relative_cov_eigenvalues)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_gram = d.matrix' * d.matrix
    design_matrix_cond_num = 0
    design_matrix_pinv_norm = 0
    try
        largest_singular_value =
            sqrt(Arpack.eigs(design_matrix_gram; nev = 1, which = :LM)[1][1])
        smallest_singular_value =
            sqrt(Arpack.eigs(design_matrix_gram; nev = 1, which = :LM, sigma = 0.0)[1][1])
        design_matrix_cond_num = largest_singular_value / smallest_singular_value
        design_matrix_pinv_norm = 1 / smallest_singular_value
    catch
        println("Pseudoinverse norm and condition number calculation failed.")
    end
    # Construct the merit object
    merit = Merit(
        d.c.circuit_param,
        d.c.noise_param,
        d.tuple_set_data,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
        d.c.gate_data.N,
        d.c.gate_data.N_marginal,
        d.c.gate_data.N_relative,
        length(d.c.total_gates),
        gls_expectation,
        gls_variance,
        gls_cov_eigenvalues,
        gls_marginal_expectation,
        gls_marginal_variance,
        gls_marginal_cov_eigenvalues,
        gls_relative_expectation,
        gls_relative_variance,
        gls_relative_cov_eigenvalues,
        wls_expectation,
        wls_variance,
        wls_cov_eigenvalues,
        wls_marginal_expectation,
        wls_marginal_variance,
        wls_marginal_cov_eigenvalues,
        wls_relative_expectation,
        wls_relative_variance,
        wls_relative_cov_eigenvalues,
        ols_expectation,
        ols_variance,
        ols_cov_eigenvalues,
        ols_marginal_expectation,
        ols_marginal_variance,
        ols_marginal_cov_eigenvalues,
        ols_relative_expectation,
        ols_relative_variance,
        ols_relative_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
    )
    return merit::Merit
end
function calc_merit(d_rand::RandDesign; warning::Bool = true)
    d = get_design(d_rand)
    merit = calc_merit(d)
    return merit::Merit
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
    nrmse_pdf(cov_eigenvalues::Vector{Float64}, nrmse_values::Vector{Float64}; epsilon::Real = 1e-5)

Returns the probability density function (PDF) for the normalised RMS error (NRMSE) of the gate eigenvalue estimator vector, which follows a generalised chi-squared distribution and whose covariance matrix has eigenvalues `cov_eigenvalues`, at the coordinates specified by `nrmse_values`.
Does not calculate values when the normal approximation to the PDF is less than a factor of `epsilon` of its maximum value.

Calculation follows Eq. 3.2 of `Computing the distribution of quadratic forms in normal variables` by J. P. Imhof (1961).
"""
function nrmse_pdf(
    cov_eigenvalues::Vector{Float64},
    nrmse_values::Vector{Float64};
    epsilon::Real = 1e-5,
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
               QuadGK.quadgk(
            u -> nrmse_pdf_integrand(u, x, norm_cov_eigenvalues),
            0,
            Inf,
        )[1] / π
    end
    nrmse_pdf = zeros(length(nrmse_values))
    for (idx, x) in pairs(nrmse_values)
        # Only perform the calculation for sufficiently large PDF values
        if nrmse_pdf_normal[idx] / nrmse_pdf_normal_max > epsilon
            # Differentiate the CDF to obtain the PDF
            # Note that Imhof's method describes the square, whereas we seek the square root
            nrmse_pdf[idx] = 2 * x * FiniteDifferences.central_fdm(5, 1)(imhof_cdf, x^2)
        end
    end
    return nrmse_pdf::Vector{Float64}
end

"""
    NoiseEstimate

Noise estimate data structure.

# Fields

  - `eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}}`: Circuit eigenvalue estimates for each tuple and circuit eigenvalue.
  - `covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}}`: Covariance circuit eigenvalue estimates for each tuple and covariance circuit eigenvalue.
  - `shot_budget::Int`: Shot budget for the (potentially randomised) experimental ensemble.
  - `gls_unproj_gate_eigenvalues::Vector{Float64}`: Generalised least squares unprojected gate eigenvalues.
  - `gls_gate_eigenvalues::Vector{Float64}`: Generalised least squares projected gate eigenvalues.
  - `gls_unproj_gate_probabilities::Dict{Gate, Vector{Float64}}`: Generalised least squares unprojected gate error probabilities.
  - `gls_gate_probabilities::Dict{Gate, Vector{Float64}}`: Generalised least squares projected gate error probabilities.
  - `gls_unproj_marginal_gate_eigenvalues::Vector{Float64}`: Generalised least squares unprojected marginal gate eigenvalues.
  - `gls_marginal_gate_eigenvalues::Vector{Float64}`: Generalised least squares projected marginal gate eigenvalues.
  - `gls_unproj_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}`: Generalised least squares unprojected marginal gate error probabilities.
  - `gls_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}`: Generalised least squares projected marginal gate error probabilities.
  - `gls_unproj_relative_gate_eigenvalues::Vector{Float64}`: Generalised least squares unprojected relative gate eigenvalues.
  - `gls_relative_gate_eigenvalues::Vector{Float64}`: Generalised least squares projected relative gate eigenvalues.
  - `wls_unproj_gate_eigenvalues::Vector{Float64}`: Weighted least squares unprojected gate eigenvalues.
  - `wls_gate_eigenvalues::Vector{Float64}`: Weighted least squares projected gate eigenvalues.
  - `wls_unproj_gate_probabilities::Dict{Gate, Vector{Float64}}`: Weighted least squares unprojected gate error probabilities.
  - `wls_gate_probabilities::Dict{Gate, Vector{Float64}}`: Weighted least squares projected gate error probabilities.
  - `wls_unproj_marginal_gate_eigenvalues::Vector{Float64}`: Weighted least squares unprojected marginal gate eigenvalues.
  - `wls_marginal_gate_eigenvalues::Vector{Float64}`: Weighted least squares projected marginal gate eigenvalues.
  - `wls_unproj_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}`: Weighted least squares unprojected marginal gate error probabilities.
  - `wls_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}`: Weighted least squares projected marginal gate error probabilities.
  - `wls_unproj_relative_gate_eigenvalues::Vector{Float64}`: Weighted least squares unprojected relative gate eigenvalues.
  - `wls_relative_gate_eigenvalues::Vector{Float64}`: Weighted least squares projected relative gate eigenvalues.
  - `ols_unproj_gate_eigenvalues::Vector{Float64}`: Ordinary least squares unprojected gate eigenvalues.
  - `ols_gate_eigenvalues::Vector{Float64}`: Ordinary least squares projected gate eigenvalues.
  - `ols_unproj_gate_probabilities::Dict{Gate, Vector{Float64}}`: Ordinary least squares unprojected gate error probabilities.
  - `ols_gate_probabilities::Dict{Gate, Vector{Float64}}`: Ordinary least squares projected gate error probabilities.
  - `ols_unproj_marginal_gate_eigenvalues::Vector{Float64}`: Ordinary least squares unprojected marginal gate eigenvalues.
  - `ols_marginal_gate_eigenvalues::Vector{Float64}`: Ordinary least squares projected marginal gate eigenvalues.
  - `ols_unproj_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}`: Ordinary least squares unprojected marginal gate error probabilities.
  - `ols_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}`: Ordinary least squares projected marginal gate error probabilities.
  - `ols_unproj_relative_gate_eigenvalues::Vector{Float64}`: Ordinary least squares unprojected relative gate eigenvalues.
  - `ols_relative_gate_eigenvalues::Vector{Float64}`: Ordinary least squares projected relative gate eigenvalues.
"""
struct NoiseEstimate
    eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}}
    covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}}
    shot_budget::Int
    gls_unproj_gate_eigenvalues::Vector{Float64}
    gls_gate_eigenvalues::Vector{Float64}
    gls_unproj_gate_probabilities::Dict{Gate, Vector{Float64}}
    gls_gate_probabilities::Dict{Gate, Vector{Float64}}
    gls_unproj_marginal_gate_eigenvalues::Vector{Float64}
    gls_marginal_gate_eigenvalues::Vector{Float64}
    gls_unproj_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
    gls_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
    gls_unproj_relative_gate_eigenvalues::Vector{Float64}
    gls_relative_gate_eigenvalues::Vector{Float64}
    wls_unproj_gate_eigenvalues::Vector{Float64}
    wls_gate_eigenvalues::Vector{Float64}
    wls_unproj_gate_probabilities::Dict{Gate, Vector{Float64}}
    wls_gate_probabilities::Dict{Gate, Vector{Float64}}
    wls_unproj_marginal_gate_eigenvalues::Vector{Float64}
    wls_marginal_gate_eigenvalues::Vector{Float64}
    wls_unproj_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
    wls_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
    wls_unproj_relative_gate_eigenvalues::Vector{Float64}
    wls_relative_gate_eigenvalues::Vector{Float64}
    ols_unproj_gate_eigenvalues::Vector{Float64}
    ols_gate_eigenvalues::Vector{Float64}
    ols_unproj_gate_probabilities::Dict{Gate, Vector{Float64}}
    ols_gate_probabilities::Dict{Gate, Vector{Float64}}
    ols_unproj_marginal_gate_eigenvalues::Vector{Float64}
    ols_marginal_gate_eigenvalues::Vector{Float64}
    ols_unproj_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
    ols_marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
    ols_unproj_relative_gate_eigenvalues::Vector{Float64}
    ols_relative_gate_eigenvalues::Vector{Float64}
end

function Base.show(io::IO, noise_est::NoiseEstimate)
    print(io, "Noise estimate data.")
    return nothing
end

@struct_hash_equal_isequal NoiseEstimate

"""
    NoiseError

Noise estimate error data structure containing normalised root-mean-square errors (NRMSE) for a [`NoiseEstimate`](@ref).

# Fields

  - `gls_nrmse::Float64`: Generalised least squares normalised root-mean-square error.
  - `gls_proj_nrmse::Float64`: Generalised least squares projected normalised root-mean-square error.
  - `gls_marginal_nrmse::Float64`: Generalised least squares marginal normalised root-mean-square error.
  - `gls_proj_marginal_nrmse::Float64`: Generalised least squares projected marginal normalised root-mean-square error.
  - `gls_relative_nrmse::Float64`: Generalised least squares relative normalised root-mean-square error.
  - `gls_proj_relative_nrmse::Float64`: Generalised least squares projected relative normalised root-mean-square error.
  - `wls_nrmse::Float64`: Weighted least squares normalised root-mean-square error.
  - `wls_proj_nrmse::Float64`: Weighted least squares projected normalised root-mean-square error.
  - `wls_marginal_nrmse::Float64`: Weighted least squares marginal normalised root-mean-square error.
  - `wls_proj_marginal_nrmse::Float64`: Weighted least squares projected marginal normalised root-mean-square error.
  - `wls_relative_nrmse::Float64`: Weighted least squares relative normalised root-mean-square error.
  - `wls_proj_relative_nrmse::Float64`: Weighted least squares projected relative normalised root-mean-square error.
  - `ols_nrmse::Float64`: Ordinary least squares normalised root-mean-square error.
  - `ols_proj_nrmse::Float64`: Ordinary least squares projected normalised root-mean-square error.
  - `ols_marginal_nrmse::Float64`: Ordinary least squares marginal normalised root-mean-square error.
  - `ols_proj_marginal_nrmse::Float64`: Ordinary least squares projected marginal normalised root-mean-square error.
  - `ols_relative_nrmse::Float64`: Ordinary least squares relative normalised root-mean-square error.
  - `ols_proj_relative_nrmse::Float64`: Ordinary least squares projected relative normalised root-mean-square error.
"""
struct NoiseError
    gls_nrmse::Float64
    gls_proj_nrmse::Float64
    gls_marginal_nrmse::Float64
    gls_proj_marginal_nrmse::Float64
    gls_relative_nrmse::Float64
    gls_proj_relative_nrmse::Float64
    wls_nrmse::Float64
    wls_proj_nrmse::Float64
    wls_marginal_nrmse::Float64
    wls_proj_marginal_nrmse::Float64
    wls_relative_nrmse::Float64
    wls_proj_relative_nrmse::Float64
    ols_nrmse::Float64
    ols_proj_nrmse::Float64
    ols_marginal_nrmse::Float64
    ols_proj_marginal_nrmse::Float64
    ols_relative_nrmse::Float64
    ols_proj_relative_nrmse::Float64
end

function Base.show(io::IO, noise_error::NoiseError)
    nrmse_table =
        round.(
            [
                noise_error.gls_nrmse noise_error.gls_marginal_nrmse noise_error.gls_relative_nrmse noise_error.gls_proj_nrmse noise_error.gls_proj_marginal_nrmse noise_error.gls_proj_relative_nrmse
                noise_error.wls_nrmse noise_error.wls_marginal_nrmse noise_error.wls_relative_nrmse noise_error.wls_proj_nrmse noise_error.wls_proj_marginal_nrmse noise_error.wls_proj_relative_nrmse
                noise_error.ols_nrmse noise_error.ols_marginal_nrmse noise_error.ols_relative_nrmse noise_error.ols_proj_nrmse noise_error.ols_proj_marginal_nrmse noise_error.ols_proj_relative_nrmse
            ],
            digits = 4,
        )
    label_nrmse_table = hcat(["GLS"; "WLS"; "OLS"], string.(nrmse_table))
    header = [
        "Type"
        "Ordinary NRMSE"
        "Marginal NRMSE"
        "Relative NRMSE"
        "Proj. Ordinary NRMSE"
        "Proj. Marginal NRMSE"
        "Proj. Relative NRMSE"
    ]
    pretty_table(io, label_nrmse_table; header = header, alignment = :l)
    return nothing
end

@struct_hash_equal_isequal NoiseError

"""
    NoiseScore

Noise estimate error z-score data structure, containing the z-scores for a [`NoiseError`](@ref) with respect to some predicted [`Merit`](@ref).

# Fields

  - `gls_z_score::Float64`: Generalised least squares z-score.
  - `gls_proj_z_score::Float64`: Generalised least squares projected z-score.
  - `gls_marginal_z_score::Float64`: Generalised least squares marginal z-score.
  - `gls_proj_marginal_z_score::Float64`: Generalised least squares projected marginal z-score.
  - `gls_relative_z_score::Float64`: Generalised least squares relative z-score.
  - `gls_proj_relative_z_score::Float64`: Generalised least squares projected relative z-score.
  - `wls_z_score::Float64`: Weighted least squares z-score.
  - `wls_proj_z_score::Float64`: Weighted least squares projected z-score.
  - `wls_marginal_z_score::Float64`: Weighted least squares marginal z-score.
  - `wls_proj_marginal_z_score::Float64`: Weighted least squares projected marginal z-score.
  - `wls_relative_z_score::Float64`: Weighted least squares relative z-score.
  - `wls_proj_relative_z_score::Float64`: Weighted least squares projected relative z-score.
  - `ols_z_score::Float64`: Ordinary least squares z-score.
  - `ols_proj_z_score::Float64`: Ordinary least squares projected z-score.
  - `ols_marginal_z_score::Float64`: Ordinary least squares marginal z-score.
  - `ols_proj_marginal_z_score::Float64`: Ordinary least squares projected marginal z-score.
  - `ols_relative_z_score::Float64`: Ordinary least squares relative z-score.
  - `ols_proj_relative_z_score::Float64`: Ordinary least squares projected relative z-score.
"""
struct NoiseScore
    gls_z_score::Float64
    gls_proj_z_score::Float64
    gls_marginal_z_score::Float64
    gls_proj_marginal_z_score::Float64
    gls_relative_z_score::Float64
    gls_proj_relative_z_score::Float64
    wls_z_score::Float64
    wls_proj_z_score::Float64
    wls_marginal_z_score::Float64
    wls_proj_marginal_z_score::Float64
    wls_relative_z_score::Float64
    wls_proj_relative_z_score::Float64
    ols_z_score::Float64
    ols_proj_z_score::Float64
    ols_marginal_z_score::Float64
    ols_proj_marginal_z_score::Float64
    ols_relative_z_score::Float64
    ols_proj_relative_z_score::Float64
end

function Base.show(io::IO, noise_score::NoiseScore)
    z_score_table =
        round.(
            [
                noise_score.gls_z_score noise_score.gls_marginal_z_score noise_score.gls_relative_z_score noise_score.gls_proj_z_score noise_score.gls_proj_marginal_z_score noise_score.gls_proj_relative_z_score
                noise_score.wls_z_score noise_score.wls_marginal_z_score noise_score.wls_relative_z_score noise_score.wls_proj_z_score noise_score.wls_proj_marginal_z_score noise_score.wls_proj_relative_z_score
                noise_score.ols_z_score noise_score.ols_marginal_z_score noise_score.ols_relative_z_score noise_score.ols_proj_z_score noise_score.ols_proj_marginal_z_score noise_score.ols_proj_relative_z_score
            ],
            digits = 4,
        )
    label_zscore_table = hcat(["GLS"; "WLS"; "OLS"], string.(z_score_table))
    header = [
        "Type"
        "Ordinary score"
        "Marginal score"
        "Relative score"
        "Proj. Ordinary score"
        "Proj. Marginal score"
        "Proj. Relative score"
    ]
    pretty_table(io, label_zscore_table; header = header, alignment = :l)
    return nothing
end

@struct_hash_equal_isequal NoiseScore

"""
    pauli_estimate(stim_counts::Matrix{UInt8}, pauli_meas_indices::Vector{Int}, shots_value::Integer)
    pauli_estimate(qiskit_counts::Dict{String,Int}, pauli_meas_indices::Vector{Int}, shots_value::Integer)

Returns the sample average from the counts `stim_counts` or `qiskit_counts` of the Pauli operator measured at the indices `pauli_meas_indices` across `shots_value` samples.

Note that [`parse_qiskit_dict`](@ref) parses Qiskit formatted counts into Stim formatted counts.
"""
function pauli_estimate(
    stim_counts::Matrix{UInt8},
    pauli_meas_indices::Vector{Int},
    shots_value::Integer,
)
    pauli_shots = 0
    for s in 1:shots_value
        pauli_shot = 0
        for m in pauli_meas_indices
            # The measurement outcomes are packed into UInt8
            pauli_shot += (stim_counts[s, 1 + fld(m - 1, 8)] >> ((m - 1) % 8)) & 1
        end
        pauli_shots += 1 - 2 * (pauli_shot % 2)
    end
    return pauli_shots::Int
end
function pauli_estimate(
    qiskit_counts::Dict{String, Int},
    pauli_meas_indices::Vector{Int},
    shots_value::Integer,
)
    pauli_shots = 0
    total_shots = 0
    for (key, value) in pairs(qiskit_counts)
        # Qiskit uses little-endian order
        key_reverse = reverse(key)
        pauli_shot = 0
        for m in pauli_meas_indices
            # The measurement outcomes are stored as strings of 0 and 1
            pauli_shot += key_reverse[m] == '1' ? 1 : 0
        end
        pauli_shots += value * (1 - 2 * (pauli_shot % 2))
        total_shots += value
    end
    if total_shots != shots_value
        @warn "The total number of Qiskit shots $(total_shots) does not match the expected value $(shots_value); something is wrong!"
    end
    return pauli_shots::Int
end

"""
    estimate_experiment(counts::Union{Matrix{UInt8}, Dict{String, Int}}, experiment_meas_indices::Vector{Vector{Int}}, experiment_pauli_signs::Vector{Bool}, shots_value::Integer)

Returns the sample average estimator contributions estimated from the counts `counts` for the circuit eigenvalues measured at `experiment_meas_indices` whose Pauli mappings have signs `experiment_pauli_signs`, across `shots_value` samples.
"""
function estimate_experiment(
    counts::Union{Matrix{UInt8}, Dict{String, Int}},
    experiment_meas_indices::Vector{Vector{Int}},
    experiment_pauli_signs::Vector{Bool},
    shots_value::Integer,
)
    # Initialise variables
    E = length(experiment_meas_indices)
    @assert E == length(experiment_pauli_signs) "The number of signs does not match the number of circuit eigenvalues."
    # Estimate the circuit eigenvalues corrected for their sign
    est_eigenvalues_experiment = zeros(Float64, E)
    @threads :static for idx in 1:E
        pauli_shots = pauli_estimate(counts, experiment_meas_indices[idx], shots_value)
        est_eigenvalues_experiment[idx] =
            ((-1)^experiment_pauli_signs[idx]) * pauli_shots / shots_value
    end
    return est_eigenvalues_experiment::Vector{Float64}
end

"""
    get_experiment_data(d::Design; qiskit_qubit_map::Union{Vector{Int}, Nothing} = nothing)

Returns the experiment data for the design `d`, including the experiment ensemble, measurement indices ensemble, and Pauli sign ensemble.
"""
function get_experiment_data(
    d::Design;
    qiskit_qubit_map::Union{Vector{Int}, Nothing} = nothing,
)
    # Initialise variables
    tuple_number = length(d.tuple_set)
    n = d.c.qubit_num
    experiment_ensemble = d.experiment_ensemble
    if qiskit_qubit_map !== nothing
        # Qiskit indexed from 0, so we must add 1
        @assert length(qiskit_qubit_map) == n "The qubit map is inconsistent with the number of qubits."
        meas_indices_ensemble = [
            [qiskit_qubit_map[get_support(m.final)] .+ 1 for m in d.mapping_ensemble[i]] for i in 1:tuple_number
        ]
    else
        meas_indices_ensemble =
            [[get_support(m.final) for m in d.mapping_ensemble[i]] for i in 1:tuple_number]
    end
    pauli_sign_ensemble =
        [[m.final.pauli[2n + 1] for m in d.mapping_ensemble[i]] for i in 1:tuple_number]
    return (
        experiment_ensemble::Vector{Vector{Vector{Int}}},
        meas_indices_ensemble::Vector{Vector{Vector{Int}}},
        pauli_sign_ensemble::Vector{Vector{Bool}},
    )
end

"""
    get_covariance_experiment_data(d::Design, covariance_keys_ensemble::Vector{Vector{CartesianIndex{2}}}, covariance_key_indices_ensemble::Vector{Dict{CartesianIndex{2}, Int}}; qiskit_qubit_map::Union{Vector{Int}, Nothing} = nothing)

Returns the covariance experiment data for the design `d`, including the covariance experiment ensemble, covariance measurement indices ensemble, and covariance Pauli sign ensemble.
"""
function get_covariance_experiment_data(
    d::Design,
    covariance_keys_ensemble::Vector{Vector{CartesianIndex{2}}},
    covariance_key_indices_ensemble::Vector{Dict{CartesianIndex{2}, Int}};
    qiskit_qubit_map::Union{Vector{Int}, Nothing} = nothing,
)
    # Initialise parameters
    tuple_number = length(d.tuple_set)
    n = d.c.qubit_num
    covariance_experiment_ensemble = Vector{Vector{Vector{Int}}}(undef, tuple_number)
    for i in 1:tuple_number
        if length(covariance_keys_ensemble[i]) == 0
            covariance_experiment_ensemble[i] = Vector{Vector{Int}}()
        else
            experiment_number = d.experiment_numbers[i]
            experiment_set = d.experiment_ensemble[i]
            covariance_dict = d.covariance_dict_ensemble[i]
            covariance_key_indices = covariance_key_indices_ensemble[i]
            covariance_experiment_ensemble[i] =
                Vector{Vector{Int}}(undef, experiment_number)
            for j in 1:experiment_number
                experiment = experiment_set[j]
                E = length(experiment)
                covariance_experiment = Vector{Int}()
                for pauli_idx_1 in 1:E
                    for pauli_idx_2 in (pauli_idx_1 + 1):E
                        pauli_idx =
                            CartesianIndex(experiment[pauli_idx_1], experiment[pauli_idx_2])
                        if haskey(covariance_dict, pauli_idx)
                            push!(covariance_experiment, covariance_key_indices[pauli_idx])
                        end
                    end
                end
                covariance_experiment_ensemble[i][j] = covariance_experiment
            end
        end
    end
    if qiskit_qubit_map !== nothing
        # Qiskit indexed from 0, so we must add 1
        @assert length(qiskit_qubit_map) == n "The qubit map is inconsistent with the number of qubits."
        covariance_meas_indices_ensemble = [
            [
                qiskit_qubit_map[get_support(
                    d.covariance_dict_ensemble[i][covariance_key][1].final,
                )] .+ 1 for covariance_key in covariance_keys_ensemble[i]
            ] for i in 1:tuple_number
        ]
    else
        covariance_meas_indices_ensemble = [
            [
                get_support(d.covariance_dict_ensemble[i][covariance_key][1].final) for
                covariance_key in covariance_keys_ensemble[i]
            ] for i in 1:tuple_number
        ]
    end
    covariance_pauli_sign_ensemble = [
        [
            d.covariance_dict_ensemble[i][covariance_key][1].final.pauli[2n + 1] for
            covariance_key in covariance_keys_ensemble[i]
        ] for i in 1:tuple_number
    ]
    return (
        covariance_experiment_ensemble::Vector{Vector{Vector{Int}}},
        covariance_meas_indices_ensemble::Vector{Vector{Vector{Int}}},
        covariance_pauli_sign_ensemble::Vector{Vector{Bool}},
    )
end
function get_covariance_experiment_data(
    d::Design;
    qiskit_qubit_map::Union{Vector{Int}, Nothing} = nothing,
)
    (covariance_keys_ensemble, covariance_key_indices_ensemble) =
        get_covariance_mapping_data(d)
    (
        covariance_experiment_ensemble,
        covariance_meas_indices_ensemble,
        covariance_pauli_sign_ensemble,
    ) = get_covariance_experiment_data(
        d,
        covariance_keys_ensemble,
        covariance_key_indices_ensemble;
        qiskit_qubit_map = qiskit_qubit_map,
    )
    return (
        covariance_experiment_ensemble::Vector{Vector{Vector{Int}}},
        covariance_meas_indices_ensemble::Vector{Vector{Vector{Int}}},
        covariance_pauli_sign_ensemble::Vector{Vector{Bool}},
    )
end

"""
    get_clipped_indices(est_eigenvalues::Vector{Float64}; min_eigenvalue::Real = 0.1, warning::Bool = true)

Returns clipped indices for the estimated circuit eigenvalues `est_eigenvalues` that not less than `min_eigenvalue`, warning if eigenvalues are clipped and `warning` is true.
"""
function get_clipped_indices(
    est_eigenvalues::Vector{Float64};
    min_eigenvalue::Real = 0.1,
    warning::Bool = true,
)
    # Find eigenvalues less than min_eigenvalue
    @assert min_eigenvalue > 0.0 "The eigenvalue clipping threshold must be positive."
    M = length(est_eigenvalues)
    clip_indices = findall(est_eigenvalues .< min_eigenvalue)
    if warning && length(clip_indices) > 0
        println(
            "$(length(clip_indices)) of $(M) estimated eigenvalues are smaller than $(min_eigenvalue) and will be omitted from the estimation procedure, the smallest being $(minimum(est_eigenvalues[clip_indices])), and the largest being $(maximum(est_eigenvalues[clip_indices])).",
        )
    end
    # Throw an error if eigenvalues are larger than 1, as this indicates a serious problem
    large_indices = findall(est_eigenvalues .> 1.0)
    if length(large_indices) > 0
        throw(
            error(
                "$(length(large_indices)) of $(M) estimated eigenvalues are larger than 1, the smallest being $(minimum(clipped_eigenvalues[large_indices])), and the largest being $(maximum(clipped_eigenvalues[large_indices])).",
            ),
        )
    end
    clipped_indices = setdiff(1:M, clip_indices)
    return clipped_indices::Vector{Int}
end

"""
    get_clipped_mapping_lengths(mapping_lengths::Vector{Int}, clipped_indices::Vector{Int})

Returns updated `mapping_lengths` given the clipped indices `clipped_indices`.
"""
function get_clipped_mapping_lengths(
    mapping_lengths::Vector{Int},
    clipped_indices::Vector{Int},
)
    # Get the clipped indices
    M = sum(mapping_lengths)
    @assert clipped_indices ⊆ collect(1:M) "The clipped indices are not a subset of the mapping indices."
    @assert clipped_indices == sort(unique(clipped_indices)) "The clipped indices are not both sorted and unique."
    # Clip the mapping lengths
    tuple_number = length(mapping_lengths)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    clipped_mapping_lengths = [
        length(intersect(mapping_lower[i]:mapping_upper[i], clipped_indices)) for
        i in 1:tuple_number
    ]
    return clipped_mapping_lengths::Vector{Int}
end

"""
    estimate_gate_eigenvalues(design_matrix::SparseMatrixCSC{Float64, Int}, est_covariance_log_inv_factor::SparseMatrixCSC{Float64, Int}, est_eigenvalues::Vector{Float64}; constrain::Bool = false)
    estimate_gate_eigenvalues(design_matrix::SparseMatrixCSC{Float64, Int}, est_eigenvalues::Vector{Float64}; constrain::Bool = false)

Returns the gate eigenvalues estimated from the estimated circuit eigenvalues `est_eigenvalues` by least squares inversion with the design matrix `design_matrix`, both weighted by the factor `est_covariance_log_inv_factor` if appropriate, constraining the gate eigenvalues to be at most 1 if `constrain` is true.
"""
function estimate_gate_eigenvalues(
    design_matrix::SparseMatrixCSC{Float64, Int},
    est_covariance_log_inv_factor::SparseMatrixCSC{Float64, Int},
    est_eigenvalues::Vector{Float64};
    constrain::Bool = false,
)
    # Estimate the gate log-eigenvalues, normalised by the factor matrix
    est_log_eigenvalues = -log.(est_eigenvalues)
    est_gate_log_eigenvalues =
        (est_covariance_log_inv_factor * design_matrix) \
        (est_covariance_log_inv_factor * est_log_eigenvalues)
    # Ensure that the (negative) gate log-eigenvalues are not less than 0, which corresponds to the gate eigenvalues being greater than 1
    if constrain
        est_gate_log_eigenvalues[(est_gate_log_eigenvalues .< 0.0)] .= 0.0
    end
    est_gate_eigenvalues = exp.(-est_gate_log_eigenvalues)
    return est_gate_eigenvalues::Vector{Float64}
end
function estimate_gate_eigenvalues(
    design_matrix::SparseMatrixCSC{Float64, Int},
    est_eigenvalues::Vector{Float64};
    constrain::Bool = false,
)
    M = length(est_eigenvalues)
    est_covariance_log_inv_factor = sparse(Diagonal(ones(Float64, M)))
    est_gate_eigenvalues = estimate_gate_eigenvalues(
        design_matrix,
        est_covariance_log_inv_factor,
        est_eigenvalues;
        constrain = constrain,
    )
    return est_gate_eigenvalues::Vector{Float64}
end

"""
    full_project_gate_eigenvalues(est_unproj_gate_eigenvalues::Vector{Float64}, est_precision_matrix::SparseMatrixCSC{Float64, Int64}, gate_data::GateData; precision::Real = 1e-8, diagnostics::Bool = false)

Returns the projection into the probability simplex of the gate eigenvalues `est_unproj_gate_eigenvalues` and the corresponding gate probabilities, projecting in the Mahalanobis distance using the precision matrix `est_precision_matrix`, with corresponding gate data `gate_data`.
The solver precision is controlled by `precision`, and solver diagnostics are printed if `diagnostics` is `true`.
"""
function full_project_gate_eigenvalues(
    est_unproj_gate_eigenvalues::Vector{Float64},
    est_precision_matrix::SparseMatrixCSC{Float64, Int64},
    gate_data::GateData;
    precision::Real = 1e-8,
    diagnostics::Bool = false,
)
    # Initialise the transform matrices
    pad_transform = get_pad_transform(gate_data)
    pad_transform_probs = get_pad_transform(gate_data; probabilities = true)
    pad_mask = get_pad_mask(gate_data)
    wht_transform = get_wht_transform(gate_data)
    wht_transform_inv = get_wht_transform(gate_data; inverse = true)
    # Generate the unprojected gate probabilities
    est_unproj_gate_probabilities_vec =
        wht_transform_inv * (pad_transform * est_unproj_gate_eigenvalues + pad_mask)
    est_unproj_unpad_gate_probabilities_vec =
        pad_transform' * est_unproj_gate_probabilities_vec
    # Generate the gate probabilities precision matrix
    probabilities_transform = pad_transform_probs' * wht_transform * pad_transform
    est_unpad_probabilities_precision_matrix =
        probabilities_transform * est_precision_matrix * probabilities_transform'
    # Project the gate probabilities
    est_unpad_gate_probabilities_vec = scs_project_nonnegative(
        est_unproj_unpad_gate_probabilities_vec,
        est_unpad_probabilities_precision_matrix;
        precision = precision,
        diagnostics = diagnostics,
    )
    # Generate the projected gate eigenvalues
    est_gate_probabilities_vec =
        (pad_transform_probs * est_unpad_gate_probabilities_vec + pad_mask)
    est_gate_eigenvalues = pad_transform' * wht_transform * est_gate_probabilities_vec
    return (
        est_gate_eigenvalues::Vector{Float64},
        est_unproj_gate_probabilities_vec::Vector{Float64},
        est_gate_probabilities_vec::Vector{Float64},
    )
end

"""
    split_project_gate_eigenvalues(est_unproj_gate_eigenvalues::Vector{Float64}, est_precision_matrix::SparseMatrixCSC{Float64, Int64}, gate_data::GateData; precision::Real = 1e-8, diagnostics::Bool = false)

Returns the projection into the probability simplex of the gate eigenvalues `est_unproj_gate_eigenvalues` and the corresponding gate probabilities, projecting separately for each gate in the Mahalanobis distance using the relevant slice of the precision matrix `est_precision_matrix`, with corresponding gate data `gate_data`.
The solver precision is controlled by `precision`, and solver diagnostics are printed if `diagnostics` is `true`.
"""
function split_project_gate_eigenvalues(
    est_unproj_gate_eigenvalues::Vector{Float64},
    est_precision_matrix::SparseMatrixCSC{Float64, Int64},
    gate_data::GateData;
    precision::Real = 1e-8,
    diagnostics::Bool = false,
)
    # Initialise the transform matrices
    pad_transform = get_pad_transform(gate_data)
    pad_transform_probs = get_pad_transform(gate_data; probabilities = true)
    pad_mask = get_pad_mask(gate_data)
    wht_transform = get_wht_transform(gate_data)
    wht_transform_inv = get_wht_transform(gate_data; inverse = true)
    # Generate the unprojected gate probabilities
    est_unproj_gate_probabilities_vec =
        wht_transform_inv * (pad_transform * est_unproj_gate_eigenvalues + pad_mask)
    est_unproj_unpad_gate_probabilities_vec =
        pad_transform' * est_unproj_gate_probabilities_vec
    # Generate the gate probabilities precision matrix
    probabilities_transform = pad_transform_probs' * wht_transform * pad_transform
    est_unpad_probabilities_precision_matrix =
        probabilities_transform * est_precision_matrix * probabilities_transform'
    # Project the gate probabilities
    est_unpad_gate_probabilities_vec =
        Vector{Float64}(undef, length(est_unproj_unpad_gate_probabilities_vec))
    for gate_idx in gate_data.gate_indices
        indices = gate_idx.indices
        vector_slice = est_unproj_unpad_gate_probabilities_vec[indices]
        precision_matrix_slice = est_unpad_probabilities_precision_matrix[indices, indices]
        proj_vector_slice = scs_project_nonnegative(
            vector_slice,
            precision_matrix_slice;
            precision = precision,
            diagnostics = diagnostics,
        )
        est_unpad_gate_probabilities_vec[indices] = proj_vector_slice
    end
    # Generate the projected gate eigenvalues
    est_gate_probabilities_vec =
        (pad_transform_probs * est_unpad_gate_probabilities_vec + pad_mask)
    est_gate_eigenvalues = pad_transform' * wht_transform * est_gate_probabilities_vec
    return (
        est_gate_eigenvalues::Vector{Float64},
        est_unproj_gate_probabilities_vec::Vector{Float64},
        est_gate_probabilities_vec::Vector{Float64},
    )
end

"""
    simple_project_gate_eigenvalues(est_unproj_gate_eigenvalues::Vector{Float64}, gate_data::GateData)

Returns the projection into the probability simplex of the gate eigenvalues `est_unproj_gate_eigenvalues` and the corresponding gate probabilities, projecting in the Euclidean distance with corresponding gate data `gate_data`.
"""
function simple_project_gate_eigenvalues(
    est_unproj_gate_eigenvalues::Vector{Float64},
    gate_data::GateData,
)
    # Initialise the transform matrices
    pad_transform = get_pad_transform(gate_data)
    pad_mask = get_pad_mask(gate_data)
    wht_transform = get_wht_transform(gate_data)
    wht_transform_inv = get_wht_transform(gate_data; inverse = true)
    # Generate the unprojected gate probabilities
    est_unproj_gate_probabilities_vec =
        wht_transform_inv * (pad_transform * est_unproj_gate_eigenvalues + pad_mask)
    # Project the gate probabilities
    est_gate_probabilities_vec =
        Vector{Float64}(undef, length(est_unproj_gate_probabilities_vec))
    for gate_idx in gate_data.gate_indices
        pad_indices = gate_idx.pad_indices
        vector_slice = est_unproj_gate_probabilities_vec[pad_indices]
        proj_vector_slice = project_simplex(vector_slice)
        est_gate_probabilities_vec[pad_indices] = proj_vector_slice
    end
    # Generate the projected gate eigenvalues
    est_gate_eigenvalues = pad_transform' * wht_transform * est_gate_probabilities_vec
    return (
        est_gate_eigenvalues::Vector{Float64},
        est_unproj_gate_probabilities_vec::Vector{Float64},
        est_gate_probabilities_vec::Vector{Float64},
    )
end

"""
    estimate_gate_noise(d::Design, noise_est::NoiseEstimate; kwargs...)
    estimate_gate_noise(d_rand::RandDesign, noise_est::NoiseEstimate; kwargs...)
    estimate_gate_noise(d::Design, est_eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}}, est_covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}}, shot_budget::Int; kwargs...)
    estimate_gate_noise(d_rand::RandDesign, est_eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}}, est_covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}}; kwargs...)

Returns a [`NoiseEstimate`](@ref) of the gate noise for the design `d` or randomised design `d_rand` derived from estimated circuit eigenvalues, either already contained in a noise estimate `noise_est`, or as estimated circuit eigenvalues `est_eigenvalues_experiment_ensemble` and covariance circuit eigenvalues `est_covariance_experiment_ensemble`, with shot budget `shot_budget`.

# Arguments

  - `d::Design`: Experimental design.
  - `d_rand::RandDesign`: Randomised experimental design.
  - `noise_est::NoiseEstimate`: Noise estimate.
  - `est_eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}}`: Estimated ensemble of circuit eigenvalues.
  - `est_covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}}`: Estimated ensemble of covariance circuit eigenvalues.
  - `shot_budget::Int`: Noise estimate shot budget.

# Keword arguments

  - `min_eigenvalue::Float64 = 0.1`: The minimum eigenvalue threshold for clipping.
  - `clip_warning::Bool = true`: Whether to warn if eigenvalues are clipped.
  - `N_warn_split::Integer = 5 * 10^3`: The number of gate eigenvalues above which to warn about splitting.
  - `split::Bool = (d.c.gate_data.N < N_warn_split ? false : true)`: Whether to split the gate eigenvalue projection across gates or simultaneously project all gate eigenvalues.
  - `split_warning::Bool = true`: Whether to warn if the gate eigenvalue projection is not split.
  - `precision::Real = 1e-8`: The solver precision for the gate eigenvalue projection.
"""
function estimate_gate_noise(
    d::Design,
    est_eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}},
    est_covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}},
    shot_budget::Int;
    min_eigenvalue::Float64 = 0.1,
    clip_warning::Bool = true,
    N_warn_split::Integer = 5 * 10^3,
    split::Bool = (d.c.gate_data.N < N_warn_split ? false : true),
    split_warning::Bool = true,
    precision::Real = 1e-8,
)
    # Initialise variables
    gate_data = d.c.gate_data
    mapping_lengths = length.(d.mapping_ensemble)
    if split_warning && gate_data.N >= N_warn_split && ~split
        @warn "This circuit has a large number of gate eigenvalues: splitting the gate eigenvalue projection across gates by setting `split` to be `true` is advised."
    end
    # Calculate estimation quantities
    tuple_number = length(d.tuple_set)
    est_eigenvalues =
        vcat([mean.(est_eigenvalues_experiment_ensemble[i]) for i in 1:tuple_number]...)
    est_covariance = calc_covariance(
        d,
        est_eigenvalues_experiment_ensemble,
        est_covariance_experiment_ensemble;
        weight_time = false,
        warning = d.full_covariance,
    )
    est_covariance_log = calc_covariance_log(est_covariance, est_eigenvalues)
    # Clip the eigenvalues, covariance matrix, and design matrix if appropriate
    clipped_indices = get_clipped_indices(
        est_eigenvalues;
        min_eigenvalue = min_eigenvalue,
        warning = clip_warning,
    )
    design_matrix = convert(SparseMatrixCSC{Float64, Int}, d.matrix[clipped_indices, :])
    est_eigenvalues = est_eigenvalues[clipped_indices]
    clipped_mapping_lengths = get_clipped_mapping_lengths(mapping_lengths, clipped_indices)
    est_covariance_log = est_covariance_log[clipped_indices, clipped_indices]
    est_diag_covariance_log = sparse(Diagonal(est_covariance_log))
    est_diag_covariance_log_inv_factor =
        sparse(Diagonal(diag(est_diag_covariance_log) .^ (-1 / 2)))
    est_diag_covariance_log_inv =
        est_diag_covariance_log_inv_factor' * est_diag_covariance_log_inv_factor
    # Do OLS and WLS first, and only do GLS if necessary
    # Estimate the gate eigenvalues
    ols_unproj_gate_eigenvalues = estimate_gate_eigenvalues(design_matrix, est_eigenvalues)
    wls_unproj_gate_eigenvalues = estimate_gate_eigenvalues(
        design_matrix,
        est_diag_covariance_log_inv_factor,
        est_eigenvalues,
    )
    # Project the gate error probabilities into the simplex
    (ols_gate_eigenvalues, ols_unproj_gate_probabilities_vec, ols_gate_probabilities_vec) =
        simple_project_gate_eigenvalues(ols_unproj_gate_eigenvalues, gate_data)
    precision_project = ~(ols_gate_eigenvalues ≈ ols_unproj_gate_eigenvalues)
    if precision_project
        wls_precision_matrix = calc_precision_matrix(
            design_matrix,
            wls_unproj_gate_eigenvalues,
            est_diag_covariance_log_inv,
        )
        if split
            (
                wls_gate_eigenvalues,
                wls_unproj_gate_probabilities_vec,
                wls_gate_probabilities_vec,
            ) = split_project_gate_eigenvalues(
                wls_unproj_gate_eigenvalues,
                wls_precision_matrix,
                gate_data;
                precision = precision,
            )
        else
            (
                wls_gate_eigenvalues,
                wls_unproj_gate_probabilities_vec,
                wls_gate_probabilities_vec,
            ) = full_project_gate_eigenvalues(
                wls_unproj_gate_eigenvalues,
                wls_precision_matrix,
                gate_data;
                precision = precision,
            )
        end
    else
        (
            wls_gate_eigenvalues,
            wls_unproj_gate_probabilities_vec,
            wls_gate_probabilities_vec,
        ) = simple_project_gate_eigenvalues(wls_unproj_gate_eigenvalues, gate_data)
    end
    # Get the gate probability dictionaries
    ols_unproj_gate_probabilities =
        get_gate_probabilities_dict(ols_unproj_gate_probabilities_vec, gate_data)
    ols_gate_probabilities =
        get_gate_probabilities_dict(ols_gate_probabilities_vec, gate_data)
    wls_unproj_gate_probabilities =
        get_gate_probabilities_dict(wls_unproj_gate_probabilities_vec, gate_data)
    wls_gate_probabilities =
        get_gate_probabilities_dict(wls_gate_probabilities_vec, gate_data)
    # Get the marginal gate eigenvalues and probabilities
    ols_unproj_marginal_gate_eigenvalues =
        get_marginal_gate_eigenvalues(ols_unproj_gate_eigenvalues, gate_data)
    ols_marginal_gate_eigenvalues =
        get_marginal_gate_eigenvalues(ols_gate_eigenvalues, gate_data)
    ols_unproj_marginal_gate_probabilities =
        get_marginal_gate_probabilities(ols_unproj_gate_probabilities)
    ols_marginal_gate_probabilities =
        get_marginal_gate_probabilities(ols_gate_probabilities)
    wls_unproj_marginal_gate_eigenvalues =
        get_marginal_gate_eigenvalues(wls_unproj_gate_eigenvalues, gate_data)
    wls_marginal_gate_eigenvalues =
        get_marginal_gate_eigenvalues(wls_gate_eigenvalues, gate_data)
    wls_unproj_marginal_gate_probabilities =
        get_marginal_gate_probabilities(wls_unproj_gate_probabilities)
    wls_marginal_gate_probabilities =
        get_marginal_gate_probabilities(wls_gate_probabilities)
    # Get the relative gate eigenvalues
    ols_unproj_relative_gate_eigenvalues =
        get_relative_gate_eigenvalues(ols_unproj_gate_eigenvalues, gate_data)
    ols_relative_gate_eigenvalues =
        get_relative_gate_eigenvalues(ols_gate_eigenvalues, gate_data)
    wls_unproj_relative_gate_eigenvalues =
        get_relative_gate_eigenvalues(wls_unproj_gate_eigenvalues, gate_data)
    wls_relative_gate_eigenvalues =
        get_relative_gate_eigenvalues(wls_gate_eigenvalues, gate_data)
    # Do GLS if necessary
    if d.full_covariance
        # Invert the covariance matrix
        est_covariance_log_inv_factor =
            sparse_covariance_inv_factor(est_covariance_log, clipped_mapping_lengths)
        est_covariance_log_inv =
            est_covariance_log_inv_factor' * est_covariance_log_inv_factor
        # Estimate the gate eigenvalues
        gls_unproj_gate_eigenvalues = estimate_gate_eigenvalues(
            design_matrix,
            est_covariance_log_inv_factor,
            est_eigenvalues,
        )
        # Project the gate error probabilities into the simplex
        if precision_project
            gls_precision_matrix = calc_precision_matrix(
                design_matrix,
                gls_unproj_gate_eigenvalues,
                est_covariance_log_inv,
            )
            if split
                (
                    gls_gate_eigenvalues,
                    gls_unproj_gate_probabilities_vec,
                    gls_gate_probabilities_vec,
                ) = split_project_gate_eigenvalues(
                    gls_unproj_gate_eigenvalues,
                    gls_precision_matrix,
                    gate_data;
                    precision = precision,
                )
            else
                (
                    gls_gate_eigenvalues,
                    gls_unproj_gate_probabilities_vec,
                    gls_gate_probabilities_vec,
                ) = full_project_gate_eigenvalues(
                    gls_unproj_gate_eigenvalues,
                    gls_precision_matrix,
                    gate_data;
                    precision = precision,
                )
            end
        else
            (
                gls_gate_eigenvalues,
                gls_unproj_gate_probabilities_vec,
                gls_gate_probabilities_vec,
            ) = simple_project_gate_eigenvalues(gls_unproj_gate_eigenvalues, gate_data)
        end
        # Get the gate probability dictionaries
        gls_unproj_gate_probabilities =
            get_gate_probabilities_dict(gls_unproj_gate_probabilities_vec, gate_data)
        gls_gate_probabilities =
            get_gate_probabilities_dict(gls_gate_probabilities_vec, gate_data)
        # Get the marginal gate eigenvalues and probabilities
        gls_unproj_marginal_gate_eigenvalues =
            get_marginal_gate_eigenvalues(gls_unproj_gate_eigenvalues, gate_data)
        gls_marginal_gate_eigenvalues =
            get_marginal_gate_eigenvalues(gls_gate_eigenvalues, gate_data)
        gls_unproj_marginal_gate_probabilities =
            get_marginal_gate_probabilities(gls_unproj_gate_probabilities)
        gls_marginal_gate_probabilities =
            get_marginal_gate_probabilities(gls_gate_probabilities)
        # Get the relative gate eigenvalues
        gls_unproj_relative_gate_eigenvalues =
            get_relative_gate_eigenvalues(gls_unproj_gate_eigenvalues, gate_data)
        gls_relative_gate_eigenvalues =
            get_relative_gate_eigenvalues(gls_gate_eigenvalues, gate_data)
    else
        @assert est_covariance_log ≈ est_diag_covariance_log "The covariance matrix is not diagonal as expected."
        # Copy the WLS results
        gls_unproj_gate_eigenvalues = wls_unproj_gate_eigenvalues
        gls_gate_eigenvalues = wls_gate_eigenvalues
        gls_unproj_gate_probabilities_vec = wls_unproj_gate_probabilities_vec
        gls_gate_probabilities_vec = wls_gate_probabilities_vec
        gls_unproj_gate_probabilities = wls_unproj_gate_probabilities
        gls_gate_probabilities = wls_gate_probabilities
        gls_unproj_marginal_gate_eigenvalues = wls_unproj_marginal_gate_eigenvalues
        gls_marginal_gate_eigenvalues = wls_marginal_gate_eigenvalues
        gls_unproj_marginal_gate_probabilities = wls_unproj_marginal_gate_probabilities
        gls_marginal_gate_probabilities = wls_marginal_gate_probabilities
        gls_unproj_relative_gate_eigenvalues = wls_unproj_relative_gate_eigenvalues
        gls_relative_gate_eigenvalues = wls_relative_gate_eigenvalues
    end
    # Return the noise estimate
    noise_est = NoiseEstimate(
        est_eigenvalues_experiment_ensemble,
        est_covariance_experiment_ensemble,
        shot_budget,
        gls_unproj_gate_eigenvalues,
        gls_gate_eigenvalues,
        gls_unproj_gate_probabilities,
        gls_gate_probabilities,
        gls_unproj_marginal_gate_eigenvalues,
        gls_marginal_gate_eigenvalues,
        gls_unproj_marginal_gate_probabilities,
        gls_marginal_gate_probabilities,
        gls_unproj_relative_gate_eigenvalues,
        gls_relative_gate_eigenvalues,
        wls_unproj_gate_eigenvalues,
        wls_gate_eigenvalues,
        wls_unproj_gate_probabilities,
        wls_gate_probabilities,
        wls_unproj_marginal_gate_eigenvalues,
        wls_marginal_gate_eigenvalues,
        wls_unproj_marginal_gate_probabilities,
        wls_marginal_gate_probabilities,
        wls_unproj_relative_gate_eigenvalues,
        wls_relative_gate_eigenvalues,
        ols_unproj_gate_eigenvalues,
        ols_gate_eigenvalues,
        ols_unproj_gate_probabilities,
        ols_gate_probabilities,
        ols_unproj_marginal_gate_eigenvalues,
        ols_marginal_gate_eigenvalues,
        ols_unproj_marginal_gate_probabilities,
        ols_marginal_gate_probabilities,
        ols_unproj_relative_gate_eigenvalues,
        ols_relative_gate_eigenvalues,
    )
    return noise_est::NoiseEstimate
end
function estimate_gate_noise(
    d_rand::RandDesign,
    est_eigenvalues_experiment_ensemble::Vector{Vector{Vector{Float64}}},
    est_covariance_experiment_ensemble::Vector{Vector{Vector{Float64}}};
    min_eigenvalue::Float64 = 0.1,
    clip_warning::Bool = true,
    N_warn_split::Integer = 5 * 10^3,
    split::Bool = (d.c.gate_data.N < N_warn_split ? false : true),
    split_warning::Bool = true,
    precision::Real = 1e-8,
)
    d = get_design(d_rand)
    shot_budget = d_rand.shot_budget
    noise_est = estimate_gate_noise(
        d,
        est_eigenvalues_experiment_ensemble,
        est_covariance_experiment_ensemble,
        shot_budget;
        min_eigenvalue = min_eigenvalue,
        clip_warning = clip_warning,
        N_warn_split = N_warn_split,
        split = split,
        split_warning = split_warning,
        precision = precision,
    )
    return noise_est::NoiseEstimate
end
function estimate_gate_noise(
    d::Design,
    noise_est::NoiseEstimate;
    min_eigenvalue::Float64 = 0.1,
    clip_warning::Bool = true,
    N_warn_split::Integer = 5 * 10^3,
    split::Bool = (d.c.gate_data.N < N_warn_split ? false : true),
    split_warning::Bool = true,
    precision::Real = 1e-8,
)
    noise_est = estimate_gate_noise(
        d,
        noise_est.eigenvalues_experiment_ensemble,
        noise_est.covariance_experiment_ensemble,
        noise_est.shot_budget;
        min_eigenvalue = min_eigenvalue,
        clip_warning = clip_warning,
        N_warn_split = N_warn_split,
        split = split,
        split_warning = split_warning,
        precision = precision,
    )
    return noise_est::NoiseEstimate
end
function estimate_gate_noise(
    d_rand::RandDesign,
    noise_est::NoiseEstimate;
    min_eigenvalue::Float64 = 0.1,
    clip_warning::Bool = true,
    N_warn_split::Integer = 5 * 10^3,
    split::Bool = (d.c.gate_data.N < N_warn_split ? false : true),
    split_warning::Bool = true,
    precision::Real = 1e-8,
)
    d = get_design(d_rand)
    @assert d_rand.shot_budget == noise_est.shot_budget "The shot budgets of the design and noise estimate do not match."
    noise_est = estimate_gate_noise(
        d,
        noise_est;
        min_eigenvalue = min_eigenvalue,
        clip_warning = clip_warning,
        N_warn_split = N_warn_split,
        split = split,
        split_warning = split_warning,
        precision = precision,
    )
    return noise_est::NoiseEstimate
end

"""
    get_eigenvalues(noise_est::NoiseEstimate)

Returns the estimated circuit eigenvalues of noise estimate `noise_est`.
"""
function get_eigenvalues(noise_est::NoiseEstimate)
    est_eigenvalues = vcat(
        [
            mean.(eigenvalues_experiment) for
            eigenvalues_experiment in noise_est.eigenvalues_experiment_ensemble
        ]...,
    )
    return est_eigenvalues::Vector{Float64}
end

"""
    calc_covariance(d::Design, noise_est::NoiseEstimate; weight_time::Bool = false)
    calc_covariance(d_rand::RandDesign, noise_est::NoiseEstimate; weight_time::Bool = false)

Returns the estimated covariance matrix for a noise estimate `noise_est` associated with the design `d` or randomised design `d_rand`, weighting the shot budget by the time factor if `weight_time` is `true`.
"""
function calc_covariance(d::Design, noise_est::NoiseEstimate; weight_time::Bool = false)
    est_covariance =
        calc_covariance(
            d,
            noise_est.eigenvalues_experiment_ensemble,
            noise_est.covariance_experiment_ensemble;
            weight_time = weight_time,
            warning = d.full_covariance,
        ) / noise_est.shot_budget
    return est_covariance::SparseMatrixCSC{Float64, Int}
end
function calc_covariance(
    d_rand::RandDesign,
    noise_est::NoiseEstimate;
    weight_time::Bool = false,
)
    d = get_design(d_rand)
    @assert d_rand.shot_budget == noise_est.shot_budget "The shot budgets of the design and noise estimate do not match."
    est_covariance = calc_covariance(d, noise_est; weight_time = weight_time)
    return est_covariance::SparseMatrixCSC{Float64, Int}
end

"""
    get_noise_error(d::Design, noise_est::NoiseEstimate)
    get_noise_error(d_rand::RandDesign, noise_est::NoiseEstimate)

Returns a [`NoiseError`](@ref) object containing the normalised root-mean-square errors (NRMSEs) for the gate noise estimates in `noise_est` and the true gate eigenvalues for the design `d`, or alternatively the randomised design `d_rand`.
"""
function get_noise_error(d::Design, noise_est::NoiseEstimate)
    # Initialise variables
    times_factor = sum(d.shot_weights .* d.tuple_times)
    meas_budget = times_factor * noise_est.shot_budget
    gate_eigenvalues = get_gate_eigenvalues(d)
    N = d.c.gate_data.N
    ordinary_scaling = sqrt(meas_budget / N)
    marginal_gate_eigenvalues = get_marginal_gate_eigenvalues(d)
    N_marginal = d.c.gate_data.N_marginal
    marginal_scaling = sqrt(meas_budget / N_marginal)
    relative_gate_eigenvalues = get_relative_gate_eigenvalues(d)
    N_relative = d.c.gate_data.N_relative
    relative_scaling = sqrt(meas_budget / N_relative)
    # GLS errors
    gls_nrmse =
        ordinary_scaling * norm(noise_est.gls_unproj_gate_eigenvalues - gate_eigenvalues, 2)
    gls_proj_nrmse =
        ordinary_scaling * norm(noise_est.gls_gate_eigenvalues - gate_eigenvalues, 2)
    gls_marginal_nrmse =
        marginal_scaling *
        norm(noise_est.gls_unproj_marginal_gate_eigenvalues - marginal_gate_eigenvalues, 2)
    gls_proj_marginal_nrmse =
        marginal_scaling *
        norm(noise_est.gls_marginal_gate_eigenvalues - marginal_gate_eigenvalues, 2)
    gls_relative_nrmse =
        relative_scaling *
        norm(noise_est.gls_unproj_relative_gate_eigenvalues - relative_gate_eigenvalues, 2)
    gls_proj_relative_nrmse =
        relative_scaling *
        norm(noise_est.gls_relative_gate_eigenvalues - relative_gate_eigenvalues, 2)
    # WLS errors
    wls_nrmse =
        ordinary_scaling * norm(noise_est.wls_unproj_gate_eigenvalues - gate_eigenvalues, 2)
    wls_proj_nrmse =
        ordinary_scaling * norm(noise_est.wls_gate_eigenvalues - gate_eigenvalues, 2)
    wls_marginal_nrmse =
        marginal_scaling *
        norm(noise_est.wls_unproj_marginal_gate_eigenvalues - marginal_gate_eigenvalues, 2)
    wls_proj_marginal_nrmse =
        marginal_scaling *
        norm(noise_est.wls_marginal_gate_eigenvalues - marginal_gate_eigenvalues, 2)
    wls_relative_nrmse =
        relative_scaling *
        norm(noise_est.wls_unproj_relative_gate_eigenvalues - relative_gate_eigenvalues, 2)
    wls_proj_relative_nrmse =
        relative_scaling *
        norm(noise_est.wls_relative_gate_eigenvalues - relative_gate_eigenvalues, 2)
    # OLS errors
    ols_nrmse =
        ordinary_scaling * norm(noise_est.ols_unproj_gate_eigenvalues - gate_eigenvalues, 2)
    ols_proj_nrmse =
        ordinary_scaling * norm(noise_est.ols_gate_eigenvalues - gate_eigenvalues, 2)
    ols_marginal_nrmse =
        marginal_scaling *
        norm(noise_est.ols_unproj_marginal_gate_eigenvalues - marginal_gate_eigenvalues, 2)
    ols_proj_marginal_nrmse =
        marginal_scaling *
        norm(noise_est.ols_marginal_gate_eigenvalues - marginal_gate_eigenvalues, 2)
    ols_relative_nrmse =
        relative_scaling *
        norm(noise_est.ols_unproj_relative_gate_eigenvalues - relative_gate_eigenvalues, 2)
    ols_proj_relative_nrmse =
        relative_scaling *
        norm(noise_est.ols_relative_gate_eigenvalues - relative_gate_eigenvalues, 2)
    # Get the noise error
    noise_error = NoiseError(
        gls_nrmse,
        gls_proj_nrmse,
        gls_marginal_nrmse,
        gls_proj_marginal_nrmse,
        gls_relative_nrmse,
        gls_proj_relative_nrmse,
        wls_nrmse,
        wls_proj_nrmse,
        wls_marginal_nrmse,
        wls_proj_marginal_nrmse,
        wls_relative_nrmse,
        wls_proj_relative_nrmse,
        ols_nrmse,
        ols_proj_nrmse,
        ols_marginal_nrmse,
        ols_proj_marginal_nrmse,
        ols_relative_nrmse,
        ols_proj_relative_nrmse,
    )
    return noise_error::NoiseError
end
function get_noise_error(d_rand::RandDesign, noise_est::NoiseEstimate)
    d = get_design(d_rand)
    times_factor = sum(d.shot_weights .* d.tuple_times)
    meas_budget = get_meas_budget(d_rand)
    @assert meas_budget ≈ times_factor * noise_est.shot_budget
    noise_error = get_noise_error(d, noise_est)
    return noise_error::NoiseError
end

"""
    get_noise_score(noise_error::NoiseError, merit::Merit)
    get_noise_score(noise_error_vector::Vector{NoiseError}, merit::Merit)
    get_noise_score(noise_error_matrix::Matrix{NoiseError}, merit::Merit)

Returns a [`NoiseScore`](@ref) object containing the z-scores for the supplied normalised root-mean-square error (NRMSE) data in `noise_error` given the merit `merit`.
"""
function get_noise_score(noise_error::NoiseError, merit::Merit)
    # Calculate the z-scores
    gls_z_score = (noise_error.gls_nrmse - merit.gls_expectation) / sqrt(merit.gls_variance)
    gls_proj_z_score =
        (noise_error.gls_proj_nrmse - merit.gls_expectation) / sqrt(merit.gls_variance)
    gls_marginal_z_score =
        (noise_error.gls_marginal_nrmse - merit.gls_marginal_expectation) /
        sqrt(merit.gls_marginal_variance)
    gls_proj_marginal_z_score =
        (noise_error.gls_proj_marginal_nrmse - merit.gls_marginal_expectation) /
        sqrt(merit.gls_marginal_variance)
    gls_relative_z_score =
        (noise_error.gls_relative_nrmse - merit.gls_relative_expectation) /
        sqrt(merit.gls_relative_variance)
    gls_proj_relative_z_score =
        (noise_error.gls_proj_relative_nrmse - merit.gls_relative_expectation) /
        sqrt(merit.gls_relative_variance)
    wls_z_score = (noise_error.wls_nrmse - merit.wls_expectation) / sqrt(merit.wls_variance)
    wls_proj_z_score =
        (noise_error.wls_proj_nrmse - merit.wls_expectation) / sqrt(merit.wls_variance)
    wls_marginal_z_score =
        (noise_error.wls_marginal_nrmse - merit.wls_marginal_expectation) /
        sqrt(merit.wls_marginal_variance)
    wls_proj_marginal_z_score =
        (noise_error.wls_proj_marginal_nrmse - merit.wls_marginal_expectation) /
        sqrt(merit.wls_marginal_variance)
    wls_relative_z_score =
        (noise_error.wls_relative_nrmse - merit.wls_relative_expectation) /
        sqrt(merit.wls_relative_variance)
    wls_proj_relative_z_score =
        (noise_error.wls_proj_relative_nrmse - merit.wls_relative_expectation) /
        sqrt(merit.wls_relative_variance)
    ols_z_score = (noise_error.ols_nrmse - merit.ols_expectation) / sqrt(merit.ols_variance)
    ols_proj_z_score =
        (noise_error.ols_proj_nrmse - merit.ols_expectation) / sqrt(merit.ols_variance)
    ols_marginal_z_score =
        (noise_error.ols_marginal_nrmse - merit.ols_marginal_expectation) /
        sqrt(merit.ols_marginal_variance)
    ols_proj_marginal_z_score =
        (noise_error.ols_proj_marginal_nrmse - merit.ols_marginal_expectation) /
        sqrt(merit.ols_marginal_variance)
    ols_relative_z_score =
        (noise_error.ols_relative_nrmse - merit.ols_relative_expectation) /
        sqrt(merit.ols_relative_variance)
    ols_proj_relative_z_score =
        (noise_error.ols_proj_relative_nrmse - merit.ols_relative_expectation) /
        sqrt(merit.ols_relative_variance)
    # Return the overall scores
    noise_score = NoiseScore(
        gls_z_score,
        gls_proj_z_score,
        gls_marginal_z_score,
        gls_proj_marginal_z_score,
        gls_relative_z_score,
        gls_proj_relative_z_score,
        wls_z_score,
        wls_proj_z_score,
        wls_marginal_z_score,
        wls_proj_marginal_z_score,
        wls_relative_z_score,
        wls_proj_relative_z_score,
        ols_z_score,
        ols_proj_z_score,
        ols_marginal_z_score,
        ols_proj_marginal_z_score,
        ols_relative_z_score,
        ols_proj_relative_z_score,
    )
    return noise_score::NoiseScore
end
function get_noise_score(noise_error_vector::Vector{NoiseError}, merit::Merit)
    noise_score_vector =
        [get_noise_score(noise_error, merit) for noise_error in noise_error_vector]
    return noise_score_vector::Vector{NoiseScore}
end
function get_noise_score(noise_error_matrix::Matrix{NoiseError}, merit::Merit)
    noise_score_matrix =
        [get_noise_score(noise_error, merit) for noise_error in noise_error_matrix]
    return noise_score_matrix::Matrix{NoiseScore}
end

"""
    get_model_score(d::Design, noise_est::NoiseEstimate; projected::Bool = false)
    get_model_score(d_rand::RandDesign, noise_est::NoiseEstimate; projected::Bool = false)
    get_model_score(d::Design, noise_est_vector::Vector{NoiseEstimate}; projected::Bool = false)
    get_model_score(d::Design, noise_est_matrix::Matrix{NoiseEstimate}; projected::Bool = false)

Returns the noise model violation z-score for the generalised residual sum of squares corresponding to the noise estimate `noise_est`, given the design `d` or alternatively the randomised design `d_rand`, calculating with the projected gate eigenvalues if `projected` is `true`.
"""
function get_model_score(d::Design, noise_est::NoiseEstimate; projected::Bool = false)
    # Initialise variables 
    mapping_lengths = length.(d.mapping_ensemble)
    est_eigenvalues = get_eigenvalues(noise_est)
    est_covariance = calc_covariance(d, noise_est)
    est_covariance_log = calc_covariance_log(est_covariance, est_eigenvalues)
    est_covariance_log_inv = sparse_covariance_inv(est_covariance_log, mapping_lengths)
    design_matrix = d.matrix
    (M, N) = size(design_matrix)
    K = M - N
    # Compute the model violation score
    if projected
        gls_log_eigenvalues_residual =
            (design_matrix * log.(noise_est.gls_gate_eigenvalues)) - log.(est_eigenvalues)
    else
        gls_log_eigenvalues_residual =
            (design_matrix * log.(noise_est.gls_unproj_gate_eigenvalues)) -
            log.(est_eigenvalues)
    end
    gls_residual_sum_squares =
        gls_log_eigenvalues_residual' *
        est_covariance_log_inv *
        gls_log_eigenvalues_residual
    model_score = (gls_residual_sum_squares - K) / sqrt(2 * K)
    return model_score::Float64
end
function get_model_score(
    d_rand::RandDesign,
    noise_est::NoiseEstimate;
    projected::Bool = false,
)
    d = get_design(d_rand)
    @assert d_rand.shot_budget == noise_est.shot_budget "The shot budgets of the design and noise estimate do not match."
    model_score = get_model_score(d, noise_est; projected = projected)
    return model_score::Float64
end
function get_model_score(
    d::Design,
    noise_est_vector::Vector{NoiseEstimate};
    projected::Bool = false,
)
    model_score_vector = [
        get_model_score(d, noise_est; projected = projected) for
        noise_est in noise_est_vector
    ]
    return model_score_vector::Vector{Float64}
end
function get_model_score(
    d::Design,
    noise_est_matrix::Matrix{NoiseEstimate};
    projected::Bool = false,
)
    model_score_matrix = [
        get_model_score(d, noise_est; projected = projected) for
        noise_est in noise_est_matrix
    ]
    return model_score_matrix::Matrix{Float64}
end

"""
    is_score_expected(noise_score::NoiseScore, z_score_cutoff::Real)
    is_score_expected(noise_score::NoiseScore, z_score_cutoff_lower::Real, z_score_cutoff_upper::Real)

Returns a Boolean indicating whether the z-scores in `noise_score` are within the specified cutoffs, and projection improves the z-scores as expected.
"""
function is_score_expected(
    noise_score::NoiseScore,
    z_score_cutoff_lower::Real,
    z_score_cutoff_upper::Real,
)
    @assert z_score_cutoff_lower < z_score_cutoff_upper "The lower z-score cutoff must be less than the upper z-score cutoff."
    # Check that the unprojected z-scores are between the cutoffs
    gls_between =
        (noise_score.gls_z_score > z_score_cutoff_lower) &&
        (noise_score.gls_z_score < z_score_cutoff_upper)
    gls_marginal_between =
        (noise_score.gls_marginal_z_score > z_score_cutoff_lower) &&
        (noise_score.gls_marginal_z_score < z_score_cutoff_upper)
    gls_relative_between =
        (noise_score.gls_relative_z_score > z_score_cutoff_lower) &&
        (noise_score.gls_relative_z_score < z_score_cutoff_upper)
    wls_between =
        (noise_score.wls_z_score > z_score_cutoff_lower) &&
        (noise_score.wls_z_score < z_score_cutoff_upper)
    wls_marginal_between =
        (noise_score.wls_marginal_z_score > z_score_cutoff_lower) &&
        (noise_score.wls_marginal_z_score < z_score_cutoff_upper)
    wls_relative_between =
        (noise_score.wls_relative_z_score > z_score_cutoff_lower) &&
        (noise_score.wls_relative_z_score < z_score_cutoff_upper)
    ols_between =
        (noise_score.ols_z_score > z_score_cutoff_lower) &&
        (noise_score.ols_z_score < z_score_cutoff_upper)
    ols_marginal_between =
        (noise_score.ols_marginal_z_score > z_score_cutoff_lower) &&
        (noise_score.ols_marginal_z_score < z_score_cutoff_upper)
    ols_relative_between =
        (noise_score.ols_relative_z_score > z_score_cutoff_lower) &&
        (noise_score.ols_relative_z_score < z_score_cutoff_upper)
    # Check that the projected z-scores are as expected
    gls_improved = noise_score.gls_proj_z_score <= noise_score.gls_z_score
    gls_marginal_under = noise_score.gls_proj_marginal_z_score < z_score_cutoff_upper
    gls_relative_under = noise_score.gls_proj_relative_z_score < z_score_cutoff_upper
    wls_improved = noise_score.wls_proj_z_score <= noise_score.wls_z_score
    wls_marginal_under = noise_score.wls_proj_marginal_z_score < z_score_cutoff_upper
    wls_relative_under = noise_score.wls_proj_relative_z_score < z_score_cutoff_upper
    ols_improved = noise_score.ols_proj_z_score <= noise_score.ols_z_score
    # Aggregate the results
    is_between =
        (gls_between && gls_marginal_between && gls_relative_between) &&
        (wls_between && wls_marginal_between && wls_relative_between) &&
        (ols_between && ols_marginal_between && ols_relative_between)
    is_projection_expected =
        (gls_improved && gls_marginal_under && gls_relative_under) &&
        (wls_improved && wls_marginal_under && wls_relative_under) &&
        ols_improved
    is_expected = is_between && is_projection_expected
    return is_expected::Bool
end
function is_score_expected(noise_score::NoiseScore, z_score_cutoff::Real)
    @assert z_score_cutoff > 0.0 "The z-score cutoff must be positive."
    is_expected = is_score_expected(noise_score, -z_score_cutoff, z_score_cutoff)
    return is_expected::Bool
end

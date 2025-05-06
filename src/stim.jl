"""
    MemoryData

Stim memory simulation results.

# Fields

  - `circuit_param::AbstractCircuitParameters`: Circuit parameters.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.
  - `rounds::Int`: Number of memory circuit rounds.
  - `shots::Int`: Number of sampled shots.
  - `seed::UInt64`: Seed for Stim simulation.
  - `reset_type::Symbol`: Reset type, which is `:meas_reset` by default but can be set to `:meas`.
  - `decoder_type::Symbol`: Decoder type, which is `:pymatching` by default but can be set to `:beliefmatching`.
  - `decoder_gate_probabilities::Vector{Dict{Gate, Vector{Float64}}}`: Vector of gate probabilities used to inform the decoder.
  - `decoder_labels::Vector{String}`: Labels for the gate probabilities in `decoder_gate_probabilities`.
  - `memory_x_observable::Vector{Bool}`: X memory experiment observable across each of the shots.
  - `memory_z_observable::Vector{Bool}`: Z memory experiment observable across each of the shots.
  - `decoder_x_predictions::Matrix{Bool}`: X memory experiment decoder predictions across each of the shots and decoders in `decoder_gate_probabilities`.
  - `decoder_z_predictions::Matrix{Bool}`: Z memory experiment decoder predictions across each of the shots and decoders in `decoder_gate_probabilities`.
"""
struct MemoryData
    circuit_param::AbstractCircuitParameters
    noise_param::AbstractNoiseParameters
    rounds::Int
    shots::Int
    seed::UInt64
    reset_type::Symbol
    decoder_type::Symbol
    decoder_gate_probabilities::Vector{Dict{Gate, Vector{Float64}}}
    decoder_labels::Vector{String}
    memory_x_observable::Vector{Bool}
    memory_z_observable::Vector{Bool}
    decoder_x_predictions::Matrix{Bool}
    decoder_z_predictions::Matrix{Bool}
end

function Base.show(io::IO, memory_data::MemoryData)
    # Process the memory data
    shots = memory_data.shots
    decoder_num = length(memory_data.decoder_labels)
    decoder_x_errors = Matrix{Bool}(undef, shots, decoder_num)
    decoder_z_errors = Matrix{Bool}(undef, shots, decoder_num)
    for i in 1:decoder_num
        decoder_x_errors[:, i] =
            memory_data.memory_x_observable .!= memory_data.decoder_x_predictions[:, i]
        decoder_z_errors[:, i] =
            memory_data.memory_z_observable .!= memory_data.decoder_z_predictions[:, i]
    end
    # Calculate means and variances
    memory_x_errors_sum = vec(sum(decoder_x_errors; dims = 1))
    memory_z_errors_sum = vec(sum(decoder_z_errors; dims = 1))
    memory_x_errors = memory_x_errors_sum / shots
    memory_z_errors = memory_z_errors_sum / shots
    memory_x_errors_cov = cov(decoder_x_errors)
    memory_z_errors_cov = cov(decoder_z_errors)
    memory_x_errors_sem = [sqrt(memory_x_errors_cov[i, i] / shots) for i in 1:decoder_num]
    memory_z_errors_sem = [sqrt(memory_z_errors_cov[i, i] / shots) for i in 1:decoder_num]
    # Print the memory summary data
    digits = 4
    memory_x_errors_string = [
        "($(round(100 * memory_x_errors[i], digits = digits)) ± $(round(100 * memory_x_errors_sem[i], digits = digits)))%"
        for i in 1:decoder_num
    ]
    memory_z_errors_string = [
        "($(round(100 * memory_z_errors[i], digits = digits)) ± $(round(100 * memory_z_errors_sem[i], digits = digits)))%"
        for i in 1:decoder_num
    ]
    pretty_table(
        io,
        hcat(memory_data.decoder_labels, memory_x_errors_string, memory_z_errors_string);
        header = ["Decoder label"; "Memory X error rate"; "Memory Z error rate"],
        alignment = :l,
    )
    return nothing
end

@struct_hash_equal_isequal MemoryData

"""
    MemorySummary

Stim memory simulation summary results.

# Fields

  - `circuit_param::AbstractCircuitParameters`: Circuit parameters.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.
  - `rounds::Int`: Number of memory circuit rounds.
  - `shots::Int`: Number of sampled shots.
  - `seed::UInt64`: Seed for Stim simulation.
  - `reset_type::Symbol`: Reset type, which is `:meas_reset` by default but can be set to `:meas`.
  - `decoder_type::Symbol`: Decoder type, which is `:pymatching` by default but can be set to `:beliefmatching`.
  - `decoder_labels::Vector{String}`: Labels for the decoder gate probabilities.
  - `memory_x_errors::Vector{Float64}`: X memory experiment error rate.
  - `memory_z_errors::Vector{Float64}`: Z memory experiment error rate.
  - `memory_x_errors_cov::Matrix{Float64}`: X memory experiment error rate covariance.
  - `memory_z_errors_cov::Matrix{Float64}`: Z memory experiment error rate covariance.
  - `change_x_errors::Matrix{Float64}`: X memory experiment change in error rate across different gate probabilities.
  - `change_z_errors::Matrix{Float64}`: Z memory experiment change in error rate across different gate probabilities.
  - `change_x_errors_sem::Matrix{Float64}`: X memory experiment change in error rate across different gate probabilities standard error of the mean.
  - `change_z_errors_sem::Matrix{Float64}`: Z memory experiment change in error rate across different gate probabilities standard error of the mean.
  - `decoder_x_confusion::Matrix{Int}`: X memory experiment decoder confusion matrix.
  - `decoder_z_confusion::Matrix{Int}`: Z memory experiment decoder confusion matrix.
"""
struct MemorySummary
    circuit_param::AbstractCircuitParameters
    noise_param::AbstractNoiseParameters
    rounds::Int
    shots::Int
    seed::UInt64
    reset_type::Symbol
    decoder_type::Symbol
    decoder_labels::Vector{String}
    memory_x_errors::Vector{Float64}
    memory_z_errors::Vector{Float64}
    memory_x_errors_cov::Matrix{Float64}
    memory_z_errors_cov::Matrix{Float64}
    change_x_errors::Matrix{Float64}
    change_z_errors::Matrix{Float64}
    change_x_errors_sem::Matrix{Float64}
    change_z_errors_sem::Matrix{Float64}
    decoder_x_confusion::Matrix{Int}
    decoder_z_confusion::Matrix{Int}
end

function Base.show(io::IO, memory_summary::MemorySummary)
    # Print the memory summary data
    digits = 4
    shots = memory_summary.shots
    decoder_num = length(memory_summary.decoder_labels)
    memory_x_errors_sem =
        [sqrt(memory_summary.memory_x_errors_cov[i, i] / shots) for i in 1:decoder_num]
    memory_z_errors_sem =
        [sqrt(memory_summary.memory_z_errors_cov[i, i] / shots) for i in 1:decoder_num]
    memory_x_errors_string = [
        "($(round(100 * memory_summary.memory_x_errors[i], digits = digits)) ± $(round(100 * memory_z_errors_sem[i], digits = digits)))%"
        for i in 1:decoder_num
    ]
    memory_z_errors_string = [
        "($(round(100 * memory_summary.memory_z_errors[i], digits = digits)) ± $(round(100 * memory_x_errors_sem[i], digits = digits)))%"
        for i in 1:decoder_num
    ]
    pretty_table(
        io,
        hcat(memory_summary.decoder_labels, memory_x_errors_string, memory_z_errors_string);
        header = ["Decoder label"; "Memory X error rate"; "Memory Z error rate"],
        alignment = :l,
    )
    return nothing
end

@struct_hash_equal_isequal MemorySummary

"""
    get_stim_qubits(qubits::Vector{Tuple{Int, Int}})
    get_stim_qubits(code_param::CodeParameters)
    get_stim_qubits(empty_code_param::EmptyCodeParameters)

Returns a Stim string representation of the qubit coordinates described by `qubits`, which may be a field of `code_param`.
"""
function get_stim_qubits(qubits::Vector{Tuple{Int, Int}})
    # Generate the qubit coordinates
    qubit_list = Vector{String}(undef, length(qubits))
    for (idx, (i, j)) in pairs(qubits)
        qubit_list[idx] = "QUBIT_COORDS($(j), $(i)) $(idx)\n"
    end
    stim_qubits = join(qubit_list)
    return stim_qubits::String
end
function get_stim_qubits(code_param::CodeParameters)
    stim_qubits = get_stim_qubits(code_param.qubits)
    return stim_qubits::String
end

"""
    get_stim_circuit(circuit::Vector{Layer}, gate_probabilities::Dict{Gate, Vector{Float64}}, noisy_prep::Bool, noisy_meas::Bool; extra_fields::Dict{Symbol, Any} = Dict{Symbol, Any}(), meas_and_reset::Bool = false)

Returns a Stim string representation of the circuit `circuit` alongside error probabilities specified by `gate_probabilities`, as well as noisy preparations if `noisy_prep` is `true` and noisy measurements if `noisy_meas` is `true`.
Qubit coordinates are specified by `extra_fields` if it contains a [`CodeParameters`](@ref) object.
Reset types are specified by `reset_type` and can be `:reset`, `:meas`, or `:meas_reset`.

Stim orders one-qubit Paulis as: X, Y, Z.
We order one-qubit Paulis as: X, Z, Y.
The indexing to transform between these orderings is: 1, 3, 2.

Stim orders two-qubit Paulis as: IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ.
We order two-qubit Paulis as: XI, IX, XX, ZI, YI, ZX, YX, IZ, XZ, IY, XY, ZZ, YZ, ZY, YY.
The indexing to transform from the Stim ordering to ours is: 4, 1, 5, 12, 8, 13, 9, 3, 7, 2, 6, 15, 11, 14, 10.
The indexing to transform from our ordering to the Stim ordering is: 2, 10, 8, 1, 3, 11, 9, 5, 7, 15, 13, 4, 6, 14, 12.
"""
function get_stim_circuit(
    circuit::Vector{Layer},
    gate_probabilities::Dict{Gate, Vector{Float64}},
    noisy_prep::Bool,
    noisy_meas::Bool;
    extra_fields::Dict{Symbol, Any} = Dict{Symbol, Any}(),
    reset_type::Symbol = :reset,
)
    # Check variables
    @assert reset_type ∈ [:reset; :meas; :meas_reset] "The reset type must be either `:reset`, `:meas`, or `:meas_reset`."
    # Stim orders its one-qubit Paulis as
    # X, Y, Z
    # We order our one-qubit Paulis as
    # X, Z, Y
    # To go from the Stim ordering to our ordering, and vice versa, we have the indexing
    # 1, 3, 2.
    # We add 1 to ignore the normalising probability of I.
    order_1 = [1; 3; 2] .+ 1
    # Stim orders its two-qubit Paulis as
    # IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ.
    # We order our two-qubit Paulis as
    # XI, IX, XX, ZI, YI, ZX, YX, IZ, XZ, IY, XY, ZZ, YZ, ZY, YY.
    # To go from the Stim ordering to our ordering, we have the indexing
    # 4, 1, 5, 12, 8, 13, 9, 3, 7, 2, 6, 15, 11, 14, 10.
    # To go from our ordering to the Stim ordering, we have the indexing
    # 2, 10, 8, 1, 3, 11, 9, 5, 7, 15, 13, 4, 6, 14, 12.
    # Here we need the latter ordering.
    # We add 1 to ignore the normalising probability of II.
    order_2 = [2; 10; 8; 1; 3; 11; 9; 5; 7; 15; 13; 4; 6; 14; 12] .+ 1
    # Measurement and preparation are both only associated with a single error probability
    # This is the second index, as the first index is the probability of no error.
    # Generate the gate list for Stim
    string_vector = Vector{String}(undef, length(circuit))
    for (i, l) in pairs(circuit)
        k = length(l.layer)
        noise_string_vector = Vector{String}(undef, k)
        gate_string_vector = Vector{String}(undef, k)
        meas_noise_string_vector = Vector{String}(undef, k)
        for (j, gate) in pairs(l.layer)
            noise_string_vector[j] = ""
            meas_noise_string_vector[j] = ""
            if is_one_qubit(gate)
                # Generate the Stim single-qubit gate string
                if is_mid_meas(gate)
                    gate_string_vector[j] = "$(gate.type)($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                    meas_noise_string_vector[j] = "X_ERROR($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                elseif is_mid_reset(gate)
                    if reset_type == :reset
                        gate_string_vector[j] = "$(gate.type) $(gate.targets[1])\n"
                    elseif reset_type == :meas || reset_type == :meas_reset
                        # Get the measurement error from the relevant SPAM error if possible
                        spam_meas_gate = Gate("MZ", 0, [gate.targets[1]])
                        spam_prep_gate = Gate("PZ", 0, [gate.targets[1]])
                        if haskey(gate_probabilities, spam_meas_gate)
                            spam_gate = spam_meas_gate
                        elseif haskey(gate_probabilities, spam_prep_gate)
                            spam_gate = spam_prep_gate
                        else
                            spam_gate = gate
                        end
                        gate_string_vector[j] = "$(reset_type == :meas ? "M" : "MR")($(gate_probabilities[spam_gate][2])) $(gate.targets[1])\n"
                    else
                        throw(error("Unsupported reset type $(reset_type)."))
                    end
                    meas_noise_string_vector[j] = "X_ERROR($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                else
                    noise_string_vector[j] = "PAULI_CHANNEL_1($(join(gate_probabilities[gate][order_1], ", "))) $(gate.targets[1])\n"
                    gate_string_vector[j] = "$((is_meas_idle(gate) ? "I" : gate.type)) $(gate.targets[1])\n"
                end
            elseif is_two_qubit(gate; stim_supported = true)
                # Generate the Stim two-qubit gate string
                noise_string_vector[j] = "PAULI_CHANNEL_2($(join(gate_probabilities[gate][order_2], ", "))) $(join(gate.targets, " "))\n"
                gate_string_vector[j] = "$(gate.type) $(join(gate.targets, " "))\n"
            elseif is_state_meas(gate)
                # Generate the Stim measurement gate string
                if noisy_meas
                    gate_string_vector[j] = "$(gate.type)($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                else
                    gate_string_vector[j] = "$(gate.type) $(gate.targets[1])\n"
                end
            elseif is_state_prep(gate)
                # Generate the Stim preparation gate string
                if noisy_prep
                    noise_string_vector[j] = "X_ERROR($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                end
                if gate.type == "PZ"
                    gate_string_vector[j] = ""
                elseif gate.type == "PX"
                    gate_string_vector[j] = "H $(gate.targets[1])\n"
                elseif gate.type == "PY"
                    gate_string_vector[j] = "H $(gate.targets[1])\nS $(gate.targets[1])\n"
                else
                    throw(error("Invalid preparation gate type $(gate.type)."))
                end
            else
                throw(error("Invalid gate $(gate)."))
            end
        end
        string_vector[i] =
            join(noise_string_vector) *
            join(gate_string_vector) *
            join(meas_noise_string_vector) *
            "TICK\n"
    end
    # Generate the qubit coordinates
    if haskey(extra_fields, :code_param)
        code_param = extra_fields[:code_param]
        insert!(string_vector, 1, get_stim_qubits(code_param))
    end
    # Generate the Stim circuit string
    stim_circuit = join(string_vector)
    return stim_circuit::String
end

"""
    get_stim_initialise(code_param::CodeParameters, memory_type::Symbol, gate_probabilities::Dict{Gate, Vector{Float64}}, noisy_prep::Bool)

Returns a Stim string representation of the initialisation for the code parameters `code_param` corresponding to the Pauli type specified by `memory_type`, which can be either `:x` or `:z`, using the gate error probabilities `gate_probabilities` if `noisy_prep` is `true`.
"""
function get_stim_initialise(
    code_param::CodeParameters,
    memory_type::Symbol,
    gate_probabilities::Dict{Gate, Vector{Float64}},
    noisy_prep::Bool,
)
    # Initialise variables
    @assert memory_type ∈ [:x; :z] "The memory type must be either `:x` or `:z`."
    data_indices = code_param.data_indices
    if memory_type == :x
        init_indices = code_param.init_x_indices
    elseif memory_type == :z
        init_indices = code_param.init_z_indices
    else
        throw(error("Unsupported memory type $(memory_type)."))
    end
    # Calculate the initialisation with noise if appropriate
    noise_list = String[]
    init_list = String[]
    for data_idx in data_indices
        if data_idx ∈ init_indices
            gate = Gate("PX", 0, [data_idx])
            push!(init_list, "H $(data_idx)\n")
        else
            gate = Gate("PZ", 0, [data_idx])
        end
        if noisy_prep
            prep_noise = "X_ERROR($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
            push!(noise_list, prep_noise)
        end
    end
    stim_initialise = join(noise_list) * join(init_list)
    if length(init_list) > 0
        stim_initialise = stim_initialise * "TICK\n"
    end
    return stim_initialise::String
end

"""
    get_stim_detectors(code_param::CodeParameters, memory_type::Symbol; previous_compare::Integer = 0, include_data::Bool = false)
    get_stim_detectors(code_param::CodeParameters, previous_compare::Integer)

Returns a Stim string representation of the detectors for the code parameters `code_param` corresponding to the Pauli type specified by `memory_type`, which can be either `:x` or `:z`, comparing the detectors against those from `previous_compare` rounds ago.
If `include_data` is true, the detector checks the relevant data qubits and is included only if it checks a nonzero number of data qubits.

If no memory type is specified, this generates detectors for both X and Z memory types but requires `previous_compare` to be specified.
"""
function get_stim_detectors(
    code_param::CodeParameters,
    memory_type::Symbol;
    previous_compare::Integer = 0,
    include_data::Bool = false,
)
    # Initialise variables
    @assert memory_type ∈ [:x; :z] "The memory type must be either `:x` or `:z`."
    @assert previous_compare ∈ [0; 1; 2] "The number of previous rounds against which to compare the detectors must be 0, 1, or 2."
    qubits = code_param.qubits
    data_indices = code_param.data_indices
    data_num = length(data_indices)
    ancilla_indices = code_param.ancilla_indices
    ancilla_num = length(ancilla_indices)
    if memory_type == :x
        check_indices = code_param.check_x_indices
    elseif memory_type == :z
        check_indices = code_param.check_z_indices
    else
        throw(error("Unsupported memory type $(memory_type)."))
    end
    if include_data
        @assert previous_compare != 0 "If the detector includes data qubits, the number of previous rounds against which to compare the detector cannot be 0."
        rec_offset = data_num + ancilla_num + 1
        prev_rec_offset = data_num + previous_compare * ancilla_num + 1
        data_rec_offset = data_num + 1
    else
        rec_offset = ancilla_num + 1
        prev_rec_offset = (1 + previous_compare) * ancilla_num + 1
    end
    # Construct the detectors
    detectors_list = String[]
    for check_idx_pair in check_indices
        # Get indices
        det_indices = [
            findfirst(ancilla_indices .== ancilla_idx) for ancilla_idx in check_idx_pair[1]
        ]
        meas_indices =
            [findfirst(data_indices .== data_idx) for data_idx in check_idx_pair[2]]
        # Construct the detector
        detector = "DETECTOR($(qubits[check_idx_pair[1][1]][2]), $(qubits[check_idx_pair[1][1]][1]), 0)"
        for det_idx in det_indices
            detector *= " rec[-$(rec_offset - det_idx)]"
        end
        if (include_data && previous_compare == 2) ||
           (~include_data && previous_compare ∈ [1; 2])
            for det_idx in det_indices
                detector *= " rec[-$(prev_rec_offset - det_idx)]"
            end
        end
        # Include data indices if appropriate, and then add the detector
        if include_data
            for meas_idx in meas_indices
                detector *= " rec[-$(data_rec_offset - meas_idx)]"
            end
            detector *= "\n"
            push!(detectors_list, detector)
        else
            detector *= "\n"
            push!(detectors_list, detector)
        end
    end
    stim_detectors = join(detectors_list)
    return stim_detectors::String
end
function get_stim_detectors(code_param::CodeParameters, previous_compare::Integer)
    # Construct the detectors
    stim_detectors_x =
        get_stim_detectors(code_param, :x; previous_compare = previous_compare)
    stim_detectors_z =
        get_stim_detectors(code_param, :z; previous_compare = previous_compare)
    stim_detectors = "SHIFT_COORDS(0, 0, 1)\n" * stim_detectors_x * stim_detectors_z
    return stim_detectors::String
end

"""
    get_stim_measure_detectors(code_param::CodeParameters, gate_probabilities::Dict{Gate, Vector{Float64}}, noisy_meas::Bool; reset_type::Symbol = :meas_reset, do_rounds::Bool = true)

Returns a Stim string representation of the measurement and detectors for the code parameters `code_param` corresponding to the Pauli type specified by `memory_type`, which can be either `:x` or `:z`, using the gate error probabilities `gate_probabilities` if `noisy_meas` is `true`, for the reset type `reset_type`, which can be either `:meas_reset` or `:meas`.
Set `do_rounds` to `false` when no rounds of the memory circuit are performed and `reset_type` is `:meas` to avoid breaking the measurement detectors.
"""
function get_stim_measure_detectors(
    code_param::CodeParameters,
    memory_type::Symbol,
    gate_probabilities::Dict{Gate, Vector{Float64}},
    noisy_meas::Bool;
    reset_type::Symbol = :meas_reset,
    do_rounds::Bool = true,
)
    # Initialise variables
    @assert memory_type ∈ [:x; :z] "The memory type must be either `:x` or `:z`."
    @assert reset_type ∈ [:meas_reset; :meas] "The reset type must be either `:meas_reset` or `:meas`."
    if reset_type == :meas_reset
        previous_compare = 1
    elseif reset_type == :meas
        if do_rounds
            previous_compare = 2
        else
            previous_compare = 1
        end
    else
        throw(error("Unsupported reset type $(reset_type)."))
    end
    data_indices = code_param.data_indices
    data_num = length(data_indices)
    if memory_type == :x
        check_indices = code_param.check_x_indices
        init_indices = code_param.init_x_indices
        logical_indices = code_param.logical_x_indices
    elseif memory_type == :z
        check_indices = code_param.check_z_indices
        init_indices = code_param.init_z_indices
        logical_indices = code_param.logical_z_indices
    else
        throw(error("Unsupported memory type $(memory_type)."))
    end
    # Measure the data qubits
    measure_list = Vector{String}(undef, data_num)
    for (idx, data_idx) in pairs(data_indices)
        if data_idx ∈ init_indices
            gate_type = "MX"
        else
            gate_type = "MZ"
        end
        gate = Gate(gate_type, 0, [data_idx])
        if noisy_meas
            measure_list[idx] = "$(gate.type)($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
        else
            measure_list[idx] = "$(gate.type) $(gate.targets[1])\n"
        end
    end
    stim_measure = join(measure_list)
    # Construct the detectors
    stim_detectors =
        "SHIFT_COORDS(0, 0, 1)\n" * get_stim_detectors(
            code_param,
            memory_type;
            previous_compare = previous_compare,
            include_data = true,
        )
    # Construct the observable
    data_rec_offset = data_num + 1
    stim_observable =
        "OBSERVABLE_INCLUDE(0)" *
        join(" rec[-$(data_rec_offset - logical_idx)]" for logical_idx in logical_indices) *
        "\n"
    # Construct the measurement and detector string
    stim_measure_detectors = stim_measure * stim_detectors * stim_observable
    return stim_measure_detectors::String
end

"""
    get_stim_memory_circuit(c::AbstractCircuit, gate_probabilities::Dict{Gate, Vector{Float64}}, memory_type::Symbol, rounds::Integer; reset_type::Symbol = :meas_reset)
    get_stim_memory_circuit(c::AbstractCircuit, memory_type::Symbol, rounds::Integer; reset_type::Symbol = :meas_reset)

Returns a Stim string representation of the memory circuit corresponding to the syndrome extraction circuit `c`, optionally using the gate error probabilities `gate_probabilities`, memory type `memory_type`, which can be either `:x` or `:z`, number of rounds `rounds`, and reset type `reset_type`, which can be either `:meas_reset` or `:meas`.

If you want to add support for custom syndrome extraction circuits `c` which do not store their code parameters in the [`CodeParameters`](@ref) object, you will need to define new methods for `get_stim_qubits`, `get_stim_initialise`, `get_stim_detectors`, and `get_stim_measure_detectors` to make this function work.
"""
function get_stim_memory_circuit(
    c::T,
    gate_probabilities::Dict{Gate, Vector{Float64}},
    memory_type::Symbol,
    rounds::Integer;
    reset_type::Symbol = :meas_reset,
) where {T <: AbstractCircuit}
    # Check variables
    @assert haskey(c.extra_fields, :code_param) "The circuit lacks a `code_param` field."
    code_param = c.extra_fields[:code_param]
    @assert typeof(code_param) == CodeParameters "The code parameters must be of type `CodeParameters`."
    circuit_params = c.circuit_param.params
    @assert haskey(circuit_params, :ancilla_measurement) &&
            circuit_params[:ancilla_measurement] "The circuit parameters lack an `ancilla_measurement` field or it is not set to `true`."
    @assert memory_type ∈ [:x; :z] "The memory type must be either `:x` or `:z`."
    @assert rounds >= 0 "The number of rounds must be at least 0."
    @assert reset_type ∈ [:meas_reset; :meas] "The reset type must be either `:meas_reset` or `:meas`."
    if rounds > 0
        do_rounds = true
    else
        do_rounds = false
    end
    if reset_type == :meas_reset
        previous_compare = 1
    elseif reset_type == :meas
        previous_compare = 2
    else
        throw(error("Unsupported reset type $(reset_type)."))
    end
    uncompared = 0
    # Initialise Stim strings
    stim_qubits = get_stim_qubits(code_param)
    stim_initialise_pauli =
        get_stim_initialise(code_param, memory_type, gate_probabilities, c.noisy_prep)
    stim_detectors_pauli = get_stim_detectors(code_param, memory_type)
    stim_detectors_uncompared = get_stim_detectors(code_param, uncompared)
    stim_detectors = get_stim_detectors(code_param, previous_compare)
    stim_circuit = get_stim_circuit(
        c.circuit,
        gate_probabilities,
        c.noisy_prep,
        c.noisy_meas;
        reset_type = reset_type,
    )
    stim_measure_pauli = get_stim_measure_detectors(
        code_param,
        memory_type,
        gate_probabilities,
        c.noisy_meas;
        reset_type = reset_type,
        do_rounds = do_rounds,
    )
    # Generate the memory circuit
    memory_circuit =
        stim_qubits * stim_initialise_pauli * stim_circuit * stim_detectors_pauli
    if reset_type == :meas_reset
        if rounds >= 1
            memory_circuit =
                memory_circuit *
                "REPEAT $(rounds) {\n" *
                stim_circuit *
                stim_detectors *
                "}\n"
        end
    elseif reset_type == :meas
        if rounds >= 1
            memory_circuit = memory_circuit * stim_circuit * stim_detectors_uncompared
        end
        if rounds >= 2
            memory_circuit =
                memory_circuit *
                "REPEAT $(rounds - 1) {\n" *
                stim_circuit *
                stim_detectors *
                "}\n"
        end
    else
        throw(error("Unsupported reset type $(reset_type)."))
    end
    memory_circuit = memory_circuit * stim_measure_pauli
    return memory_circuit::String
end
function get_stim_memory_circuit(
    c::T,
    memory_type::Symbol,
    rounds::Integer;
    reset_type::Symbol = :meas_reset,
) where {T <: AbstractCircuit}
    memory_circuit = get_stim_memory_circuit(
        c,
        c.gate_probabilities,
        memory_type,
        rounds;
        reset_type = reset_type,
    )
    return memory_circuit::String
end

"""
    batch_simulate_memory(c::AbstractCircuit, memory_type::Symbol, rounds::Integer, shots::Integer; kwargs...)

Returns the memory and decoding simulation results for a memory experiment conducted with the given syndrome extraction circuit.

WARNING: Seeding has the same features as in Stim.
The behaviour of the same random seed will differ across different versions of Stim.
Also, when measurement shots are sampled in batches, which occurs when `max_samples` is exceeded, the results will differ from when all shots are sampled at once.

# Arguments

  - `c::AbstractCircuit`: Syndrome extraction circuit.
  - `memory_type::Symbol`: Memory type, which can be either `:x` or `:z`.
  - `rounds::Integer`: Number of memory circuit rounds, which must be at least 0.
  - `shots::Integer`: Number of shots, which must be at least 1.

# Keyword arguments

  - `seed::Union{UInt64, Nothing} = nothing`: Seed to use for random number generation.
  - `reset_type::Symbol = :meas_reset`: Reset type, which can be either `:meas_reset` or `:meas`.
  - `decoder_type::Symbol = :pymatching`: Decoder type, which can be either `:pymatching` or `:beliefmatching`.
  - `decoder_gate_probabilities::Vector{Dict{Gate, Vector{Float64}}} = [c.gate_probabilities]`: Gate probabilities used for decoding, the first of which must be the gate probabilities of the supplied circuit `c`.
  - `max_samples::Integer = 10^9`: Maximum number of Stim detector samples collected in a single simulation.
  - `diagnostics::Bool = false`: Whether to print diagnostics.
"""
function batch_simulate_memory(
    c::T,
    memory_type::Symbol,
    rounds::Integer,
    shots::Integer;
    seed::Union{UInt64, Nothing} = nothing,
    reset_type::Symbol = :meas_reset,
    decoder_type::Symbol = :pymatching,
    decoder_gate_probabilities::Vector{Dict{Gate, Vector{Float64}}} = [
        c.gate_probabilities,
    ],
    max_samples::Integer = 10^9,
    diagnostics::Bool = false,
) where {T <: AbstractCircuit}
    # Check variables
    @assert haskey(c.extra_fields, :code_param) "The circuit lacks a `code_param` field."
    code_param = c.extra_fields[:code_param]
    @assert typeof(code_param) == CodeParameters "The code parameters must be of type `CodeParameters`."
    circuit_params = c.circuit_param.params
    @assert haskey(circuit_params, :ancilla_measurement) &&
            circuit_params[:ancilla_measurement] "The circuit parameters lack an `ancilla_measurement` field or it is not set to `true`."
    @assert memory_type ∈ [:x; :z] "The memory type must be either `:x` or `:z`."
    @assert rounds >= 0 "The number of rounds must be at least 0."
    @assert shots >= 1 "The number of shots must be at least 1."
    @assert reset_type ∈ [:meas_reset; :meas] "The reset type must be either `:meas_reset` or `:meas`."
    @assert decoder_type ∈ [:pymatching; :beliefmatching] "The decoder type must be either `:pymatching` or `:beliefmatching`."
    @assert decoder_gate_probabilities[1] == c.gate_probabilities "The first decoder gate probabilities must be the same as the circuit gate probabilities."
    @assert all(
        sort(collect(keys(gate_probabilities))) == sort(c.total_gates) for
        gate_probabilities in decoder_gate_probabilities
    ) "All of the decoder gate probabilities must be consistent with the circuit."
    # Generate the memory circuit
    memory_circuit =
        get_stim_memory_circuit(c, memory_type, rounds; reset_type = reset_type)
    stim_memory_circuit = stim.Circuit(memory_circuit)
    # Generate the decoders
    decoder_num = length(decoder_gate_probabilities)
    detector_count = 0
    decoders = Vector{Py}(undef, decoder_num)
    for (idx, gate_probabilities) in pairs(decoder_gate_probabilities)
        decoder_circuit = get_stim_memory_circuit(
            c,
            gate_probabilities,
            memory_type,
            rounds;
            reset_type = reset_type,
        )
        decoder_dem = stim.Circuit(decoder_circuit).detector_error_model(;
            decompose_errors = true,
            approximate_disjoint_errors = true,
        )
        if idx == 1
            detector_count = pyconvert(Int, decoder_dem.num_detectors)
        else
            @assert detector_count == pyconvert(Int, decoder_dem.num_detectors) "The number of detectors in the decoders does not match the number of detectors in the first decoder."
        end
        if decoder_type == :pymatching
            decoder = Matching(decoder_dem)
        elseif decoder_type == :beliefmatching
            try
                decoder = BeliefMatching(decoder_dem)
            catch
                @warn "Failed to load the BeliefMatching decoder. Falling back to PyMatching."
                decoder = Matching(decoder_dem)
            end
        else
            throw(error("Unsupported decoder type $(decoder_type)."))
        end
        decoders[idx] = decoder
    end
    # Generate the shot batches and seeds
    shot_batches = batch_shots(shots, detector_count, max_samples)
    batch_num = length(shot_batches)
    if seed !== nothing
        Random.seed!(seed)
    end
    batch_seeds = rand(UInt64, batch_num)
    if seed !== nothing
        Random.seed!()
    end
    # Sample and decode the memory circuit batches
    start_time = time()
    if diagnostics
        println(
            "Simulating memory type $(memory_type), dividing the shots into $(batch_num) batches.",
        )
    end
    memory_observable = Vector{Bool}()
    decoder_predictions = Matrix{Bool}(undef, 0, decoder_num)
    for batch_idx in 1:batch_num
        # Initialise batch parameters
        batch_seed = batch_seeds[batch_idx]
        batch_shots = shot_batches[batch_idx]
        # Sample the memory circuits
        sample_start = time()
        memory_sampler = stim_memory_circuit.compile_detector_sampler(; seed = batch_seed)
        (batch_memory_detection, batch_observable) =
            memory_sampler.sample(batch_shots; separate_observables = true)
        append!(memory_observable, vec(pyconvert(Matrix{Bool}, batch_observable)))
        if diagnostics
            println(
                "Sampled batch $(batch_idx) of $(batch_num) in $(round(time() - sample_start, digits = 2)) s. The time elapsed since simulation started is $(round(time() - start_time, digits = 2)) s.",
            )
        end
        # Decode the memory circuits
        batch_decoder_predictions = Matrix{Bool}(undef, batch_shots, decoder_num)
        for (idx, decoder) in pairs(decoders)
            decode_start = time()
            batch_decoder_prediction =
                vec(pyconvert(Matrix{Bool}, decoder.decode_batch(batch_memory_detection)))
            batch_decoder_predictions[:, idx] = batch_decoder_prediction
            if diagnostics
                println(
                    "Decoded batch $(batch_idx) with decoder $(idx) in $(round(time() - decode_start, digits = 2)) s.",
                )
            end
        end
        decoder_predictions = vcat(decoder_predictions, batch_decoder_predictions)
        if batch_num > 1
            GC.gc()
        end
    end
    if diagnostics
        println(
            "Simulated memory type $(memory_type) in $(round(time() - start_time, digits = 2)) s.",
        )
    end
    @assert length(memory_observable) == shots "The number of memory observables does not match the number of batches."
    @assert size(decoder_predictions, 1) == shots "The number of decoder predictions does not match the number of shots."
    return (memory_observable::Vector{Bool}, decoder_predictions::Matrix{Bool})
end

"""
    simulate_memory(c::AbstractCircuit, rounds::Integer, shots::Integer; kwargs...)

Returns the memory and decoding simulation results for a memory experiment conducted with the given syndrome extraction circuit.

WARNING: Seeding has the same features as in Stim.
The behaviour of the same random seed will differ across different versions of Stim.
Also, when measurement shots are sampled in batches, which occurs when `max_samples` is exceeded, the results will differ from when all shots are sampled at once.

# Arguments

  - `c::AbstractCircuit`: Syndrome extraction circuit.
  - `rounds::Integer`: Number of memory circuit rounds, which must be at least 0.
  - `shots::Integer`: Number of shots, which must be at least 1.

# Keyword arguments

  - `seed::Union{UInt64, Nothing} = nothing`: Seed to use for random number generation.
  - `reset_type::Symbol = :meas_reset`: Reset type, which can be either `:meas_reset` or `:meas`.
  - `decoder_type::Symbol = :pymatching`: Decoder type, which can be either `:pymatching` or `:beliefmatching`.
  - `decoder_gate_probabilities::Vector{Dict{Gate, Vector{Float64}}} = [c.gate_probabilities]`: Gate probabilities used for decoding, the first of which must be the gate probabilities of the supplied circuit `c`.
  - `decoder_labels::Vector{String} = [string(idx) for idx in 1:length(decoder_gate_probabilities)]`: Labels for the decoder gate probabilities.
  - `max_samples::Integer = 10^9`: Maximum number of Stim detector samples collected in a single simulation.
  - `diagnostics::Bool = false`: Whether to print diagnostics.
"""
function simulate_memory(
    c::T,
    rounds::Integer,
    shots::Integer;
    seed::Union{UInt64, Nothing} = nothing,
    reset_type::Symbol = :meas_reset,
    decoder_type::Symbol = :pymatching,
    decoder_gate_probabilities::Vector{Dict{Gate, Vector{Float64}}} = [
        c.gate_probabilities,
    ],
    decoder_labels::Vector{String} = [
        "Decoder " * string(idx) for idx in 1:length(decoder_gate_probabilities)
    ],
    max_samples::Integer = 10^9,
    diagnostics::Bool = false,
) where {T <: AbstractCircuit}
    # Initialise the random seeds
    if seed !== nothing
        Random.seed!(seed)
    end
    (seed_x, seed_z) = rand(UInt64, 2)
    if seed !== nothing
        Random.seed!()
    end
    # Simulate the memory circuits
    (memory_x_observable, decoder_x_predictions) = batch_simulate_memory(
        c,
        :x,
        rounds,
        shots;
        seed = seed_x,
        reset_type = reset_type,
        decoder_type = decoder_type,
        decoder_gate_probabilities = decoder_gate_probabilities,
        max_samples = max_samples,
        diagnostics = diagnostics,
    )
    (memory_z_observable, decoder_z_predictions) = batch_simulate_memory(
        c,
        :z,
        rounds,
        shots;
        seed = seed_z,
        reset_type = reset_type,
        decoder_type = decoder_type,
        decoder_gate_probabilities = decoder_gate_probabilities,
        max_samples = max_samples,
        diagnostics = diagnostics,
    )
    # Return the memory data
    memory_data = MemoryData(
        c.circuit_param,
        c.noise_param,
        rounds,
        shots,
        seed,
        reset_type,
        decoder_type,
        decoder_gate_probabilities,
        decoder_labels,
        memory_x_observable,
        memory_z_observable,
        decoder_x_predictions,
        decoder_z_predictions,
    )
    return memory_data::MemoryData
end

"""
    process_memory_data(memory_observable::Vector{Bool}, decoder_predictions::Matrix{Bool})

Returns summary data for the memory observables `memory_observable` and decoder predictions `decoder_predictions`.
"""
function process_memory_data(
    memory_observable::Vector{Bool},
    decoder_predictions::Matrix{Bool},
)
    # Initialise variables
    (shots, decoder_num) = size(decoder_predictions)
    @assert length(memory_observable) == shots "The number of memory observables does not match the number of shots."
    # Process the memory data
    decoder_errors = Matrix{Bool}(undef, shots, decoder_num)
    for i in 1:decoder_num
        decoder_errors[:, i] = memory_observable .!= decoder_predictions[:, i]
    end
    # Calculate means and variances
    memory_errors_sum = vec(sum(decoder_errors; dims = 1))
    memory_errors = memory_errors_sum / shots
    memory_errors_cov = cov(decoder_errors)
    # Calculate the change in the logical error rate from using different decoder gate probabilities
    change_errors = [
        memory_errors_sum[i] / memory_errors_sum[j] - 1 for i in 1:decoder_num,
        j in 1:decoder_num
    ]
    change_errors_sem = [
        (memory_errors_sum[i] / memory_errors_sum[j]) * sqrt(
            shots * (
                memory_errors_cov[i, i] / memory_errors_sum[i]^2 +
                memory_errors_cov[j, j] / memory_errors_sum[j]^2 -
                2 * memory_errors_cov[i, j] / (memory_errors_sum[i] * memory_errors_sum[j])
            ),
        ) for i in 1:decoder_num, j in 1:decoder_num
    ]
    # Calculate the decoder confusion matrix
    decoder_confusion = zeros(Int, decoder_num, decoder_num)
    for i in 1:decoder_num
        for j in 1:decoder_num
            if j == i
                decoder_confusion[i, i] = sum(decoder_errors[:, i])
            else
                decoder_confusion[i, j] =
                    sum((.~decoder_errors[:, i]) .& decoder_errors[:, j])
            end
        end
    end
    # Return the summary data
    return (
        memory_errors::Vector{Float64},
        memory_errors_cov::Matrix{Float64},
        change_errors::Matrix{Float64},
        change_errors_sem::Matrix{Float64},
        decoder_confusion::Matrix{Int},
    )
end

"""
    get_memory_summary(memory_data::MemoryData)

Returns memory summary data as a [`MemorySummary`](@ref) object from the memory data `memory_data`.
"""
function get_memory_summary(memory_data::MemoryData)
    # Process the memory data
    (
        memory_x_errors,
        memory_x_errors_cov,
        change_x_errors,
        change_x_errors_sem,
        decoder_x_confusion,
    ) = process_memory_data(
        memory_data.memory_x_observable,
        memory_data.decoder_x_predictions,
    )
    (
        memory_z_errors,
        memory_z_errors_cov,
        change_z_errors,
        change_z_errors_sem,
        decoder_z_confusion,
    ) = process_memory_data(
        memory_data.memory_z_observable,
        memory_data.decoder_z_predictions,
    )
    # Store the memory data
    memory_summary = MemorySummary(
        memory_data.circuit_param,
        memory_data.noise_param,
        memory_data.rounds,
        memory_data.shots,
        memory_data.seed,
        memory_data.reset_type,
        memory_data.decoder_type,
        memory_data.decoder_labels,
        memory_x_errors,
        memory_z_errors,
        memory_x_errors_cov,
        memory_z_errors_cov,
        change_x_errors,
        change_z_errors,
        change_x_errors_sem,
        change_z_errors_sem,
        decoder_x_confusion,
        decoder_z_confusion,
    )
    return memory_summary::MemorySummary
end

"""
    round_exponential_model(rounds, params)

Exponential decay model as a function of the number of rounds `rounds` with parameters `params`.
"""
function round_exponential_model(rounds, params)
    return 0.5 * (1 .- params[2] * exp.(-params[1] .* rounds))
end

"""
    fit_round_error(rounds_set::Vector{Int}, rounds_memory_errors::Vector{Float64}, rounds_memory_errors_sem::Vector{Float64}; return_cov::Bool = false)

Returns parameters for [`round_exponential_model`](@ref) fitted with weighted least squares from the memory errors `rounds_memory_errors` and their standard errors of the mean `rounds_memory_errors_sem` across the rounds in `rounds_set`.
"""
function fit_round_error(
    rounds_set::Vector{Int},
    rounds_memory_errors::Vector{Float64},
    rounds_memory_errors_sem::Vector{Float64};
    return_cov::Bool = false,
)
    # Fit the error scaling with rounds with weighted least squares
    rounds_fit = LsqFit.curve_fit(
        round_exponential_model,
        rounds_set,
        rounds_memory_errors,
        rounds_memory_errors_sem .^ (-2),
        [1.0; 1.0],
    )
    round_params = rounds_fit.param
    if return_cov
        round_params_cov = LsqFit.estimate_covar(rounds_fit)
        return (round_params::Vector{Float64}, round_params_cov::Matrix{Float64})
    else
        return round_params::Vector{Float64}
    end
end

"""
    get_round_error(round_params::Vector{Float64})
    get_round_error(round_params::Vector{Float64}, round_params_cov::Matrix{Float64})

Returns the error per round determined from `round_params`, and its standard error if its covariance matrix `round_params_cov` is supplied.
"""
function get_round_error(round_params::Vector{Float64})
    round_error = 0.5 * (1 - exp(-round_params[1]))
    return round_error::Float64
end
function get_round_error(round_params::Vector{Float64}, round_params_cov::Matrix{Float64})
    round_error = get_round_error(round_params)
    round_error_se = 0.5 * exp(-round_params[1]) * sqrt(round_params_cov[1, 1])
    return (round_error::Float64, round_error_se::Float64)
end

"""
    dist_linear_model(dist, params)

Linear model as a function of the distance `dist` with parameters `params`.
"""
function dist_linear_model(dist, params)
    return params[1] * dist .+ params[2]
end

"""
    fit_dist_error(dist_set::Vector{Int}, dist_memory_errors::Vector{Float64}, dist_memory_errors_sem::Vector{Float64}; return_cov::Bool = false)

Returns parameters for [`dist_linear_model`](@ref) fitted with weighted least squares from the memory errors `dist_memory_errors` and their standard errors of the mean `dist_memory_errors_sem` across the distances in `dist_set`.
Note that the fitting is performed in log space: to fit the original parameters, exponentiate the output of [`dist_linear_model`](@ref).
"""
function fit_dist_error(
    dist_set::Vector{Int},
    dist_memory_errors::Vector{Float64},
    dist_memory_errors_sem::Vector{Float64};
    return_cov::Bool = false,
)
    # Fit in log-space the error scaling with distance with weighted least squares
    log_dist_memory_errors = log.(dist_memory_errors)
    log_dist_memory_errors_sem = dist_memory_errors_sem ./ dist_memory_errors
    log_dist_fit = LsqFit.curve_fit(
        dist_linear_model,
        dist_set,
        log_dist_memory_errors,
        log_dist_memory_errors_sem .^ (-2),
        [1.0; 0.0],
    )
    dist_params = log_dist_fit.param
    if return_cov
        dist_params_cov = LsqFit.estimate_covar(log_dist_fit)
        return (dist_params::Vector{Float64}, dist_params_cov::Matrix{Float64})
    else
        return dist_params::Vector{Float64}
    end
end

"""
    get_dist_error(dist_params::Vector{Float64})
    get_dist_error(dist_params::Vector{Float64}, dist_params_cov::Matrix{Float64})

Returns the error suppression factor, the change in error rate in a code as the distance is increased by 2, determined from `dist_params`, and its standard error if its covariance matrix `dist_params_cov` is supplied.
"""
function get_dist_error(dist_params::Vector{Float64})
    error_suppression = exp(-2 * dist_params[1])
    return error_suppression::Float64
end
function get_dist_error(dist_params::Vector{Float64}, dist_params_cov::Matrix{Float64})
    error_suppression = get_dist_error(dist_params)
    error_suppression_se = 2 * error_suppression * sqrt(dist_params_cov[1, 1])
    return (error_suppression::Float64, error_suppression_se::Float64)
end

"""
    calc_memory_distances(c::AbstractCircuit; reset_type::Symbol = :meas_reset, max_event_set_size::Integer = 4, max_edge_degree::Integer = 2, explore_increasing_degree::Bool = false)

Calculates the logical Z and X error distances corresponding to the X and Z memory circuits corresponding to the syndrome extraction circuit `c`, which are conventionally the vertical and horizontal distances, respectively.
Reset types are specified by `reset_type` and can be `:meas_reset` or `:meas`.
Search parameters used by Stim are `max_event_set_size`, `max_edge_degree`, and `explore_increasing_degree`.
"""
function calc_memory_distances(
    c::T;
    reset_type::Symbol = :meas_reset,
    max_event_set_size::Integer = 4,
    max_edge_degree::Integer = 2,
    explore_increasing_degree::Bool = false,
) where {T <: AbstractCircuit}
    # Check variables
    @assert haskey(c.extra_fields, :code_param) "The circuit lacks a `code_param` field."
    code_param = c.extra_fields[:code_param]
    @assert typeof(code_param) == CodeParameters "The code parameters must be of type `CodeParameters`."
    circuit_params = c.circuit_param.params
    @assert haskey(circuit_params, :vertical_dist) "The circuit parameters lack a `vertical_dist` field."
    @assert haskey(circuit_params, :horizontal_dist) "The circuit parameters lack  a `horizontal_dist` field."
    @assert haskey(circuit_params, :ancilla_measurement) &&
            circuit_params[:ancilla_measurement] "The circuit parameters lack an `ancilla_measurement` field or it is not set to `true`."
    @assert reset_type ∈ [:meas_reset; :meas] "The reset type must be either `:meas_reset` or `:meas`."
    # Generate the memory circuits
    rounds = max(circuit_params[:vertical_dist], circuit_params[:horizontal_dist])
    memory_x_circuit = get_stim_memory_circuit(c, :x, rounds; reset_type = reset_type)
    memory_z_circuit = get_stim_memory_circuit(c, :z, rounds; reset_type = reset_type)
    # Calculate the logical error distances
    z_dist = length(
        stim.Circuit(memory_x_circuit).search_for_undetectable_logical_errors(;
            dont_explore_detection_event_sets_with_size_above = max_event_set_size,
            dont_explore_edges_with_degree_above = max_edge_degree,
            dont_explore_edges_increasing_symptom_degree = ~explore_increasing_degree,
        ),
    )
    x_dist = length(
        stim.Circuit(memory_z_circuit).search_for_undetectable_logical_errors(;
            dont_explore_detection_event_sets_with_size_above = max_event_set_size,
            dont_explore_edges_with_degree_above = max_edge_degree,
            dont_explore_edges_increasing_symptom_degree = ~explore_increasing_degree,
        ),
    )
    return (z_dist::Int, x_dist::Int)
end

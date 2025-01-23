"""
    Circuit

Circuit information, including noise parameters.

# Fields

  - `circuit_param::AbstractCircuitParameters`: Circuit parameters.
  - `circuit::Vector{Layer}`: Circuit.
  - `circuit_tuple::Vector{Int}`: Tuple which arranges the order of the circuit layers; this is initialised as trivial.
  - `qubit_num::Int`: Number of qubits in the circuit.
  - `unique_layer_indices::Vector{Int}`: Unique non-measurement gate layer indices of the circuit.
  - `layer_types::Vector{Symbol}`: Types of the layers in the circuit, used for layer times and dynamical decoupling.
  - `layer_times::Vector{Float64}`: Times taken to implement each layer in the circuit, as well as measurement and reset.
  - `gates::Vector{Gate}`: Gates in the circuit arranged by the tuple.
  - `total_gates::Vector{Gate}`: Total gates in the circuit, including preparations if `noisy_prep` and measurements if `noisy_meas`.
  - `gate_data::GateData`: Gate data for the circuit.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.
  - `gate_probabilities::Dict{Gate, Vector{Float64}}`: Pauli error probabilities for each gate, stored as a dictionary.
  - `gate_eigenvalues::Vector{Float64}`: Eigenvalues for each gate, stored as a vector whose order is determined by the gate order in `total_gates` and indexed by `gate_data`.
  - `noisy_prep::Bool`: Whether to treat preparations as noisy and characterise the associated noise, which should default to `false`; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `noisy_meas::Bool`: Whether to treat measurements as noisy and characterise the associated noise, which should default to `true`; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `extra_fields::Dict{Symbol, Any}`: Extra data for the circuit, including code parameters for syndrome extraction circuits stored as a `:code_param` field which is a [`CodeParameters`](@ref) object.
"""
struct Circuit <: AbstractCircuit
    circuit_param::AbstractCircuitParameters
    circuit::Vector{Layer}
    circuit_tuple::Vector{Int}
    qubit_num::Int
    unique_layer_indices::Vector{Int}
    layer_types::Vector{Symbol}
    layer_times::Vector{Float64}
    gates::Vector{Gate}
    total_gates::Vector{Gate}
    gate_data::GateData
    noise_param::AbstractNoiseParameters
    gate_probabilities::Dict{Gate, Vector{Float64}}
    gate_eigenvalues::Vector{Float64}
    noisy_prep::Bool
    noisy_meas::Bool
    extra_fields::Dict{Symbol, Any}
end

function Base.show(io::IO, c::Circuit)
    println(
        io,
        "Circuit of type $(c.circuit_param.circuit_name) with $(c.qubit_num) qubits.",
    )
    if haskey(c.extra_fields, :code_param)
        code_param = c.extra_fields[:code_param]
        show(io, code_param)
    end
    return nothing
end

@struct_hash_equal_isequal Circuit

"""
    EmptyCircuitParameters

Empty circuit parameters.
"""
struct EmptyCircuitParameters <: AbstractCircuitParameters
    params::Dict{Symbol, Any}
    circuit_name::String
    # Constructor
    function EmptyCircuitParameters()
        params = Dict{Symbol, Any}()
        circuit_name = "circuit"
        return new(params, circuit_name)::EmptyCircuitParameters
    end
end

@struct_hash_equal_isequal EmptyCircuitParameters

"""
    EmptyNoiseParameters

Empty noise parameters.
"""
struct EmptyNoiseParameters <: AbstractNoiseParameters
    params::Dict{Symbol, Any}
    noise_name::String
    # Constructor
    function EmptyNoiseParameters()
        params = Dict{Symbol, Any}()
        noise_name = "noise"
        return new(params, noise_name)::EmptyNoiseParameters
    end
end

@struct_hash_equal_isequal EmptyNoiseParameters

"""
    CodeParameters

Code parameters for syndrome extraction circuits.

Note that these expect the qubits to be arranged such that the data qubits are first, followed by the ancilla qubits, and these parameters should be checked with [`check_code_parameters`](@ref).

# Fields

  - `qubits::Vector{Tuple{Int, Int}}`: Code qubit lattice locations.
  - `qubit_layout::Matrix{String}`: Diagram of the layout of the code qubits.
  - `inverse_indices::Dict{Tuple{Int, Int}, Int}`: Inverse mapping from the qubit lattice locations to their indices.
  - `data_indices::Vector{Int}`: Data qubit indices.
  - `ancilla_indices::Vector{Int}`: Ancilla qubit indices.
  - `ancilla_x_indices::Vector{Int}`: Ancilla X check qubit indices.
  - `ancilla_z_indices::Vector{Int}`: Ancilla Z check qubit indices.
  - `check_x_indices::Vector{Tuple{Vector{Int}, Vector{Int}}}`: Pairs of (ancilla, data) qubit indices for each of the X checks.
  - `check_z_indices::Vector{Tuple{Vector{Int}, Vector{Int}}}`: Pairs of (ancilla, data) qubit indices for each of the Z checks.
  - `init_x_indices::Vector{Int}`: Logical X initialisation qubit indices for which to initialise in the X basis.
  - `init_z_indices::Vector{Int}`: Logical Z initialisation qubit indices for which to initialise in the X basis.
  - `logical_x_indices::Vector{Int}`: Logical X operator qubit indices.
  - `logical_z_indices::Vector{Int}`: Logical Z operator qubit indices.
"""
struct CodeParameters <: AbstractCodeParameters
    qubits::Vector{Tuple{Int, Int}}
    qubit_layout::Matrix{String}
    inverse_indices::Dict{Tuple{Int, Int}, Int}
    data_indices::Vector{Int}
    ancilla_indices::Vector{Int}
    ancilla_x_indices::Vector{Int}
    ancilla_z_indices::Vector{Int}
    check_x_indices::Vector{Tuple{Vector{Int}, Vector{Int}}}
    check_z_indices::Vector{Tuple{Vector{Int}, Vector{Int}}}
    init_x_indices::Vector{Int}
    init_z_indices::Vector{Int}
    logical_x_indices::Vector{Int}
    logical_z_indices::Vector{Int}
end

function Base.show(io::IO, code_param::CodeParameters)
    show(io, MIME("text/plain"), code_param.qubit_layout)
    return nothing
end

@struct_hash_equal_isequal CodeParameters

"""
    apply_tuple(c::AbstractCircuit, circuit_tuple::Vector{Int})

Returns a copy of the circuit `c` arranged by the tuple `circuit_tuple`.
"""
function apply_tuple(c::T, circuit_tuple::Vector{Int}) where {T <: AbstractCircuit}
    tuple_circuit = c.circuit[circuit_tuple]
    tuple_gates = get_gates(tuple_circuit)
    # Update the parameters
    # Applying a tuple to the circuit makes the unique layer and meas indices meaningless
    c_tuple = deepcopy(c)
    @reset c_tuple.circuit_tuple = circuit_tuple
    @reset c_tuple.gates = tuple_gates
    return c_tuple::T
end

"""
    update_noise(c::AbstractCircuit, noise_param::AbstractNoiseParameters)
    update_noise(c::AbstractCircuit, gate_probabilities::Dict{Gate, Vector{Float64}})

Returns a copy of `c` where the circuit has been updated with noise generated according to `noise_param`, or the gate probabilities `gate_probabilities`.
"""
function update_noise(
    c::T,
    noise_param::U,
) where {T <: AbstractCircuit, U <: AbstractNoiseParameters}
    # Generate the noise
    @assert ~(c.noisy_prep && c.noisy_meas) "Preparation and measurement noise cannot be characterised simultaneously."
    c_update = deepcopy(c)
    gate_probabilities = init_gate_probabilities(c.total_gates, noise_param)
    if haskey(noise_param.params, :combined)
        if noise_param.params[:combined]
            if c.gate_data.combined
                gate_data = c.gate_data
            else
                gate_data = get_gate_data(
                    c.total_gates;
                    combined = true,
                    strict = c.gate_data.strict,
                )
            end
            combined_gate_probabilities =
                get_combined_gate_probabilities(gate_probabilities, gate_data)
            @assert all(
                gate_probabilities[gate] ≈ combined_gate_probabilities[gate] for
                gate in c.total_gates
            ) "The gate probabilities are not consistent with the combined gate probabilities."
        else
            if c.gate_data.combined
                gate_data = get_gate_data(
                    c.total_gates;
                    combined = false,
                    strict = c.gate_data.strict,
                )
            else
                gate_data = c.gate_data
            end
        end
    elseif c.gate_data.combined
        gate_data = c.gate_data
        gate_probabilities = get_combined_gate_probabilities(gate_probabilities, gate_data)
    else
        gate_data = c.gate_data
    end
    gate_eigenvalues = get_gate_eigenvalues(gate_probabilities, gate_data)
    # Update the circuit
    @reset c_update.noise_param = noise_param
    @reset c_update.gate_data = gate_data
    @reset c_update.gate_probabilities = gate_probabilities
    @reset c_update.gate_eigenvalues = gate_eigenvalues
    return c_update::T
end
function update_noise(
    c::T,
    gate_probabilities::Dict{Gate, Vector{Float64}},
) where {T <: AbstractCircuit}
    # Generate the noise
    @assert ~(c.noisy_prep && c.noisy_meas) "Preparation and measurement noise cannot be characterised simultaneously."
    @assert sort(c.total_gates) == sort(collect(keys(gate_probabilities))) "The gate probabilities must correspond to the gates in the circuit."
    gate_data = c.gate_data
    c_update = deepcopy(c)
    if c.gate_data.combined
        gate_probabilities = get_combined_gate_probabilities(gate_probabilities, gate_data)
    end
    gate_eigenvalues = get_gate_eigenvalues(gate_probabilities, gate_data)
    # Update the circuit
    @reset c_update.noise_param = EmptyNoiseParameters()
    @reset c_update.gate_probabilities = gate_probabilities
    @reset c_update.gate_eigenvalues = gate_eigenvalues
    return c_update::T
end

"""
    get_layer_times(layer_types::Vector{Symbol}, layer_time_dict::Dict{Symbol, Float64})

Returns the times taken to implement each layer in the circuit based on their types in `layer_types` and the times specified in `layer_time_dict, including the time for final measurement and reset.
"""
function get_layer_times(
    layer_types::Vector{Symbol},
    layer_time_dict::Dict{Symbol, Float64},
)
    # Append the layer times
    layer_times = Vector{Float64}(undef, length(layer_types) + 1)
    for (idx, type) in pairs(layer_types)
        if haskey(layer_time_dict, type)
            layer_times[idx] = layer_time_dict[type]
        else
            throw(error("Unsupported layer type $(type)."))
        end
    end
    if haskey(layer_time_dict, :meas_reset)
        layer_times[end] = layer_time_dict[:meas_reset]
    else
        throw(error("Unsupported layer time $(:meas_reset)."))
    end
    @assert all(layer_times .> 0.0) "The layer times must be positive."
    return layer_times::Vector{Float64}
end

"""
    unwrap_circuit(circuit::Vector{Layer})

Returns a copy of the circuit `circuit` where each layer has been unwrapped into a vector of gates.
"""
function unwrap_circuit(circuit::Vector{Layer})
    @assert all(l.qubit_num == circuit[1].qubit_num for l in circuit) "All layers in the circuit must act on the same number of qubits."
    unwrapped_circuit = [[gate for gate in l.layer] for l in circuit]
    return unwrapped_circuit::Vector{Vector{Gate}}
end

"""
    get_gates(circuit::Vector{Vector{Gate}})
    get_gates(circuit::Vector{Layer})

Returns the unique gates in the circuit `circuit`.
"""
function get_gates(circuit::Vector{Vector{Gate}})
    gates = unique([gate for layer in circuit for gate in layer])
    return gates::Vector{Gate}
end
function get_gates(circuit::Vector{Layer})
    @assert all(l.qubit_num == circuit[1].qubit_num for l in circuit) "All layers in the circuit must act on the same number of qubits."
    gates = unique([gate for l in circuit for gate in l.layer])
    return gates::Vector{Gate}
end

"""
    label_circuit(circuit::Vector{Vector{Gate}})
    label_circuit(circuit::Vector{Layer})

Returns a labelled copy of the circuit `circuit`, with indices indicating the unique layer appearances of each gate in the circuit.
This function should only be applied to the original circuit before it is arranged by a tuple.
"""
function label_circuit(circuit::Vector{Vector{Gate}})
    # Remove any previous labels
    unlabelled_circuit =
        [sort!([Gate(gate.type, 0, gate.targets) for gate in layer]) for layer in circuit]
    labelled_circuit = deepcopy(unlabelled_circuit)
    # Determine the gates
    gates = get_gates(unlabelled_circuit)
    # Construct a version of the circuit pruned of all non-unique layers
    L = length(circuit)
    unique_circuit = Vector{Vector{Gate}}()
    unique_layer_indices = Vector{Int}()
    for l in 1:L
        unlabelled_layer = unlabelled_circuit[l]
        if unlabelled_layer ∉ unique_circuit
            push!(unique_circuit, unlabelled_layer)
            push!(unique_layer_indices, l)
        end
    end
    # For each of the unique layers in the circuit, determine which layers they correspond to in the original circuit
    U = length(unique_circuit)
    layer_indices = Vector{Vector{Int}}(undef, U)
    for u in 1:U
        indices = Vector{Int}()
        for l in 1:L
            if unique_circuit[u] == unlabelled_circuit[l]
                push!(indices, l)
            end
        end
        layer_indices[u] = indices
    end
    # Find the unique layers in which each gate appears and relabel if appropriate
    for gate in gates
        appears_index = findall(gate ∈ unique_circuit[u] for u in 1:U)
        for a in eachindex(appears_index)
            for i in layer_indices[appears_index[a]]
                for j in eachindex(labelled_circuit[i])
                    if gate == labelled_circuit[i][j]
                        # Append the appearance number to the gate label vector
                        labelled_circuit[i][j] = Gate(gate.type, a, gate.targets)
                    end
                end
            end
        end
    end
    return (labelled_circuit::Vector{Vector{Gate}}, unique_layer_indices::Vector{Int})
end
function label_circuit(circuit::Vector{Layer})
    n = circuit[1].qubit_num
    @assert all(l.qubit_num == n for l in circuit) "All layers in the circuit must act on the same number of qubits."
    unwrapped_circuit = unwrap_circuit(circuit)
    (unwrapped_labelled_circuit, unique_layer_indices) = label_circuit(unwrapped_circuit)
    labelled_circuit = [Layer(layer, n) for layer in unwrapped_labelled_circuit]
    return (labelled_circuit::Vector{Layer}, unique_layer_indices::Vector{Int})
end

"""
    get_total_gates(circuit::Vector{Layer}, noisy_prep::Bool, noisy_meas::Bool)

Returns the total gates in the circuit `circuit`, including preparations if `noisy_prep` and measurements if `noisy_meas` each on all of the qubits, as well as a separate list of the added measurement gates.
"""
function get_total_gates(circuit::Vector{Layer}, noisy_prep::Bool, noisy_meas::Bool)
    # Get the gates
    n = convert(Int, circuit[1].qubit_num)
    @assert all(l.qubit_num == n for l in circuit) "All layers in the circuit must act on the same number of qubits."
    total_gates = get_gates(circuit)
    # Append preparations to the gate list if appropriate
    if noisy_prep
        for i in 1:n
            append!(
                total_gates,
                [
                    Gate("PZ", 0, [i])
                    Gate("PX", 0, [i])
                    Gate("PY", 0, [i])
                ],
            )
        end
    end
    # Append measurements to the gate list if appropriate
    if noisy_meas
        for i in 1:n
            append!(
                total_gates,
                [
                    Gate("MZ", 0, [i])
                    Gate("MX", 0, [i])
                    Gate("MY", 0, [i])
                ],
            )
        end
    end
    @assert total_gates == unique(total_gates) "The total gates must be unique."
    return total_gates::Vector{Gate}
end

"""
    check_circuit(circuit::Vector{Layer}, circuit_tuple::Vector{Int}, qubit_num::Integer, layer_types::Vector{Symbol}, layer_times::Vector{Float64})

Checks the provided circuit parameters are consistent with each other.
"""
function check_circuit(
    circuit::Vector{Layer},
    circuit_tuple::Vector{Int},
    qubit_num::Integer,
    layer_types::Vector{Symbol},
    layer_times::Vector{Float64},
)
    # Check the parameters
    L = length(circuit)
    @assert all(l.qubit_num == qubit_num for l in circuit) "All layers in the circuit must act on the same number of qubits."
    @assert circuit_tuple == collect(1:L) "Unexpected circuit tuple."
    @assert length(layer_times) == L + 1 "The layer times must correspond to the times taken for the circuit layers, alongside measurement and reset at the end."
    @assert all(layer_times .> 0.0) "The layer times must be positive."
    @assert length(layer_types) == L "The layer types must correspond to the circuit layers."
    for (idx, type) in pairs(layer_types)
        l = circuit[idx]
        # Check that only supported gates are present
        @assert ~any(is_state_prep(gate) for gate in l.layer) "The layer $(l) contains mid-circuit state preparations, which are not supported."
        @assert ~any(is_state_meas(gate) for gate in l.layer) "The layer $(l) contains mid-circuit state measurements, which are not supported."
        @assert all(
            is_one_qubit(gate) || is_two_qubit(gate; stim_supported = true) for
            gate in l.layer
        ) "The layer $(l) contains unsupported gates."
        @assert ~any(is_mid_meas(gate) for gate in l.layer) "The layer $(l) contains mid-circuit measurements, which are not currently supported."
        # Check the layer types are marked correctly
        if any(is_two_qubit(gate; stim_supported = true) for gate in l.layer)
            @assert type == :two_qubit "Layers $(l) with two-qubit gates should be marked :two_qubit."
        end
        if any(is_mid_reset(gate) || is_meas_idle(gate) for gate in l.layer)
            @assert type == :mid_reset "Layers $(l) with mid-circuit resets or measurement idles should be marked :mid_reset."
        end
        # Check the gates are consistent with the marked layer types
        if type == :single_qubit
            # Check single-qubit gate layers
            @assert all(is_one_qubit(g) for g in l.layer) "The single-qubit gate layer $(l) contains non-single-qubit gates."
        elseif type == :dynamical
            # Check dynamical decoupling gate layers
            @assert all(is_pauli(gate) for gate in l.layer) "The dynamical decoupling gate layer $(l) contains non-Pauli gates."
        elseif type == :two_qubit
            # Check two-qubit gate layers
            @assert any(is_two_qubit(gate; stim_supported = true) for gate in l.layer) "The two-qubit gate layer $(l) does not contain at least one supported two-qubit gate."
        elseif type == :mid_reset
            # Check mid-circuit reset layers
            @assert any(is_mid_reset(gate) for gate in l.layer) "The mid-circuit reset layer $(l) does not contain any resets."
            @assert all(is_mid_reset(gate) || is_meas_idle(gate) for gate in l.layer) "The mid-circuit reset layer $(l) contains gates that are not R (computational basis reset) or IM (idle during measurement or reset)."
        else
            @warn "Unknown or unsupported layer type $(type)."
        end
    end
    return nothing
end

"""
    check_code_parameters(code_param::CodeParameters, qubit_num::Integer)

Checks the provided code parameters are consistent with each other.
"""
function check_code_parameters(code_param::CodeParameters, qubit_num::Integer)
    # Check the parameters
    qubits = code_param.qubits
    @assert qubit_num == length(qubits) "The number of qubits must match the number of qubits in the code parameters."
    inverse_indices = code_param.inverse_indices
    @assert inverse_indices == Dict(qubits[i] => i for i in 1:qubit_num) "The inverse indices must invert the qubits."
    data_indices = code_param.data_indices
    data_num = length(data_indices)
    ancilla_indices = code_param.ancilla_indices
    ancilla_num = length(ancilla_indices)
    @assert data_num + ancilla_num == qubit_num "The number of data and ancilla qubits must sum to the total number of qubits."
    @assert data_indices == collect(1:data_num) "The data indices must come before the ancilla qubits."
    @assert ancilla_indices == collect((data_num + 1):qubit_num) "The ancilla indices must come after the data qubits."
    ancilla_x_indices = code_param.ancilla_x_indices
    ancilla_z_indices = code_param.ancilla_z_indices
    @assert ancilla_indices == sort(vcat(ancilla_x_indices, ancilla_z_indices)) "The ancilla X and Z check qubit indices must be a subset of the ancilla qubit indices."
    check_x_indices = code_param.check_x_indices
    check_z_indices = code_param.check_z_indices
    @assert all(
        check_x_idx_pair[1] ⊆ ancilla_indices && check_x_idx_pair[2] ⊆ data_indices for
        check_x_idx_pair in check_x_indices
    ) "The X check qubit indices must be ordered as (ancilla, data)."
    @assert all(
        check_z_idx_pair[1] ⊆ ancilla_indices && check_z_idx_pair[2] ⊆ data_indices for
        check_z_idx_pair in check_z_indices
    ) "The Z check qubit indices must be a subset of the ancilla Z check qubit indices and the data qubit indices."
    init_x_indices = code_param.init_x_indices
    init_z_indices = code_param.init_z_indices
    @assert init_x_indices ⊆ data_indices "The logical X initialisation qubit indices must be a subset of the data qubit indices."
    @assert init_z_indices ⊆ data_indices "The logical Z initialisation qubit indices must be a subset of the data qubit indices."
    logical_x_indices = code_param.logical_x_indices
    logical_z_indices = code_param.logical_z_indices
    @assert logical_x_indices ⊆ data_indices "The logical X operator qubit indices must be a subset of the data qubit indices."
    @assert logical_z_indices ⊆ data_indices "The logical Z operator qubit indices must be a subset of the data qubit indices."
    return nothing
end

"""
    prepare_circuit(circuit::Vector{Layer}, noise_param::AbstractNoiseParameters; noisy_prep::Bool = false, noisy_meas::Bool = true)

Returns a labelled copy of the circuit as well as a number of required fields for subtypes `T <: AbstractCircuit`.

# Arguments

  - `circuit::Vector{Layer}`: Circuit.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.

# Keyword arguments

  - `noisy_prep::Bool = false`: Whether to treat preparations as noisy and aim to characterise them.
  - `noisy_meas::Bool = true`: Whether to treat measurements as noisy and aim to characterise them.
  - `combined::Bool`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
  - `strict::Bool = false`: Whether to be strict about which gates count as estimable to relative precision.

# Returns

  - `circuit::Vector{Layer}`: Circuit with labelled gates.
  - `unique_layer_indices::Vector{Int}`: Unique non-measurement gate layer indices of the circuit.
  - `gates::Vector{Gate}`: Gates in the circuit.
  - `total_gates::Vector{Gate}`: Total gates in the circuit, including preparations if `noisy_prep` and measurements if `noisy_meas`.
  - `gate_data::GateData`: Gate data for the circuit.
  - `gate_probabilities::Dict{Gate, Vector{Float64}}`: Pauli error probabilities for each gate, stored as a dictionary.
  - `gate_eigenvalues::Vector{Float64}`: Eigenvalues for each gate, stored as a vector whose order is determined by the gate order in `total_gates` and indexed by `gate_data`.
"""
function prepare_circuit(
    circuit::Vector{Layer},
    noise_param::T;
    noisy_prep::Bool = false,
    noisy_meas::Bool = true,
    combined::Bool = false,
    strict::Bool = false,
) where {T <: AbstractNoiseParameters}
    # Check preparation and measurement are not simultaneously characterised
    @assert ~(noisy_prep && noisy_meas) "Preparation and measurement noise cannot be characterised simultaneously."
    # Label the circuit
    (labelled_circuit, unique_layer_indices) = label_circuit(circuit)
    # Generate the gates and total gates
    gates = get_gates(labelled_circuit)
    total_gates = get_total_gates(labelled_circuit, noisy_prep, noisy_meas)
    # Calculate the gate data
    gate_data = get_gate_data(total_gates; combined = combined, strict = strict)
    # Generate the gate noise
    gate_probabilities = init_gate_probabilities(total_gates, noise_param)
    @assert all(
        all(gate_probs .>= 0.0) && sum(gate_probs) ≈ 1.0 for
        gate_probs in values(gate_probabilities)
    ) "The gate probabilities are not all valid probability distributions."
    gate_eigenvalues = get_gate_eigenvalues(gate_probabilities, gate_data)
    # Return the data
    return (
        labelled_circuit::Vector{Layer},
        unique_layer_indices::Vector{Int},
        gates::Vector{Gate},
        total_gates::Vector{Gate},
        gate_data::GateData,
        gate_probabilities::Dict{Gate, Vector{Float64}},
        gate_eigenvalues::Vector{Float64},
    )
end

"""
    get_circuit(circuit::Vector{Layer}, layer_types::Vector{Symbol}, layer_times::Vector{Float64}, noise_param::AbstractNoiseParameters; kwargs...)

Returns a [`Circuit`](@ref) object given the supplied circuit and noise parameters.

# Arguments

  - `circuit::Vector{Layer}`: Circuit.
  - `layer_types::Vector{Symbol}`: Types of the layers in the circuit.
  - `layer_times::Vector{Float64}`: Times taken to implement each layer in the circuit, as well as measurement and reset.
  - `noise_param::AbstractNoiseParameters`: Noise parameters for the circuit.

# Keyword arguments

  - `circuit_param::AbstractCircuitParameters = EmptyCircuitParameters()`: Circuit parameters.
  - `extra_fields::Dict{Symbol, Any} = Dict{Symbol, Any}()`: Extra fields, including for [`CodeParameters`](@ref).
  - `noisy_prep::Bool = false`: Whether to treat preparations as noisy and characterise the associated noise; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `noisy_meas::Bool = true`: Whether to treat measurements as noisy and characterise the associated noise; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] : false,`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
  - `strict::Bool = false`: Whether to be strict about which gates count as estimable to relative precision.
"""
function get_circuit(
    circuit::Vector{Layer},
    layer_types::Vector{Symbol},
    layer_times::Vector{Float64},
    noise_param::T;
    circuit_param::U = EmptyCircuitParameters(),
    extra_fields::Dict{Symbol, Any} = Dict{Symbol, Any}(),
    noisy_prep::Bool = false,
    noisy_meas::Bool = true,
    combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] :
                     false,
    strict::Bool = false,
) where {T <: AbstractNoiseParameters, U <: AbstractCircuitParameters}
    # Check the circuit
    circuit_tuple = collect(1:length(circuit))
    qubit_num = convert(Int, circuit[1].qubit_num)
    check_circuit(circuit, circuit_tuple, qubit_num, layer_types, layer_times)
    if haskey(extra_fields, :code_param)
        code_param = extra_fields[:code_param]
        @assert typeof(code_param) == CodeParameters "The code parameters must be of type `CodeParameters`."
        check_code_parameters(code_param, qubit_num)
    end
    # Prepare the circuit by generating additional parameters
    (
        labelled_circuit,
        unique_layer_indices,
        gates,
        total_gates,
        gate_data,
        gate_probabilities,
        gate_eigenvalues,
    ) = prepare_circuit(
        circuit,
        noise_param;
        noisy_prep = noisy_prep,
        noisy_meas = noisy_meas,
        combined = combined,
        strict = strict,
    )
    # Return the circuit
    c = Circuit(
        circuit_param,
        labelled_circuit,
        circuit_tuple,
        qubit_num,
        unique_layer_indices,
        layer_types,
        layer_times,
        gates,
        total_gates,
        gate_data,
        noise_param,
        gate_probabilities,
        gate_eigenvalues,
        noisy_prep,
        noisy_meas,
        extra_fields,
    )
    return c::Circuit
end

"""
    get_combined_circuit(c::AbstractCircuit)

Returns a copy of the circuit `c` where the three parameters describing Pauli X, Y, and Z basis measurements have been combined into a single parameter for each qubit.
"""
function get_combined_circuit(c::T) where {T <: AbstractCircuit}
    if c.gate_data.combined
        return c::T
    else
        # Combine the preparation and measurement noise for different Pauli types
        combined_gate_data =
            get_gate_data(c.total_gates; combined = true, strict = c.gate_data.strict)
        @assert ~(c.noisy_prep && c.noisy_meas) "Preparation and measurement noise cannot be characterised simultaneously."
        combined_gate_probabilities =
            get_combined_gate_probabilities(c.gate_probabilities, combined_gate_data)
        combined_gate_eigenvalues =
            get_gate_eigenvalues(combined_gate_probabilities, combined_gate_data)
        # Create the combined circuit
        c_combined = deepcopy(c)
        @reset c_combined.gate_data = combined_gate_data
        @reset c_combined.gate_probabilities = combined_gate_probabilities
        @reset c_combined.gate_eigenvalues = combined_gate_eigenvalues
        return c_combined::T
    end
end

"""
    For an AbstractCircuitParameters type, we only particularly expect a `circuit_name` field.

The circuit name should implicitly describe all of the usually toggled parameters.

"""

"""
    For an AbstractCircuit type, we expect a large list of subfields

circuit_param
circuit
circuit_tuple
qubit_num
unique_layer_indices
layer_types
layer_times
gates
total_gates
gate_index
N
noise_param
gate_probabilities
gate_eigenvalues
add_prep
add_meas

To ensure that sign configurations are handled properly, we also need a `partition` field.
"""

struct Circuit <: AbstractCircuit
    # Circuit parameters
    circuit_param::AbstractCircuitParameters
    # Circuit
    circuit::Vector{Layer}
    # Tuple indexing the order of the circuit layers
    circuit_tuple::Vector{Int}
    # Qubit number
    qubit_num::Int
    # Indices of the unique layers in the original circuit
    # These become meaningless and are removed when a tuple is applied
    unique_layer_indices::Vector{Int}
    # Type of each layer
    layer_types::Vector{Symbol}
    # Time taken to perform each layer, including measurement and reset at the end
    layer_times::Vector{Float64}
    # Gates in the circuit tuple
    gates::Vector{Gate}
    # Total gates in the original circuit
    # Includes preparations if add_prep and measurements if add_meas
    total_gates::Vector{Gate}
    # Gate index labelling the ordering of the gate eigenvalues
    gate_index::Dict{Gate, Int}
    # Total number of gate eigenvalues
    N::Int
    # Noise parameters
    noise_param::AbstractNoiseParameters
    # Gate probabilities
    gate_probabilities::Dict{Gate, Vector{Float64}}
    # Gate eigenvalues
    gate_eigenvalues::Vector{Float64}
    # Whether to treat preparations as noisy and aim to characterise them
    add_prep::Bool
    # Whether to treat preparations as noisy and aim to characterise them
    add_meas::Bool
end

@struct_hash_equal_isequal Circuit

#
function update_noise(
    c::T,
    noise_param::U,
) where {T <: AbstractCircuit, U <: AbstractNoiseParameters}
    # Generate the noise
    gate_probabilities = get_gate_probabilities(c.total_gates, noise_param)
    gate_eigenvalues =
        get_gate_eigenvalues(gate_probabilities, c.total_gates, c.gate_index, c.N)
    # Update the circuit
    c_update = deepcopy(c)
    @reset c_update.noise_param = noise_param
    @reset c_update.gate_probabilities = gate_probabilities
    @reset c_update.gate_eigenvalues = gate_eigenvalues
    return c_update::T
end

#
function get_layer_times(
    layer_types::Vector{Symbol};
    single_qubit_time::Float64 = 0.0,
    two_qubit_time::Float64 = 0.0,
    dynamical_decoupling_time::Float64 = 0.0,
    meas_reset_time::Float64 = 0.0,
)
    # Append the layer times
    layer_times = Vector{Float64}(undef, length(layer_types) + 1)
    for (idx, type) in enumerate(layer_types)
        if type == :single_qubit
            layer_times[idx] = single_qubit_time
        elseif type == :two_qubit
            layer_times[idx] = two_qubit_time
        elseif type == :dynamical
            layer_times[idx] = dynamical_decoupling_time
        else
            throw(error("Unsupported layer type $(type)."))
        end
    end
    layer_times[end] = meas_reset_time
    @assert all(layer_times .> 0.0) "The layer times must be positive; check that layer times have been supplied for all of the relevant layer types."
    return layer_times::Vector{Float64}
end

#
function prepare_circuit(
    circuit_param::T,
    circuit::Vector{Layer},
    qubit_num::Int,
    layer_types::Vector{Symbol},
    layer_times::Vector{Float64},
    noise_param::U;
    add_prep::Bool = false,
    add_meas::Bool = true,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    # Check the parameters
    @assert length(layer_times) == length(circuit) + 1 "The layer times must correspond to the times taken for the circuit layers, alongside measurement and reset at the end."
    for (idx, type) in enumerate(layer_types)
        if type == :single_qubit
            @assert layer_times[idx] == circuit_param.single_qubit_time "The layer time $(layer_times[idx]) does not match the single-qubit layer gate time $(circuit_param.single_qubit_time)."
            @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 1 "The single-qubit layer $(circuit[idx]) does not contain only single-qubit gates."
        elseif type == :two_qubit
            @assert layer_times[idx] == circuit_param.two_qubit_time "The layer time $(layer_times[idx]) does not match the two-qubit layer gate time $(circuit_param.two_qubit_time)."
            @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 2 "The two-qubit layer $(circuit[idx]) does not contain two-qubit gates."
        elseif type == :dynamical
            @assert layer_times[idx] == circuit_param.dynamical_decoupling_time "The layer time $(layer_times[idx]) does not match the dynamical decoupling layer gate time $(circuit_param.dynamical_decoupling_time)."
            @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 1 "The dynamical decoupling layer $(circuit[idx]) does not contain only single-qubit gates."
        else
            throw(error("Unsupported layer type $(type)."))
        end
    end
    @assert layer_times[end] == circuit_param.meas_reset_time "The layer time $(layer_times[end]) does not match the measurement and reset time $(circuit_param.meas_reset_time)."
    # Label the circuit
    (circuit, unique_layer_indices) = label_circuit(circuit, qubit_num)
    # Generate the gates, total gates, and noise
    gates = get_gates(circuit)
    (total_gates, gate_index, N) = index_gates(gates, qubit_num, add_prep, add_meas)
    gate_probabilities = get_gate_probabilities(total_gates, noise_param)
    gate_eigenvalues = get_gate_eigenvalues(gate_probabilities, total_gates, gate_index, N)
    # Return the data
    return (
        circuit::Vector{Layer},
        unique_layer_indices::Vector{Int},
        gates::Vector{Gate},
        total_gates::Vector{Gate},
        gate_index::Dict{Gate, Int},
        N::Int,
        gate_probabilities::Dict{Gate, Vector{Float64}},
        gate_eigenvalues::Vector{Float64},
    )
end

struct RotatedPlanarParameters <: AbstractCircuitParameters
    # Vertical and Z distance of the code
    vertical_dist::Int
    # Horizontal and X distance of the code
    horizontal_dist::Int
    # Type of stabiliser used in the circuit
    check_type::Symbol
    # Type of two-qubit gate used in the circuit
    gate_type::Symbol
    # Whether to dynamically decouple the circuit
    dynamically_decouple::Bool
    # Whether to pad layers with single-qubit identity gates
    pad_identity::Bool
    # Single-qubit gate layer time
    single_qubit_time::Float64
    # Two-qubit gate layer time
    two_qubit_time::Float64
    # Dynamical decoupling gate layer time
    dynamical_decoupling_time::Float64
    # Measurement and reset time
    meas_reset_time::Float64
    # Name of the circuit for saving data
    circuit_name::String
    # Default constructor
    function RotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int,
        check_type::Symbol,
        gate_type::Symbol,
        dynamically_decouple::Bool,
        pad_identity::Bool,
        single_qubit_time::Float64,
        two_qubit_time::Float64,
        dynamical_decoupling_time::Float64,
        meas_reset_time::Float64,
        circuit_name::String,
    )
        # Check some conditions
        @assert (vertical_dist >= 3 && horizontal_dist >= 3) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 3 x 3."
        @assert (check_type == :xzzx || check_type == :standard) "Invalid check type $(check_type). Must be either :xzzx or :standard."
        @assert (gate_type == :cx || gate_type == :cz) "Invalid gate type $(gate_type). Must be either :cx or :cz."
        @assert (check_type == :xzzx && gate_type == :cz) ||
                (check_type == :standard && gate_type == :cx) "Unsupported pairing of check type $(check_type) and gate type $(gate_type)."
        if dynamically_decouple && ~(check_type == :xzzx && gate_type == :cz)
            @warn "Dynamical decoupling is only supported for check type :xzzx and gate type :cz."
        end
        test_circuit_name = "rotated_planar_$(vertical_dist)_$(horizontal_dist)"
        if check_type != :xzzx
            test_circuit_name *= "_check_type_$(check_type)"
        end
        if gate_type != :cz
            test_circuit_name *= "_gate_type_$(gate_type)"
        end
        if dynamically_decouple != true
            test_circuit_name *= "_dynamically_decouple_$(dynamically_decouple)"
        end
        if pad_identity != true
            test_circuit_name *= "_pad_identity_$(pad_identity)"
        end
        @assert circuit_name == test_circuit_name "The circuit name $(circuit_name) does not match the circuit name generated by the supplied parameters $(test_circuit_name)."
        # Return parameters
        return new(
            vertical_dist,
            horizontal_dist,
            check_type,
            gate_type,
            dynamically_decouple,
            pad_identity,
            single_qubit_time,
            two_qubit_time,
            dynamical_decoupling_time,
            meas_reset_time,
            circuit_name,
        )::RotatedPlanarParameters
    end
    # Constructor
    # The default gate times are specified in nanoseconds (though units ultimately don't matter) and estimated from Google device data
    # In `Suppressing quantum errors by scaling a surface code logical qubit`, they specify measurement takes 500 ns and reset takes 160 ns
    # They also specify that the overall circuit, including measurement and reset, takes 921 ns
    # They say they achieve similar or improved results as `Exponential suppression of bit or phase errors with cyclic error correction`
    # This specifies 26 ns CZ gates, and 80 ns for two layers of H and two layers of CZ, implying 14 nz H gates
    # This implies Hadamard gates take 14 ns, but the single-qubit gate layers in the original paper would then take on average 31.4 ns
    # If we assume the dynamical decoupling X gates are decomposed into 2 H gates and a Z rotation, then we can imagine the single-qubit gate layers taking 31.4 ns
    # However, there's sufficient ambiguity that we'll simply treat all layers as taking the same amount of time, namely 29 ns
    function RotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int;
        check_type::Symbol = :xzzx,
        gate_type::Symbol = :cz,
        dynamically_decouple::Bool = true,
        pad_identity::Bool = true,
        single_qubit_time::Float64 = 29.0,
        two_qubit_time::Float64 = 29.0,
        dynamical_decoupling_time::Float64 = 29.0,
        meas_reset_time::Float64 = 660.0,
    )
        # Check some conditions
        @assert (vertical_dist >= 3 && horizontal_dist >= 3) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 3 x 3."
        @assert (check_type == :xzzx || check_type == :standard) "Invalid check type $(check_type). Must be either :xzzx or :standard."
        @assert (gate_type == :cx || gate_type == :cz) "Invalid gate type $(gate_type). Must be either :cx or :cz."
        @assert (check_type == :xzzx && gate_type == :cz) ||
                (check_type == :standard && gate_type == :cx) "Unsupported pairing of check type $(check_type) and gate type $(gate_type)."
        if dynamically_decouple && ~(check_type == :xzzx && gate_type == :cz)
            @warn "Dynamical decoupling is only supported for check type :xzzx and gate type :cz."
        end
        # Generate the circuit name
        circuit_name = "rotated_planar_$(vertical_dist)_$(horizontal_dist)"
        if check_type != :xzzx
            circuit_name *= "_check_type_$(check_type)"
        end
        if gate_type != :cz
            circuit_name *= "_gate_type_$(gate_type)"
        end
        if dynamically_decouple != true
            circuit_name *= "_dynamically_decouple_$(dynamically_decouple)"
        end
        if pad_identity != true
            circuit_name *= "_pad_identity_$(pad_identity)"
        end
        # Return parameters
        return new(
            vertical_dist,
            horizontal_dist,
            check_type,
            gate_type,
            dynamically_decouple,
            pad_identity,
            single_qubit_time,
            two_qubit_time,
            dynamical_decoupling_time,
            meas_reset_time,
            circuit_name,
        )::RotatedPlanarParameters
    end
    # Square constructor
    function RotatedPlanarParameters(
        dist::Int;
        check_type::Symbol = :xzzx,
        gate_type::Symbol = :cz,
        dynamically_decouple::Bool = true,
        pad_identity::Bool = true,
        single_qubit_time::Float64 = 29.0,
        two_qubit_time::Float64 = 29.0,
        dynamical_decoupling_time::Float64 = 29.0,
        meas_reset_time::Float64 = 660.0,
    )
        # Return parameters
        return RotatedPlanarParameters(
            dist,
            dist;
            check_type = check_type,
            gate_type = gate_type,
            dynamically_decouple = dynamically_decouple,
            pad_identity = pad_identity,
            single_qubit_time = single_qubit_time,
            two_qubit_time = two_qubit_time,
            dynamical_decoupling_time = dynamical_decoupling_time,
            meas_reset_time = meas_reset_time,
        )::RotatedPlanarParameters
    end
end

@struct_hash_equal_isequal RotatedPlanarParameters

struct RotatedPlanarCircuit <: AbstractCircuit
    # Circuit parameters
    circuit_param::AbstractCircuitParameters
    # Circuit
    circuit::Vector{Layer}
    # Tuple indexing the order of the circuit layers
    circuit_tuple::Vector{Int}
    # Qubit number
    qubit_num::Int
    # Indices of the unique layers in the original circuit
    # These become meaningless and are removed when a tuple is applied
    unique_layer_indices::Vector{Int}
    # Type of each layer
    layer_types::Vector{Symbol}
    # Time taken to perform each layer, including measurement and reset at the end
    layer_times::Vector{Float64}
    # Gates in the circuit tuple
    gates::Vector{Gate}
    # Total gates in the original circuit
    # Includes preparations if add_prep and measurements if add_meas
    total_gates::Vector{Gate}
    # Gate index labelling the ordering of the gate eigenvalues
    gate_index::Dict{Gate, Int}
    # Total number of gate eigenvalues
    N::Int
    # Noise parameters
    noise_param::AbstractNoiseParameters
    # Gate probabilities
    gate_probabilities::Dict{Gate, Vector{Float64}}
    # Gate eigenvalues
    gate_eigenvalues::Vector{Float64}
    # Whether to treat preparations as noisy and aim to characterise them
    add_prep::Bool
    # Whether to treat preparations as noisy and aim to characterise them
    add_meas::Bool
    # Partition of the qubits (data, ancilla) for preparation sign configurations
    partition::Tuple{Vector{Int}, Vector{Int}}
    # Code qubits
    qubits::Vector{Tuple{Int, Int}}
    # Inverse of the code qubit indices
    inverse_indices::Dict{Tuple{Int, Int}, Int}
    # Data qubit indices
    data_indices::Vector{Int}
    # Ancilla qubit indices
    ancilla_indices::Vector{Int}
    # Ancilla X-check qubit indices
    ancilla_X_indices::Vector{Int}
    # Ancilla Z-check qubit indices
    ancilla_Z_indices::Vector{Int}
    # Code qubit layout
    qubit_layout::Matrix{String}
end

Base.show(io::IO, c::RotatedPlanarCircuit) = show(io, MIME("text/plain"), c.qubit_layout)

@struct_hash_equal_isequal RotatedPlanarCircuit

"""
    rotated_planar_circuit(rotated_param::RotatedPlanarParameters)

Generate the syndrome extraction circuit for a rotated planar surface code with the specified vertical (Z) and horizontal (X) distances. The checks can be either "xzzx" or "standard" and the gates can be either "cz" or "cx". Dynamical decoupling currently only works for "xzzx" and "cz". Note that the most natural pairings are "xzzx" and "cz", and "standard" and "cx".
"""
function rotated_planar_circuit(rotated_param::RotatedPlanarParameters)
    # Set up variables
    v = rotated_param.vertical_dist
    h = rotated_param.horizontal_dist
    check_type = rotated_param.check_type
    gate_type = rotated_param.gate_type
    dynamically_decouple = rotated_param.dynamically_decouple
    pad_identity = rotated_param.pad_identity
    single_qubit_time = rotated_param.single_qubit_time
    two_qubit_time = rotated_param.two_qubit_time
    dynamical_decoupling_time = rotated_param.dynamical_decoupling_time
    meas_reset_time = rotated_param.meas_reset_time
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    dynamical_decoupling_type = :dynamical
    # Generate the qubit indices
    data_qubits = vcat([(2i, 2j) for i in 1:v, j in 1:h]...)
    data_num = length(data_qubits)
    @assert data_num == v * h
    inner_ancilla_qubits = vcat([(2i + 1, 2j + 1) for i in 1:(v - 1), j in 1:(h - 1)]...)
    inner_num = length(inner_ancilla_qubits)
    @assert inner_num == (v - 1) * (h - 1)
    boundary_locations = vcat(
        [(2i + 1, 1) for i in 1:(v - 1)],
        [(2v + 1, 2j + 1) for j in 1:(h - 1)],
        [(2v + 1 - 2i, 2h + 1) for i in 1:(v - 1)],
        [(1, 2h + 1 - 2j) for j in 1:(h - 1)],
    )
    boundary_ancilla_qubits = boundary_locations[1:2:end]
    boundary_num = length(boundary_ancilla_qubits)
    @assert boundary_num == v + h - 2
    ancilla_qubits = vcat(inner_ancilla_qubits, boundary_ancilla_qubits)
    ancilla_num = length(ancilla_qubits)
    @assert ancilla_num == inner_num + boundary_num
    qubits = vcat(data_qubits, ancilla_qubits)
    qubit_num = length(qubits)
    @assert qubit_num == data_num + ancilla_num
    @assert qubit_num == 2v * h - 1
    data_indices = collect(1:data_num)
    ancilla_indices = collect((data_num + 1):qubit_num)
    ancilla_X_indices = ancilla_indices[(sum.(ancilla_qubits) .% 4) .== 0]
    ancilla_Z_indices = ancilla_indices[(sum.(ancilla_qubits) .% 4) .== 2]
    inverse_indices = Dict(qubits[i] => i for i in 1:qubit_num)
    # Generate the qubit layout
    qubit_layout = [" " for i in 1:(2v + 1), j in 1:(2h + 1)]
    for (i, j) in qubits[data_indices]
        qubit_layout[i, j] = "o"
    end
    for (i, j) in qubits[ancilla_X_indices]
        qubit_layout[i, j] = "x"
    end
    for (i, j) in qubits[ancilla_Z_indices]
        qubit_layout[i, j] = "z"
    end
    # Generate the X-type ancilla gate layers
    layers_X = Vector{Vector{Int}}[[], [], [], []]
    for ancilla_X_qubit in qubits[ancilla_X_indices]
        # Get indices of adjacent qubits in a vertically-mirrored N order
        adj_N_indices = [
            (ancilla_X_qubit[1] - 1, ancilla_X_qubit[2] - 1),
            (ancilla_X_qubit[1] + 1, ancilla_X_qubit[2] - 1),
            (ancilla_X_qubit[1] - 1, ancilla_X_qubit[2] + 1),
            (ancilla_X_qubit[1] + 1, ancilla_X_qubit[2] + 1),
        ]
        for (idx, adj_idx) in enumerate(adj_N_indices)
            if haskey(inverse_indices, adj_idx)
                push!(
                    layers_X[idx],
                    [inverse_indices[ancilla_X_qubit], inverse_indices[adj_idx]],
                )
            end
        end
    end
    # Generate the Z-type ancilla gate layers
    layers_Z = Vector{Vector{Int}}[[], [], [], []]
    for ancilla_Z_qubit in qubits[ancilla_Z_indices]
        # Get indices of adjacent qubits in a Z order
        adj_Z_indices = [
            (ancilla_Z_qubit[1] - 1, ancilla_Z_qubit[2] - 1),
            (ancilla_Z_qubit[1] - 1, ancilla_Z_qubit[2] + 1),
            (ancilla_Z_qubit[1] + 1, ancilla_Z_qubit[2] - 1),
            (ancilla_Z_qubit[1] + 1, ancilla_Z_qubit[2] + 1),
        ]
        for (idx, adj_idx) in enumerate(adj_Z_indices)
            if haskey(inverse_indices, adj_idx)
                push!(
                    layers_Z[idx],
                    [inverse_indices[adj_idx], inverse_indices[ancilla_Z_qubit]],
                )
            end
        end
    end
    # Construct the circuit
    if check_type == :xzzx && gate_type == :cz
        if dynamically_decouple
            circuit = [
                make_layer(["X", "H"], [data_indices, ancilla_indices], qubit_num),
                make_layer("CZ", vcat(layers_Z[1], layers_X[1]), qubit_num),
                make_layer(["H", "X"], [data_indices, ancilla_indices], qubit_num),
                make_layer("CZ", vcat(layers_Z[2], layers_X[2]), qubit_num),
                make_layer(["X", "X"], [data_indices, ancilla_indices], qubit_num),
                make_layer("CZ", vcat(layers_Z[3], layers_X[3]), qubit_num),
                make_layer(["H", "X"], [data_indices, ancilla_indices], qubit_num),
                make_layer("CZ", vcat(layers_Z[4], layers_X[4]), qubit_num),
                make_layer(["X", "H"], [data_indices, ancilla_indices], qubit_num),
            ]
            layer_types = [
                single_qubit_type,
                two_qubit_type,
                single_qubit_type,
                two_qubit_type,
                dynamical_decoupling_type,
                two_qubit_type,
                single_qubit_type,
                two_qubit_type,
                single_qubit_type,
            ]
        else
            circuit = [
                make_layer("H", ancilla_indices, qubit_num),
                make_layer("CZ", vcat(layers_Z[1], layers_X[1]), qubit_num),
                make_layer("H", data_indices, qubit_num),
                make_layer("CZ", vcat(layers_Z[2], layers_X[2]), qubit_num),
                make_layer("CZ", vcat(layers_Z[3], layers_X[3]), qubit_num),
                make_layer("H", data_indices, qubit_num),
                make_layer("CZ", vcat(layers_Z[4], layers_X[4]), qubit_num),
                make_layer("H", ancilla_indices, qubit_num),
            ]
            layer_types = [
                single_qubit_type,
                two_qubit_type,
                single_qubit_type,
                two_qubit_type,
                two_qubit_type,
                single_qubit_type,
                two_qubit_type,
                single_qubit_type,
            ]
        end
    elseif check_type == :standard && gate_type == :cx
        circuit = [
            make_layer("H", ancilla_X_indices, qubit_num),
            make_layer("CX", vcat(layers_Z[1], layers_X[1]), qubit_num),
            make_layer("CX", vcat(layers_Z[2], layers_X[2]), qubit_num),
            make_layer("CX", vcat(layers_Z[3], layers_X[3]), qubit_num),
            make_layer("CX", vcat(layers_Z[4], layers_X[4]), qubit_num),
            make_layer("H", ancilla_X_indices, qubit_num),
        ]
        layer_types = [
            single_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            single_qubit_type,
        ]
    end
    layer_times = get_layer_times(
        layer_types;
        single_qubit_time = single_qubit_time,
        two_qubit_time = two_qubit_time,
        dynamical_decoupling_time = dynamical_decoupling_time,
        meas_reset_time = meas_reset_time,
    )
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [pad_layer(l) for l in circuit]
    end
    # Return the code data
    return (
        circuit::Vector{Layer},
        qubit_num::Int,
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
        qubits::Vector{Tuple{Int, Int}},
        inverse_indices::Dict{Tuple{Int, Int}, Int},
        data_indices::Vector{Int},
        ancilla_indices::Vector{Int},
        ancilla_X_indices::Vector{Int},
        ancilla_Z_indices::Vector{Int},
        qubit_layout::Matrix{String},
    )
end

function get_circuit(
    rotated_param::RotatedPlanarParameters,
    noise_param::T;
    add_prep::Bool = false,
    add_meas::Bool = true,
) where {T <: AbstractNoiseParameters}
    # Construct the circuit
    (
        circuit,
        qubit_num,
        layer_types,
        layer_times,
        qubits,
        inverse_indices,
        data_indices,
        ancilla_indices,
        ancilla_X_indices,
        ancilla_Z_indices,
        qubit_layout,
    ) = rotated_planar_circuit(rotated_param)
    circuit_tuple = collect(1:length(circuit))
    partition = (data_indices, ancilla_indices)
    # Prepare the circuit and generate additional parameters
    (
        circuit,
        unique_layer_indices,
        gates,
        total_gates,
        gate_index,
        N,
        gate_probabilities,
        gate_eigenvalues,
    ) = prepare_circuit(
        rotated_param,
        circuit,
        qubit_num,
        layer_types,
        layer_times,
        noise_param;
        add_prep = add_prep,
        add_meas = add_meas,
    )
    # Return the circuit
    c = RotatedPlanarCircuit(
        rotated_param,
        circuit,
        circuit_tuple,
        qubit_num,
        unique_layer_indices,
        layer_types,
        layer_times,
        gates,
        total_gates,
        gate_index,
        N,
        noise_param,
        gate_probabilities,
        gate_eigenvalues,
        add_prep,
        add_meas,
        partition,
        qubits,
        inverse_indices,
        data_indices,
        ancilla_indices,
        ancilla_X_indices,
        ancilla_Z_indices,
        qubit_layout,
    )
    return c::RotatedPlanarCircuit
end

struct UnrotatedPlanarParameters <: AbstractCircuitParameters
    # Vertical and Z distance of the code
    vertical_dist::Int
    # Horizontal and X distance of the code
    horizontal_dist::Int
    # Type of two-qubit gate used in the circuit
    gate_type::Symbol
    # Whether to pad layers with single-qubit identity gates
    pad_identity::Bool
    # Single-qubit gate layer time
    single_qubit_time::Float64
    # Two-qubit gate layer time
    two_qubit_time::Float64
    # Measurement and reset time
    meas_reset_time::Float64
    # Name of the circuit for saving data
    circuit_name::String
    # Default constructor
    function UnrotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int,
        gate_type::Symbol,
        pad_identity::Bool,
        single_qubit_time::Float64,
        two_qubit_time::Float64,
        meas_reset_time::Float64,
        circuit_name::String,
    )
        # Check some conditions
        @assert (vertical_dist >= 3 && horizontal_dist >= 3) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 3 x 3."
        @assert gate_type == :cx "Invalid gate type $(gate_type). Must be :cx."
        test_circuit_name = "unrotated_planar_$(vertical_dist)_$(horizontal_dist)"
        if gate_type != :cx
            test_circuit_name *= "_gate_type_$(gate_type)"
        end
        if pad_identity != true
            test_circuit_name *= "_pad_identity_$(pad_identity)"
        end
        @assert circuit_name == test_circuit_name "The circuit name $(circuit_name) does not match the circuit name generated by the supplied parameters $(test_circuit_name)."
        # Return parameters
        return new(
            vertical_dist,
            horizontal_dist,
            gate_type,
            pad_identity,
            single_qubit_time,
            two_qubit_time,
            meas_reset_time,
            circuit_name,
        )::UnrotatedPlanarParameters
    end
    # Constructor
    # The default gate times are specified in nanoseconds (though units ultimately don't matter) and estimated from Google device data
    # In `Suppressing quantum errors by scaling a surface code logical qubit`, they specify measurement takes 500 ns and reset takes 160 ns
    # They also specify that the overall circuit, including measurement and reset, takes 921 ns
    # They say they achieve similar or improved results as `Exponential suppression of bit or phase errors with cyclic error correction`
    # This specifies 26 ns CZ gates, and 80 ns for two layers of H and two layers of CZ, implying 14 nz H gates
    # This implies Hadamard gates take 14 ns, but the single-qubit gate layers in the original paper would then take on average 31.4 ns
    # If we assume the dynamical decoupling X gates are decomposed into 2 H gates and a Z rotation, then we can imagine the single-qubit gate layers taking 31.4 ns
    # However, there's sufficient ambiguity that we'll simply treat all layers as taking the same amount of time, namely 29 ns
    function UnrotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int;
        gate_type::Symbol = :cx,
        pad_identity::Bool = true,
        single_qubit_time::Float64 = 29.0,
        two_qubit_time::Float64 = 29.0,
        meas_reset_time::Float64 = 660.0,
    )
        # Check some conditions
        @assert (vertical_dist >= 3 && horizontal_dist >= 3) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 3 x 3."
        @assert gate_type == :cx "Invalid gate type $(gate_type). Must be :cx."
        # Generate the circuit name
        circuit_name = "unrotated_planar_$(vertical_dist)_$(horizontal_dist)"
        if gate_type != :cx
            circuit_name *= "_gate_type_$(gate_type)"
        end
        if pad_identity != true
            circuit_name *= "_pad_identity_$(pad_identity)"
        end
        # Return parameters
        return new(
            vertical_dist,
            horizontal_dist,
            gate_type,
            pad_identity,
            single_qubit_time,
            two_qubit_time,
            meas_reset_time,
            circuit_name,
        )::UnrotatedPlanarParameters
    end
    # Square constructor
    function UnrotatedPlanarParameters(
        dist::Int;
        gate_type::Symbol = :cx,
        pad_identity::Bool = true,
        single_qubit_time::Float64 = 29.0,
        two_qubit_time::Float64 = 29.0,
        meas_reset_time::Float64 = 660.0,
    )
        # Return parameters
        return UnrotatedPlanarParameters(
            dist,
            dist;
            gate_type = gate_type,
            pad_identity = pad_identity,
            single_qubit_time = single_qubit_time,
            two_qubit_time = two_qubit_time,
            meas_reset_time = meas_reset_time,
        )::UnrotatedPlanarParameters
    end
end

@struct_hash_equal_isequal UnrotatedPlanarParameters

struct UnrotatedPlanarCircuit <: AbstractCircuit
    # Circuit parameters
    circuit_param::AbstractCircuitParameters
    # Circuit
    circuit::Vector{Layer}
    # Tuple indexing the order of the circuit layers
    circuit_tuple::Vector{Int}
    # Qubit number
    qubit_num::Int
    # Indices of the unique layers in the original circuit
    # These become meaningless and are removed when a tuple is applied
    unique_layer_indices::Vector{Int}
    # Type of each layer
    layer_types::Vector{Symbol}
    # Time taken to perform each layer, including measurement and reset at the end
    layer_times::Vector{Float64}
    # Gates in the circuit tuple
    gates::Vector{Gate}
    # Total gates in the original circuit
    # Includes preparations if add_prep and measurements if add_meas
    total_gates::Vector{Gate}
    # Gate index labelling the ordering of the gate eigenvalues
    gate_index::Dict{Gate, Int}
    # Total number of gate eigenvalues
    N::Int
    # Noise parameters
    noise_param::AbstractNoiseParameters
    # Gate probabilities
    gate_probabilities::Dict{Gate, Vector{Float64}}
    # Gate eigenvalues
    gate_eigenvalues::Vector{Float64}
    # Whether to treat preparations as noisy and aim to characterise them
    add_prep::Bool
    # Whether to treat preparations as noisy and aim to characterise them
    add_meas::Bool
    # Partition of the qubits (data, ancilla) for preparation sign configurations
    partition::Tuple{Vector{Int}, Vector{Int}}
    # Code qubits
    qubits::Vector{Tuple{Int, Int}}
    # Inverse of the code qubit indices
    inverse_indices::Dict{Tuple{Int, Int}, Int}
    # Data qubit indices
    data_indices::Vector{Int}
    # Ancilla qubit indices
    ancilla_indices::Vector{Int}
    # Ancilla X-check qubit indices
    ancilla_X_indices::Vector{Int}
    # Ancilla Z-check qubit indices
    ancilla_Z_indices::Vector{Int}
    # Code qubit layout
    qubit_layout::Matrix{String}
end

Base.show(io::IO, c::UnrotatedPlanarCircuit) = show(io, MIME("text/plain"), c.qubit_layout)

@struct_hash_equal_isequal UnrotatedPlanarCircuit

"""
    unrotated_planar_circuit(vertical_dist::Int, horizontal_dist::Int, gate_type::Symbol, pad_identity::Bool, single_qubit_time::Float64, two_qubit_time::Float64, meas_reset_time::Float64)

Generate the syndrome extraction circuit for an unrotated planar surface code with the specified vertical and horizontal distances. The gates can be either "cx" or "cz".
"""
function unrotated_planar_circuit(unrotated_param::UnrotatedPlanarParameters)
    # Set up variables
    v = unrotated_param.vertical_dist
    h = unrotated_param.horizontal_dist
    gate_type = unrotated_param.gate_type
    pad_identity = unrotated_param.pad_identity
    single_qubit_time = unrotated_param.single_qubit_time
    two_qubit_time = unrotated_param.two_qubit_time
    meas_reset_time = unrotated_param.meas_reset_time
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    # Generate the qubit indices
    qubits = vcat([(i, j) for i in 1:(2v - 1), j in 1:(2h - 1)]...)
    qubit_num = length(qubits)
    qubit_indices = collect(1:qubit_num)
    @assert qubit_num == (2v - 1) * (2h - 1)
    data_indices = qubit_indices[(sum.(qubits) .% 2) .== 0]
    data_num = length(data_indices)
    @assert data_num == v * h + (v - 1) * (h - 1)
    ancilla_indices = qubit_indices[(sum.(qubits) .% 2) .== 1]
    ancilla_num = length(ancilla_indices)
    @assert ancilla_num == v * (h - 1) + (v - 1) * h
    ancilla_X_indices =
        ancilla_indices[([qubits[ancilla_indices][i][1] for i in 1:ancilla_num] .% 2) .== 0]
    ancilla_Z_indices =
        ancilla_indices[([qubits[ancilla_indices][i][2] for i in 1:ancilla_num] .% 2) .== 0]
    inverse_indices = Dict(qubits[i] => i for i in 1:qubit_num)
    # Generate the qubit layout
    qubit_layout = [" " for i in 1:(2v - 1), j in 1:(2h - 1)]
    for (i, j) in qubits[data_indices]
        qubit_layout[i, j] = "o"
    end
    for (i, j) in qubits[ancilla_X_indices]
        qubit_layout[i, j] = "x"
    end
    for (i, j) in qubits[ancilla_Z_indices]
        qubit_layout[i, j] = "z"
    end
    # Generate the X-type ancilla gate layers
    layers_X = Vector{Vector{Int}}[[], [], [], []]
    for ancilla_X_qubit in qubits[ancilla_X_indices]
        # Get indices of adjacent qubits in a reverse N order, viewed after rotating the code 45 degrees clockwise
        adj_N_indices = [
            (ancilla_X_qubit[1] - 1, ancilla_X_qubit[2]),
            (ancilla_X_qubit[1], ancilla_X_qubit[2] + 1),
            (ancilla_X_qubit[1], ancilla_X_qubit[2] - 1),
            (ancilla_X_qubit[1] + 1, ancilla_X_qubit[2]),
        ]
        for (idx, adj_idx) in enumerate(adj_N_indices)
            if haskey(inverse_indices, adj_idx)
                push!(
                    layers_X[idx],
                    [inverse_indices[ancilla_X_qubit], inverse_indices[adj_idx]],
                )
            end
        end
    end
    # Generate the Z-type ancilla gate layers
    layers_Z = Vector{Vector{Int}}[[], [], [], []]
    for ancilla_Z_qubit in qubits[ancilla_Z_indices]
        # Get indices of adjacent qubits in a horizontally-mirrored Z order, viewed after rotating the code 45 degrees clockwise
        adj_Z_indices = [
            (ancilla_Z_qubit[1] - 1, ancilla_Z_qubit[2]),
            (ancilla_Z_qubit[1], ancilla_Z_qubit[2] - 1),
            (ancilla_Z_qubit[1], ancilla_Z_qubit[2] + 1),
            (ancilla_Z_qubit[1] + 1, ancilla_Z_qubit[2]),
        ]
        for (idx, adj_idx) in enumerate(adj_Z_indices)
            if haskey(inverse_indices, adj_idx)
                push!(
                    layers_Z[idx],
                    [inverse_indices[adj_idx], inverse_indices[ancilla_Z_qubit]],
                )
            end
        end
    end
    # Construct the circuit
    if gate_type == :cx
        circuit = [
            make_layer("H", ancilla_X_indices, qubit_num),
            make_layer("CX", vcat(layers_Z[1], layers_X[1]), qubit_num),
            make_layer("CX", vcat(layers_Z[2], layers_X[2]), qubit_num),
            make_layer("CX", vcat(layers_Z[3], layers_X[3]), qubit_num),
            make_layer("CX", vcat(layers_Z[4], layers_X[4]), qubit_num),
            make_layer("H", ancilla_X_indices, qubit_num),
        ]
        layer_types = [
            single_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            single_qubit_type,
        ]
    end
    layer_times = get_layer_times(
        layer_types;
        single_qubit_time = single_qubit_time,
        two_qubit_time = two_qubit_time,
        meas_reset_time = meas_reset_time,
    )
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [pad_layer(l) for l in circuit]
    end
    # Return the code data
    return (
        circuit::Vector{Layer},
        qubit_num::Int,
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
        qubits::Vector{Tuple{Int, Int}},
        inverse_indices::Dict{Tuple{Int, Int}, Int},
        data_indices::Vector{Int},
        ancilla_indices::Vector{Int},
        ancilla_X_indices::Vector{Int},
        ancilla_Z_indices::Vector{Int},
        qubit_layout::Matrix{String},
    )
end

function get_circuit(
    unrotated_param::UnrotatedPlanarParameters,
    noise_param::T;
    add_prep::Bool = false,
    add_meas::Bool = true,
) where {T <: AbstractNoiseParameters}
    # Construct the circuit
    (
        circuit,
        qubit_num,
        layer_types,
        layer_times,
        qubits,
        inverse_indices,
        data_indices,
        ancilla_indices,
        ancilla_X_indices,
        ancilla_Z_indices,
        qubit_layout,
    ) = unrotated_planar_circuit(unrotated_param)
    circuit_tuple = collect(1:length(circuit))
    partition = (data_indices, ancilla_indices)
    # Prepare the circuit and generate additional parameters
    (
        circuit,
        unique_layer_indices,
        gates,
        total_gates,
        gate_index,
        N,
        gate_probabilities,
        gate_eigenvalues,
    ) = prepare_circuit(
        unrotated_param,
        circuit,
        qubit_num,
        layer_types,
        layer_times,
        noise_param;
        add_prep = add_prep,
        add_meas = add_meas,
    )
    # Return the circuit
    c = UnrotatedPlanarCircuit(
        unrotated_param,
        circuit,
        circuit_tuple,
        qubit_num,
        unique_layer_indices,
        layer_types,
        layer_times,
        gates,
        total_gates,
        gate_index,
        N,
        noise_param,
        gate_probabilities,
        gate_eigenvalues,
        add_prep,
        add_meas,
        partition,
        qubits,
        inverse_indices,
        data_indices,
        ancilla_indices,
        ancilla_X_indices,
        ancilla_Z_indices,
        qubit_layout,
    )
    return c::UnrotatedPlanarCircuit
end

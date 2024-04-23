"""
    Circuit

Circuit information, including noise parameters.

# Fields

  - `circuit_param::AbstractCircuitParameters`: Circuit parameters.
  - `circuit::Vector{Layer}`: Circuit arranged by the tuple.
  - `circuit_tuple::Vector{Int}`: Tuple which arranges the order of the circuit layers; this is initialised as trivial.
  - `qubit_num::Int`: Number of qubits in the circuit.
  - `unique_layer_indices::Vector{Int}`: Unique layer indices of the circuit, which become meaningless and are removed the circuit is arranged by the tuple.
  - `layer_types::Vector{Symbol}`: Types of the layers in the circuit, used for layer times and dynamical decoupling.
  - `layer_times::Vector{Float64}`: Times taken to implement each layer in the circuit, as well as measurement and reset.
  - `gates::Vector{Gate}`: Gates in the circuit arranged by the tuple.
  - `total_gates::Vector{Gate}`: Gates in the original circuit, which includes noisy preparations if `add_prep` and noisy measurements if `add_meas`.
  - `gate_index::Dict{Gate, Int}`: Index of the gate eigenvalues for each gates in the original circuit.
  - `N::Int`: Number of gate eigenvalues.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.
  - `gate_probabilities::Dict{Gate, Vector{Float64}}`: Pauli error probabilities for each gate, stored as a dictionary.
  - `gate_eigenvalues::Vector{Float64}`: Eigenvalues for each gate, stored as a vector whose order is determined by `gate_index`.
  - `add_prep::Bool`: Whether to treat preparations as noisy and characterise the associated noise, defaulting to `false`; a full-rank design cannot be produced if both `add_prep` and `add_meas` are `true`.
  - `add_meas::Bool`: Whether to treat measurements as noisy and characterise the associated noise, defaulting to `true`; a full-rank design cannot be produced if both `add_prep` and `add_meas` are `true`.
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
    gate_index::Dict{Gate, Int}
    N::Int
    noise_param::AbstractNoiseParameters
    gate_probabilities::Dict{Gate, Vector{Float64}}
    gate_eigenvalues::Vector{Float64}
    add_prep::Bool
    add_meas::Bool
end

@struct_hash_equal_isequal Circuit

"""
    apply_tuple(c::AbstractCircuit, circuit_tuple::Vector{Int})

Returns a copy of the circuit `c` arranged by the tuple `circuit_tuple`.
"""
function apply_tuple(c::T, circuit_tuple::Vector{Int}) where {T <: AbstractCircuit}
    tuple_circuit = c.circuit[circuit_tuple]
    tuple_gates = get_gates(tuple_circuit)
    # Applying a tuple to the circuit makes the unique layer indices meaningless
    # Accordingly, we get rid of them to avoid confusion
    tuple_unique_layer_indices = Int[]
    # Update the parameters
    c_tuple = deepcopy(c)
    @reset c_tuple.circuit = tuple_circuit
    @reset c_tuple.circuit_tuple = circuit_tuple
    @reset c_tuple.unique_layer_indices = tuple_unique_layer_indices
    @reset c_tuple.gates = tuple_gates
    return c_tuple::T
end

"""
    update_noise(c::AbstractCircuit, noise_param::AbstractNoiseParameters)

Returns a copy of `c` where the circuit has been updated with noise generated according to `noise_param`.
"""
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

"""
    get_layer_times(layer_types::Vector{Symbol}, layer_time_dict::Dict{Symbol, Float64})

Returns the times taken to implement each layer in the circuit based on their types in `layer_types` and the times specified in `layer_time_dict, including the time for measurement and reset at the end.
"""
function get_layer_times(
    layer_types::Vector{Symbol},
    layer_time_dict::Dict{Symbol, Float64},
)
    # Append the layer times
    layer_times = Vector{Float64}(undef, length(layer_types) + 1)
    for (idx, type) in enumerate(layer_types)
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
    gates = sort(unique([gate for layer in circuit for gate in layer]))
    return gates::Vector{Gate}
end
function get_gates(circuit::Vector{Layer})
    @assert all(l.qubit_num == circuit[1].qubit_num for l in circuit) "All layers in the circuit must act on the same number of qubits."
    gates = sort(unique([gate for l in circuit for gate in l.layer]))
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
    unique_circuit = Vector{Vector{Gate}}(undef, 0)
    unique_layer_indices = Vector{Int}(undef, 0)
    for l in 1:L
        if unlabelled_circuit[l] ∉ unique_circuit
            push!(unique_circuit, unlabelled_circuit[l])
            push!(unique_layer_indices, l)
        end
    end
    # For each of the unique layers in the circuit, determine which layers they correspond to in the original circuit
    U = length(unique_circuit)
    layer_indices = Vector{Vector{Int}}(undef, U)
    for u in 1:U
        indices = Vector{Int}(undef, 0)
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
    index_gates(gates::Vector{Gate}, n::Integer, add_prep::Bool, add_meas::Bool)

Returns a dictionary indexing the gates in `gates`, adding preparations and measurements on all `n` qubits to the gates `add_prep` and `add_meas` are `true`, respectively, as well as the total collection of gates once these have been added and the number of gate eigenvalues.
"""
function index_gates(gates::Vector{Gate}, n::Integer, add_prep::Bool, add_meas::Bool)
    # Append preparations to the gate list if appropriate
    total_gates = deepcopy(gates)
    if add_prep
        append!(total_gates, make_layer("PZ-", collect(1:n), n).layer)
        append!(total_gates, make_layer("PX+", collect(1:n), n).layer)
        append!(total_gates, make_layer("PZ+", collect(1:n), n).layer)
        append!(total_gates, make_layer("PX-", collect(1:n), n).layer)
        append!(total_gates, make_layer("PY+", collect(1:n), n).layer)
        append!(total_gates, make_layer("PY-", collect(1:n), n).layer)
    end
    # Append measurements to the gate list if appropriate
    if add_meas
        append!(total_gates, make_layer("MZ", collect(1:n), n).layer)
        append!(total_gates, make_layer("MX", collect(1:n), n).layer)
        append!(total_gates, make_layer("MY", collect(1:n), n).layer)
    end
    # Determine the gate indices
    gate_index = Dict{Gate, Int}()
    N = 0
    for gate in total_gates
        not_prep = (gate.type ∉ ["PZ+", "PZ-", "PX+", "PX-", "PY+", "PY-"])
        not_meas = (gate.type ∉ ["MZ", "MX", "MY"])
        if not_prep && not_meas
            gate_index[gate] = N
            if length(gate.targets) == 1
                N += 3
            elseif length(gate.targets) == 2
                N += 15
            else
                throw(error("The gate $(gate) does not operate on either 1 or 2 qubits."))
            end
        end
        if (add_prep && (~not_prep))
            gate_index[gate] = N
            N += 1
        end
        if (add_meas && (~not_meas))
            gate_index[gate] = N
            N += 1
        end
    end
    return (total_gates::Vector{Gate}, gate_index::Dict{Gate, Int}, N::Int)
end

"""
    prepare_circuit(circuit::Vector{Layer}, qubit_num::Int, layer_types::Vector{Symbol}, layer_times::Vector{Float64}, noise_param::AbstractNoiseParameters; add_prep::Bool = false, add_meas::Bool = true)

Returns a labelled copy of the circuit as well as a number of required fields for subtypes `T <: AbstractCircuit`.

# Arguments

  - `circuit::Vector{Layer}`: Circuit.
  - `qubit_num::Int`: Number of qubits in the circuit.
  - `layer_types::Vector{Symbol}`: Types of the layers in the circuit.
  - `layer_times::Vector{Float64}`: Times taken to implement each layer in the circuit, including measurement and reset at the end.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.

# Keyword arguments

  - `add_prep::Bool = false`: Whether to treat preparations as noisy and aim to characterise them.
  - `add_meas::Bool = true`: Whether to treat measurements as noisy and aim to characterise them.

# Returns

  - `circuit::Vector{Layer}`: Circuit with labelled gates.
  - `unique_layer_indices::Vector{Int}`: Indices of the unique layers in the circuit.
  - `gates::Vector{Gate}`: Gates in the circuit.
  - `total_gates::Vector{Gate}`: Total gates in the circuit, including preparations if `add_prep` and measurements if `add_meas`.
  - `gate_index::Dict{Gate, Int}`: Index of the gate eigenvalues for each gate in the original circuit.
  - `N::Int`: Number of gate eigenvalues.
  - `gate_probabilities::Dict{Gate, Vector{Float64}}`: Pauli error probabilities for each gate, stored as a dictionary.
  - `gate_eigenvalues::Vector{Float64}`: Eigenvalues for each gate, stored as a vector whose order is determined by `gate_index`.
"""
function prepare_circuit(
    circuit::Vector{Layer},
    qubit_num::Int,
    layer_types::Vector{Symbol},
    layer_times::Vector{Float64},
    noise_param::T;
    add_prep::Bool = false,
    add_meas::Bool = true,
) where {T <: AbstractNoiseParameters}
    # Check the parameters
    @assert length(layer_times) == length(circuit) + 1 "The layer times must correspond to the times taken for the circuit layers, alongside measurement and reset at the end."
    for (idx, type) in enumerate(layer_types)
        if type == :single_qubit
            @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 1 "The single-qubit layer $(circuit[idx]) does not contain only single-qubit gates."
        elseif type == :two_qubit
            @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 2 "The two-qubit layer $(circuit[idx]) does not contain two-qubit gates."
        elseif type == :dynamical
            @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 1 "The dynamical decoupling layer $(circuit[idx]) does not contain only single-qubit gates."
        end
    end
    # Label the circuit
    (labelled_circuit, unique_layer_indices) = label_circuit(circuit)
    # Generate the gates, total gates, and noise
    gates = get_gates(labelled_circuit)
    (total_gates, gate_index, N) = index_gates(gates, qubit_num, add_prep, add_meas)
    gate_probabilities = get_gate_probabilities(total_gates, noise_param)
    gate_eigenvalues = get_gate_eigenvalues(gate_probabilities, total_gates, gate_index, N)
    # Return the data
    return (
        labelled_circuit::Vector{Layer},
        unique_layer_indices::Vector{Int},
        gates::Vector{Gate},
        total_gates::Vector{Gate},
        gate_index::Dict{Gate, Int},
        N::Int,
        gate_probabilities::Dict{Gate, Vector{Float64}},
        gate_eigenvalues::Vector{Float64},
    )
end

"""
    RotatedPlanarParameters

Parameters for the syndrome extraction circuit of a rotated surface code.

# Fields

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `check_type::Symbol`: Type of stabiliser used in the circuit, either `:xzzx` or `:standard`.
  - `gate_type::Symbol`: Type of two-qubit gate used in the circuit, either `:cx` or `:cz`.
  - `dynamically_decouple::Bool`: Whether to dynamically decouple the circuit; `true` is currently only supported for `:xzzx` and `:cz`.
  - `pad_identity::Bool`: Whether to pad layers with single-qubit identity gates.
  - `layer_time_dict::Dict{Symbol, Float64}`: Dictionary of layer times.
  - `circuit_name::String`: Name of the circuit used for saving data.
"""
struct RotatedPlanarParameters <: AbstractCircuitParameters
    vertical_dist::Int
    horizontal_dist::Int
    check_type::Symbol
    gate_type::Symbol
    dynamically_decouple::Bool
    pad_identity::Bool
    layer_time_dict::Dict{Symbol, Float64}
    circuit_name::String
    # Constructor
    function RotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int,
        check_type::Symbol,
        gate_type::Symbol,
        dynamically_decouple::Bool,
        pad_identity::Bool,
        layer_time_dict::Dict{Symbol, Float64},
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
        @assert haskey(layer_time_dict, :single_qubit) "The layer time dictionary must contain the key :single_qubit."
        @assert haskey(layer_time_dict, :two_qubit) "The layer time dictionary must contain the key :two_qubit."
        @assert haskey(layer_time_dict, :meas_reset) "The layer time dictionary must contain the key :meas_reset."
        @assert haskey(layer_time_dict, :dynamical) "The layer time dictionary must contain the key :dynamical."
        # Test the circuit name
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
            layer_time_dict,
            circuit_name,
        )::RotatedPlanarParameters
    end
end

@struct_hash_equal_isequal RotatedPlanarParameters

"""
    get_rotated_param(vertical_dist::Int, horizontal_dist::Int; kwargs...)
    get_rotated_param(dist::Int; kwargs...)

Returns a [`RotatedPlanarParameters`](@ref) object that parameterises the syndrome extraction circuit of a rotated surface code.

Default gate layer times are estimated from `Suppressing quantum errors by scaling a surface code logical qubit` by Google Quantum AI.

# Arguments

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `dist::Int`: Distance of the code; this is equivalent to setting `vertical_dist = dist` and `horizontal_dist = dist`.

# Keyword arguments

  - `check_type::Symbol = :xzzx`: Type of stabiliser used in the circuit, either `:xzzx` or `:standard`.
  - `gate_type::Symbol = :cz`: Type of two-qubit gate used in the circuit, either `:cx` or `:cz`.
  - `dynamically_decouple::Bool = true`: Whether to dynamically decouple the circuit; `true` is currently only supported for `:xzzx` and `:cz`.
  - `pad_identity::Bool = true`: Whether to pad layers with single-qubit identity gates.
  - `single_qubit_time::Float64 = 29.0`: Time taken to implement a single-qubit gate in nanoseconds.
  - `two_qubit_time::Float64 = 29.0`: Time taken to implement a two-qubit gate in nanoseconds.
  - `dynamical_decoupling_time::Float64 = 29.0`: Time taken to implement a dynamical decoupling layer in nanoseconds.
  - `meas_reset_time::Float64 = 660.0`: Time taken to perform measurement and reset at the end of the circuit in nanoseconds.
"""
function get_rotated_param(
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
    # Construct the layer time dictionary
    layer_time_dict = Dict(
        :single_qubit => single_qubit_time,
        :two_qubit => two_qubit_time,
        :dynamical => dynamical_decoupling_time,
        :meas_reset => meas_reset_time,
    )
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
    rotated_param = RotatedPlanarParameters(
        vertical_dist,
        horizontal_dist,
        check_type,
        gate_type,
        dynamically_decouple,
        pad_identity,
        layer_time_dict,
        circuit_name,
    )
    return rotated_param::RotatedPlanarParameters
end
function get_rotated_param(
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
    rotated_param = get_rotated_param(
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
    )
    return rotated_param::RotatedPlanarParameters
end

"""
    RotatedPlanarCircuit

Syndrome extraction circuit for a rotated surface code.

# Fields

  - `circuit_param::RotatedPlanarParameters`: Circuit parameters.
  - `circuit::Vector{Layer}`: Circuit arranged by the tuple.
  - `circuit_tuple::Vector{Int}`: Tuple which arranges the order of the circuit layers; this is initialised as trivial.
  - `qubit_num::Int`: Number of qubits in the circuit.
  - `unique_layer_indices::Vector{Int}`: Unique layer indices of the circuit, which become meaningless and are removed the circuit is arranged by the tuple.
  - `layer_types::Vector{Symbol}`: Types of the layers in the circuit, used for layer times and dynamical decoupling.
  - `layer_times::Vector{Float64}`: Times taken to implement each layer in the circuit, as well as measurement and reset.
  - `gates::Vector{Gate}`: Gates in the circuit arranged by the tuple.
  - `total_gates::Vector{Gate}`: Gates in the original circuit, which includes noisy preparations if `add_prep` and noisy measurements if `add_meas`.
  - `gate_index::Dict{Gate, Int}`: Index of the gate eigenvalues for each gates in the original circuit.
  - `N::Int`: Number of gate eigenvalues.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.
  - `gate_probabilities::Dict{Gate, Vector{Float64}}`: Pauli error probabilities for each gate, stored as a dictionary.
  - `gate_eigenvalues::Vector{Float64}`: Eigenvalues for each gate, stored as a vector whose order is determined by `gate_index`.
  - `add_prep::Bool`: Whether to treat preparations as noisy and characterise the associated noise, defaulting to `false`; a full-rank design cannot be produced if both `add_prep` and `add_meas` are `true`.
  - `add_meas::Bool`: Whether to treat measurements as noisy and characterise the associated noise, defaulting to `true`; a full-rank design cannot be produced if both `add_prep` and `add_meas` are `true`.
  - `partition::Tuple{Vector{Int}, Vector{Int}}`: Partition of the qubits (data, ancilla), allowing for easy preparation of sign configurations for Pauli eigenstates.
  - `qubits::Vector{Tuple{Int, Int}}`: Code qubit lattice locations.
  - `inverse_indices::Dict{Tuple{Int, Int}, Int}`: Inverse mapping from the qubit lattice locations to their indices.
  - `data_indices::Vector{Int}`: Data qubit indices.
  - `ancilla_indices::Vector{Int}`: Ancilla qubit indices.
  - `ancilla_X_indices::Vector{Int}`: Ancilla X-check qubit indices.
  - `ancilla_Z_indices::Vector{Int}`: Ancilla Z-check qubit indices.
  - `qubit_layout::Matrix{String}`: Diagram of the layout of the code qubits.
"""
struct RotatedPlanarCircuit <: AbstractCircuit
    circuit_param::RotatedPlanarParameters
    circuit::Vector{Layer}
    circuit_tuple::Vector{Int}
    qubit_num::Int
    unique_layer_indices::Vector{Int}
    layer_types::Vector{Symbol}
    layer_times::Vector{Float64}
    gates::Vector{Gate}
    total_gates::Vector{Gate}
    gate_index::Dict{Gate, Int}
    N::Int
    noise_param::AbstractNoiseParameters
    gate_probabilities::Dict{Gate, Vector{Float64}}
    gate_eigenvalues::Vector{Float64}
    add_prep::Bool
    add_meas::Bool
    partition::Tuple{Vector{Int}, Vector{Int}}
    qubits::Vector{Tuple{Int, Int}}
    inverse_indices::Dict{Tuple{Int, Int}, Int}
    data_indices::Vector{Int}
    ancilla_indices::Vector{Int}
    ancilla_X_indices::Vector{Int}
    ancilla_Z_indices::Vector{Int}
    qubit_layout::Matrix{String}
end

Base.show(io::IO, c::RotatedPlanarCircuit) = show(io, MIME("text/plain"), c.qubit_layout)

@struct_hash_equal_isequal RotatedPlanarCircuit

"""
    rotated_planar_circuit(rotated_param::RotatedPlanarParameters)

Returns fields used to construct the syndrome extraction circuit of a rotated surface code in the form of a [`RotatedPlanarCircuit`](@ref) object, based on the supplied parameters `rotated_param`.
"""
function rotated_planar_circuit(rotated_param::RotatedPlanarParameters)
    # Set up variables
    v = rotated_param.vertical_dist
    h = rotated_param.horizontal_dist
    check_type = rotated_param.check_type
    gate_type = rotated_param.gate_type
    dynamically_decouple = rotated_param.dynamically_decouple
    pad_identity = rotated_param.pad_identity
    layer_time_dict = rotated_param.layer_time_dict
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
    layer_times = get_layer_times(layer_types, layer_time_dict)
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [pad_layer(l) for l in circuit]
    end
    # Return circuit data
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

"""
    UnrotatedPlanarParameters

Parameters for the syndrome extraction circuit of an unrotated surface code.

# Fields

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `gate_type::Symbol`: Type of two-qubit gate used in the circuit, which must be `:cx`.
  - `pad_identity::Bool`: Whether to pad layers with single-qubit identity gates.
  - `layer_time_dict::Dict{Symbol, Float64}`: Dictionary of layer times.
  - `circuit_name::String`: Name of the circuit used for saving data.
"""
struct UnrotatedPlanarParameters <: AbstractCircuitParameters
    vertical_dist::Int
    horizontal_dist::Int
    gate_type::Symbol
    pad_identity::Bool
    layer_time_dict::Dict{Symbol, Float64}
    circuit_name::String
    # Constructor
    function UnrotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int,
        gate_type::Symbol,
        pad_identity::Bool,
        layer_time_dict::Dict{Symbol, Float64},
        circuit_name::String,
    )
        # Check some conditions
        @assert (vertical_dist >= 3 && horizontal_dist >= 3) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 3 x 3."
        @assert gate_type == :cx "Invalid gate type $(gate_type). Must be :cx."
        @assert haskey(layer_time_dict, :single_qubit) "The layer time dictionary must contain the key :single_qubit."
        @assert haskey(layer_time_dict, :two_qubit) "The layer time dictionary must contain the key :two_qubit."
        @assert haskey(layer_time_dict, :meas_reset) "The layer time dictionary must contain the key :meas_reset."
        # Test the circuit name
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
            layer_time_dict,
            circuit_name,
        )::UnrotatedPlanarParameters
    end
end

@struct_hash_equal_isequal UnrotatedPlanarParameters

"""
    get_unrotated_param(vertical_dist::Int, horizontal_dist::Int; kwargs...)
    get_unrotated_param(dist::Int; kwargs...)

Returns an [`UnrotatedPlanarParameters`](@ref) object that parameterises the syndrome extraction circuit of a unrotated surface code.

Default gate layer times are estimated from `Suppressing quantum errors by scaling a surface code logical qubit` by Google Quantum AI.

# Arguments

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `dist::Int`: Distance of the code; this is equivalent to setting `vertical_dist = dist` and `horizontal_dist = dist`.

# Keyword arguments

  - `gate_type::Symbol = :cx`: Type of two-qubit gate used in the circuit, which must be `:cx`.
  - `pad_identity::Bool = true`: Whether to pad layers with single-qubit identity gates.
  - `single_qubit_time::Float64 = 29.0`: Time taken to implement a single-qubit gate in nanoseconds.
  - `two_qubit_time::Float64 = 29.0`: Time taken to implement a two-qubit gate in nanoseconds.
  - `dynamical_decoupling_time::Float64 = 29.0`: Time taken to implement a dynamical decoupling layer in nanoseconds.
  - `meas_reset_time::Float64 = 660.0`: Time taken to perform measurement and reset at the end of the circuit in nanoseconds.
"""
function get_unrotated_param(
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
    # Construct the layer time dictionary
    layer_time_dict = Dict(
        :single_qubit => single_qubit_time,
        :two_qubit => two_qubit_time,
        :meas_reset => meas_reset_time,
    )
    # Generate the circuit name
    circuit_name = "unrotated_planar_$(vertical_dist)_$(horizontal_dist)"
    if gate_type != :cx
        circuit_name *= "_gate_type_$(gate_type)"
    end
    if pad_identity != true
        circuit_name *= "_pad_identity_$(pad_identity)"
    end
    # Return parameters
    unrotated_param = UnrotatedPlanarParameters(
        vertical_dist,
        horizontal_dist,
        gate_type,
        pad_identity,
        layer_time_dict,
        circuit_name,
    )
    return unrotated_param::UnrotatedPlanarParameters
end
function get_unrotated_param(
    dist::Int;
    gate_type::Symbol = :cx,
    pad_identity::Bool = true,
    single_qubit_time::Float64 = 29.0,
    two_qubit_time::Float64 = 29.0,
    meas_reset_time::Float64 = 660.0,
)
    # Return parameters
    unrotated_param = get_unrotated_param(
        dist,
        dist;
        gate_type = gate_type,
        pad_identity = pad_identity,
        single_qubit_time = single_qubit_time,
        two_qubit_time = two_qubit_time,
        meas_reset_time = meas_reset_time,
    )
    return unrotated_param::UnrotatedPlanarParameters
end

"""
    UnrotatedPlanarCircuit

Syndrome extraction circuit for a unrotated surface code.

# Fields

  - `circuit_param::UnrotatedPlanarParameters`: Circuit parameters.
  - `circuit::Vector{Layer}`: Circuit arranged by the tuple.
  - `circuit_tuple::Vector{Int}`: Tuple which arranges the order of the circuit layers; this is initialised as trivial.
  - `qubit_num::Int`: Number of qubits in the circuit.
  - `unique_layer_indices::Vector{Int}`: Unique layer indices of the circuit, which become meaningless and are removed the circuit is arranged by the tuple.
  - `layer_types::Vector{Symbol}`: Types of the layers in the circuit, used for layer times and dynamical decoupling.
  - `layer_times::Vector{Float64}`: Times taken to implement each layer in the circuit, as well as measurement and reset.
  - `gates::Vector{Gate}`: Gates in the circuit arranged by the tuple.
  - `total_gates::Vector{Gate}`: Gates in the original circuit, which includes noisy preparations if `add_prep` and noisy measurements if `add_meas`.
  - `gate_index::Dict{Gate, Int}`: Index of the gate eigenvalues for each gates in the original circuit.
  - `N::Int`: Number of gate eigenvalues.
  - `noise_param::AbstractNoiseParameters`: Noise parameters.
  - `gate_probabilities::Dict{Gate, Vector{Float64}}`: Pauli error probabilities for each gate, stored as a dictionary.
  - `gate_eigenvalues::Vector{Float64}`: Eigenvalues for each gate, stored as a vector whose order is determined by `gate_index`.
  - `add_prep::Bool`: Whether to treat preparations as noisy and characterise the associated noise, defaulting to `false`; a full-rank design cannot be produced if both `add_prep` and `add_meas` are `true`.
  - `add_meas::Bool`: Whether to treat measurements as noisy and characterise the associated noise, defaulting to `true`; a full-rank design cannot be produced if both `add_prep` and `add_meas` are `true`.
  - `partition::Tuple{Vector{Int}, Vector{Int}}`: Partition of the qubits (data, ancilla), allowing for easy preparation of sign configurations for Pauli eigenstates.
  - `qubits::Vector{Tuple{Int, Int}}`: Code qubit lattice locations.
  - `inverse_indices::Dict{Tuple{Int, Int}, Int}`: Inverse mapping from the qubit lattice locations to their indices.
  - `data_indices::Vector{Int}`: Data qubit indices.
  - `ancilla_indices::Vector{Int}`: Ancilla qubit indices.
  - `ancilla_X_indices::Vector{Int}`: Ancilla X-check qubit indices.
  - `ancilla_Z_indices::Vector{Int}`: Ancilla Z-check qubit indices.
  - `qubit_layout::Matrix{String}`: Diagram of the layout of the code qubits.
"""
struct UnrotatedPlanarCircuit <: AbstractCircuit
    circuit_param::UnrotatedPlanarParameters
    circuit::Vector{Layer}
    circuit_tuple::Vector{Int}
    qubit_num::Int
    unique_layer_indices::Vector{Int}
    layer_types::Vector{Symbol}
    layer_times::Vector{Float64}
    gates::Vector{Gate}
    total_gates::Vector{Gate}
    gate_index::Dict{Gate, Int}
    N::Int
    noise_param::AbstractNoiseParameters
    gate_probabilities::Dict{Gate, Vector{Float64}}
    gate_eigenvalues::Vector{Float64}
    add_prep::Bool
    add_meas::Bool
    partition::Tuple{Vector{Int}, Vector{Int}}
    qubits::Vector{Tuple{Int, Int}}
    inverse_indices::Dict{Tuple{Int, Int}, Int}
    data_indices::Vector{Int}
    ancilla_indices::Vector{Int}
    ancilla_X_indices::Vector{Int}
    ancilla_Z_indices::Vector{Int}
    qubit_layout::Matrix{String}
end

Base.show(io::IO, c::UnrotatedPlanarCircuit) = show(io, MIME("text/plain"), c.qubit_layout)

@struct_hash_equal_isequal UnrotatedPlanarCircuit

"""
    unrotated_planar_circuit(unrotated_param::UnrotatedPlanarParameters)

Returns fields used to construct the syndrome extraction circuit of an unrotated surface code in the form of a [`UnrotatedPlanarCircuit`](@ref) object, based on the supplied parameters `unrotated_param`.
"""
function unrotated_planar_circuit(unrotated_param::UnrotatedPlanarParameters)
    # Set up variables
    v = unrotated_param.vertical_dist
    h = unrotated_param.horizontal_dist
    gate_type = unrotated_param.gate_type
    pad_identity = unrotated_param.pad_identity
    layer_time_dict = unrotated_param.layer_time_dict
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
    layer_times = get_layer_times(layer_types, layer_time_dict)
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [pad_layer(l) for l in circuit]
    end
    # Return circuit data
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

"""
    get_circuit(rotated_param::RotatedPlanarParameters, noise_param::AbstractNoiseParameters; kwargs...)
    get_circuit(unrotated_param::UnrotatedPlanarParameters, noise_param::AbstractNoiseParameters; kwargs...)

Returns a circuit object, a subtype `T <: AbstractCircuit`, parameterised by the supplied circuit and noise parameters.

# Arguments

  - `rotated_param::RotatedPlanarParameters`: Parameters for a rotated surface code.
  - `unrotated_param::UnrotatedPlanarParameters`: Parameters for an unrotated surface code.
  - `noise_param::AbstractNoiseParameters`: Noise parameters for the circuit.

# Keyword arguments

  - `add_prep::Bool = false`: Whether to treat preparations as noisy and characterise the associated noise, defaulting to `false`; a full-rank design cannot be produced if both `add_prep` and `add_meas` are `true`.
  - `add_meas::Bool = true`: Whether to treat measurements as noisy and characterise the associated noise, defaulting to `true`; a full-rank design cannot be produced if both `add_prep` and `add_meas` are `true`.
"""
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
        labelled_circuit,
        unique_layer_indices,
        gates,
        total_gates,
        gate_index,
        N,
        gate_probabilities,
        gate_eigenvalues,
    ) = prepare_circuit(
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
        labelled_circuit,
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
        labelled_circuit,
        unique_layer_indices,
        gates,
        total_gates,
        gate_index,
        N,
        gate_probabilities,
        gate_eigenvalues,
    ) = prepare_circuit(
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
        labelled_circuit,
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

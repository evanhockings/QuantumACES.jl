"""
    Layer(gate_type::String, range_set::Vector{Vector{Int}}, n::Int)

Creates a layer of potentially multi-qubit gates, each acting on the qubits in the range set.
"""
function Layer(gate_type::String, range_set::Vector{Vector{Int}}, n::Int)
    l = Layer([Gate(gate_type, 0, range) for range in range_set], n)
    return l::Layer
end

"""
    Layer(gate_type::String, range::Vector{Int}, n::Int)

Creates a layer of single-qubit gates acting on each of the qubits in the range.
"""
function Layer(gate_type::String, range::Vector{Int}, n::Int)
    l = Layer([Gate(gate_type, 0, [qubit]) for qubit in range], n)
    return l::Layer
end

"""
    Layer(gate_types::Vector{String}, ranges::Vector{Vector{Int}}, n::Int)

Creates a layer of single-qubit gates of different types, with each acting on each of the qubits in the corresponding range.
"""
function Layer(gate_types::Vector{String}, ranges::Vector{Vector{Int}}, n::Int)
    @assert length(gate_types) == length(ranges) "The number of gate types and ranges must be the same."
    @assert length(vcat(ranges...)) == length(unique(vcat(ranges...))) "All elements in the ranges must be unique."
    gate_num = length(gate_types)
    l = Layer(
        vcat(
            [
                [Gate(gate_types[i], 0, [range]) for range in ranges[i]] for i in 1:gate_num
            ]...,
        ),
        n,
    )
    return l::Layer
end

#
function PadIdentity(l::Layer)
    target_set = sort(vcat([gate.targets for gate in l.layer]...))
    complement_set = setdiff(collect(1:(l.qubit_num)), target_set)
    layer_padded = l.layer
    for qubit in complement_set
        push!(layer_padded, Gate("I", 0, [qubit]))
    end
    l_padded = Layer(layer_padded, l.qubit_num)
    return l_padded::Layer
end

"""
    Unwrap(c::Vector{Layer})

Unwraps a circuit into the vector of vector of gates format.
"""
function Unwrap(c::Vector{Layer})
    circuit = [[gate for gate in l.layer] for l in c]
    return circuit::Vector{Vector{Gate}}
end

"""
    Gates(circuit::Vector{Vector{Gate}})

Returns the unique gates in the circuit.
"""
function Gates(circuit::Vector{Vector{Gate}})
    gates = sort(unique([gate for layer in circuit for gate in layer]))
    return gates::Vector{Gate}
end

"""
    Gates(circuit::Vector{Vector{Gate}})

Returns the unique gates in the circuit.
"""
function Gates(circuit::Vector{Layer})
    gates = sort(unique([gate for l in circuit for gate in l.layer]))
    return gates::Vector{Gate}
end

"""
    Label(circuit::Vector{Vector{Gate}})

When a gate on a particular set of qubits appears multiple in different environments, that is, non-identical layers, the different occurrences of the gate are numerically labelled in order to distinguish them.
"""
function Label(circuit::Vector{Vector{Gate}})
    # Remove any previous labels
    unlabelled_circuit =
        [sort!([Gate(gate.type, 0, gate.targets) for gate in layer]) for layer in circuit]
    # Determine the gates
    gates = Gates(unlabelled_circuit)
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
                for j in 1:length(circuit[i])
                    if gate == circuit[i][j]
                        # Append the appearance number to the gate label vector
                        circuit[i][j] = Gate(gate.type, a, gate.targets)
                    end
                end
            end
        end
    end
    return (circuit::Vector{Vector{Gate}}, unique_layer_indices::Vector{Int})
end

"""
    Label(circuit::Vector{Layer}, n::Int)

When a gate on a particular set of qubits appears multiple in different environments, that is, non-identical layers, the different occurrences of the gate are numerically labelled in order to distinguish them.
"""
function Label(circuit::Vector{Layer}, n::Int)
    unwrapped_circuit = Unwrap(circuit)
    (unwrapped_circuit, unique_layer_indices) = Label(unwrapped_circuit)
    circuit = [Layer(layer, n) for layer in unwrapped_circuit]
    return (circuit::Vector{Layer}, unique_layer_indices::Vector{Int})
end

"""
    Index(gates::Vector{Gate},n::Integer,add_prep::Bool,add_meas::Bool)

Returns the index for the gates and the total number of Paulis for the gates, and modifies the gates to include any additional preparation or measurement gates.
"""
function Index(gates::Vector{Gate}, n::Integer, add_prep::Bool, add_meas::Bool)
    # Append preparations to the gate list if appropriate
    total_gates = deepcopy(gates)
    if add_prep
        append!(total_gates, Layer("PZ-", collect(1:n), n).layer)
        append!(total_gates, Layer("PX+", collect(1:n), n).layer)
        append!(total_gates, Layer("PZ+", collect(1:n), n).layer)
        append!(total_gates, Layer("PX-", collect(1:n), n).layer)
        append!(total_gates, Layer("PY+", collect(1:n), n).layer)
        append!(total_gates, Layer("PY-", collect(1:n), n).layer)
    end
    # Append measurements to the gate list if appropriate
    if add_meas
        append!(total_gates, Layer("MZ", collect(1:n), n).layer)
        append!(total_gates, Layer("MX", collect(1:n), n).layer)
        append!(total_gates, Layer("MY", collect(1:n), n).layer)
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
    ApplyTuple(c::AbstractCircuit, circuit_tuple::Vector{Int})

Returns a new circuit struct whose circuit layers are the original circuit layers rearranged by a tuple, or arrangement with repetition.
"""
function ApplyTuple(c::AbstractCircuit, circuit_tuple::Vector{Int})
    tuple_circuit = c.circuit[circuit_tuple]
    tuple_gates = Gates(tuple_circuit)
    # Applying a tuple to the circuit makes the unique layer indices meaningless
    # Accordingly, we get rid of them to avoid confusion
    tuple_unique_layer_indices = Int[]
    # Update the parameters
    c_tuple = deepcopy(c)
    @reset c_tuple.circuit = tuple_circuit
    @reset c_tuple.circuit_tuple = circuit_tuple
    @reset c_tuple.unique_layer_indices = tuple_unique_layer_indices
    @reset c_tuple.gates = tuple_gates
    return c_tuple::AbstractCircuit
end

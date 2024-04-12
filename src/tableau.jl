# Stabiliser circuit simulation following `Improved Simulation of Stabilizer Circuits` by S. Aaronson and D. Gottesman (2004)
# arXiv:quant-ph/0406196

mutable struct Tableau
    # Tableau
    tableau::Matrix{Bool}
    # Qubit number
    qubit_num::Int16
    # Defaunt constructor
    function Tableau(n::Integer)
        tableau = zeros(Bool, 2n + 1, 2n + 1)
        for i in 1:(2n)
            tableau[i, i] = 1
        end
        return new(tableau, n)::Tableau
    end
end

function Base.show(io::IO, t::Tableau)
    for i in 1:(t.qubit_num)
        println(io, get_pauli_string(Pauli(t.tableau[t.qubit_num + i, :], t.qubit_num)))
    end
end

@struct_hash_equal_isequal Tableau

struct Gate
    # Gate type
    type::String
    # Index labelling the unique occurrences of the gate in a circuit
    # Two gates are considered unique if they occur in non-identical layers
    index::Int32
    # Gate targets
    targets::Vector{Int16}
end

Base.show(io::IO, g::Gate) = print(io, "($(g.type)-$(Int(g.index)):$(Int.(g.targets)))")

function Base.isless(g₁::Gate, g₂::Gate)
    return isless([g₁.type; g₁.index; g₁.targets], [g₂.type; g₂.index; g₂.targets])
end

@struct_hash_equal_isequal Gate

struct Layer
    # Layer of gates
    layer::Vector{Gate}
    # Qubit number
    qubit_num::Int16
    # Default constructor
    function Layer(layer::Vector{Gate}, qubit_num::Integer)
        # Parameter checks
        target_set = sort(vcat([gate.targets for gate in layer]...))
        @assert all(target_set .>= 1) && all(target_set .<= qubit_num) "The layer $(layer) acts on a qubit that is out of bounds on $(qubit_num) qubits."
        @assert target_set == unique(target_set) "The layer $(layer) acts on the same target more than once."
        return new(layer, qubit_num)::Layer
    end
end

Base.show(io::IO, l::Layer) = print(io, [gate for gate in l.layer])

@struct_hash_equal_isequal Layer

"""
    cx!(t::Tableau,control::Integer,target::Integer)

Perform a CNOT on the tableau with the specified control and target qubits.
"""
function cx!(t::Tableau, control::Integer, target::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert (target >= 1) && (target <= n) "The target $(target) on $(n) qubits is out of bounds."
    @assert (control >= 1) && (control <= n) "The control $(control) on $(n) qubits is out of bounds."
    @assert target != control "The control and target cannot both be the same qubit $(target)."
    for i in 1:(2n)
        # Set the phase variables
        tableau[i, 2n + 1] =
            (
                tableau[i, 2n + 1] +
                (tableau[i, control] * tableau[i, n + target]) *
                ((tableau[i, target] + tableau[i, n + control] + 1) % 2)
            ) % 2
        # Set the x variables
        tableau[i, target] = (tableau[i, target] + tableau[i, control]) % 2
        # Set the z variables
        tableau[i, n + control] = (tableau[i, n + control] + tableau[i, n + target]) % 2
    end
    return nothing
end

"""
    hadamard!(t::Tableau,target::Integer)

Perform a Hadamard gate on the tableau with the specified target qubit.
"""
function hadamard!(t::Tableau, target::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert (target >= 1) && (target <= n) "The target $(target) on $(n) qubits is out of bounds."
    for i in 1:(2n)
        # Set the phase variables
        tableau[i, 2n + 1] =
            (tableau[i, 2n + 1] + (tableau[i, target] * tableau[i, n + target])) % 2
        # Swap the x and z variables
        temp_tableau_var = tableau[i, target]
        tableau[i, target] = tableau[i, n + target]
        tableau[i, n + target] = temp_tableau_var
    end
    return nothing
end

"""
    phase!(t::Tableau,target::Integer)

Perform a phase gate on the tableau with the specified target qubit.
"""
function phase!(t::Tableau, target::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert (target >= 1) && (target <= n) "The target $(target) on $(n) qubits is out of bounds."
    for i in 1:(2n)
        # Set the phase variables
        tableau[i, 2n + 1] =
            (tableau[i, 2n + 1] + (tableau[i, target] * tableau[i, n + target])) % 2
        # Set the z variables
        tableau[i, n + target] = (tableau[i, n + target] + tableau[i, target]) % 2
    end
    return nothing
end

"""
    x!(t::Tableau,target::Integer)

Perform a Pauli X gate on the tableau with the specified target qubit, noting that ``X=HS^2H``.
"""
function x!(t::Tableau, target::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert (target >= 1) && (target <= n) "The target $(target) on $(n) qubits is out of bounds."
    for i in 1:(2n)
        # Set the phase variables
        tableau[i, 2n + 1] = (tableau[i, 2n + 1] + tableau[i, n + target]) % 2
    end
    return nothing
end

"""
    z!(t::Tableau,target::Integer)

Perform a Pauli Z gate on the tableau with the specified target qubit, noting that ``Z=S^2``.
"""
function z!(t::Tableau, target::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert (target >= 1) && (target <= n) "The target $(target) on $(n) qubits is out of bounds."
    for i in 1:(2n)
        # Set the phase variables
        tableau[i, 2n + 1] = (tableau[i, 2n + 1] + tableau[i, target]) % 2
    end
    return nothing
end
"""
    y!(t::Tableau,target::Integer)

Perform a Pauli Y gate on the tableau with the specified target qubit, noting that ``Y=SHS^2HS^3``.
"""
function y!(t::Tableau, target::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert (target >= 1) && (target <= n) "The target $(target) on $(n) qubits is out of bounds."
    for i in 1:(2n)
        # Set the phase variables
        tableau[i, 2n + 1] =
            (tableau[i, 2n + 1] + tableau[i, target] + tableau[i, n + target]) % 2
    end
    return nothing
end

#
function cz!(t::Tableau, control::Integer, target::Integer)
    hadamard!(t, target)
    cx!(t, control, target)
    hadamard!(t, target)
    return nothing
end

#
function sqrt_zz!(t::Tableau, control::Integer, target::Integer)
    cz!(t, control, target)
    phase!(t, control)
    phase!(t, target)
    return nothing
end

#
function sqrt_zz_dag!(t::Tableau, control::Integer, target::Integer)
    sqrt_zz!(t, control, target)
    z!(t, control)
    z!(t, target)
    return nothing
end

"""
    row_phase(x₁::Bool,z₁::Bool,x₂::Bool,z₂::Bool)

Calculate a phase for RowSum.
"""
function row_phase(x₁::Bool, z₁::Bool, x₂::Bool, z₂::Bool)
    # Returns the exponent to which i is raised if we multiply x₁z₁ by x₂z₂
    if x₁ == 0 && z₁ == 0
        g = 0
    elseif x₁ == 0 && z₁ == 1
        g = x₂ * (1 - 2z₂)
    elseif x₁ == 1 && z₁ == 0
        g = z₂ * (2x₂ - 1)
    elseif x₁ == 1 && z₁ == 1
        g = z₂ - x₂
    else
        throw(error("How are the Booleans $(x₁), $(z₁), $(x₂), $(z₂) so messed up?"))
    end
    return g::Integer
end

"""
    row_sum!(t::Tableau,target::Integer,control::Integer)

Add the control row to the target row, keeping track of the phase bit of the target row.
"""
function row_sum!(t::Tableau, target::Integer, control::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert target != control "The control and target cannot both be the same qubit $(target)."
    # Adds the row `control` (i in Aaronson 2004) to the row `target` (h in Aaronson 2004), keeping track of the phase bit of the row `target`
    row_sum =
        2 * tableau[target, 2n + 1] +
        2 * tableau[control, 2n + 1] +
        sum(
            row_phase(
                tableau[control, i],
                tableau[control, n + i],
                tableau[target, i],
                tableau[target, n + i],
            ) for i in 1:n
        )
    # Determine the phase
    if (row_sum + 4n) % 4 == 0
        tableau[target, 2n + 1] = 0
    elseif (row_sum + 4n) % 4 == 2
        tableau[target, 2n + 1] = 1
    else
        throw(error("The row sum $(row_sum) is odd when it should be even."))
    end
    # Add the control to the target
    for i in 1:n
        tableau[target, i] = (tableau[target, i] + tableau[control, i]) % 2
        tableau[target, n + i] = (tableau[target, n + i] + tableau[control, n + i]) % 2
    end
    return nothing
end

"""
    measure!(t::Tableau,target::Integer)

Perform a measurement on the tableau with the specified target qubit, and return the measurement outcome.
"""
function measure!(t::Tableau, target::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert (target >= 1) && (target <= n) "The target $(target) on $(n) qubits is out of bounds."
    p = findfirst(tableau[(n + 1):(2n), target])
    if isnothing(p)
        # The measurement outcome is determined
        # Zero the last row
        tableau[2n + 1, :] = zeros(Bool, 2n + 1)
        # Perform rowsum operations
        for i in 1:n
            if tableau[i, target] == 1
                row_sum!(t, 2n + 1, n + i)
            end
        end
        # Return the measurement outcome
        measurement = tableau[2n + 1, 2n + 1]
    else
        # The measurement outcome is random
        # Set the true value for p
        p += n
        # Perform rowsum operations
        for i in 1:(2n)
            if tableau[i, target] == 1 && i != p && i != p - n
                row_sum!(t, i, p)
            end
        end
        # Set the pth row as a destabiliser
        tableau[p - n, :] = tableau[p, :]
        # Set the pth row
        tableau[p, :] = zeros(Bool, 2n + 1)
        tableau[p, n + target] = 1
        tableau[p, 2n + 1] = rand(Bool)
        measurement = tableau[p, 2n + 1]
    end
    return measurement::Bool
end

"""
    reset!(t::Tableau,target::Integer)

Reset the specified target qubit into the ``|0\\rangle`` state by measuring in the Z basis and flipping it if it is in the ``|1\\rangle`` state.
"""
function reset!(t::Tableau, target::Integer)
    measurement = measure!(t, target)
    if measurement == true
        x!(t, target)
    end
    return nothing
end

"""
    apply!(t::Tableau, l::Layer; return_measurements::Bool = false)

Perform all the gates in the layer on the tableau and return the list of measurement outcomes if return_measurements.
"""
function apply!(t::Tableau, l::Layer; return_measurements::Bool = false)
    sqrt_rot = [
        "SQRT_XX"
        "SQRT_XZ"
        "SQRT_XY"
        "SQRT_ZX"
        "SQRT_ZZ"
        "SQRT_ZY"
        "SQRT_YX"
        "SQRT_YZ"
        "SQRT_YY"
    ]
    sqrt_rot_dag = [
        "SQRT_XX_DAG"
        "SQRT_XZ_DAG"
        "SQRT_XY_DAG"
        "SQRT_ZX_DAG"
        "SQRT_ZZ_DAG"
        "SQRT_ZY_DAG"
        "SQRT_YX_DAG"
        "SQRT_YZ_DAG"
        "SQRT_YY_DAG"
    ]
    measurements = Vector{Tuple{Int, Int}}(undef, 0)
    for g in l.layer
        # Perform the relevant operation
        if g.type == "CX" || g.type == "CNOT"
            cx!(t, g.targets[1], g.targets[2])
        elseif g.type == "CZ"
            cz!(t, g.targets[1], g.targets[2])
        elseif g.type == "H"
            hadamard!(t, g.targets[1])
        elseif g.type == "S"
            phase!(t, g.targets[1])
        elseif g.type == "S_DAG"
            phase!(t, g.targets[1])
            z!(t, g.targets[1])
        elseif g.type == "I"
            continue
        elseif g.type == "X"
            x!(t, g.targets[1])
        elseif g.type == "Z"
            z!(t, g.targets[1])
        elseif g.type == "Y"
            y!(t, g.targets[1])
        elseif g.type == "PZ+"
        elseif g.type == "PZ-"
            x!(t, g.targets[1])
        elseif g.type == "PX+"
            hadamard!(t, g.targets[1])
        elseif g.type == "PX-"
            x!(t, g.targets[1])
            hadamard!(t, g.targets[1])
        elseif g.type == "PY+"
            hadamard!(t, g.targets[1])
            phase!(t, g.targets[1])
        elseif g.type == "PY-"
            x!(t, g.targets[1])
            hadamard!(t, g.targets[1])
            phase!(t, g.targets[1])
        elseif g.type == "II"
            continue
        elseif g.type ∈ sqrt_rot || g.type ∈ sqrt_rot_dag
            # Set up variables
            pauli_1 = g.type[6]
            pauli_2 = g.type[7]
            # Rotate into the correct Pauli basis
            if pauli_1 == 'X'
                hadamard!(t, g.targets[1])
            elseif pauli_1 == 'Y'
                phase!(t, g.targets[1])
                z!(t, g.targets[1])
                hadamard!(t, g.targets[1])
            elseif pauli_1 == 'Z'
            else
                throw(error("There's a problem with $(g)."))
            end
            if pauli_2 == 'X'
                hadamard!(t, g.targets[2])
            elseif pauli_2 == 'Y'
                phase!(t, g.targets[2])
                z!(t, g.targets[2])
                hadamard!(t, g.targets[2])
            elseif pauli_2 == 'Z'
            else
                throw(error("There's a problem with $(g)."))
            end
            # Apply the Pauli rotation
            if g.type ∈ sqrt_rot && g.type ∉ sqrt_rot_dag
                sqrt_zz!(t, g.targets[1], g.targets[2])
            elseif g.type ∈ sqrt_rot_dag && g.type ∉ sqrt_rot
                sqrt_zz_dag!(t, g.targets[1], g.targets[2])
            else
                throw(error("There's a problem with $(g)."))
            end
            # Rotate back into the correct Pauli basis
            if pauli_1 == 'X'
                hadamard!(t, g.targets[1])
            elseif pauli_1 == 'Y'
                hadamard!(t, g.targets[1])
                phase!(t, g.targets[1])
            elseif pauli_1 == 'Z'
            else
                throw(error("There's a problem with $(g)."))
            end
            if pauli_2 == 'X'
                hadamard!(t, g.targets[2])
            elseif pauli_2 == 'Y'
                hadamard!(t, g.targets[2])
                phase!(t, g.targets[2])
            elseif pauli_2 == 'Z'
            else
                throw(error("There's a problem with $(g)."))
            end
        elseif g.type == "XX"
            x!(t, g.targets[1])
            x!(t, g.targets[2])
        elseif g.type == "XZ"
            x!(t, g.targets[1])
            z!(t, g.targets[2])
        elseif g.type == "XY"
            x!(t, g.targets[1])
            y!(t, g.targets[2])
        elseif g.type == "ZX"
            z!(t, g.targets[1])
            x!(t, g.targets[2])
        elseif g.type == "ZZ"
            z!(t, g.targets[1])
            z!(t, g.targets[2])
        elseif g.type == "ZY"
            z!(t, g.targets[1])
            y!(t, g.targets[2])
        elseif g.type == "YX"
            y!(t, g.targets[1])
            x!(t, g.targets[2])
        elseif g.type == "YZ"
            y!(t, g.targets[1])
            z!(t, g.targets[2])
        elseif g.type == "YY"
            y!(t, g.targets[1])
            y!(t, g.targets[2])
        elseif g.type == "MZ" || g.type == "M"
            measurement = measure!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "MX"
            hadamard!(t, g.targets[1])
            measurement = measure!(t, g.targets[1])
            hadamard!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "MY"
            phase!(t, g.targets[1])
            z!(t, g.targets[1])
            hadamard!(t, g.targets[1])
            measurement = measure!(t, g.targets[1])
            hadamard!(t, g.targets[1])
            phase!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "R"
            reset!(t, g.targets[1])
        else
            throw(error("The layer has an unsupported operation $(g)."))
        end
    end
    if return_measurements
        return measurements::Vector{Tuple{Int, Int}}
    else
        return nothing
    end
end

"""
    make_layer(gate_type::String, range_set::Vector{Vector{Int}}, n::Int)

Creates a layer of potentially multi-qubit gates, each acting on the qubits in the range set.
"""
function make_layer(gate_type::String, range_set::Vector{Vector{Int}}, n::Int)
    l = Layer([Gate(gate_type, 0, range) for range in range_set], n)
    return l::Layer
end

"""
    make_layer(gate_type::String, range::Vector{Int}, n::Int)

Creates a layer of single-qubit gates acting on each of the qubits in the range.
"""
function make_layer(gate_type::String, range::Vector{Int}, n::Int)
    l = Layer([Gate(gate_type, 0, [qubit]) for qubit in range], n)
    return l::Layer
end

"""
    make_layer(gate_types::Vector{String}, ranges::Vector{Vector{Int}}, n::Int)

Creates a layer of single-qubit gates of different types, with each acting on each of the qubits in the corresponding range.
"""
function make_layer(gate_types::Vector{String}, ranges::Vector{Vector{Int}}, n::Int)
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
function pad_layer(l::Layer)
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
    unwrap_circuit(c::Vector{Layer})

unwrap_circuits a circuit into the vector of vector of gates format.
"""
function unwrap_circuit(c::Vector{Layer})
    circuit = [[gate for gate in l.layer] for l in c]
    return circuit::Vector{Vector{Gate}}
end

"""
    get_gates(circuit::Vector{Vector{Gate}})

Returns the unique gates in the circuit.
"""
function get_gates(circuit::Vector{Vector{Gate}})
    gates = sort(unique([gate for layer in circuit for gate in layer]))
    return gates::Vector{Gate}
end

"""
    get_gates(circuit::Vector{Vector{Gate}})

Returns the unique gates in the circuit.
"""
function get_gates(circuit::Vector{Layer})
    gates = sort(unique([gate for l in circuit for gate in l.layer]))
    return gates::Vector{Gate}
end

"""
    label_circuit(circuit::Vector{Vector{Gate}})

When a gate on a particular set of qubits appears multiple in different environments, that is, non-identical layers, the different occurrences of the gate are numerically labelled in order to distinguish them.
"""
function label_circuit(circuit::Vector{Vector{Gate}})
    # Remove any previous labels
    unlabelled_circuit =
        [sort!([Gate(gate.type, 0, gate.targets) for gate in layer]) for layer in circuit]
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
                for j in eachindex(circuit[i])
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
    label_circuit(circuit::Vector{Layer}, n::Int)

When a gate on a particular set of qubits appears multiple in different environments, that is, non-identical layers, the different occurrences of the gate are numerically labelled in order to distinguish them.
"""
function label_circuit(circuit::Vector{Layer}, n::Int)
    unwrapped_circuit = unwrap_circuit(circuit)
    (unwrapped_circuit, unique_layer_indices) = label_circuit(unwrapped_circuit)
    circuit = [Layer(layer, n) for layer in unwrapped_circuit]
    return (circuit::Vector{Layer}, unique_layer_indices::Vector{Int})
end

"""
    index_gates(gates::Vector{Gate},n::Integer,add_prep::Bool,add_meas::Bool)

Returns the index for the gates and the total number of Paulis for the gates, and modifies the gates to include any additional preparation or measurement gates.
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
    apply_tuple(c::T, circuit_tuple::Vector{Int}) where {T<:AbstractCircuit}

Returns a new circuit struct whose circuit layers are the original circuit layers rearranged by a tuple, or arrangement with repetition.
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

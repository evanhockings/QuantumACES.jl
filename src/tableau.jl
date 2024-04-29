"""
    Tableau

A tableau representation of a stabiliser state.

Stabiliser circuit simulations follow `Improved simulation of stabilizer circuits` by S. Aaronson and D. Gottesman (2004).

# Fields

  - `tableau::Matrix{Bool}`: The tableau representation of the stabiliser state.
  - `qubit_num::Int16`: The number of qubits in the stabiliser state.
"""
mutable struct Tableau
    tableau::Matrix{Bool}
    qubit_num::Int16
    # Constructor
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

"""
    Gate

A gate in a stabiliser circuit.

# Fields

  - `type::String`: String describing the gate type.
  - `index::Int32`: The index labelling the unique layer occurrences of the gate in a circuit.
  - `targets::Vector{Int16}`: The qubit targets of the gate.

# Supported gates

  - `H`: Hadamard gate.
  - `S`: Phase gate.
  - `CX` or `CNOT`: Controlled-X gate.
  - `CZ`: Controlled-Z gate.
  - `I`: Identity gate.
  - `Z`: Pauli Z gate.
  - `X`: Pauli X gate.
  - `Y`: Pauli Y gate.
  - `II`: Two-qubit identity gate.
  - `AB`: Two-qubit Pauli gate, where `A` and `B` are Paulis `Z`, `X`, or `Y`.
  - `SQRT_AB`: Two-qubit Pauli rotation, where `A` and `B` are Paulis `Z`, `X`, or `Y`.
  - `SQRT_AB_DAG` : Two-qubit Pauli rotation, where `A` and `B` are Paulis `Z`, `X`, or `Y`.
  - `PZ+`: Prepare the Pauli +Z eigenstate.
  - `PZ-`: Prepare the Pauli -Z eigenstate.
  - `PX+`: Prepare the Pauli +X eigenstate.
  - `PX-`: Prepare the Pauli -X eigenstate.
  - `PY+`: Prepare the Pauli +Y eigenstate.
  - `PY-`: Prepare the Pauli -Y eigenstate.
  - `M` or `MZ`: Measure in the computational Pauli Z basis.
  - `MX`: Measure in the Pauli X basis.
  - `MY`: Measure in the Pauli Y basis.
  - `R`: Reset to the computational Z basis.
"""
struct Gate
    type::String
    index::Int32
    targets::Vector{Int16}
end

Base.show(io::IO, g::Gate) = print(io, "($(g.type):$(Int(g.index)):$(Int.(g.targets)))")

function Base.isless(g₁::Gate, g₂::Gate)
    return isless([g₁.type; g₁.index; g₁.targets], [g₂.type; g₂.index; g₂.targets])
end

@struct_hash_equal_isequal Gate

"""
    Layer

A layer of gates in a stabiliser circuit.
Gates in a layer are simultaneously implemented by the device, and act on disjoint sets of qubits such that they trivially commute with each other.

# Fields

  - `layer::Vector{Gate}`: The gates in the layer.
  - `qubit_num::Int16`: The number of qubits in the circuit.
"""
struct Layer
    layer::Vector{Gate}
    qubit_num::Int16
    # Constructor
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
    hadamard!(t::Tableau,target::Integer)

Perform a Hadamard gate on the tableau `t` with target qubit `target`.
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

Perform a phase gate on the tableau `t` with target qubit `target`.
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
    cx!(t::Tableau, control::Integer, target::Integer)

Perform a controlled-X gate on the tableau `t` with control qubit `control` and target qubit `target`.
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
    cz!(t::Tableau, control::Integer, target::Integer)

Perform a controlled-Z gate on the tableau `t` with control qubit `control` and target qubit `target`.
The gate is symmetric, so the control and target qubits can be swapped.
"""
function cz!(t::Tableau, control::Integer, target::Integer)
    hadamard!(t, target)
    cx!(t, control, target)
    hadamard!(t, target)
    return nothing
end

"""
    x!(t::Tableau,target::Integer)

Perform a Pauli X gate on the tableau `t` with target qubit `target`, noting that ``X = H S^2 H``.
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

Perform a Pauli Z gate on the tableau `t` with target qubit `target`, noting that ``Z = S^2``.
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

Perform a Pauli Z gate on the tableau `t` with target qubit `target`, noting that ``Y = S H S^2 H S^3``.
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

"""
    sqrt_zz!(t::Tableau, control::Integer, target::Integer)

Perform a π/2 ZZ Pauli rotation gate on the tableau `t` with control qubit `control` and target qubit `target`.
The gate is symmetric, so the control and target qubits can be swapped.
"""
function sqrt_zz!(t::Tableau, control::Integer, target::Integer)
    cz!(t, control, target)
    phase!(t, control)
    phase!(t, target)
    return nothing
end

"""
    sqrt_zz_dag!(t::Tableau, control::Integer, target::Integer)

Perform a -π/2 ZZ Pauli rotation gate on the tableau `t` with control qubit `control` and target qubit `target`.
The gate is symmetric, so the control and target qubits can be swapped.
"""
function sqrt_zz_dag!(t::Tableau, control::Integer, target::Integer)
    sqrt_zz!(t, control, target)
    z!(t, control)
    z!(t, target)
    return nothing
end

"""
    row_phase(x₁::Bool, z₁::Bool, x₂::Bool, z₂::Bool)

Calculate a phase for [`row_sum!`](@ref).
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
    row_sum!(t::Tableau, target::Integer, control::Integer)

In the tableau `t`, add the control row `control` to the target row `target`, while tracking the phase bit of the target row.
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
    measure!(t::Tableau, target::Integer)

Measure the tableau `t` at the target qubit `target`, and return the measurement outcome.
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
    reset!(t::Tableau, target::Integer)

Reset the tableau `t` at the target qubit `target` by measuring in the computational basis and flipping the phase if the measurement outcome is -1.
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

Perform on the tableau `t` all gates in the layer `l`, and return the list of measurement outcomes if `return_measurements` is `true`.
"""
function apply!(t::Tableau, l::Layer; return_measurements::Bool = false)
    pauli_rot = ["XX"; "XZ"; "XY"; "ZX"; "ZZ"; "ZY"; "YX"; "YZ"; "YY"]
    sqrt_rot = ["SQRT_" * pauli for pauli in pauli_rot]
    sqrt_rot_dag = ["SQRT_" * pauli * "_DAG" for pauli in pauli_rot]
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
        elseif g.type == "II"
            continue
        elseif g.type ∈ pauli_rot
            # Set up variables
            pauli_1 = g.type[1]
            pauli_2 = g.type[2]
            if pauli_1 == 'X'
                x!(t, g.targets[1])
            elseif pauli_1 == 'Y'
                y!(t, g.targets[1])
            elseif pauli_1 == 'Z'
                z!(t, g.targets[1])
            else
                throw(error("There's a problem with $(g)."))
            end
            if pauli_2 == 'X'
                x!(t, g.targets[2])
            elseif pauli_2 == 'Y'
                y!(t, g.targets[2])
            elseif pauli_2 == 'Z'
                z!(t, g.targets[2])
            else
                throw(error("There's a problem with $(g)."))
            end
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
                continue
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
                continue
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
                continue
            else
                throw(error("There's a problem with $(g)."))
            end
            if pauli_2 == 'X'
                hadamard!(t, g.targets[2])
            elseif pauli_2 == 'Y'
                hadamard!(t, g.targets[2])
                phase!(t, g.targets[2])
            elseif pauli_2 == 'Z'
                continue
            else
                throw(error("There's a problem with $(g)."))
            end
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

Returns a layer of `gate_type` gates, each acting on the qubits in `range_set`, where the layer acts on `n` qubits.
"""
function make_layer(gate_type::String, range_set::Vector{Vector{Int}}, n::Int)
    l = Layer([Gate(gate_type, 0, range) for range in range_set], n)
    return l::Layer
end

"""
    make_layer(gate_type::String, range::Vector{Int}, n::Int)

Returns a layer of single-qubit `gate_type` gates acting on the qubits in `range`, where the layer acts on `n` qubits.
"""
function make_layer(gate_type::String, range::Vector{Int}, n::Int)
    l = Layer([Gate(gate_type, 0, [qubit]) for qubit in range], n)
    return l::Layer
end

"""
    make_layer(gate_types::Vector{String}, ranges::Vector{Vector{Int}}, n::Int)

Returns a layer of single-qubit gates, with gate types specified by `gate_types` and the qubits upon which they act specified by `ranges`, where the layer acts on `n` qubits.
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

"""
    pad_layer(l::Layer)

Returns a copy of the layer `l` padded by single-qubit identity gates that act on each of the qubits not already acted upon by some gate in the layer.
"""
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

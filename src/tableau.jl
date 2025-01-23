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
        print(
            io,
            get_pauli_string(Pauli(t.tableau[t.qubit_num + i, :], t.qubit_num)) *
            (i == t.qubitnum ? "" : "\n"),
        )
    end
    return nothing
end

@struct_hash_equal_isequal Tableau

"""
    Gate

A gate in a stabiliser circuit.

# Fields

  - `type::String`: String describing the gate type.
  - `index::Int32`: The index labelling the unique layer occurrences of the gate in a circuit.
  - `targets::Vector{Int16}`: The qubit target or targets of the gate.

# Supported gates

  - `H`: Hadamard gate.
  - `S`: Phase gate.
  - `S_DAG`: Conjugate phase gate.
  - `CX` or `CNOT`: Controlled-X gate; the first qubit is the control and the second qubit is the target.
  - `CZ`: Controlled-Z gate.
  - `I`: Identity gate.
  - `Z`: Pauli Z gate.
  - `X`: Pauli X gate.
  - `Y`: Pauli Y gate.
  - `II`: Two-qubit identity gate.
  - `AB`: Two-qubit Pauli gate, where `A` and `B` are Paulis `Z`, `X`, or `Y`.
  - `SQRT_AB`: Two-qubit Pauli rotation, where `A` and `B` are Paulis `Z`, `X`, or `Y`.
  - `SQRT_AB_DAG` : Two-qubit Pauli rotation, where `A` and `B` are Paulis `Z`, `X`, or `Y`.
  - `PZ`: Prepare the Pauli Z eigenstate.
  - `PX`: Prepare the Pauli X eigenstate.
  - `PY`: Prepare the Pauli Y eigenstate.
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

function Base.show(io::IO, g::Gate)
    print(io, "($(g.type):$(Int(g.index)):$(Int.(g.targets)))")
    return nothing
end

function Base.isless(g_1::Gate, g_2::Gate)
    return isless([g_1.type; g_1.index; g_1.targets], [g_2.type; g_2.index; g_2.targets])
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

function Base.show(io::IO, l::Layer)
    show(io, [gate for gate in l.layer])
    return nothing
end

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
    row_phase(x_1::Bool, z_1::Bool, x_2::Bool, z_2::Bool)

Calculate a phase for [`row_sum!`](@ref).
"""
function row_phase(x_1::Bool, z_1::Bool, x_2::Bool, z_2::Bool)
    # Returns the exponent to which i is raised if we multiply x_1z_1 by x_2z_2
    if x_1 == 0 && z_1 == 0
        g = 0
    elseif x_1 == 0 && z_1 == 1
        g = x_2 * (1 - 2z_2)
    elseif x_1 == 1 && z_1 == 0
        g = z_2 * (2x_2 - 1)
    elseif x_1 == 1 && z_1 == 1
        g = z_2 - x_2
    else
        throw(error("How are the Booleans $(x_1), $(z_1), $(x_2), $(z_2) so messed up?"))
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
    if measurement
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
    measurements = Vector{Tuple{Int, Int}}()
    for g in l.layer
        # Perform the relevant operation
        if g.type == "CX" || g.type == "CNOT"
            @assert length(g.targets) == 2 "The gate $(g) must act on two qubits."
            cx!(t, g.targets[1], g.targets[2])
        elseif g.type == "CZ"
            @assert length(g.targets) == 2 "The gate $(g) must act on two qubits."
            cz!(t, g.targets[1], g.targets[2])
        elseif g.type == "H"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            hadamard!(t, g.targets[1])
        elseif g.type == "S"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            phase!(t, g.targets[1])
        elseif g.type == "S_DAG"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            phase!(t, g.targets[1])
            z!(t, g.targets[1])
        elseif g.type == "I" || g.type == "IM"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
        elseif g.type == "X"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            x!(t, g.targets[1])
        elseif g.type == "Z"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            z!(t, g.targets[1])
        elseif g.type == "Y"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            y!(t, g.targets[1])
        elseif g.type == "PZ"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
        elseif g.type == "PX"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            hadamard!(t, g.targets[1])
        elseif g.type == "PY"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            hadamard!(t, g.targets[1])
            phase!(t, g.targets[1])
        elseif g.type == "MZ" || g.type == "M"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            measurement = measure!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "MX"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            hadamard!(t, g.targets[1])
            measurement = measure!(t, g.targets[1])
            hadamard!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "MY"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            phase!(t, g.targets[1])
            z!(t, g.targets[1])
            hadamard!(t, g.targets[1])
            measurement = measure!(t, g.targets[1])
            hadamard!(t, g.targets[1])
            phase!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "R"
            @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
            reset!(t, g.targets[1])
        elseif g.type == "II"
            @assert length(g.targets) == 2 "The gate $(g) must act on two qubits."
        elseif g.type ∈ pauli_rot
            @assert length(g.targets) == 2 "The gate $(g) must act on two qubits."
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
            @assert length(g.targets) == 2 "The gate $(g) must act on two qubits."
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
    make_layer(gate_type::String, range_set::Vector{Vector{Int}}, n::Integer; index::Integer = 0)

Returns a layer of `gate_type` gates, each acting on the qubits in `range_set`, where the layer acts on `n` qubits, optionally specifying the gate index `index`.
"""
function make_layer(
    gate_type::String,
    range_set::Vector{Vector{Int}},
    n::Integer;
    index::Integer = 0,
)
    l = Layer([Gate(gate_type, index, range) for range in range_set], n)
    return l::Layer
end

"""
    make_layer(gate_type::String, range::Vector{Int}, n::Integer; index::Integer = 0)

Returns a layer of single-qubit `gate_type` gates acting on the qubits in `range`, where the layer acts on `n` qubits, optionally specifying the gate index `index`.
"""
function make_layer(gate_type::String, range::Vector{Int}, n::Integer; index::Integer = 0)
    l = Layer([Gate(gate_type, index, [qubit]) for qubit in range], n)
    return l::Layer
end

"""
    make_layer(gate_types::Vector{String}, ranges::Vector{Vector{Int}}, n::Integer; index::Integer = 0)

Returns a layer of single-qubit gates, with gate types specified by `gate_types` and the qubits upon which they act specified by `ranges`, where the layer acts on `n` qubits, optionally specifying the gate index `index`.
"""
function make_layer(
    gate_types::Vector{String},
    ranges::Vector{Vector{Int}},
    n::Integer;
    index::Integer = 0,
)
    @assert length(gate_types) == length(ranges) "The number of gate types and ranges must be the same."
    @assert length(vcat(ranges...)) == length(unique(vcat(ranges...))) "All elements in the ranges must be unique."
    gate_num = length(gate_types)
    l = Layer(
        vcat(
            [
                [Gate(gate_types[i], index, [range]) for range in ranges[i]] for
                i in 1:gate_num
            ]...,
        ),
        n,
    )
    return l::Layer
end

"""
    get_one_qubit_gates()

Returns a list of the supported single-qubit gate types.
"""
function get_one_qubit_gates()
    one_qubit_gates = ["I"; "X"; "Y"; "Z"; "H"; "S"; "S_DAG"; "IM"; "M"; "R"]
    return one_qubit_gates::Vector{String}
end

"""
    get_two_qubit_gates(; stim_supported::Bool = false)

Returns a list of the supported two-qubit gate types, restricting to those also supported by Stim if `stim_supported` is `true`.
"""
function get_two_qubit_gates(; stim_supported::Bool = false)
    if stim_supported
        pauli_rot = ["XX"; "ZZ"; "YY"]
    else
        pauli_rot = ["XX"; "XZ"; "XY"; "ZX"; "ZZ"; "ZY"; "YX"; "YZ"; "YY"]
    end
    sqrt_rot = ["SQRT_" * pauli for pauli in pauli_rot]
    sqrt_rot_dag = ["SQRT_" * pauli * "_DAG" for pauli in pauli_rot]
    if stim_supported
        two_qubit_gates = [pauli_rot...; "CX"; "CNOT"; "CZ"; sqrt_rot...; sqrt_rot_dag...]
    else
        two_qubit_gates =
            ["II"; pauli_rot...; "CX"; "CNOT"; "CZ"; sqrt_rot...; sqrt_rot_dag...]
    end
    return two_qubit_gates::Vector{String}
end

"""
    is_state_prep(g::Gate)

Returns `true` if the gate `g` is a state preparation gate, and `false` otherwise.
"""
function is_state_prep(g::Gate)
    is_state_prep_bool = g.type ∈ ["PZ"; "PX"; "PY"]
    if is_state_prep_bool
        @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
    end
    return is_state_prep_bool::Bool
end

"""
    is_state_meas(g::Gate)

Returns `true` if the gate `g` is a state measurement gate, and `false` otherwise.
"""
function is_state_meas(g::Gate)
    is_state_meas_bool = g.type ∈ ["MZ"; "MX"; "MY"]
    if is_state_meas_bool
        @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
    end
    return is_state_meas_bool::Bool
end

"""
    is_spam(g::Gate)

Returns `true` if the gate `g` is a state preparation or measurement gate, and `false` otherwise.
"""
function is_spam(g::Gate)
    return is_state_prep(g) || is_state_meas(g)
end

"""
    is_mid_meas(g::Gate)

Returns `true` if the gate `g` is a mid-circuit measurement gate, and `false` otherwise.
"""
function is_mid_meas(g::Gate)
    is_mid_meas_bool = g.type == "M"
    if is_mid_meas_bool
        @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
    end
    return is_mid_meas_bool::Bool
end

"""
    is_mid_reset(g::Gate)

Returns `true` if the gate `g` is a mid-circuit reset gate, and `false` otherwise.
"""
function is_mid_reset(g::Gate)
    is_mid_reset_bool = g.type == "R"
    if is_mid_reset_bool
        @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
    end
    return is_mid_reset_bool::Bool
end

"""
    is_mid_meas_reset(g::Gate)

Returns `true` if the gate `g` is a mid-circuit measurement or reset gate, and `false` otherwise.
"""
function is_mid_meas_reset(g::Gate)
    is_mid_meas_reset_bool = is_mid_meas(g) || is_mid_reset(g)
    return is_mid_meas_reset_bool::Bool
end

"""
    is_meas_idle(g::Gate)

Returns `true` if the gate `g` is a measurement idle gate, and `false` otherwise.
"""
function is_meas_idle(g::Gate)
    is_meas_idle_bool = g.type == "IM"
    if is_meas_idle_bool
        @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
    end
    return is_meas_idle_bool::Bool
end

"""
    is_additive(g::Gate; strict::Bool = false)

Returns `true` if the noise associated with the gate `g` can only be estimated to additive precision, and `false` otherwise.
If `strict` is `true`, this is `true` only for state preparation and measurement noise, as is proper, otherwise it is also `true` for mid-circuit measurement and reset noise, and measurement idle noise.
"""
function is_additive(g::Gate; strict::Bool = false)
    if strict
        is_additive_bool = is_spam(g) || is_mid_reset(g)
    else
        is_additive_bool = is_spam(g) || is_mid_meas_reset(g) || is_meas_idle(g)
    end
    return is_additive_bool::Bool
end

"""
    is_pauli(g::Gate)

Returns `true` if the gate `g` is a Pauli gate, and `false` otherwise.
"""
function is_pauli(g::Gate)
    is_pauli_bool = g.type ∈ ["I"; "X"; "Y"; "Z"]
    if is_pauli_bool
        @assert length(g.targets) == 1 "The Pauli gate $(g) must act on one qubit."
    end
    return is_pauli_bool::Bool
end

"""
    is_one_qubit(g::Gate)

Returns `true` if the gate `g` is a supported single-qubit gate, and `false` otherwise.
"""
function is_one_qubit(g::Gate)
    is_one_qubit_bool = g.type ∈ get_one_qubit_gates()
    if is_one_qubit_bool
        @assert length(g.targets) == 1 "The gate $(g) must act on one qubit."
    end
    return is_one_qubit_bool::Bool
end

"""
    is_two_qubit(g::Gate; stim_supported::Bool = false)

Returns `true` if the gate `g` is a supported two-qubit gate, and `false` otherwise.
If `stim_supported` is `true`, restrict to those gates also supported by Stim.
"""
function is_two_qubit(g::Gate; stim_supported::Bool = false)
    is_two_qubit_bool = g.type ∈ get_two_qubit_gates(; stim_supported = stim_supported)
    if is_two_qubit_bool
        @assert length(g.targets) == 2 "The gate $(g) must act on two qubits."
    end
    return is_two_qubit_bool::Bool
end

"""
    pad_layer(l::Layer)

Returns a copy of the layer `l` padded by single-qubit identity gates that act on each of the qubits not already acted upon by some gate in the layer.
"""
function pad_layer(l::Layer)
    # If mid-circuit measurement or reset, pad with IM gates rather than I gates
    any_spam = any(is_mid_meas_reset(gate) for gate in l.layer)
    if any_spam
        id_type = "IM"
    else
        id_type = "I"
    end
    # Construct the padded layer
    target_set = sort(vcat([gate.targets for gate in l.layer]...))
    complement_set = setdiff(collect(1:(l.qubit_num)), target_set)
    layer_padded = l.layer
    for qubit in complement_set
        push!(layer_padded, Gate(id_type, 0, [qubit]))
    end
    l_padded = Layer(layer_padded, l.qubit_num)
    return l_padded::Layer
end

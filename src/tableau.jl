# Stabiliser circuit simulation following `Improved Simulation of Stabilizer Circuits` by S. Aaronson and D. Gottesman (2004)
# arXiv:quant-ph/0406196

"""
    CX!(t::Tableau,control::Integer,target::Integer)

Perform a CNOT on the tableau with the specified control and target qubits.
"""
function CX!(t::Tableau, control::Integer, target::Integer)
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
    Hadamard!(t::Tableau,target::Integer)

Perform a Hadamard gate on the tableau with the specified target qubit.
"""
function Hadamard!(t::Tableau, target::Integer)
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
    Phase!(t::Tableau,target::Integer)

Perform a phase gate on the tableau with the specified target qubit.
"""
function Phase!(t::Tableau, target::Integer)
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
    X!(t::Tableau,target::Integer)

Perform a Pauli X gate on the tableau with the specified target qubit, noting that ``X=HS^2H``.
"""
function X!(t::Tableau, target::Integer)
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
    Z!(t::Tableau,target::Integer)

Perform a Pauli Z gate on the tableau with the specified target qubit, noting that ``Z=S^2``.
"""
function Z!(t::Tableau, target::Integer)
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
    Y!(t::Tableau,target::Integer)

Perform a Pauli Y gate on the tableau with the specified target qubit, noting that ``Y=SHS^2HS^3``.
"""
function Y!(t::Tableau, target::Integer)
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
function CZ!(t::Tableau, control::Integer, target::Integer)
    Hadamard!(t, target)
    CX!(t, control, target)
    Hadamard!(t, target)
    return nothing
end

#
function SQRT_ZZ!(t::Tableau, control::Integer, target::Integer)
    CZ!(t, control, target)
    Phase!(t, control)
    Phase!(t, target)
    return nothing
end

#
function SQRT_ZZ_DAG!(t::Tableau, control::Integer, target::Integer)
    SQRT_ZZ!(t, control, target)
    Z!(t, control)
    Z!(t, target)
    return nothing
end

"""
    RowPhase(x₁::Bool,z₁::Bool,x₂::Bool,z₂::Bool)

Calculate a phase for RowSum.
"""
function RowPhase(x₁::Bool, z₁::Bool, x₂::Bool, z₂::Bool)
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
    RowSum!(t::Tableau,target::Integer,control::Integer)

Add the control row to the target row, keeping track of the phase bit of the target row.
"""
function RowSum!(t::Tableau, target::Integer, control::Integer)
    tableau = t.tableau
    n = t.qubit_num
    @assert target != control "The control and target cannot both be the same qubit $(target)."
    # Adds the row `control` (i in Aaronson 2004) to the row `target` (h in Aaronson 2004), keeping track of the phase bit of the row `target`
    row_sum =
        2 * tableau[target, 2n + 1] +
        2 * tableau[control, 2n + 1] +
        sum(
            RowPhase(
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
    Measure!(t::Tableau,target::Integer)

Perform a measurement on the tableau with the specified target qubit, and return the measurement outcome.
"""
function Measure!(t::Tableau, target::Integer)
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
                RowSum!(t, 2n + 1, n + i)
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
                RowSum!(t, i, p)
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
    Reset!(t::Tableau,target::Integer)

Reset the specified target qubit into the ``|0\\rangle`` state by measuring in the Z basis and flipping it if it is in the ``|1\\rangle`` state.
"""
function Reset!(t::Tableau, target::Integer)
    measurement = Measure!(t, target)
    if measurement == true
        X!(t, target)
    end
    return nothing
end

"""
    Apply!(t::Tableau, l::Layer; return_measurements::Bool = false)

Perform all the gates in the layer on the tableau and return the list of measurement outcomes if return_measurements.
"""
function Apply!(t::Tableau, l::Layer; return_measurements::Bool = false)
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
            CX!(t, g.targets[1], g.targets[2])
        elseif g.type == "CZ"
            CZ!(t, g.targets[1], g.targets[2])
        elseif g.type == "H"
            Hadamard!(t, g.targets[1])
        elseif g.type == "S"
            Phase!(t, g.targets[1])
        elseif g.type == "S_DAG"
            Phase!(t, g.targets[1])
            Z!(t, g.targets[1])
        elseif g.type == "I"
            continue
        elseif g.type == "X"
            X!(t, g.targets[1])
        elseif g.type == "Z"
            Z!(t, g.targets[1])
        elseif g.type == "Y"
            Y!(t, g.targets[1])
        elseif g.type == "PZ+"
        elseif g.type == "PZ-"
            X!(t, g.targets[1])
        elseif g.type == "PX+"
            Hadamard!(t, g.targets[1])
        elseif g.type == "PX-"
            X!(t, g.targets[1])
            Hadamard!(t, g.targets[1])
        elseif g.type == "PY+"
            Hadamard!(t, g.targets[1])
            Phase!(t, g.targets[1])
        elseif g.type == "PY-"
            X!(t, g.targets[1])
            Hadamard!(t, g.targets[1])
            Phase!(t, g.targets[1])
        elseif g.type == "II"
            continue
        elseif g.type ∈ sqrt_rot || g.type ∈ sqrt_rot_dag
            # Set up variables
            pauli_1 = g.type[6]
            pauli_2 = g.type[7]
            # Rotate into the correct Pauli basis
            if pauli_1 == 'X'
                Hadamard!(t, g.targets[1])
            elseif pauli_1 == 'Y'
                Phase!(t, g.targets[1])
                Z!(t, g.targets[1])
                Hadamard!(t, g.targets[1])
            elseif pauli_1 == 'Z'
            else
                throw(error("There's a problem with $(g)."))
            end
            if pauli_2 == 'X'
                Hadamard!(t, g.targets[2])
            elseif pauli_2 == 'Y'
                Phase!(t, g.targets[2])
                Z!(t, g.targets[2])
                Hadamard!(t, g.targets[2])
            elseif pauli_2 == 'Z'
            else
                throw(error("There's a problem with $(g)."))
            end
            # Apply the Pauli rotation
            if g.type ∈ sqrt_rot && g.type ∉ sqrt_rot_dag
                SQRT_ZZ!(t, g.targets[1], g.targets[2])
            elseif g.type ∈ sqrt_rot_dag && g.type ∉ sqrt_rot
                SQRT_ZZ_DAG!(t, g.targets[1], g.targets[2])
            else
                throw(error("There's a problem with $(g)."))
            end
            # Rotate back into the correct Pauli basis
            if pauli_1 == 'X'
                Hadamard!(t, g.targets[1])
            elseif pauli_1 == 'Y'
                Hadamard!(t, g.targets[1])
                Phase!(t, g.targets[1])
            elseif pauli_1 == 'Z'
            else
                throw(error("There's a problem with $(g)."))
            end
            if pauli_2 == 'X'
                Hadamard!(t, g.targets[2])
            elseif pauli_2 == 'Y'
                Hadamard!(t, g.targets[2])
                Phase!(t, g.targets[2])
            elseif pauli_2 == 'Z'
            else
                throw(error("There's a problem with $(g)."))
            end
        elseif g.type == "XX"
            X!(t, g.targets[1])
            X!(t, g.targets[2])
        elseif g.type == "XZ"
            X!(t, g.targets[1])
            Z!(t, g.targets[2])
        elseif g.type == "XY"
            X!(t, g.targets[1])
            Y!(t, g.targets[2])
        elseif g.type == "ZX"
            Z!(t, g.targets[1])
            X!(t, g.targets[2])
        elseif g.type == "ZZ"
            Z!(t, g.targets[1])
            Z!(t, g.targets[2])
        elseif g.type == "ZY"
            Z!(t, g.targets[1])
            Y!(t, g.targets[2])
        elseif g.type == "YX"
            Y!(t, g.targets[1])
            X!(t, g.targets[2])
        elseif g.type == "YZ"
            Y!(t, g.targets[1])
            Z!(t, g.targets[2])
        elseif g.type == "YY"
            Y!(t, g.targets[1])
            Y!(t, g.targets[2])
        elseif g.type == "MZ" || g.type == "M"
            measurement = Measure!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "MX"
            Hadamard!(t, g.targets[1])
            measurement = Measure!(t, g.targets[1])
            Hadamard!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "MY"
            Phase!(t, g.targets[1])
            Z!(t, g.targets[1])
            Hadamard!(t, g.targets[1])
            measurement = Measure!(t, g.targets[1])
            Hadamard!(t, g.targets[1])
            Phase!(t, g.targets[1])
            push!(measurements, ((-1)^measurement, g.targets[1]))
        elseif g.type == "R"
            Reset!(t, g.targets[1])
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

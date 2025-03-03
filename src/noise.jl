"""
    GateIndex

An object that stores all of the different indices of a gate.

# Fields

  - `gate::Gate`: The gate.
  - `indices::Vector{Int32}`: Gate eigenvalue or probability indices.
  - `pad_indices::Vector{Int32}`: Padded gate eigenvalue or probability indices.
  - `marg_indices::Vector{Int32}`: Marginal gate eigenvalue or probability indices.
  - `pad_marg_indices::Vector{Int32}`: Padded marginal gate eigenvalue or probability indices.
  - `rel_indices::Vector{Int32}`: Relative gate eigenvalue or probability indices.
  - `pad_rel_indices::Vector{Int32}`: Padded relative gate eigenvalue or probability indices.
"""
struct GateIndex
    gate::Gate
    indices::Vector{Int32}
    pad_indices::Vector{Int32}
    marg_indices::Vector{Int32}
    pad_marg_indices::Vector{Int32}
    rel_indices::Vector{Int32}
    pad_rel_indices::Vector{Int32}
end

@struct_hash_equal_isequal GateIndex

"""
    GateData

An object that stores all of the gate indices and data.

# Fields

  - `gate_indices::Vector{GateIndex}`: Indices for all gates.
  - `G::Int32`: Number of gates.
  - `N::Int32`: Number of indices.
  - `N_pad::Int32`: Number of padded indices.
  - `N_marginal::Int32`: Number of marginal indices.
  - `N_pad_marginal::Int32`: Number of padded marginal indices.
  - `N_relative::Int32`: Number of relative indices.
  - `N_pad_relative::Int32`: Number of padded relative indices.
  - `combined::Bool`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
  - `strict::Bool`: Whether to be strict about whether gates are deemed to be estimable to additive or relative precision.
"""
struct GateData
    gate_indices::Vector{GateIndex}
    G::Int32
    N::Int32
    N_pad::Int32
    N_marginal::Int32
    N_pad_marginal::Int32
    N_relative::Int32
    N_pad_relative::Int32
    combined::Bool
    strict::Bool
end

@struct_hash_equal_isequal GateData

"""
    get_orbit_indices_dict()

Returns a hard-coded dictionary of the Pauli orbit indices for each gate type, tested by comparing the results to [`get_orbit_indices`](@ref).
"""
function get_orbit_indices_dict()
    orbit_indices_dict = Dict{String, Vector{Vector{Int}}}()
    # Single-qubit gates
    one_qubit_trivial = [[1], [2], [3]]
    orbit_indices_dict["I"] = one_qubit_trivial
    orbit_indices_dict["IM"] = one_qubit_trivial
    orbit_indices_dict["X"] = one_qubit_trivial
    orbit_indices_dict["Y"] = one_qubit_trivial
    orbit_indices_dict["Z"] = one_qubit_trivial
    orbit_indices_dict["H"] = [[1; 2], [3]]
    orbit_indices_dict["S"] = [[1; 3], [2]]
    orbit_indices_dict["S_DAG"] = orbit_indices_dict["S"]
    # Two-qubit gates
    two_qubit_trivial =
        [[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13], [14], [15]]
    orbit_indices_dict["II"] = two_qubit_trivial
    orbit_indices_dict["XX"] = two_qubit_trivial
    orbit_indices_dict["XY"] = two_qubit_trivial
    orbit_indices_dict["XZ"] = two_qubit_trivial
    orbit_indices_dict["YX"] = two_qubit_trivial
    orbit_indices_dict["YY"] = two_qubit_trivial
    orbit_indices_dict["YZ"] = two_qubit_trivial
    orbit_indices_dict["ZX"] = two_qubit_trivial
    orbit_indices_dict["ZY"] = two_qubit_trivial
    orbit_indices_dict["ZZ"] = two_qubit_trivial
    orbit_indices_dict["CX"] =
        [[1; 3], [2], [4], [5; 7], [6], [8; 12], [9; 15], [10; 14], [11; 13]]
    orbit_indices_dict["CNOT"] = orbit_indices_dict["CX"]
    orbit_indices_dict["CZ"] =
        [[1; 9], [2; 6], [3; 15], [4], [5; 13], [7; 11], [8], [10; 14], [12]]
    orbit_indices_dict["SQRT_XX"] =
        [[1], [2], [3], [4; 7], [5; 6], [8; 11], [9; 10], [12], [13], [14], [15]]
    orbit_indices_dict["SQRT_XX_DAG"] = orbit_indices_dict["SQRT_XX"]
    orbit_indices_dict["SQRT_XY"] =
        [[1], [2; 9], [3; 8], [4; 15], [5; 14], [6], [7], [10], [11], [12], [13]]
    orbit_indices_dict["SQRT_XY_DAG"] = orbit_indices_dict["SQRT_XY"]
    orbit_indices_dict["SQRT_XZ"] =
        [[1], [2; 11], [3; 10], [4; 13], [5; 12], [6], [7], [8], [9], [14], [15]]
    orbit_indices_dict["SQRT_XZ_DAG"] = orbit_indices_dict["SQRT_XZ"]
    orbit_indices_dict["SQRT_YX"] =
        [[1; 6], [2], [3; 4], [5], [7], [8; 15], [9], [10; 13], [11], [12], [14]]
    orbit_indices_dict["SQRT_YX_DAG"] = orbit_indices_dict["SQRT_YX"]
    orbit_indices_dict["SQRT_YY"] =
        [[1; 14], [2; 13], [3], [4; 11], [5], [6], [7; 8], [9], [10], [12], [15]]
    orbit_indices_dict["SQRT_YY_DAG"] = orbit_indices_dict["SQRT_YY"]
    orbit_indices_dict["SQRT_YZ"] =
        [[1; 12], [2; 15], [3], [4; 9], [5], [6], [7; 10], [8], [11], [13], [14]]
    orbit_indices_dict["SQRT_YZ_DAG"] = orbit_indices_dict["SQRT_YZ"]
    orbit_indices_dict["SQRT_ZX"] =
        [[1; 7], [2], [3; 5], [4], [6], [8; 14], [9], [10; 12], [11], [13], [15]]
    orbit_indices_dict["SQRT_ZX_DAG"] = orbit_indices_dict["SQRT_ZX"]
    orbit_indices_dict["SQRT_ZY"] =
        [[1; 15], [2; 12], [3], [4], [5; 11], [6; 8], [7], [9], [10], [13], [14]]
    orbit_indices_dict["SQRT_ZY_DAG"] = orbit_indices_dict["SQRT_ZY"]
    orbit_indices_dict["SQRT_ZZ"] =
        [[1; 13], [2; 14], [3], [4], [5; 9], [6; 10], [7], [8], [11], [12], [15]]
    orbit_indices_dict["SQRT_ZZ_DAG"] = orbit_indices_dict["SQRT_ZZ"]
    return orbit_indices_dict::Dict{String, Vector{Vector{Int}}}
end

"""
    get_gate_data(total_gates::Vector{Gate}; combined::Bool = false, strict::Bool = false)

Returns the gate data for the gates `total_gates` in the form of a [`GateData`](@ref) object, combining Pauli X, Y, and Z basis SPAM noise if `combined` is `true`, and being strict about which gates count as estimable to additive or relative precision if `strict` is `true`.
"""
function get_gate_data(
    total_gates::Vector{Gate};
    combined::Bool = false,
    strict::Bool = false,
)
    # Determine the gate indices
    orbit_indices_dict = get_orbit_indices_dict()
    G = length(total_gates)
    N = 0
    N_pad = 0
    N_marginal = 0
    N_pad_marginal = 0
    N_relative = 0
    N_pad_relative = 0
    if combined
        spam_gates = Gate[]
    end
    # Calculate the gate indices
    gate_indices = Vector{GateIndex}(undef, G)
    for (idx, gate) in pairs(total_gates)
        if is_spam(gate) && combined
            if gate ∉ spam_gates
                # Generate the SPAM gates
                spam_z = Gate(gate.type[1] * "Z", 0, gate.targets)
                spam_x = Gate(gate.type[1] * "X", 0, gate.targets)
                spam_y = Gate(gate.type[1] * "Y", 0, gate.targets)
                # Check the ordering is as expected
                @assert gate == spam_z &&
                        total_gates[idx + 1] == spam_x &&
                        total_gates[idx + 2] == spam_y "The SPAM gates $(spam_z), $(spam_x), and $(spam_y) are not in the expected order."
                # Calculate the indices up to the N offsets
                indices = [1]
                pad_indices = [1; 2]
                marg_indices = [1]
                pad_marg_indices = [1; 2]
                rel_indices = []
                pad_rel_indices = []
                # Create the gate indices
                spam_z_idx = GateIndex(
                    spam_z,
                    N .+ indices,
                    N_pad .+ pad_indices,
                    N_marginal .+ marg_indices,
                    N_pad_marginal .+ pad_marg_indices,
                    N_relative .+ rel_indices,
                    N_pad_relative .+ pad_rel_indices,
                )
                spam_x_idx = GateIndex(
                    spam_x,
                    N .+ indices,
                    N_pad .+ pad_indices,
                    N_marginal .+ marg_indices,
                    N_pad_marginal .+ pad_marg_indices,
                    N_relative .+ rel_indices,
                    N_pad_relative .+ pad_rel_indices,
                )
                spam_y_idx = GateIndex(
                    spam_y,
                    N .+ indices,
                    N_pad .+ pad_indices,
                    N_marginal .+ marg_indices,
                    N_pad_marginal .+ pad_marg_indices,
                    N_relative .+ rel_indices,
                    N_pad_relative .+ pad_rel_indices,
                )
                # Store the gate indices and update the N offsets
                gate_indices[idx] = spam_z_idx
                gate_indices[idx + 1] = spam_x_idx
                gate_indices[idx + 2] = spam_y_idx
                append!(spam_gates, [spam_z; spam_x; spam_y])
                N += length(indices)
                N_pad += length(pad_indices)
                N_marginal += length(marg_indices)
                N_pad_marginal += length(pad_marg_indices)
                N_relative += length(rel_indices)
                N_pad_relative += length(pad_rel_indices)
            end
        else
            # Calculate the indices up to the N offsets
            if is_spam(gate) || is_mid_meas_reset(gate)
                indices = [1]
                pad_indices = [1; 2]
                marg_indices = [1]
                pad_marg_indices = [1; 2]
            else
                gate_support_size = length(gate.targets)
                gate_orbits = orbit_indices_dict[gate.type]
                indices = collect(1:(4^gate_support_size - 1))
                pad_indices = collect(1:(4^gate_support_size))
                marg_indices = collect(1:length(gate_orbits))
                pad_marg_indices = collect(1:(length(gate_orbits) + 1))
            end
            if ~is_additive(gate; strict = strict)
                rel_indices = marg_indices
                pad_rel_indices = pad_marg_indices
            else
                rel_indices = []
                pad_rel_indices = []
            end
            # Create the gate index
            gate_idx = GateIndex(
                gate,
                N .+ indices,
                N_pad .+ pad_indices,
                N_marginal .+ marg_indices,
                N_pad_marginal .+ pad_marg_indices,
                N_relative .+ rel_indices,
                N_pad_relative .+ pad_rel_indices,
            )
            # Store the gate index and update the N offsets
            gate_indices[idx] = gate_idx
            N += length(indices)
            N_pad += length(pad_indices)
            N_marginal += length(marg_indices)
            N_pad_marginal += length(pad_marg_indices)
            N_relative += length(rel_indices)
            N_pad_relative += length(pad_rel_indices)
        end
    end
    gate_data = GateData(
        gate_indices,
        G,
        N,
        N_pad,
        N_marginal,
        N_pad_marginal,
        N_relative,
        N_pad_relative,
        combined,
        strict,
    )
    @assert [gate_idx.gate for gate_idx in gate_indices] == total_gates
    return gate_data::GateData
end

"""
    get_gate_index_dict(gate_data::GateData)

Returns a dictionary that maps each gate to its index.
"""
function get_gate_index_dict(gate_data::GateData)
    gate_index_dict = Dict{Gate, Int}()
    for gate_idx in gate_data.gate_indices
        gate_index_dict[gate_idx.gate] = gate_idx.indices[1] - 1
    end
    return gate_index_dict::Dict{Gate, Int}
end

"""
    get_wht_matrix(n::Integer; inverse::Bool = false)
    get_wht_matrix(gate::Gate; inverse::Bool = false)

Returns the symplectically ordered Walsh-Hadamard transform matrix of order `n`, the number of qubits on which the gate `gate` acts, which maps an n-qubit Pauli error probability distribution to its eigenvalues, or the inverse transform if `inverse` is `true`.
"""
function get_wht_matrix(n::Integer; inverse::Bool = false)
    a = BitArray(undef, 2n)
    b = BitArray(undef, 2n)
    wht_matrix = Matrix{Float64}(undef, 4^n, 4^n)
    for a.chunks[1] in 0:(4^n - 1)
        for b.chunks[1] in 0:(4^n - 1)
            symplectic_form =
                convert(Bool, (a[1:n]' * b[(n + 1):(2n)] + a[(n + 1):(2n)]' * b[1:n]) % 2)
            wht_matrix[a.chunks[1] + 1, b.chunks[1] + 1] = (-1)^symplectic_form
        end
    end
    if inverse
        wht_matrix = wht_matrix ./ 4^n
    end
    return wht_matrix::Matrix{Float64}
end
function get_wht_matrix(gate::Gate; inverse::Bool = false)
    if is_spam(gate) || is_mid_meas_reset(gate)
        wht_matrix = [1.0 1.0; 1.0 -1.0]
        if inverse
            wht_matrix = wht_matrix ./ 2
        end
    else
        gate_support_size = length(gate.targets)
        wht_matrix = get_wht_matrix(gate_support_size; inverse = inverse)
    end
    return wht_matrix::Matrix{Float64}
end

"""
    get_marginal_wht_matrix(gate::Gate; inverse::Bool = false)

Returns the symplectically ordered Walsh-Hadamard transform matrix for the gate `gate`, marginalised over gate orbits, which maps the marginal Pauli error probability distribution to its marginal eigenvalues, or the inverse transform if `inverse` is `true`.
"""
function get_marginal_wht_matrix(gate::Gate; inverse::Bool = false)
    # Initialise variables
    orbit_indices_dict = get_orbit_indices_dict()
    if is_spam(gate) || is_mid_meas_reset(gate)
        wht_matrix = [1.0 1.0; 1.0 -1.0]
        if inverse
            wht_matrix = wht_matrix ./ 2
        end
        marg_wht_matrix = wht_matrix
    else
        # Calculate the marginal Walsh-Hadamard transform matrix
        gate_support_size = length(gate.targets)
        gate_orbits = orbit_indices_dict[gate.type]
        gate_orbit_num = length(gate_orbits) + 1
        wht_matrix = get_wht_matrix(gate_support_size; inverse = inverse)
        # Generate the marginalisation matrix
        marginal_transform = zeros(Float64, gate_orbit_num, 4^gate_support_size)
        marginal_transform[1, 1] = 1
        for (idx, orbit) in pairs(gate_orbits)
            if inverse
                marginal_transform[idx + 1, orbit .+ 1] .= 1
            else
                marginal_transform[idx + 1, orbit .+ 1] .= 1 // length(orbit)
            end
        end
        # Marginalise the Walsh-Hadamard transform matrix
        marg_wht_matrix = marginal_transform * wht_matrix * marginal_transform'
    end
    return marg_wht_matrix::Matrix{Float64}
end

"""
    get_gate_eigenvalues(gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the gate eigenvalues for the gate probabilities `gate_probabilities` calculated using the gate data `gate_data`.
"""
function get_gate_eigenvalues(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    # Calculate the gate eigenvalues
    gate_eigenvalues = zeros(Float64, gate_data.N)
    for idx in 1:(gate_data.G)
        gate = gate_data.gate_indices[idx].gate
        gate_eig_indices = gate_data.gate_indices[idx].indices
        gate_eigs = get_wht_matrix(gate) * gate_probabilities[gate]
        gate_eigenvalues[gate_eig_indices] = gate_eigs[2:end]
    end
    return gate_eigenvalues::Vector{Float64}
end

"""
    get_gate_probabilities(gate_eigenvalues::Vector{Float64}, gate_data::GateData)

Returns the gate probabilities for the gate eigenvalues `gate_eigenvalues` calculated using the gate data `gate_data`.
"""
function get_gate_probabilities(gate_eigenvalues::Vector{Float64}, gate_data::GateData)
    # Calculate the gate probabilities
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for idx in 1:(gate_data.G)
        gate = gate_data.gate_indices[idx].gate
        gate_eigs = [1; gate_eigenvalues[gate_data.gate_indices[idx].indices]]
        gate_probs = get_wht_matrix(gate; inverse = true) * gate_eigs
        gate_probabilities[gate] = gate_probs
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_combined_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the combined gate probabilities obtained from the gate probabilities `gate_probabilities` by averaging SPAM noise parameters on each qubit, calculated using combined gate data `gate_data`.
"""
function get_combined_gate_probabilities(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    # Initialise variables
    @assert gate_data.combined "This function must be supplied with combined gate data."
    total_gates = [gate_idx.gate for gate_idx in gate_data.gate_indices]
    spam_gates = Gate[]
    combined_gate_probabilities = Dict{Gate, Vector{Float64}}()
    for (idx, gate) in pairs(total_gates)
        if is_spam(gate)
            if gate ∉ spam_gates
                # Generate the SPAM gates
                spam_z = Gate(gate.type[1] * "Z", 0, gate.targets)
                spam_x = Gate(gate.type[1] * "X", 0, gate.targets)
                spam_y = Gate(gate.type[1] * "Y", 0, gate.targets)
                # Check the ordering is as expected
                @assert gate == spam_z &&
                        total_gates[idx + 1] == spam_x &&
                        total_gates[idx + 2] == spam_y "The SPAM gates $(spam_z), $(spam_x), and $(spam_y) are not in the expected order."
                # Set the SPAM error probabilities by the computational basis
                spam_probability = gate_probabilities[spam_z]
                # Set the combined SPAM error probabilities
                combined_gate_probabilities[spam_z] = spam_probability
                combined_gate_probabilities[spam_x] = spam_probability
                combined_gate_probabilities[spam_y] = spam_probability
                # Append the SPAM gates
                append!(spam_gates, [spam_z; spam_x; spam_y])
            end
        else
            combined_gate_probabilities[gate] = gate_probabilities[gate]
        end
    end
    return combined_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_average_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}})
    get_average_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the averaged gate probabilities obtained from the gate probabilities `gate_probabilities` by averaging error probabilities within the orbits of each of the gates.
"""
function get_average_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}})
    # Initialise variables
    orbit_indices_dict = get_orbit_indices_dict()
    average_gate_probabilities = Dict{Gate, Vector{Float64}}()
    for (gate, gate_probs) in pairs(gate_probabilities)
        if is_spam(gate) || is_mid_meas_reset(gate)
            # Preparations and measurements are trivial
            average_gate_probabilities[gate] = gate_probabilities[gate]
        else
            # Gates require averaging over the probabilities within each orbit
            gate_orbits = orbit_indices_dict[gate.type]
            avg_gate_probs = Vector{Float64}(undef, length(gate_probs))
            avg_gate_probs[1] = gate_probs[1]
            for orbit in gate_orbits
                avg_gate_probs[orbit .+ 1] .= sum(gate_probs[orbit .+ 1]) / length(orbit)
            end
            average_gate_probabilities[gate] = avg_gate_probs
        end
    end
    return average_gate_probabilities::Dict{Gate, Vector{Float64}}
end
function get_average_gate_probabilities(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    average_gate_probabilities = get_average_gate_probabilities(gate_probabilities)
    return average_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_full_average_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}})
    get_full_average_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the fully averaged gate probabilities obtained from the gate probabilities `gate_probabilities` by averaging all error probabilities within each of the gates.
"""
function get_full_average_gate_probabilities(
    gate_probabilities::Dict{Gate, Vector{Float64}},
)
    # Initialise variables
    full_average_gate_probabilities = Dict{Gate, Vector{Float64}}()
    for (gate, gate_probs) in pairs(gate_probabilities)
        if is_spam(gate) || is_mid_meas_reset(gate)
            # Preparations and measurements are trivial
            full_average_gate_probabilities[gate] = gate_probabilities[gate]
        else
            # Average gates over all error probabilities
            avg_gate_probs = Vector{Float64}(undef, length(gate_probs))
            avg_gate_probs[1] = gate_probs[1]
            avg_gate_probs[2:end] .= sum(gate_probs[2:end]) / length(gate_probs[2:end])
            full_average_gate_probabilities[gate] = avg_gate_probs
        end
    end
    return full_average_gate_probabilities::Dict{Gate, Vector{Float64}}
end
function get_full_average_gate_probabilities(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    full_average_gate_probabilities =
        get_full_average_gate_probabilities(gate_probabilities)
    return full_average_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_marginal_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}})
    get_marginal_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the marginal gate probabilities obtained from the gate probabilities `gate_probabilities` by marginalising error probabilities across the orbits of each of the gates.
"""
function get_marginal_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}})
    # Initialise variables
    orbit_indices_dict = get_orbit_indices_dict()
    marginal_gate_probabilities = Dict{Gate, Vector{Float64}}()
    for (gate, gate_probs) in pairs(gate_probabilities)
        if is_spam(gate) || is_mid_meas_reset(gate)
            # Preparations and measurements are trivial
            marginal_gate_probabilities[gate] = gate_probabilities[gate]
        else
            # Gates require marginalisation of the probabilities within each orbit
            gate_orbits = orbit_indices_dict[gate.type]
            marg_gate_probs = Vector{Float64}(undef, length(gate_orbits) + 1)
            marg_gate_probs[1] = gate_probs[1]
            for (idx, orbit) in pairs(gate_orbits)
                marg_gate_probs[idx + 1] = sum(gate_probs[orbit .+ 1])
            end
            marginal_gate_probabilities[gate] = marg_gate_probs
        end
    end
    return marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
end
function get_marginal_gate_probabilities(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    marginal_gate_probabilities = get_marginal_gate_probabilities(gate_probabilities)
    return marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_relative_gate_probabilities(gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the marginal gate probabilities obtained from the gate probabilities `gate_probabilities` by marginalising error probabilities across the orbits of each of the gates which are estimable to relative precision.
"""
function get_relative_gate_probabilities(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    # Initialise variables
    orbit_indices_dict = get_orbit_indices_dict()
    relative_gate_probabilities = Dict{Gate, Vector{Float64}}()
    for (gate, gate_probs) in pairs(gate_probabilities)
        if ~is_additive(gate; strict = gate_data.strict)
            if is_spam(gate) || is_mid_meas_reset(gate)
                # Preparations and measurements are trivial
                relative_gate_probabilities[gate] = gate_probabilities[gate]
            else
                # Gates require marginalisation of the probabilities within each orbit
                gate_orbits = orbit_indices_dict[gate.type]
                marg_gate_probs = Vector{Float64}(undef, length(gate_orbits) + 1)
                marg_gate_probs[1] = gate_probs[1]
                for (idx, orbit) in pairs(gate_orbits)
                    marg_gate_probs[idx + 1] = sum(gate_probs[orbit .+ 1])
                end
                relative_gate_probabilities[gate] = marg_gate_probs
            end
        end
    end
    return relative_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_gate_probabilities_vec(gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the gate probabilities vector obtained from the corresponding dictionary `gate_probabilities` calculated using the gate data `gate_data`.
"""
function get_gate_probabilities_vec(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    # Convert the gate probabilities
    gate_probabilities_vec = zeros(Float64, gate_data.N_pad)
    for gate_idx in gate_data.gate_indices
        gate_probabilities_vec[gate_idx.pad_indices] = gate_probabilities[gate_idx.gate]
    end
    return gate_probabilities_vec::Vector{Float64}
end

"""
    get_gate_probabilities_dict(gate_probabilities_vec::Vector{Float64}, gate_data::GateData)

Returns the gate probabilities dictionary obtained from the corresponding vector `gate_probabilities_vec` calculated using the gate data `gate_data`.
"""
function get_gate_probabilities_dict(
    gate_probabilities_vec::Vector{Float64},
    gate_data::GateData,
)
    # Convert the gate probabilities
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate_idx in gate_data.gate_indices
        gate_probabilities[gate_idx.gate] = gate_probabilities_vec[gate_idx.pad_indices]
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_marginal_gate_probabilities_vec(marginal_gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the marginal gate probabilities vector obtained from the corresponding dictionary `marginal_gate_probabilities` calculated using the gate data `gate_data`.
"""
function get_marginal_gate_probabilities_vec(
    marginal_gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    # Convert the marginal gate probabilities
    marginal_gate_probabilities_vec = zeros(Float64, gate_data.N_pad_marginal)
    for gate_idx in gate_data.gate_indices
        marginal_gate_probabilities_vec[gate_idx.pad_marg_indices] =
            marginal_gate_probabilities[gate_idx.gate]
    end
    return marginal_gate_probabilities_vec::Vector{Float64}
end

"""
    get_marginal_gate_probabilities_dict(marginal_gate_probabilities_vec::Vector{Float64}, gate_data::GateData)

Returns the marginal gate probabilities dictionary obtained from the corresponding vector `marginal_gate_probabilities_vec` calculated using the gate data `gate_data`.
"""
function get_marginal_gate_probabilities_dict(
    marginal_gate_probabilities_vec::Vector{Float64},
    gate_data::GateData,
)
    # Convert the marginal gate probabilities
    marginal_gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate_idx in gate_data.gate_indices
        marginal_gate_probabilities[gate_idx.gate] =
            marginal_gate_probabilities_vec[gate_idx.pad_marg_indices]
    end
    return marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_relative_gate_probabilities_vec(relative_gate_probabilities::Dict{Gate, Vector{Float64}}, gate_data::GateData)

Returns the marginal gate probabilities vector for those gates estimable to relative precision obtained from the corresponding dictionary `relative_gate_probabilities` calculated using the gate data `gate_data`.
"""
function get_relative_gate_probabilities_vec(
    relative_gate_probabilities::Dict{Gate, Vector{Float64}},
    gate_data::GateData,
)
    # Convert the relative gate probabilities
    relative_gate_probabilities_vec = zeros(Float64, gate_data.N_pad_relative)
    for gate_idx in gate_data.gate_indices
        if length(gate_idx.pad_rel_indices) > 0
            relative_gate_probabilities_vec[gate_idx.pad_rel_indices] =
                relative_gate_probabilities[gate_idx.gate]
        end
    end
    return relative_gate_probabilities_vec::Vector{Float64}
end

"""
    get_relative_gate_probabilities_dict(relative_gate_probabilities_vec::Vector{Float64}, gate_data::GateData)

Returns the marginal gate probabilities dictionary for those gates estimable to relative precision obtained from the corresponding vector `relative_gate_probabilities_vec` calculated using the gate data `gate_data`.
"""
function get_relative_gate_probabilities_dict(
    relative_gate_probabilities_vec::Vector{Float64},
    gate_data::GateData,
)
    # Convert the relative gate probabilities
    relative_gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate_idx in gate_data.gate_indices
        if length(gate_idx.pad_rel_indices) > 0
            relative_gate_probabilities[gate_idx.gate] =
                relative_gate_probabilities_vec[gate_idx.pad_rel_indices]
        end
    end
    return relative_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_ordinary_transform(gate_data::GateData)

Returns an identity transform matrix that maps gate eigenvalues to themselves, calculated using the gate data `gate_data`.
"""
function get_ordinary_transform(gate_data::GateData)
    ordinary_transform = sparse(Diagonal(ones(Float64, gate_data.N)))
    return ordinary_transform::SparseMatrixCSC{Float64, Int}
end

"""
    get_marginal_transform(gate_data::GateData)

Returns a transform matrix that maps gate eigenvalues to marginal gate eigenvalues, calculated using the gate data `gate_data`.
"""
function get_marginal_transform(gate_data::GateData)
    # Initialise variables
    orbit_indices_dict = get_orbit_indices_dict()
    marginal_transform = spzeros(Float64, gate_data.N_marginal, gate_data.N)
    for gate_idx in gate_data.gate_indices
        gate = gate_idx.gate
        idx = gate_idx.indices[1] - 1
        marg_idx = gate_idx.marg_indices[1] - 1
        if is_spam(gate) || is_mid_meas_reset(gate)
            # Preparations and measurements are trivial
            marginal_transform[marg_idx + 1, idx + 1] = 1
        else
            # Gates require marginalisation of the probabilities within each orbit
            gate_orbits = orbit_indices_dict[gate.type]
            for (orbit_idx, orbit) in pairs(gate_orbits)
                marginal_transform[marg_idx + orbit_idx, idx .+ orbit] .= 1 / length(orbit)
            end
        end
    end
    return marginal_transform::SparseMatrixCSC{Float64, Int}
end

"""
    get_relative_transform(gate_data::GateData)

Returns a transform matrix that maps gate eigenvalues to the marginal gate eigenvalues for gates estimable to relative precision, calculated using the gate data `gate_data`.
"""
function get_relative_transform(gate_data::GateData)
    # Initialise variables
    orbit_indices_dict = get_orbit_indices_dict()
    relative_transform = spzeros(Float64, gate_data.N_relative, gate_data.N)
    for gate_idx in gate_data.gate_indices
        if length(gate_idx.rel_indices) > 0
            gate = gate_idx.gate
            idx = gate_idx.indices[1] - 1
            rel_idx = gate_idx.rel_indices[1] - 1
            if is_spam(gate) || is_mid_meas_reset(gate)
                # Preparations and measurements are trivial
                relative_transform[rel_idx + 1, idx + 1] = 1
            else
                # Gates require marginalisation of the probabilities within each orbit
                gate_orbits = orbit_indices_dict[gate.type]
                for (orbit_idx, orbit) in pairs(gate_orbits)
                    relative_transform[rel_idx + orbit_idx, idx .+ orbit] .=
                        1 / length(orbit)
                end
            end
        end
    end
    return relative_transform::SparseMatrixCSC{Float64, Int}
end

"""
    get_marginal_gate_eigenvalues(gate_eigenvalues::Vector{Float64}, gate_data::GateData)

Returns the marginal gate eigenvalues corresponding to the gate eigenvalues `gate_eigenvalues` calculated using the gate data `gate_data`.
"""
function get_marginal_gate_eigenvalues(
    gate_eigenvalues::Vector{Float64},
    gate_data::GateData,
)
    marginal_gate_eigenvalues = get_marginal_transform(gate_data) * gate_eigenvalues
    return marginal_gate_eigenvalues::Vector{Float64}
end

"""
    get_relative_gate_eigenvalues(gate_eigenvalues::Vector{Float64}, gate_data::GateData)

Returns the marginal gate eigenvalues for gates estimable to relative precision corresponding to the gate eigenvalues `gate_eigenvalues` calculated using the gate data `gate_data`.
"""
function get_relative_gate_eigenvalues(
    gate_eigenvalues::Vector{Float64},
    gate_data::GateData,
)
    relative_gate_eigenvalues = get_relative_transform(gate_data) * gate_eigenvalues
    return relative_gate_eigenvalues::Vector{Float64}
end

"""
    get_pad_transform(gate_data::GateData; probabilities::Bool = false)

Returns a transform matrix that pads gate eigenvalues, or gate error probabilities if `probabilities` is `true`, with identity eigenvaleus or error probabilities respectively, up to a constant given by [`get_pad_mask`](@ref), calculated using the gate data `gate_data`.
"""
function get_pad_transform(gate_data::GateData; probabilities::Bool = false)
    # Generate the transform
    pad_transform = spzeros(Float64, gate_data.N_pad, gate_data.N)
    for gate_idx in gate_data.gate_indices
        pad_indices = gate_idx.pad_indices
        indices = gate_idx.indices
        for idx in eachindex(indices)
            pad_transform[pad_indices[idx + 1], indices[idx]] = 1
            if probabilities
                pad_transform[pad_indices[1], indices[idx]] = -1
            end
        end
    end
    return pad_transform::SparseMatrixCSC{Float64, Int}
end

"""
    get_pad_mask(gate_data::GateData)

Returns a mask vector that sets the values of the identity gate eigenvalues or gate error probabilites, calculated using the gate data `gate_data`.
"""
function get_pad_mask(gate_data::GateData)
    # Generate the mask
    pad_mask = zeros(Float64, gate_data.N_pad)
    for gate_idx in gate_data.gate_indices
        pad_mask[gate_idx.pad_indices[1]] = 1
    end
    return pad_mask::Vector{Float64}
end

"""
    get_marginal_pad_transform(gate_data::GateData; probabilities::Bool = false)

Returns a transform matrix that pads marginal gate eigenvalues, or marginal gate error probabilities if `probabilities` is `true`, with identity eigenvaleus or error probabilities respectively, up to a constant given by [`get_marginal_pad_mask`](@ref), calculated using the gate data `gate_data`.
"""
function get_marginal_pad_transform(gate_data::GateData; probabilities::Bool = false)
    # Generate the transform
    marg_pad_transform = spzeros(Float64, gate_data.N_pad_marginal, gate_data.N_marginal)
    for gate_idx in gate_data.gate_indices
        pad_marg_indices = gate_idx.pad_marg_indices
        marg_indices = gate_idx.marg_indices
        for idx in eachindex(marg_indices)
            marg_pad_transform[pad_marg_indices[idx + 1], marg_indices[idx]] = 1
            if probabilities
                marg_pad_transform[pad_marg_indices[1], marg_indices[idx]] = -1
            end
        end
    end
    return marg_pad_transform::SparseMatrixCSC{Float64, Int}
end

"""
    get_marginal_pad_mask(gate_data::GateData)

Returns a mask vector that sets the values of the identity marginal gate eigenvalues or marginal gate error probabilites, calculated using the gate data `gate_data`.
"""
function get_marginal_pad_mask(gate_data::GateData)
    # Generate the mask
    marginal_pad_mask = zeros(Float64, gate_data.N_pad_marginal)
    for gate_idx in gate_data.gate_indices
        marginal_pad_mask[gate_idx.pad_marg_indices[1]] = 1
    end
    return marginal_pad_mask::Vector{Float64}
end

"""
    get_relative_pad_transform(gate_data::GateData; probabilities::Bool = false)

Returns a transform matrix that pads marginal gate eigenvalues, or marginal gate error probabilities if `probabilities` is `true`, for gates estimable to relative precision, with identity eigenvaleus or error probabilities respectively, up to a constant given by [`get_relative_pad_mask`](@ref), calculated using the gate data `gate_data`.
"""
function get_relative_pad_transform(gate_data::GateData; probabilities::Bool = false)
    # Generate the transform
    rel_pad_transform = spzeros(Float64, gate_data.N_pad_relative, gate_data.N_relative)
    for gate_idx in gate_data.gate_indices
        pad_rel_indices = gate_idx.pad_rel_indices
        rel_indices = gate_idx.rel_indices
        if length(pad_rel_indices) > 0
            for idx in eachindex(rel_indices)
                rel_pad_transform[pad_rel_indices[idx + 1], rel_indices[idx]] = 1
                if probabilities
                    rel_pad_transform[pad_rel_indices[1], rel_indices[idx]] = -1
                end
            end
        end
    end
    return rel_pad_transform::SparseMatrixCSC{Float64, Int}
end

"""
    get_relative_pad_mask(gate_data::GateData)

Returns a mask vector that sets the values of the identity marginal gate eigenvalues for gates estimable to relative precision, calculated using the gate data `gate_data`.
"""
function get_relative_pad_mask(gate_data::GateData)
    # Generate the mask
    relative_pad_mask = zeros(Float64, gate_data.N_pad_relative)
    for gate_idx in gate_data.gate_indices
        if length(gate_idx.pad_rel_indices) > 0
            relative_pad_mask[gate_idx.pad_rel_indices[1]] = 1
        end
    end
    return relative_pad_mask::Vector{Float64}
end

"""
    get_wht_transform(gate_data::GateData; inverse::Bool = false)

Returns a transform matrix that maps padded gate error probabilities to padded gate eigenvalues, or the inverse transform if `inverse` is `true`, calculated using the gate data `gate_data`.
"""
function get_wht_transform(gate_data::GateData; inverse::Bool = false)
    # Generate the indices, checking for a combined design
    G = gate_data.G
    covered_indices = Vector{Vector{Int}}()
    reduced_indices = Vector{Int}()
    for idx in 1:G
        pad_indices = gate_data.gate_indices[idx].pad_indices
        if length(pad_indices) > 0 && pad_indices ∉ covered_indices
            push!(reduced_indices, idx)
            push!(covered_indices, pad_indices)
        end
    end
    # Generate the transform
    G_red = length(reduced_indices)
    wht_transform_list = Vector{SparseMatrixCSC{Float64, Int}}(undef, G_red)
    @threads :static for idx in 1:G_red
        wht_transform_list[idx] = sparse(
            get_wht_matrix(
                gate_data.gate_indices[reduced_indices[idx]].gate;
                inverse = inverse,
            ),
        )
    end
    wht_transform = blockdiag(wht_transform_list...)
    return wht_transform::SparseMatrixCSC{Float64, Int}
end

"""
    get_marginal_wht_transform(gate_data::GateData; inverse::Bool = false)

Returns a transform matrix that maps padded marginal gate error probabilities to padded marginal gate eigenvalues, or the inverse transform if `inverse` is `true`, calculated using the gate data `gate_data`.
"""
function get_marginal_wht_transform(gate_data::GateData; inverse::Bool = false)
    # Generate the indices, checking for a combined design
    G = gate_data.G
    covered_indices = Vector{Vector{Int}}()
    reduced_indices = Vector{Int}()
    for idx in 1:G
        pad_marg_indices = gate_data.gate_indices[idx].pad_marg_indices
        if length(pad_marg_indices) > 0 && pad_marg_indices ∉ covered_indices
            push!(reduced_indices, idx)
            push!(covered_indices, pad_marg_indices)
        end
    end
    # Generate the transform
    G_red = length(reduced_indices)
    marginal_wht_transform_list = Vector{SparseMatrixCSC{Float64, Int}}(undef, G_red)
    @threads :static for idx in 1:G_red
        marginal_wht_transform_list[idx] = sparse(
            get_marginal_wht_matrix(
                gate_data.gate_indices[reduced_indices[idx]].gate;
                inverse = inverse,
            ),
        )
    end
    marginal_wht_transform = blockdiag(marginal_wht_transform_list...)
    return marginal_wht_transform::SparseMatrixCSC{Float64, Int}
end

"""
    get_relative_wht_transform(gate_data::GateData; inverse::Bool = false)

Returns a transform matrix that maps padded marginal gate error probabilities to padded marginal gate eigenvalues for gates estimable to relative precision, or the inverse transform if `inverse` is `true`, calculated using the gate data `gate_data`.
"""
function get_relative_wht_transform(gate_data::GateData; inverse::Bool = false)
    # Generate the indices, checking for a combined design
    G = gate_data.G
    covered_indices = Vector{Vector{Int}}()
    reduced_indices = Vector{Int}()
    for idx in 1:G
        pad_rel_indices = gate_data.gate_indices[idx].pad_rel_indices
        if length(pad_rel_indices) > 0 && pad_rel_indices ∉ covered_indices
            push!(reduced_indices, idx)
            push!(covered_indices, pad_rel_indices)
        end
    end
    # Generate the transform
    G_red = length(reduced_indices)
    relative_wht_transform_list = Vector{SparseMatrixCSC{Float64, Int}}(undef, G_red)
    @threads :static for idx in 1:G_red
        relative_wht_transform_list[idx] = sparse(
            get_marginal_wht_matrix(
                gate_data.gate_indices[reduced_indices[idx]].gate;
                inverse = inverse,
            ),
        )
    end
    relative_wht_transform = blockdiag(relative_wht_transform_list...)
    return relative_wht_transform::SparseMatrixCSC{Float64, Int}
end

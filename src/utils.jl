"""
    batch_shots(shots::Integer, measurements::Integer, max_samples::Integer)

Returns the shots divided into batches for sampling from Stim.
"""
function batch_shots(shots::Integer, measurements::Integer, max_samples::Integer)
    # Divide the shots into batches
    # Stim samples in batches of 256, so we want batch sizes that are multiples of 256
    base = 256
    batch_num = cld(measurements * shots, max_samples)
    base_mult = cld(fld(shots, base) + 1, batch_num)
    shot_batches = base_mult * base * ones(Int, batch_num)
    # Remove the excess shots
    excess_shots = base_mult * base * batch_num - shots
    @assert excess_shots >= 0 "The number of excess shots is negative."
    excess_batch_num = fld(excess_shots, base)
    remainder = excess_shots % base
    if remainder > 0
        shot_batches[end - excess_batch_num] -= remainder
    end
    for idx in 0:(excess_batch_num - 1)
        shot_batches[end - idx] -= base
    end
    # Filter and check the shot batches
    shot_batches = filter(x -> x > 0, shot_batches)
    @assert all(shot_batches .> 0) "The batch sizes are not all positive."
    @assert sum(shot_batches) == shots "The number of shots in the batches does not sum to the total number of shots."
    return shot_batches::Vector{Int}
end

"""
    project_simplex(probabilities::Vector{Float64})

Returns a copy of the probability distribution `probabilities` projected into the probability simplex according to the Euclidean norm.
"""
function project_simplex(probabilities::Vector{Float64})
    @assert ~any(isnan.(probabilities)) "The supplied vector contains NaN."
    sorted_probabilities = sort(probabilities; rev = true)
    sum_probabilities =
        [1 / i for i in eachindex(probabilities)] .* (1 .- cumsum(sorted_probabilities))
    projected_probabilities =
        max.(
            probabilities .+
            sum_probabilities[findlast(sorted_probabilities + sum_probabilities .> 0)],
            0,
        )
    @assert all(projected_probabilities .>= 0) "The projected probabilities contain negative values."
    @assert sum(projected_probabilities) ≈ 1.0 "The projected probabilities do not sum to 1."
    return projected_probabilities::Vector{Float64}
end

"""
    scs_project_nonnegative(vector::Vector{Float64}, precision_matrix::SparseMatrixCSC{Float64, Int}; precision::Real = 1e-8, diagnostics::Bool = false)

Returns a copy of the vector `vector` projected into the nonnegative orthant according to the quadratic form defined by the precision matrix `precision_matrix` using the solver SCS, with the specified precision `precision`, and printing diagonstics if `diagnostics` is `true`.
"""
function scs_project_nonnegative(
    vector::Vector{Float64},
    precision_matrix::SparseMatrixCSC{Float64, Int};
    precision::Real = 1e-8,
    diagnostics::Bool = false,
)
    # Prepare variables for projection
    N = length(vector)
    magnitude_scaling = 1 / median(abs.(vector))
    A = sparse(Diagonal(ones(Float64, N)))
    b = zeros(Float64, N)
    P = sparse(UpperTriangular(precision_matrix))
    c = magnitude_scaling * precision_matrix * vector
    # Use the low-level SCS interface, avoiding JuMP overhead
    # SCS with the DirectSolver appears to be most performant solver for this problem
    solution = SCS.scs_solve(
        SCS.DirectSolver,           # Solver type, DirectSolver or IndirectSolver
        N,                          # Number of constraints
        N,                          # Number of variables
        A,                          # Constraint matrix
        P,                          # Quadratic form matrix
        b,                          # Constraint vector
        c,                          # Linear objective vector
        0,                          # Primal equality constraint number
        N,                          # Linear cone number
        Float64[],                  # Box constraint upper bounds
        Float64[],                  # Box constraint lower bounds
        Int[],                      # Second order cone sizes
        Int[],                      # Semidefinite cone sizes
        0,                          # Primal exponential cone number
        0,                          # Dual exponential cone number
        Float64[];                  # Power cone parameters
        eps_abs = precision,        # Absolute solver precision
        eps_rel = precision,        # Relative solver precision
        verbose = diagnostics,      # Print solver output
    )
    # Remove numerical noise
    projected_vector = -solution.x / magnitude_scaling
    projected_vector[projected_vector .< precision] .= 0.0
    return projected_vector::Vector{Float64}
end

"""
    get_paulis(n::Integer)

Returns a list of all n-qubit Paulis ordered lexicographically according to their bit string representation described in [`Pauli`](@ref).
For single-qubit gates, the Pauli error probabilities are ordered as `I`, `X`, `Z`, `Y`.
For two-qubit gates, the Pauli error probabilities are ordered as `II`, `XI`, `IX`, `XX`, `ZI`, `YI`, `ZX`, `YX`, `IZ`, `XZ`, `IY`, `XY`, `ZZ`, `YZ`, `ZY`, `YY`.
"""
function get_paulis(n::Integer)
    bit_pauli = BitArray(undef, 2n + 1)
    pauli_list = Vector{Pauli}()
    for bit_pauli.chunks[1] in 0:(4^n - 1)
        pauli = Pauli(convert(Vector{Bool}, bit_pauli), n)
        push!(pauli_list, pauli)
    end
    return pauli_list::Vector{Pauli}
end

"""
    get_support(p::Pauli)

Returns the support of the Pauli `p`.
"""
function get_support(p::Pauli)
    n = p.qubit_num
    pauli = p.pauli
    support = sort(findall(pauli[1:n] + pauli[(n + 1):(2n)] .> 0))
    return support::Vector{Int}
end

"""
    get_pauli_string(p::Pauli)

Returns the string representation of the Pauli `p`.
"""
function get_pauli_string(p::Pauli)
    pauli = p.pauli
    n = p.qubit_num
    if pauli[2n + 1] == 0
        pauli_string = "+"
    elseif pauli[2n + 1] == 1
        pauli_string = "-"
    else
        throw(error("The Pauli $(pauli) has an invalid sign."))
    end
    for j in 1:n
        signature = Bool[pauli[j]; pauli[n + j]]
        if signature == Bool[0; 0]
        elseif signature == Bool[0; 1]
            pauli_string *= " Z$(j)"
        elseif signature == Bool[1; 0]
            pauli_string *= " X$(j)"
        elseif signature == Bool[1; 1]
            pauli_string *= " Y$(j)"
        else
            throw(
                error(
                    "The Pauli $(pauli) has an invalid signature $(signature) on qubit $(j).",
                ),
            )
        end
    end
    return pauli_string::String
end

"""
    get_mapping_string(m::Mapping, c::AbstractCircuit; two_qubit_only::Bool = false)

Returns the string representation of the mapping `m` for the circuit `c`, including eigenvalues.
"""
function get_mapping_string(
    m::Mapping,
    c::T;
    two_qubit_only::Bool = false,
) where {T <: AbstractCircuit}
    # Construct the Pauli string for the mapping
    initial_string = get_pauli_string(m.initial)
    final_string = get_pauli_string(m.final)
    mapping_string = initial_string * " => " * final_string * " : "
    # We order our one-qubit Paulis as
    # X, Z, Y
    one_qubit_strings = ["X", "Z", "Y"]
    # We order our two-qubit Paulis as
    # XI, IX, XX, ZI, YI, ZX, YX, IZ, XZ, IY, XY, ZZ, YZ, ZY, YY.
    two_qubit_strings = [
        "XI",
        "IX",
        "XX",
        "ZI",
        "YI",
        "ZX",
        "YX",
        "IZ",
        "XZ",
        "IY",
        "XY",
        "ZZ",
        "YZ",
        "ZY",
        "YY",
    ]
    # Determine the Pauli string for each gate eigenvalue in the design row
    non_zero_indices = findall(!iszero, m.design_row)
    total_gates = c.total_gates
    index_values = [gate_idx.indices[1] - 1 for gate_idx in c.gate_data.gate_indices]
    gate_eigenvalue_strings = Vector{Tuple{String, String}}()
    for nz_idx in non_zero_indices
        g_idx = findmax(index_values[index_values .< nz_idx])[1]
        g = total_gates[findfirst(index_values .== g_idx)]
        pauli_idx = nz_idx - g_idx
        if length(g.targets) == 1
            pauli_string = one_qubit_strings[pauli_idx]
        elseif length(g.targets) == 2
            pauli_string = two_qubit_strings[pauli_idx]
        else
            throw(error("The gate $(g) does not operate on either 1 or 2 qubits."))
        end
        eigenvalue_power = m.design_row[nz_idx]
        gate_string = "($(g.type):$(Int(g.index)):$(Int.(g.targets)))"
        eigenvalue_string = pauli_string * " ($(Int(eigenvalue_power)))"
        if ~two_qubit_only || length(g.targets) == 2
            push!(gate_eigenvalue_strings, (gate_string, eigenvalue_string))
        end
    end
    # Consolidate the eigenvalues for each gate
    gate_strings = unique([gate_string[1] for gate_string in gate_eigenvalue_strings])
    for gate_string in gate_strings
        mapping_string *= "\n    " * gate_string * " : "
        gate_indices = [
            gate_eigenvalue_string[1] == gate_string for
            gate_eigenvalue_string in gate_eigenvalue_strings
        ]
        for gate_eigenvalue_string in gate_eigenvalue_strings[gate_indices]
            eigenvalue_string = gate_eigenvalue_string[2]
            mapping_string *= eigenvalue_string * ", "
        end
        mapping_string = mapping_string[1:(end - 2)]
    end
    return mapping_string::String
end

"""
    get_random_pauli(n::Integer)

Returns the string representation of a random n-qubit Pauli.
"""
function get_random_pauli(n::Integer)
    # Generate a random n-qubit Pauli
    random_pauli = join(rand(["I", "X", "Z", "Y"]) for i in 1:n)
    return random_pauli::String
end

"""
    pauli_to_string(p::Pauli)

Returns the string representation of the Pauli `p`.
"""
function pauli_to_string(p::Pauli)
    # Convert the Pauli into a string
    pauli = p.pauli
    n = p.qubit_num
    if ~pauli[2n + 1]
        pauli_string = "+"
    elseif pauli[2n + 1]
        pauli_string = "-"
    else
        throw(error("The Pauli $(pauli) has an invalid sign."))
    end
    for j in 1:n
        signature = [pauli[j]; pauli[n + j]]
        if signature == [false; false]
            pauli_string *= "I"
        elseif signature == [true; false]
            pauli_string *= "X"
        elseif signature == [false; true]
            pauli_string *= "Z"
        elseif signature == [true; true]
            pauli_string *= "Y"
        else
            throw(
                error(
                    "The Pauli $(pauli) has an invalid signature $(signature) on qubit $(j).",
                ),
            )
        end
    end
    return pauli_string::String
end

"""
    string_to_pauli(pauli_string::String)

Returns the Pauli corresponding to the string representation `pauli_string`.
"""
function string_to_pauli(pauli_string::String)
    # Convert the string into a Pauli
    n = length(pauli_string) - 1
    pauli = zeros(Bool, 2n + 1)
    offset = 0
    if pauli_string[1] == '-'
        pauli[2n + 1] = true
        offset = 1
    elseif pauli_string[1] == '+'
        pauli[2n + 1] = false
        offset = 1
    end
    for j in 1:n
        if pauli_string[j + offset] == 'I'
        elseif pauli_string[j + offset] == 'X'
            pauli[j] = true
        elseif pauli_string[j + offset] == 'Z'
            pauli[n + j] = true
        elseif pauli_string[j + offset] == 'Y'
            pauli[j] = true
            pauli[n + j] = true
        else
            throw(
                error(
                    "The Pauli string $(pauli_string) has an invalid Pauli on qubit $(j).",
                ),
            )
        end
    end
    p = Pauli(pauli, n)
    return p::Pauli
end

"""
    calc_pauli(initial::Pauli, circuit::Vector{Layer})

Returns the Pauli to which the initial Pauli `initial` is mapped after the circuit `circuit` is applied.
"""
function calc_pauli(initial::Pauli, circuit::Vector{Layer})
    # Set up the tableau
    n = initial.qubit_num
    t = Tableau(n)
    # Generate the layer which prepares the appropriate Pauli
    @assert initial.pauli[2n + 1] == 0 "This function requires an initial Pauli with positive sign."
    initial_support = get_support(initial)
    @assert length(initial_support) > 0 "The initial Pauli is not supported on any qubits."
    prep_layer = get_prep_layer(initial, initial_support)
    # Prepare the tableau with the initial Pauli
    apply!(t, prep_layer)
    # Use row_sum! to assemble the Pauli
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    @assert t.tableau[n + initial_support[1], :] == initial.pauli "The initial Pauli has not been appropriately initialised in the tableau."
    # Use row_sum! to disassemble the Pauli
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    # Apply the main circuit
    for l in circuit
        apply!(t, l)
    end
    # Use row_sum! to assemble the Pauli
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    # Determine the final Pauli to which the initial Pauli is mapped
    final = Pauli(t.tableau[n + initial_support[1], :], n)
    final_support = get_support(final)
    # Use row_sum! to disassemble the Pauli
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    # Generate the layer which measures the appropriate Pauli
    meas_layer = get_meas_layer(final, final_support)
    # Check that the measurement outcome is appropriate
    measurements = apply!(t, meas_layer; return_measurements = true)
    @assert prod(meas[1] for meas in measurements) == (-1)^final.pauli[2n + 1] "The measurement does not match the final Pauli."
    return final::Pauli
end

"""
    calc_orbit(initial::Pauli, gate::Gate; ignore_sign::Bool = true)

Returns the Pauli orbit generated by the gate `gate` acting on the Pauli `initial`, ignoring the sign of the Pauli operators if `ignore_sign` is `true`.
"""
function calc_orbit(initial::Pauli, gate::Gate; ignore_sign::Bool = true)
    # Create the tableau and gate layer
    n = length(gate.targets)
    @assert gate.targets == convert(Vector{Int16}, collect(1:n)) "This function only supports n-qubit gates that act on qubits 1 to n."
    t = Tableau(n)
    l = Layer([gate], n)
    # Generate the layer which prepares the appropriate Pauli
    @assert initial.pauli[2n + 1] == 0 "This function requires an initial Pauli with positive sign."
    initial_support = get_support(initial)
    @assert length(initial_support) > 0 "The initial Pauli is not supported on any qubits."
    prep_layer = get_prep_layer(initial, initial_support)
    # Prepare the tableau with the initial Pauli
    apply!(t, prep_layer)
    # Store the Pauli in the row corresponding to the first element of the initial support
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    @assert t.tableau[n + initial_support[1], :] == initial.pauli "The initial Pauli has not been appropriately initialised in the tableau."
    # Apply the gate and determine the Paulis in the orbit
    pauli_orbit = Pauli[initial]
    apply!(t, l)
    orbit = Pauli(t.tableau[n + initial_support[1], :], n)
    if ignore_sign
        orbit.pauli[2n + 1] = 0
    end
    while orbit != initial
        push!(pauli_orbit, orbit)
        apply!(t, l)
        orbit = Pauli(t.tableau[n + initial_support[1], :], n)
        if ignore_sign
            orbit.pauli[2n + 1] = 0
        end
    end
    return pauli_orbit::Vector{Pauli}
end

"""
    calc_gate_orbits(gate::Gate; ignore_sign::Bool = true)

Returns the Pauli orbits generated by the gate `gate`, for all Paulis supported on that gate, ignoring the sign of the Pauli operators if `ignore_sign` is `true`.
"""
function calc_gate_orbits(gate::Gate; ignore_sign::Bool = true)
    # Calculate the orbits for each Pauli acted upon by the gate
    n = length(gate.targets)
    @assert gate.targets == convert(Vector{Int16}, collect(1:n)) "This function only supports n-qubit gates that act on qubits 1 to n."
    gate_paulis = get_paulis(n)[2:end]
    gate_orbits = Vector{Vector{Pauli}}(undef, length(gate_paulis))
    for (idx, initial) in pairs(gate_paulis)
        gate_orbits[idx] = calc_orbit(initial, gate; ignore_sign = ignore_sign)
    end
    return gate_orbits::Vector{Vector{Pauli}}
end

"""
    get_orbit_indices(gate_orbits::Vector{Vector{Pauli}})

Returns the Pauli orbit indices for each of the Pauli orbits in `gate_orbits`.
"""
function get_orbit_indices(gate_orbits::Vector{Vector{Pauli}})
    # Index the Paulis occurring in each orbit
    orbit_num = length(gate_orbits)
    repeat_orbit_indices = Vector{Vector{Int}}(undef, orbit_num)
    for idx in 1:orbit_num
        repeat_orbit_indices[idx] = sort(
            findall([
                issetequal(gate_orbits[idx], pauli_orbit) for pauli_orbit in gate_orbits
            ]),
        )
    end
    orbit_indices = unique(repeat_orbit_indices)
    @assert sort(unique(vcat(orbit_indices...))) == collect(1:orbit_num) "The orbit indices are not unique."
    return orbit_indices::Vector{Vector{Int}}
end

"""
    display_circuit(stim_circuit_string::String; without_noise::Bool = true)
    display_circuit(qiskit_circuit::Py)

Displays the supplied Stim or Qiskit circuit `stim_circuit_string` or `qiskit_circuit`, and if `without_noise` is `true`, displays the Stim circuit without Pauli noise channels.
Be careful: Python arrays stored in PythonCall index from 0, so it is easy to index into the wrong Qiskit circuit.
"""
function display_circuit(stim_circuit_string::String; without_noise::Bool = true)
    if without_noise
        display(stim.Circuit(stim_circuit_string).without_noise().diagram("timeslice-svg"))
    else
        display(stim.Circuit(stim_circuit_string).diagram("timeslice-svg"))
    end
    return nothing
end
function display_circuit(qiskit_circuit::Py)
    display(qiskit_circuit.draw(; output = "mpl", fold = -1))
    return nothing
end

"""
    parse_uint8_vector(bitstring::String)

Returns the bitstring `bitstring` parsed into a `UInt8` vector.
"""
function parse_uint8_vector(bitstring::String)
    bit_num = length(bitstring)
    uint8_num = ceil(Int, bit_num / 8)
    uint8_vector = Vector{UInt8}(undef, uint8_num)
    for idx in 1:uint8_num
        start_idx = (idx - 1) * 8 + 1
        if idx == uint8_num
            end_idx = bit_num
        else
            end_idx = idx * 8
        end
        uint8_vector[idx] = parse(UInt8, bitstring[end_idx:-1:start_idx]; base = 2)
    end
    return uint8_vector::Vector{UInt8}
end

"""
    parse_bitstring(uint8_vector::Vector{UInt8}, bit_num::Integer)

Returns the `UInt8` vector `uint8_vector` parsed into a bitstring of length `bit_num`.
"""
function parse_bitstring(uint8_vector::Vector{UInt8}, bit_num::Integer)
    @assert length(uint8_vector) == ceil(Int, bit_num / 8)
    bitstring =
        join(((uint8_vector[1 + fld(m - 1, 8)] >> ((m - 1) % 8)) & 1) for m in 1:bit_num)
    return bitstring::String
end

"""
    parse_qiskit_dict(counts::Dict{String, Int}, shots::Integer, qiskit_qubit_num::Integer; qiskit_qubit_map::Union{Vector{Int}, Nothing} = nothing)

Returns a parsed Stim-compatible matrix of `UInt8` circuit outcomes determined from the dictionary `counts` output by Qiskit with `shots` shots, where the circuit acts on `qiskit_qubit_num` qubits, and if `qiskit_qubit_map` is not `nothing`, checks that the measurement outcomes for the qubits not included in the map are 0.
"""
function parse_qiskit_dict(
    counts::Dict{String, Int},
    shots::Integer,
    qiskit_qubit_num::Integer;
    qiskit_qubit_map::Union{Vector{Int}, Nothing} = nothing,
)
    # Check the input parameters.
    if qiskit_qubit_map !== nothing
        @assert qiskit_qubit_num >= maximum(qiskit_qubit_map) + 1 "The supplied number of qubits is inconsistent with the qubit map."
        qiskit_qubit_complement =
            setdiff(collect(0:(qiskit_qubit_num - 1)), qiskit_qubit_map)
    end
    # Parse the dictionary
    uint8_num = ceil(Int, qiskit_qubit_num / 8)
    parsed_counts = Matrix{UInt8}(undef, shots, uint8_num)
    idx_start = 0
    for (key, value) in pairs(counts)
        # Qiskit uses little-endian order
        reverse_key = reverse(key)
        if qiskit_qubit_map !== nothing
            @assert all(
                outcome == '0' for outcome in reverse_key[qiskit_qubit_complement .+ 1]
            )
        end
        parsed_key = parse_uint8_vector(reverse_key)
        for idx in (idx_start + 1):(idx_start + value)
            parsed_counts[idx, :] = parsed_key
        end
        idx_start += value
    end
    parsed_counts = parsed_counts[randperm(shots), :]
    @assert idx_start == shots "Expected $(shots) shots, but instead found $(idx_start)."
    return parsed_counts::Matrix{UInt8}
end

"""
    get_stim_circuit(c::AbstractCircuit)
    get_stim_circuit(c::AbstractCircuit, circuit_tuple::Vector{Int})
    get_stim_circuit(c::AbstractCircuit, gate_probabilities::Dict{Gate, Vector{Float64}})
    get_stim_circuit(c::AbstractCircuit, circuit_tuple::Vector{Int}, gate_probabilities::Dict{Gate, Vector{Float64}})
    get_stim_circuit(d::Design)
    get_stim_circuit(d::Design, gate_probabilities::Dict{Gate, Vector{Float64}})
    get_stim_circuit(d_rand::RandDesign)
    get_stim_circuit(d_rand::RandDesign, gate_probabilities::Dict{Gate, Vector{Float64}})

Returns the Stim circuit for the circuit `c`, stored in the design `d`, or randomised design `d_rand`, optionally applying the tuple `circuit_tuple` with the gate probabilities `gate_probabilities`.
"""
function get_stim_circuit(
    c::AbstractCircuit,
    circuit_tuple::Vector{Int},
    gate_probabilities::Dict{Gate, Vector{Float64}},
)
    stim_circuit_string = get_stim_circuit(
        c.circuit[circuit_tuple],
        gate_probabilities,
        c.noisy_prep,
        c.noisy_meas;
        extra_fields = c.extra_fields,
    )
    return stim_circuit_string::String
end
function get_stim_circuit(c::T, circuit_tuple::Vector{Int}) where {T <: AbstractCircuit}
    stim_circuit_string = get_stim_circuit(c, circuit_tuple, c.gate_probabilities)
    return stim_circuit_string::String
end
function get_stim_circuit(
    c::T,
    gate_probabilities::Dict{Gate, Vector{Float64}},
) where {T <: AbstractCircuit}
    stim_circuit_string = get_stim_circuit(c, c.circuit_tuple, gate_probabilities)
    return stim_circuit_string::String
end
function get_stim_circuit(c::T) where {T <: AbstractCircuit}
    stim_circuit_string = get_stim_circuit(c, c.circuit_tuple, c.gate_probabilities)
    return stim_circuit_string::String
end
function get_stim_circuit(d::Design)
    stim_circuit_string = get_stim_circuit(d.c)
    return stim_circuit_string::String
end
function get_stim_circuit(d::Design, gate_probabilities::Dict{Gate, Vector{Float64}})
    stim_circuit_string = get_stim_circuit(d.c, gate_probabilities)
    return stim_circuit_string::String
end
function get_stim_circuit(d_rand::RandDesign)
    d = get_design(d_rand)
    stim_circuit_string = get_stim_circuit(d)
    return stim_circuit_string::String
end
function get_stim_circuit(
    d_rand::RandDesign,
    gate_probabilities::Dict{Gate, Vector{Float64}},
)
    d = get_design(d_rand)
    stim_circuit_string = get_stim_circuit(d, gate_probabilities)
    return stim_circuit_string::String
end

"""
    pretty_print(c::AbstractCircuit)
    pretty_print(c::AbstractCircuit, circuit_tuple::Vector{Int})

Prints the circuit `c` using Stim in a readable format, optionally rearranged by the tuple `circuit_tuple`.
"""
function pretty_print(c::T, circuit_tuple::Vector{Int}) where {T <: AbstractCircuit}
    stim_circuit_string = get_stim_circuit(c, circuit_tuple)
    display_circuit(stim_circuit_string)
    return nothing
end
function pretty_print(c::T) where {T <: AbstractCircuit}
    pretty_print(c, c.circuit_tuple)
    return nothing
end

"""
    pretty_print(io::IO, d::Design)
    pretty_print(d::Design)

Prints the tuple set and shot weight data of the design `d` in a readable format.
"""
function pretty_print(io::IO, d::Design)
    # Initialise data
    tuple_set_data = d.tuple_set_data
    repeat_tuple_set = tuple_set_data.repeat_tuple_set
    repeat_numbers = tuple_set_data.repeat_numbers
    repeat_indices = tuple_set_data.repeat_indices
    nonrepeat_tuple_set = tuple_set_data.tuple_set
    repeat_tuple_num = length(repeat_indices)
    nonrepeat_tuple_num = length(nonrepeat_tuple_set)
    @assert length(d.tuple_set) == nonrepeat_tuple_num + repeat_tuple_num
    # Set up the tuple set data
    tuple_set_data = hcat(
        round.(d.shot_weights[1:nonrepeat_tuple_num], digits = 7),
        nonrepeat_tuple_set,
        ones(Int, nonrepeat_tuple_num),
    )
    formatted_repeat_tuple_set = [
        repeat_tuple_set[tuple_idx] for
        (tuple_idx, repeat_idx) in repeat_indices if repeat_numbers[repeat_idx] > 0
    ]
    formatted_repeat_numbers = [
        repeat_numbers[repeat_idx] for
        (tuple_idx, repeat_idx) in repeat_indices if repeat_numbers[repeat_idx] > 0
    ]
    repeat_tuple_set_data = hcat(
        round.(d.shot_weights[(nonrepeat_tuple_num + 1):end], digits = 7),
        formatted_repeat_tuple_set,
        formatted_repeat_numbers,
    )
    header = [
        "Shot weight"
        "Tuple"
        "Repeat number"
    ]
    pretty_table(
        io,
        vcat(tuple_set_data, repeat_tuple_set_data);
        header = header,
        alignment = :l,
        formatters = ft_printf("%.6f", 1),
    )
    return nothing
end
function pretty_print(d::Design)
    pretty_print(stdout, d)
    return nothing
end

"""
    pretty_print(io::IO, aces_data::ACESData, merit::Merit; projected::Bool = false)
    pretty_print(aces_data::ACESData, merit::Merit; projected::Bool = false)

Prints the z-scores of the normalised RMS errors of the gate eigenvalue estimator vector for the GLS, WLS, and OLS estimators in `aces_data` using the predictions in the merit data `merit`, or for the projected gate eigenvalue estimator vector if `projected` is `true`.
"""
function pretty_print(io::IO, aces_data::ACESData, merit::Merit; projected::Bool = false)
    # Set up parameters
    repetition_count = aces_data.repetitions
    budget_count = length(aces_data.budget_set)
    noise_score_coll = get_noise_score(aces_data, merit)
    if projected
        gls_z_scores = [noise_score.gls_proj_z_score for noise_score in noise_score_coll]
        wls_z_scores = [noise_score.wls_proj_z_score for noise_score in noise_score_coll]
        ols_z_scores = [noise_score.ols_proj_z_score for noise_score in noise_score_coll]
    else
        gls_z_scores = [noise_score.gls_z_score for noise_score in noise_score_coll]
        wls_z_scores = [noise_score.wls_z_score for noise_score in noise_score_coll]
        ols_z_scores = [noise_score.ols_z_score for noise_score in noise_score_coll]
    end
    # Print the data
    repetitions = convert(Vector{Int}, collect(1:repetition_count))
    shot_number =
        ["10^$(round(log10(aces_data.budget_set[i]), digits = 3))" for i in 1:budget_count]
    header = [
        "Repetition"
        "GLS S = " .* shot_number
        "WLS S = " .* shot_number
        "OLS S = " .* shot_number
    ]
    pretty_table(
        io,
        hcat(repetitions, gls_z_scores, wls_z_scores, ols_z_scores);
        header = header,
        alignment = :l,
        formatters = (ft_printf("%i", 1), ft_printf("%.4f")),
    )
    return nothing
end
function pretty_print(aces_data::ACESData, merit::Merit; projected::Bool = false)
    pretty_print(stdout, aces_data, merit; projected = projected)
    return nothing
end

"""
    eig_to_pair_idx(d::Design, eigenvalue_idx::Integer)
    eig_to_pair_idx(d::Design, eigenvalue_indices::Vector{Int})

Returns the `(tuple_idx, mapping_idx)` pair of indices for each eigenvalue index `eigenvalue_idx` corresponding to the design `d`.
"""
function eig_to_pair_idx(d::Design, eigenvalue_indices::Vector{Int})
    # Initialise variables
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([0; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    M = sum(mapping_lengths)
    # Determine the tuple and mapping indices
    pair_indices = Vector{Tuple{Int, Int}}(undef, length(eigenvalue_indices))
    for (idx, eigenvalue_idx) in pairs(eigenvalue_indices)
        @assert eigenvalue_idx >= 1 && eigenvalue_idx <= M "Invalid eigenvalue index $(eigenvalue_idx)."
        tuple_idx = findfirst(mapping_upper .>= eigenvalue_idx)
        mapping_idx = eigenvalue_idx - mapping_lower[tuple_idx]
        pair_indices[idx] = (tuple_idx, mapping_idx)
    end
    return pair_indices::Vector{Tuple{Int, Int}}
end
function eig_to_pair_idx(d::Design, eigenvalue_idx::Integer)
    pair_idx = eig_to_pair_idx(d, [eigenvalue_idx])[1]
    return pair_idx::Tuple{Int, Int}
end

"""
    pair_to_eig_idx(d::Design, pair_idx::Tuple{Int, Int})
    pair_to_eig_idx(d::Design, pair_indices::Vector{Tuple{Int, Int}})

Returns the eigenvalue index `eigenvalue_idx` corresponding to the pair index `(tuple_idx, mapping_idx)` corresponding to the design `d`.
"""
function pair_to_eig_idx(d::Design, pair_indices::Vector{Tuple{Int, Int}})
    # Initialise variables
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([0; mapping_lengths[1:(end - 1)]])
    # Determine the eigenvalue indices
    eigenvalue_indices = Vector{Int}(undef, length(pair_indices))
    for (idx, pair_idx) in pairs(pair_indices)
        eigenvalue_idx = mapping_lower[pair_idx[1]] + pair_idx[2]
        eigenvalue_indices[idx] = eigenvalue_idx
    end
    return eigenvalue_indices::Vector{Int}
end
function pair_to_eig_idx(d::Design, pair_idx::Tuple{Int, Int})
    eigenvalue_idx = pair_to_eig_idx(d, [pair_idx])[1]
    return eigenvalue_idx::Int
end

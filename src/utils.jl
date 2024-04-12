"""
    project_simplex(probabilities::Vector{Float64})

Projects a probability distribution onto the probability simplex in the Euclidean norm.
"""
function project_simplex(probabilities::Vector{Float64})
    @assert any(isnan.(probabilities)) == false "The supplied vector contains NaN."
    sorted_probabilities = sort(probabilities; rev = true)
    sum_probabilities =
        [1 / i for i in eachindex(probabilities)] .* (1 .- cumsum(sorted_probabilities))
    projected_probabilities =
        max.(
            probabilities .+
            sum_probabilities[findlast(sorted_probabilities + sum_probabilities .> 0)],
            0,
        )
    return projected_probabilities::Vector{Float64}
end

"""
    get_support(p::Pauli)

Returns the support of the Pauli.
"""
function get_support(p::Pauli)
    n = p.qubit_num
    pauli = p.pauli
    support = sort(findall(pauli[1:n] + pauli[(n + 1):(2n)] .> 0))
    return support::Vector{Int}
end

"""
    wht_matrix(n::Int)

Returns the symplectically ordered Walsh-Hadamard transform matrix of order n; to perform the inverse transform, divide the transform matrix by a scalar factor 4^n.
"""
function wht_matrix(n::Int)
    a = BitArray(undef, 2n)
    b = BitArray(undef, 2n)
    WHT_matrix = Matrix{Int}(undef, 4^n, 4^n)
    for a.chunks[1] in 0:(4^n - 1)
        for b.chunks[1] in 0:(4^n - 1)
            symplectic_form =
                convert(Bool, (a[1:n]' * b[(n + 1):(2n)] + a[(n + 1):(2n)]' * b[1:n]) % 2)
            WHT_matrix[a.chunks[1] + 1, b.chunks[1] + 1] = (-1)^symplectic_form
        end
    end
    return WHT_matrix::Matrix{Int}
end

#
function pretty_print(merit_array::Matrix{Float64})
    # Check the input parameters
    @assert size(merit_array) == (3, 6) "The input array should be of size (3, 6) and formatted as produced by the function compare_ls_optimise_weights."
    if merit_array[1, 1] > merit_array[2, 1]
        @warn "The GLS merit is better with shot weights optimised for WLS than for GLS."
    end
    if merit_array[1, 1] > merit_array[3, 1]
        @warn "The GLS merit is better with shot weights optimised for OLS than for GLS."
    end
    if merit_array[2, 2] > merit_array[1, 2]
        @warn "The WLS merit is better with shot weights optimised for GLS than for WLS."
    end
    if merit_array[2, 2] > merit_array[3, 2]
        @warn "The WLS merit is better with shot weights optimised for OLS than for WLS."
    end
    if merit_array[3, 3] > merit_array[1, 3]
        @warn "The OLS merit is better with shot weights optimised for GLS than for OLS."
    end
    if merit_array[3, 3] > merit_array[2, 3]
        @warn "The OLS merit is better with shot weights optimised for WLS than for OLS."
    end
    # Print the data
    shot_types = ["GLS"; "WLS"; "OLS"]
    header = [
        "Shot optim."
        "GLS expectation"
        "WLS expectation"
        "OLS expectation"
        "GLS std. dev."
        "WLS std. dev."
        "OLS std. dev."
    ]
    pretty_table(hcat(shot_types, merit_array); header = header, alignment = :l)
    return nothing
end

#
function pretty_print(d::Design)
    # Initialise data
    tuple_set_data = d.tuple_set_data
    repeat_tuple_set = tuple_set_data.repeat_tuple_set
    repeat_numbers = tuple_set_data.repeat_numbers
    repeat_indices = tuple_set_data.repeat_indices
    nonrepeat_tuple_set = tuple_set_data.tuple_set
    repeat_tuple_num = length(repeat_tuple_set)
    nonrepeat_tuple_num = length(nonrepeat_tuple_set)
    @assert length(d.tuple_set) == nonrepeat_tuple_num + repeat_tuple_num
    # Set up the repeat tuple set data
    formatted_repeat_numbers = [
        repeat_numbers[repeat_indices[idx]] for
        idx in eachindex(repeat_tuple_set) if repeat_numbers[repeat_indices[idx]] > 0
    ]
    repeat_tuple_set_data = hcat(
        round.(d.shot_weights[1:repeat_tuple_num], digits = 7),
        repeat_tuple_set,
        formatted_repeat_numbers,
    )
    tuple_set_data = hcat(
        round.(d.shot_weights[(repeat_tuple_num + 1):end], digits = 7),
        nonrepeat_tuple_set,
        ones(Int, nonrepeat_tuple_num),
    )
    header = [
        "Shot weight"
        "Tuple"
        "Repeat number"
    ]
    pretty_table(
        vcat(repeat_tuple_set_data, tuple_set_data);
        header = header,
        alignment = :l,
        formatters = ft_printf("%.6f", 1),
    )
    return nothing
end

#
function pretty_print(aces_data::ACESData, merit_set::Tuple{Merit, Merit, Merit})
    # Set up parameters
    repetition_count = aces_data.repetitions
    shot_count = length(aces_data.shots_set)
    gls_gate_norm_coll = aces_data.fgls_gate_norm_coll
    wls_gate_norm_coll = aces_data.wls_gate_norm_coll
    ols_gate_norm_coll = aces_data.ols_gate_norm_coll
    @assert size(gls_gate_norm_coll) == (repetition_count, shot_count) "The GLS gate norm collection has the wrong size $(size(gls_gate_norm_coll))."
    @assert size(wls_gate_norm_coll) == (repetition_count, shot_count) "The WLS gate norm collection has the wrong size $(size(wls_gate_norm_coll))."
    @assert size(ols_gate_norm_coll) == (repetition_count, shot_count) "The OLS gate norm collection has the wrong size $(size(ols_gate_norm_coll))."
    (gls_merit, wls_merit, ols_merit) = merit_set
    @assert gls_merit.ls_type == :gls "The GLS merit has the wrong type $(gls_merit.ls_type)."
    @assert wls_merit.ls_type == :wls "The WLS merit has the wrong type $(wls_merit.ls_type)."
    @assert ols_merit.ls_type == :ols "The OLS merit has the wrong type $(ols_merit.ls_type)."
    # Calculate the z-scores
    gls_z_scores = (gls_gate_norm_coll .- gls_merit.expectation) ./ sqrt(gls_merit.variance)
    wls_z_scores = (wls_gate_norm_coll .- wls_merit.expectation) ./ sqrt(wls_merit.variance)
    ols_z_scores = (ols_gate_norm_coll .- ols_merit.expectation) ./ sqrt(ols_merit.variance)
    # Print the data
    repetitions = convert(Vector{Int}, collect(1:repetition_count))
    shot_number =
        ["10^$(round(log10(aces_data.shots_set[i]), digits = 3))" for i in 1:shot_count]
    header = [
        "Repetition"
        "GLS S = " .* shot_number
        "WLS S = " .* shot_number
        "OLS S = " .* shot_number
    ]
    pretty_table(
        hcat(repetitions, gls_z_scores, wls_z_scores, ols_z_scores);
        header = header,
        alignment = :l,
        formatters = (ft_printf("%i", 1), ft_printf("%.4f")),
    )
    return nothing
end

"""
    get_pauli_string(p::Pauli)

Returns the string corresponding to the input Pauli.
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

#
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
    index_keys = collect(keys(c.gate_index))
    index_values = collect(values(c.gate_index))
    gate_eigenvalue_strings = Vector{Tuple{String, String}}(undef, 0)
    for nz_idx in non_zero_indices
        g_idx = findmax(index_values[index_values .< nz_idx])[1]
        g = index_keys[findfirst(index_values .== g_idx)]
        pauli_idx = nz_idx - g_idx
        if length(g.targets) == 1
            pauli_string = one_qubit_strings[pauli_idx]
        elseif length(g.targets) == 2
            pauli_string = two_qubit_strings[pauli_idx]
        else
            throw(error("The gate $(g) does not operate on either 1 or 2 qubits."))
        end
        eigenvalue_power = m.design_row[nz_idx]
        gate_string = "($(g.type)-$(Int(g.index)):$(Int.(g.targets)))"
        eigenvalue_string = pauli_string * " ^$(Int(eigenvalue_power))"
        if !two_qubit_only || length(g.targets) == 2
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
        mapping_string *= ";"
    end
    return mapping_string::String
end

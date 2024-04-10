"""
PauliPreparationSet(gates::Vector{Gate}, n::Int)

Returns all weight-1 and weight-2 Paulis whose support coincides with the support of a gate in the provided circuit.
"""
function PauliPreparationSet(gates::Vector{Gate}, n::Int)
    # Hard-coded weight-1 and weight-2 Pauli signatures
    # X, Z, Y
    weight_1_signature = [[1; 0], [0; 1], [1; 1]]
    # XX, ZX, YX, XZ, ZZ, YZ, XY, ZY, YY
    weight_2_signature = [
        [1; 1; 0; 0],
        [0; 1; 1; 0],
        [1; 1; 1; 0],
        [1; 0; 0; 1],
        [0; 0; 1; 1],
        [1; 0; 1; 1],
        [1; 1; 0; 1],
        [0; 1; 1; 1],
        [1; 1; 1; 1],
    ]
    # Find the supports of all the two-qubit gates
    weight_2_supports = Vector{Vector{Int}}(undef, 0)
    for gate in gates
        support = gate.targets
        if length(support) == 1
        elseif length(support) == 2
            push!(weight_2_supports, support)
        else
            throw(
                error(
                    "The gate $(gate) has a $(support) which is not size 1 or 2, which is not currently supported by this circuit.",
                ),
            )
        end
    end
    # Sort and remove duplicates
    sort!.(weight_2_supports)
    sort!(weight_2_supports)
    unique!(weight_2_supports)
    # Generate all combinations of Paulis
    l_2 = length(weight_2_supports)
    pauli_preparation_set = Vector{Pauli}(undef, 3n + 9l_2)
    for i in 1:n
        for a in 1:3
            pauli = zeros(Bool, 2n + 1)
            pauli[i] = weight_1_signature[a][1]
            pauli[n + i] = weight_1_signature[a][2]
            pauli_preparation_set[3 * (i - 1) + a] = Pauli(pauli, n)
        end
    end
    for i in 1:l_2
        for a in 1:9
            pauli = zeros(Bool, 2n + 1)
            pauli[weight_2_supports[i]] = weight_2_signature[a][1:2]
            pauli[n .+ weight_2_supports[i]] = weight_2_signature[a][3:4]
            pauli_preparation_set[3n + 9 * (i - 1) + a] = Pauli(pauli, n)
        end
    end
    return pauli_preparation_set::Vector{Pauli}
end

"""
    PrepLayer(initial::Pauli, initial_support::Vector{Int})

Returns a layer which prepares the provided Pauli on a tableau, flipping the sign if appropriate.
"""
function PrepLayer(initial::Pauli, initial_support::Vector{Int})
    n = initial.qubit_num
    prep = Vector{Gate}(undef, 0)
    for i in initial_support
        signature = Bool[initial.pauli[i]; initial.pauli[n + i]]
        if signature == Bool[0; 1]
            push!(prep, Gate("PZ+", 0, [i]))
        elseif signature == Bool[1; 0]
            push!(prep, Gate("PX+", 0, [i]))
        elseif signature == Bool[1; 1]
            push!(prep, Gate("PY+", 0, [i]))
        else
            throw(
                error(
                    "The Pauli $(initial) is meant to be supported on qubit $(i) and yet has signature $(signature) there.",
                ),
            )
        end
    end
    prep_layer = Layer(prep, n)
    return prep_layer::Layer
end

"""
    MeasLayer(final::Pauli, final_support::Vector{Int})

Returns a circuit which measures the provided Pauli on a tableau.
"""
function MeasLayer(final::Pauli, final_support::Vector{Int})
    n = final.qubit_num
    meas = Vector{Gate}(undef, 0)
    for i in final_support
        signature = Bool[final.pauli[i]; final.pauli[n + i]]
        if signature == Bool[0; 1]
            push!(meas, Gate("MZ", 0, [i]))
        elseif signature == Bool[1; 0]
            push!(meas, Gate("MX", 0, [i]))
        elseif signature == Bool[1; 1]
            push!(meas, Gate("MY", 0, [i]))
        else
            throw(
                error(
                    "The Pauli $(final) is meant to be supported on qubit $(i) and yet has signature $(signature) there.",
                ),
            )
        end
    end
    meas_layer = Layer(meas, n)
    return meas_layer::Layer
end

"""
    Update!(design_row::Vector{Int}, pauli::Pauli, l::Layer, gate_index::Dict{Gate, Int})

Updates the row of the design matrix according to the supplied Pauli for each of the gates in the layer.
"""
function Update!(
    design_row::Vector{Int},
    pauli::Pauli,
    l::Layer,
    gate_index::Dict{Gate, Int},
)
    n = pauli.qubit_num
    support = Support(pauli)
    # Calculate the relevant additions to the design matrix row for each gate in the layer
    for gate in l.layer
        if length(intersect(support, gate.targets)) > 0
            prep_gate = (gate.type ∈ ["PZ+", "PZ-", "PX+", "PX-", "PY+", "PY-"])
            meas_gate = (gate.type ∈ ["MZ", "MX", "MY"])
            if prep_gate || meas_gate
                # Check that the appropriate preparation or measurement is being performed
                target_1 = gate.targets[1]
                signature = Bool[pauli.pauli[target_1]; pauli.pauli[n + target_1]]
                wrong_Z = (signature == Bool[0; 1] && gate.type ∉ ["PZ+", "PZ-", "MZ"])
                wrong_X = (signature == Bool[1; 0] && gate.type ∉ ["PX+", "PX-", "MX"])
                wrong_Y = (signature == Bool[1; 1] && gate.type ∉ ["PY+", "PY-", "MY"])
                @assert !wrong_Z && !wrong_X && !wrong_Y "The Pauli $(pauli) is not prepared or measured by the gate $(gate)."
                @assert signature != Bool[0; 0] "The Pauli $(pauli) is meant to be supported on qubit $(target_1) and yet has signature $(signature) there."
                # Preparations and measurements get unique variables
                pauli_index = 1
            elseif length(gate.targets) == 2
                target_1 = gate.targets[1]
                target_2 = gate.targets[2]
                # The index orders 2-qubit Paulis as:
                # (II=0), XI=1,   IX=2,   XX=3
                # ZI=4,   YI=5,   ZX=6,   YX=7
                # IZ=8,   XZ=9,   IY=10,  XY=11
                # ZZ=12,  YZ=13,  ZY=14,  YY=15
                # This is the natural bit string ordering
                pauli_index =
                    pauli.pauli[target_1] +
                    2 * pauli.pauli[target_2] +
                    4 * pauli.pauli[n + target_1] +
                    8 * pauli.pauli[n + target_2]
            elseif length(gate.targets) == 1
                target_1 = gate.targets[1]
                # The index orders 1-qubit Paulis as:
                # (I=0),  X=1,    Z=2,    Y=3
                # This is the natural bit string ordering
                pauli_index = pauli.pauli[target_1] + 2 * pauli.pauli[n + target_1]
            else
                throw(error("The gate $(gate) does not operate on either 1 or 2 qubits."))
            end
            design_row[gate_index[gate] + pauli_index] += 1
        end
    end
    return nothing
end

"""
    Mapping(initial::Pauli, c::AbstractCircuit)

Returns a `Mapping` object describing how the initial Pauli is mapped by the circuit, and its corresponding row in the design matrix.
"""
function Mapping(initial::Pauli, c::AbstractCircuit)
    # Initialise variables
    n = c.qubit_num
    gate_index = c.gate_index
    N = c.N
    add_prep = c.add_prep
    add_meas = c.add_meas
    # Set up the tableau
    t = Tableau(n)
    # Set up the row of the design matrix
    design_row = zeros(Int, N)
    # Generate the layer which prepares the appropriate Pauli
    @assert initial.pauli[2n + 1] == 0 "This function requires an initial Pauli with positive sign."
    initial_support = Support(initial)
    prep_layer = PrepLayer(initial, initial_support)
    # Update the design matrix row if appropriate
    if add_prep
        Update!(design_row, initial, prep_layer, gate_index)
    end
    # Prepare the tableau with the initial Pauli
    Apply!(t, prep_layer)
    # Store the Pauli in the row corresponding to the first element of the initial support
    for i in eachindex(initial_support[2:end])
        RowSum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    @assert t.tableau[n + initial_support[1], :] == initial.pauli "The initial Pauli has not been appropriately initialised in the tableau."
    # Apply the main circuit
    spread_track = Vector{Vector{Int16}}(undef, 0)
    for layer in c.circuit
        # Track the spread of the Pauli
        layer_pauli = Pauli(t.tableau[n + initial_support[1], :], n)
        push!(spread_track, Support(layer_pauli))
        # Update the design matrix row
        Update!(design_row, layer_pauli, layer, gate_index)
        # Apply the circuit layer
        Apply!(t, layer)
    end
    # Determine the final Pauli to which the initial Pauli is mapped
    final = Pauli(t.tableau[n + initial_support[1], :], n)
    final_support = Support(final)
    push!(spread_track, final_support)
    # Undo the initial RowSums
    for i in eachindex(initial_support[2:end])
        RowSum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    # Generate the layer which measures the appropriate Pauli
    meas_layer = MeasLayer(final, final_support)
    if add_meas
        Update!(design_row, final, meas_layer, gate_index)
    end
    # Check that the measurement outcome is appropriate
    measurements = Apply!(t, meas_layer; return_measurements = true)
    @assert prod(meas[1] for meas in measurements) == (-1)^final.pauli[2n + 1] "The measurement does not match the final Pauli."
    # Generate a Mapping object
    design_row = convert(SparseVector{Int32, Int32}, dropzeros(sparse(design_row)))
    m = Mapping(initial, final, design_row, spread_track)
    return m::Mapping
end

"""
    MappingSet(c::AbstractCircuit)

Returns a vector of `Mapping` objects for each Pauli supported on some gate in the circuit which describes how that Pauli is mapped by the circuit, and its corresponding row in the design matrix.
"""
function MappingSet(c::AbstractCircuit)
    # Determine the Paulis whose mapping set we will calculate
    pauli_preparation_set = PauliPreparationSet(c.gates, c.qubit_num)
    pauli_num = length(pauli_preparation_set)
    # Calculate the mapping set
    mapping_set = Vector{Mapping}(undef, pauli_num)
    @threads :static for idx in 1:pauli_num
        mapping_set[idx] = Mapping(pauli_preparation_set[idx], c)
    end
    # Design matrix row corresponding to the mapping set
    mapping_matrix = convert(SparseMatrixCSC{Int32, Int32}, spzeros(Int32, pauli_num, c.N))
    for (idx, m) in enumerate(mapping_set)
        mapping_matrix[idx, :] = m.design_row
    end
    return (mapping_set::Vector{Mapping}, mapping_matrix::SparseMatrixCSC{Int32, Int32})
end

"""
    ConsistencySet(mapping_set::Vector{Mapping})

Returns a list, for each element of the mapping set, of the indices of all the other elements of the mapping set with which it is compatible.
"""
function ConsistencySet(mapping_set::Vector{Mapping})
    L = length(mapping_set)
    n = mapping_set[1].initial.qubit_num
    initial_paulis = [m.initial for m in mapping_set]
    initial_supports = [Support(initial) for initial in initial_paulis]
    final_paulis = [m.final for m in mapping_set]
    final_supports = [Support(final) for final in final_paulis]
    # Determine all other simultaneously preparable and measurement mappings for each mapping in the set
    consistency_set = Vector{Vector{Int}}(undef, L)
    for i in 1:L
        consistency_set[i] = [i]
    end
    reentrant_lock = ReentrantLock()
    @threads :static for i in 1:L
        for j in (i + 1):L
            # Initial Paulis are compatible if both Paulis are the same on each qubit in the support of both Paulis
            is_initial_compatible = true
            initial_int = intersect(initial_supports[i], initial_supports[j])
            for int in initial_int
                signature_1 = Bool[
                    initial_paulis[i].pauli[int]
                    initial_paulis[i].pauli[n + int]
                ]
                signature_2 = Bool[
                    initial_paulis[j].pauli[int]
                    initial_paulis[j].pauli[n + int]
                ]
                if signature_1 != signature_2
                    is_initial_compatible = false
                end
            end
            # Final Paulis are compatible if both Paulis are the same on each qubit in the support of both Paulis
            is_final_compatible = true
            final_int = intersect(final_supports[i], final_supports[j])
            for int in final_int
                signature_1 = Bool[
                    final_paulis[i].pauli[int]
                    final_paulis[i].pauli[n + int]
                ]
                signature_2 = Bool[
                    final_paulis[j].pauli[int]
                    final_paulis[j].pauli[n + int]
                ]
                if signature_1 != signature_2
                    is_initial_compatible = false
                end
            end
            # Both the initial and final Paulis must be compatible
            if is_initial_compatible && is_final_compatible
                lock(reentrant_lock) do
                    push!(consistency_set[i], j)
                    push!(consistency_set[j], i)
                end
            end
        end
    end
    return consistency_set::Vector{Vector{Int}}
end

"""
    PackPaulis(mapping_set::Vector{Mapping}, consistency_set::Vector{Vector{Int}})

Returns a set of sets of simultaneously preparable and measurable initial-final Pauli pairs in the mapping set.
"""
function PackPaulis(mapping_set::Vector{Mapping}, consistency_set::Vector{Vector{Int}})
    # Final Paulis
    final_paulis = [m.final for m in mapping_set]
    final_supports = [Support(final) for final in final_paulis]
    # Greedily add the Paulis into sets of compatible Paulis
    experiment_set = Vector{Vector{Int}}(undef, 0)
    # Set up the Paulis and sort them by size of the final Paulis
    total_paulis = sortperm(final_supports; by = x -> length(x), rev = true)
    unadded_paulis = deepcopy(total_paulis)
    while length(unadded_paulis) > 0
        # Add the first element of the unadded two-qubit Paulis to the set
        # The addable Paulis are the Paulis that have not yet been added which are compatible with all Paulis currently added to the set
        experiment = [unadded_paulis[1]]
        set_support = final_supports[experiment[1]]
        deleteat!(unadded_paulis, 1)
        addable_paulis = intersect(unadded_paulis, consistency_set[experiment[1]])
        # Add addable Paulis until you cannot
        while length(addable_paulis) > 0
            # Add the Pauli with the most overlap with the final Paulis of the Pauli set
            # Default to the first, which has the largest final Pauli support
            overlap_set = Vector{Int}(undef, length(addable_paulis))
            @threads :static for idx in eachindex(addable_paulis)
                overlap_set[idx] =
                    length(intersect(final_supports[addable_paulis[idx]], set_support))
            end
            max_overlap = maximum(overlap_set)
            best_add_idx = findfirst(x -> x == max_overlap, overlap_set)
            best_pauli = addable_paulis[best_add_idx]
            # Update all relevant sets
            best_unadd_idx = findfirst(x -> x == best_pauli, unadded_paulis)
            push!(experiment, best_pauli)
            union!(set_support, final_supports[best_pauli])
            deleteat!(addable_paulis, best_add_idx)
            intersect!(addable_paulis, consistency_set[best_pauli])
            deleteat!(unadded_paulis, best_unadd_idx)
        end
        # The compatible Paulis are the Paulis compatible with all the Paulis currently added to the set
        compatible_paulis = setdiff(total_paulis, experiment)
        for pauli in experiment
            intersect!(compatible_paulis, consistency_set[pauli])
        end
        # Add compatible Paulis until you cannot
        while length(compatible_paulis) > 0
            # Add the Pauli with the most overlap with the final Paulis of the Pauli set
            # Default to the first, which has the largest final Pauli support
            overlap_set = Vector{Int}(undef, length(compatible_paulis))
            @threads :static for idx in eachindex(compatible_paulis)
                overlap_set[idx] =
                    length(intersect(final_supports[compatible_paulis[idx]], set_support))
            end
            max_overlap = maximum(overlap_set)
            best_comp_idx = findfirst(x -> x == max_overlap, overlap_set)
            best_pauli = compatible_paulis[best_comp_idx]
            # Update all relevant sets
            push!(experiment, best_pauli)
            union!(set_support, final_supports[best_pauli])
            deleteat!(compatible_paulis, best_comp_idx)
            intersect!(compatible_paulis, consistency_set[best_pauli])
        end
        # Add the set of simultaneously measurable Paulis to the set of sets
        sort!(experiment)
        push!(experiment_set, experiment)
    end
    return experiment_set::Vector{Vector{Int}}
end

"""
    CovarianceDict(c::AbstractCircuit, mapping_set::Vector{Mapping}, experiment_set::Vector{Vector{Int}}, full_covariance::Bool)

Generates a dictionary with the requisite information to calculate the covariance matrix for the circuit eigenvalue estimators given the `experiment_set`, which describe the sets of simultaneously measured circuit eigenvalue estimators, corresponding to the `mapping_set`, under the specified circuit.
"""
function CovarianceDict(
    c::AbstractCircuit,
    mapping_set::Vector{Mapping},
    experiment_set::Vector{Vector{Int}},
    full_covariance::Bool,
)
    # Initialise variables
    tuple_circuit_number = length(experiment_set)
    trivial_pauli = Pauli(Bool[0], 0)
    zero_row = SparseVector{Int32, Int32}(0, Int32[], Int32[])
    trivial_row = SparseVector{Int32, Int32}(1, Int32[], Int32[])
    trivial_track = Vector{Int16}[]
    zero_mapping = Mapping(trivial_pauli, trivial_pauli, zero_row, trivial_track)
    trivial_mapping = Mapping(trivial_pauli, trivial_pauli, trivial_row, trivial_track)
    initial_paulis = [m.initial for m in mapping_set]
    # Determine diagonal terms in the covariance matrix
    covariance_dict = Dict{CartesianIndex{2}, Tuple{Mapping, Int}}()
    for j in 1:tuple_circuit_number
        experiment = experiment_set[j]
        # Determine the number of sign configurations as generated by PackCircuit
        experiment_prep_supports = [Support(initial_paulis[pauli]) for pauli in experiment]
        max_prep_support = maximum(length.(experiment_prep_supports))
        @assert max_prep_support ∈ [1; 2] "Currently only supports up to two-qubit Pauli preparations."
        experiment_sign_factor = 2^max_prep_support
        for index in experiment
            # Calculate the diagonal term
            diag_index = CartesianIndex(index, index)
            if haskey(covariance_dict, diag_index)
                covariance_dict[diag_index] = (
                    covariance_dict[diag_index][1],
                    covariance_dict[diag_index][2] + experiment_sign_factor,
                )
            else
                covariance_dict[diag_index] = (trivial_mapping, experiment_sign_factor)
            end
        end
    end
    # Determine off-diagonal terms in the covariance matrix if required
    # Only store for the upper diagonal as the covariance matrix is symmetric
    if full_covariance
        # Initialise variables
        L = length(c.circuit)
        reentrant_lock = ReentrantLock()
        design_rows = [m.design_row for m in mapping_set]
        spread_tracks = [m.spread_track for m in mapping_set]
        # Calculate the terms
        for j in 1:tuple_circuit_number
            experiment = experiment_set[j]
            S = length(experiment)
            # Determine the number of sign configurations as generated by PackCircuit
            experiment_prep_supports =
                [Support(initial_paulis[pauli]) for pauli in experiment]
            max_prep_support = maximum(length.(experiment_prep_supports))
            @assert max_prep_support ∈ [1; 2] "Currently only supports up to two-qubit Pauli preparations."
            experiment_sign_factor = 2^max_prep_support
            @threads :static for s_1 in 1:S
                # Get the first Pauli
                index_1 = experiment[s_1]
                pauli_1 = initial_paulis[index_1]
                design_row_1 = design_rows[index_1]
                spread_track_1 = spread_tracks[index_1]
                for s_2 in (s_1 + 1):S
                    # Get the second Pauli
                    index_2 = experiment[s_2]
                    pauli_2 = initial_paulis[index_2]
                    design_row_2 = design_rows[index_2]
                    spread_track_2 = spread_tracks[index_2]
                    # Check if the design matrix row has previously been computed
                    index = CartesianIndex(index_1, index_2)
                    if haskey(covariance_dict, index)
                        # If the design matrix row is non-trivial, increment the shots counter
                        lock(reentrant_lock) do
                            if covariance_dict[index][1] != zero_mapping
                                covariance_dict[index] = (
                                    covariance_dict[index][1],
                                    covariance_dict[index][2] + experiment_sign_factor,
                                )
                            end
                        end
                    else
                        # If the Paulis intersect, compute the mapping and add it to the dictionary if it's non-trivial
                        tracks_intersect = any(
                            (length(intersect(spread_track_1[i], spread_track_2[i])) > 0) for i in 1:(L + 1)
                        )
                        if tracks_intersect
                            sum_pauli = pauli_1 + pauli_2
                            # Compute the mapping if necessary
                            # We get thread safety issues when trying to check the dictionary
                            sum_index = findfirst(
                                sum_pauli == initial for initial in initial_paulis
                            )
                            if sum_index === nothing
                                sum_mapping = Mapping(sum_pauli, c)
                            else
                                sum_mapping = mapping_set[sum_index]
                            end
                            # Add the mapping if non-trivial
                            lock(reentrant_lock) do
                                if sum_mapping.design_row != design_row_1 + design_row_2
                                    covariance_dict[index] =
                                        (sum_mapping, experiment_sign_factor)
                                else
                                    covariance_dict[index] = (zero_mapping, 0)
                                end
                            end
                        else
                            lock(reentrant_lock) do
                                covariance_dict[index] = (zero_mapping, 0)
                            end
                        end
                    end
                end
            end
        end
    end
    # Remove the zero elements of the covariance matrix dictionary
    for key in keys(covariance_dict)
        if covariance_dict[key] == (zero_mapping, 0)
            delete!(covariance_dict, key)
        end
    end
    return covariance_dict::Dict{CartesianIndex{2}, Tuple{Mapping, Int}}
end

"""
    PackCircuits(code::Code, mapping_set::Vector{Mapping}, experiment_set::Vector{Vector{Int}})

Returns a set of circuits which prepares and measures the initial-final Pauli pairs, as dictated by the Pauli sets.
"""
function PackCircuits(
    code::Code,
    mapping_set::Vector{Mapping},
    experiment_set::Vector{Vector{Int}},
)
    # Generate an appropriate set of circuits based on the experiment set
    n = code.qubit_num
    tuple_circuit_number = length(experiment_set)
    prep_set = Vector{Vector{Layer}}(undef, tuple_circuit_number)
    meas_set = Vector{Layer}(undef, tuple_circuit_number)
    for j in 1:tuple_circuit_number
        # Set up some variables
        experiment = experiment_set[j]
        S = length(experiment)
        initial_set = [mapping_set[pauli].initial for pauli in experiment]
        initial_support_set = [Support(initial) for initial in initial_set]
        final_set = [mapping_set[pauli].final for pauli in experiment]
        final_support_set = [Support(final) for final in final_set]
        # Collate the appropriate preparation and measurement circuits
        prep = Vector{Gate}(undef, 0)
        meas = Vector{Gate}(undef, 0)
        for s in 1:S
            append!(prep, PrepLayer(initial_set[s], initial_support_set[s]).layer)
            append!(meas, MeasLayer(final_set[s], final_support_set[s]).layer)
        end
        # Remove repeated gates
        unique!(prep)
        unique!(meas)
        # Order the gates by the qubits on which they act
        sort!(prep)
        sort!(meas)
        # Store the measurement layer
        meas_set[j] = Layer(meas, n)
        # Set up the eigenstate sign combinations
        max_prep_support = maximum(length.(initial_support_set))
        if max_prep_support == 1
            prep_2 = Vector{Gate}(undef, 0)
            for gate in prep
                negative_type = replace(gate.type, "+" => "-")
                negative_gate = Gate(negative_type, gate.index, gate.targets)
                @assert gate.type != negative_type "The preparation gate $(gate) is unsigned."
                push!(prep_2, negative_gate)
            end
            @assert length(prep) == length(prep_2) "The constructed preparation layers do not have the same length as the original layer $(prep)."
            # Store the preparation layers
            prep_set[j] = [Layer(prep, n), Layer(prep_2, n)]
        elseif max_prep_support == 2
            prep_2 = Vector{Gate}(undef, 0)
            prep_3 = Vector{Gate}(undef, 0)
            prep_4 = Vector{Gate}(undef, 0)
            for gate in prep
                negative_type = replace(gate.type, "+" => "-")
                negative_gate = Gate(negative_type, gate.index, gate.targets)
                @assert gate.type != negative_type "The preparation gate $(gate) is unsigned."
                if gate.targets[1] ∈ code.ancilla_indices
                    push!(prep_2, negative_gate)
                    push!(prep_3, gate)
                    push!(prep_4, negative_gate)
                elseif gate.targets[1] ∈ code.data_indices
                    push!(prep_2, gate)
                    push!(prep_3, negative_gate)
                    push!(prep_4, negative_gate)
                else
                    throw(
                        error(
                            "The preparation gate $(gate) targets neither an ancilla qubit nor a data qubit.",
                        ),
                    )
                end
            end
            @assert (
                (length(prep) == length(prep_2)) &&
                (length(prep) == length(prep_3)) &&
                (length(prep) == length(prep_4))
            ) "The constructed preparation layers do not have the same length as the original layer $(prep)."
            # Store the preparation layers
            prep_set[j] =
                [Layer(prep, n), Layer(prep_2, n), Layer(prep_3, n), Layer(prep_4, n)]
        else
            throw(error("Currently only supports up to two-qubit Pauli preparations."))
        end
    end
    return (prep_set::Vector{Vector{Layer}}, meas_set::Vector{Layer})
end

"""
    GenerateDesign(code::Code, tuple_set::Vector{Vector{Int}}; shot_weights::Union{Vector{Float64}, Nothing} = nothing, full_covariance::Bool = true, save_data::Bool = false, diagnostics::Bool = false)

Generates a design matrix as well as the collection of circuits required to measure all the requisite initial-final Pauli pairs using the supplied tuple set.
"""
function GenerateDesign(
    code::Code,
    tuple_set::Vector{Vector{Int}};
    shot_weights::Union{Vector{Float64}, Nothing} = nothing,
    full_covariance::Bool = true,
    diagnostics::Bool = false,
    save_data::Bool = false,
    suppress_warnings::Bool = false,
)
    # Set some parameters
    start_time = time()
    T = length(tuple_set)
    N = code.N
    @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
    if shot_weights !== nothing
        @assert length(shot_weights) == T "The number of shot weights does not match the tuple set."
        @assert sum(shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
        @assert all(shot_weights .> 0.0) "The shot weights are not all positive."
    end
    # Warn the user if they have unadvisable settings for a large circuit
    if N >= 10^4 && !suppress_warnings
        if full_covariance
            @warn "This design is for a very large circuit: generating the full covariance matrix is unadvised."
        end
        if !diagnostics
            @warn "This design is for a very large circuit: turning on diagnostics is advised."
        end
        if !save_data
            @warn "This design is for a very large circuit: saving the data is advised."
        end
    end
    # Initialise the variables
    design_matrix = convert(SparseMatrixCSC{Int32, Int32}, spzeros(Int32, 0, N))
    mapping_ensemble = Vector{Vector{Mapping}}(undef, T)
    experiment_ensemble = Vector{Vector{Vector{Int}}}(undef, T)
    covariance_dict_ensemble =
        Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}(undef, T)
    prep_ensemble = Vector{Vector{Vector{Layer}}}(undef, T)
    meas_ensemble = Vector{Vector{Layer}}(undef, T)
    calculation_times = Matrix{Float64}(undef, T, 5)
    for t in 1:T
        # Determine the design data for the tuple
        circuit_tuple = tuple_set[t]
        code_tuple = ApplyTuple(code, circuit_tuple)
        time_1 = time()
        (mapping_set, mapping_matrix) = MappingSet(code_tuple)
        time_2 = time()
        consistency_set = ConsistencySet(mapping_set)
        time_3 = time()
        experiment_set = PackPaulis(mapping_set, consistency_set)
        time_4 = time()
        covariance_dict =
            CovarianceDict(code_tuple, mapping_set, experiment_set, full_covariance)
        time_5 = time()
        (prep_set, meas_set) = PackCircuits(code_tuple, mapping_set, experiment_set)
        time_6 = time()
        # Track the times taken
        mapping_time = time_2 - time_1
        consistency_time = time_3 - time_2
        pauli_time = time_4 - time_3
        covariance_time = time_5 - time_4
        circuit_time = time_6 - time_5
        if diagnostics
            println(
                "For tuple $(t), the mappings took $(round(mapping_time, digits = 3)) s, the consistency sets took $(round(consistency_time, digits = 3)) s, packing the Paulis took $(round(pauli_time, digits = 3)) s, the covariance dictionary took $(round(covariance_time, digits = 3)) s, and the circuits took $(round(circuit_time, digits = 3)) s. Since starting, $(round(time_5 - start_time, digits = 3)) s have elapsed.",
            )
        end
        # Store the data
        tuple_set[t] = circuit_tuple
        design_matrix = vcat(design_matrix, mapping_matrix)
        mapping_ensemble[t] = mapping_set
        experiment_ensemble[t] = experiment_set
        covariance_dict_ensemble[t] = covariance_dict
        prep_ensemble[t] = prep_set
        meas_ensemble[t] = meas_set
        calculation_times[t, :] =
            [mapping_time; consistency_time; pauli_time; covariance_time; circuit_time]
    end
    # Calculate the tuple times and shot weights
    experiment_numbers = length.([vcat(prep_set...) for prep_set in prep_ensemble])
    (tuple_times, tuple_shot_weights) =
        TupleSetParameters(code, tuple_set, experiment_numbers)
    if shot_weights === nothing
        shot_weights = tuple_shot_weights
    end
    # Save and return the results
    overall_time = time() - start_time
    d = Design(
        code,
        full_covariance,
        design_matrix,
        tuple_set,
        mapping_ensemble,
        experiment_ensemble,
        covariance_dict_ensemble,
        prep_ensemble,
        meas_ensemble,
        tuple_times,
        shot_weights,
        calculation_times,
        overall_time,
    )
    if save_data
        save_design(d)
    end
    if diagnostics
        println(
            "Generated the design for all $(T) tuples in $(round(overall_time, digits = 3)) s.",
        )
    end
    return d::Design
end

#
function GenerateDesign(
    code::Code,
    tuple_set_data::TupleSetData;
    shot_weights::Union{Vector{Float64}, Nothing} = nothing,
    full_covariance::Bool = true,
    diagnostics::Bool = false,
    save_data::Bool = false,
)
    # Save the results
    suppress_warnings = false
    if code.N >= 10^4
        suppress_warnings = true
        if full_covariance
            @warn "This design is for a very large circuit: generating the full covariance matrix is unadvised."
        end
        if !diagnostics
            @warn "This design is for a very large circuit: turning on diagnostics is advised."
        end
        if !save_data
            @warn "This design is for a very large circuit: saving the data is advised."
        end
    end
    # Generate the design
    tuple_set = TupleSet(tuple_set_data)
    d = GenerateDesign(
        code,
        tuple_set;
        shot_weights = shot_weights,
        full_covariance = full_covariance,
        diagnostics = diagnostics,
        suppress_warnings = suppress_warnings,
    )
    @reset d.tuple_set_data = tuple_set_data
    # Save and return the results
    if save_data
        save_design(d)
    end
    return d::Design
end

"""
    CompleteDesign(d::Design; diagnostics::Bool = false)

Completes the covariance dictionary for the design if it does not already have a full covariance matrix.
"""
function CompleteDesign(d::Design; diagnostics::Bool = false)
    if d.full_covariance
        if diagnostics
            println("The supplied design matrix already has a full covariance matrix.")
        end
        return d::Design
    else
        # Set some parameters
        start_time = time()
        full_covariance = true
        T = length(d.tuple_set)
        # Regenerate the covariance dictionary set with a full covariance matrix
        covariance_dict_ensemble =
            Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}(undef, T)
        covariance_times = Vector{Float64}(undef, T)
        for t in 1:T
            circuit_tuple = d.tuple_set[t]
            code_tuple = ApplyTuple(d.code, circuit_tuple)
            mapping_set_tuple = d.mapping_ensemble[t]
            experiment_set_tuple = d.experiment_ensemble[t]
            # Calculate the covariance dictionary
            time_4 = time()
            covariance_dict_tuple = CovarianceDict(
                code_tuple,
                mapping_set_tuple,
                experiment_set_tuple,
                full_covariance,
            )
            time_5 = time()
            # Track the times taken
            covariance_time = time_5 - time_4
            if diagnostics
                println(
                    "For tuple $(t), generating the covariance dictionary took $(round(covariance_time, digits = 3)) s. Since starting, $(round(time_5 - start_time, digits = 3)) s have elapsed.",
                )
            end
            # Store the data
            covariance_dict_ensemble[t] = covariance_dict_tuple
            covariance_times[t] = covariance_time
        end
        overall_time = time() - start_time
        # Set the new variables 
        d_complete = deepcopy(d)
        @reset d_complete.full_covariance = full_covariance
        @reset d_complete.covariance_dict_ensemble = covariance_dict_ensemble
        @reset d_complete.overall_time =
            d_complete.overall_time - sum(d_complete.calculation_times[:, 4]) + overall_time
        @reset d_complete.calculation_times[:, 4] = covariance_times
        if diagnostics
            println(
                "Generated the full covariance dictionaries for all $(T) tuples in $(round(overall_time, digits = 3)) s.",
            )
        end
        return d_complete::Design
    end
end

"""
    Pauli

Boolean representation of a Pauli operator.

# Fields

  - `pauli::Vector{Bool}`: The Pauli operator stored as a Boolean vector. The first `qubit_num` elements represent Pauli X on each qubit, the next `qubit_num` elements represent Pauli Z on each qubit, and the final element represents the sign.
  - `qubit_num::Int16`: The number of qubits on which the Pauli operator acts; the length of the vector is `2 * qubit_num + 1`.
"""
struct Pauli
    pauli::Vector{Bool}
    qubit_num::Int16
    # Constructor
    function Pauli(pauli::Vector{Bool}, n::Integer)
        @assert length(pauli) == 2n + 1 "The Pauli $(pauli) does not have a length commensurate with the qubit number $(n)."
        return new(pauli, n)::Pauli
    end
end

Base.show(io::IO, p::Pauli) = print(io, get_pauli_string(p))

function Base.:(+)(p₁::Pauli, p₂::Pauli)
    @assert p₁.qubit_num == p₂.qubit_num "The Paulis $(p₁) and $(p₂) do not have the same number of qubits."
    return Pauli(convert(Vector{Bool}, (p₁.pauli .+ p₂.pauli) .% 2), p₁.qubit_num)
end

@struct_hash_equal_isequal Pauli

"""
    Mapping

Mapping of a Pauli operator by some circuit.

# Fields

  - `initial::Pauli`: Initial Pauli operator before the action of the circuit.
  - `final::Pauli`: Final Pauli operator after the action of the circuit.
  - `design_row::SparseVector{Int32, Int32}`: Design matrix row for the circuit eigenvalue corresponding to the initial Pauli and the circuit used for the mapping.
  - `spread_track::Vector{Vector{Int16}}`: Track of the support of the Pauli as it is acted upon by the layers of the circuit.
"""
struct Mapping
    initial::Pauli
    final::Pauli
    design_row::SparseVector{Int32, Int32}
    spread_track::Vector{Vector{Int16}}
    # Constructor
    function Mapping(
        initial::Pauli,
        final::Pauli,
        design_row::SparseVector{Int32, Int32},
        spread_track::Vector{Vector{Int16}},
    )
        @assert initial.qubit_num == final.qubit_num "The initial Pauli $(initial) and final Pauli $(final) do not have the same number of qubits."
        return new(initial, final, design_row, spread_track)::Mapping
    end
end

function Base.show(io::IO, m::Mapping)
    return print(io, get_pauli_string(m.initial) * " => " * get_pauli_string(m.final))
end

@struct_hash_equal_isequal Mapping

"""
    Design

Experimental design for a noise characterisation experiment for a circuit.

# Fields

  - `c::AbstractCircuit`: Circuit characterised by the design.
  - `full_covariance::Bool`: If `true`, generates parameters to construct the full covariance matrix in `covariance_dict_ensemble`, else if `false`, only generates parameters to construct the terms on the diagonal.
  - `matrix::SparseMatrixCSC{Int32, Int32}`: Sparse M x N design matrix, corresponding to M circuit eigenvalues and N gate eigenvalues.
  - `tuple_set::Vector{Vector{Int}}`: Set of tuples which arrange the circuit layers.
  - `tuple_set_data::TupleSetData`: [`TupleSetData`](@ref) object that generates the tuple set.
  - `mapping_ensemble::Vector{Vector{Mapping}}`: Vector of the [`Mapping`](@ref) objects for each of the circuit eigenvalue for the Paulis corresponding to that tuple, for each tuple in the set.
  - `experiment_ensemble::Vector{Vector{Vector{Int}}}`: Vector of the experiments that index [`Mapping`](@ref) objects, which correspond to simultaneously preparable and measurable circuit eigenvalues, for each tuple in the set.
  - `covariance_dict_ensemble::Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}`: Dictionary of [`Mapping`](@ref) objects describing the non-zero entries of the sparse circuit eigenvalue estimator covariance matrix, alongside the number of times the entry is estimated by the experiment set, for each tuple in the set.
  - `prep_ensemble::Vector{Vector{Vector{Layer}}}`: Vector of [`Layer`](@ref) objects that prepare qubits in Pauli eigenstates for each experiment in the set, indeed a vector preparing the necessary sign configurations, for each tuple in the set.
  - `meas_ensemble::Vector{Vector{Layer}}`: Vector of [`Layer`](@ref) objects that measure qubits in Pauli bases for each experiment in the set, for each tuple in the set.
  - `tuple_times::Vector{Float64}`: Time taken to implement the circuit arranged by each tuple in the set, normalised according to the time factor for the basic tuple set.
  - `shot_weights::Vector{Float64}`: Shot weights for each tuple in the set, which add to 1.
  - `experiment_numbers::Vector{Int}`: Number of experiments for each tuple in the set.
  - `experiment_number::Int`: Total number of experiments.
  - `calculation_times::Matrix{Float64}`: Time taken to generate components of the design for each tuple, which correspond to generating: the mappings, the sets of tuple-consistent Pauli preparations, the experiment sets, the covariance matrix dictionary, and the circuits.
  - `overall_time::Float64`: Overall time taken to generate the design.
  - `optimisation_time::Float64`: Time taken to optimise the design.
  - `ls_type::Symbol`: Type of least squares for which the shot weights were optimised.
"""
struct Design
    c::AbstractCircuit
    full_covariance::Bool
    matrix::SparseMatrixCSC{Int32, Int32}
    tuple_set::Vector{Vector{Int}}
    tuple_set_data::TupleSetData
    mapping_ensemble::Vector{Vector{Mapping}}
    experiment_ensemble::Vector{Vector{Vector{Int}}}
    covariance_dict_ensemble::Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}
    prep_ensemble::Vector{Vector{Vector{Layer}}}
    meas_ensemble::Vector{Vector{Layer}}
    tuple_times::Vector{Float64}
    shot_weights::Vector{Float64}
    experiment_numbers::Vector{Int}
    experiment_number::Int
    calculation_times::Matrix{Float64}
    overall_time::Float64
    optimisation_time::Float64
    ls_type::Symbol
    # Default constructor
    function Design(
        c::T,
        full_covariance::Bool,
        matrix::SparseMatrixCSC{Int32, Int32},
        tuple_set::Vector{Vector{Int}},
        tuple_set_data::TupleSetData,
        mapping_ensemble::Vector{Vector{Mapping}},
        experiment_ensemble::Vector{Vector{Vector{Int}}},
        covariance_dict_ensemble::Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}},
        prep_ensemble::Vector{Vector{Vector{Layer}}},
        meas_ensemble::Vector{Vector{Layer}},
        tuple_times::Vector{Float64},
        shot_weights::Vector{Float64},
        experiment_numbers::Vector{Int},
        experiment_number::Int,
        calculation_times::Matrix{Float64},
        overall_time::Float64,
        optimisation_time::Float64,
        ls_type::Symbol,
    ) where {T <: AbstractCircuit}
        # Check parameters
        tuple_number = length(tuple_set)
        @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
        @assert tuple_set == get_tuple_set(tuple_set_data) "The tuple set doesn't align with the tuple set data."
        @assert length(mapping_ensemble) == tuple_number "The size of the mapping ensemble does not match the tuple set."
        @assert length(experiment_ensemble) == tuple_number "The size of the experiment ensemble does not match the tuple set."
        @assert length(covariance_dict_ensemble) == tuple_number "The size of the covariance dictionary ensemble does not match the tuple set."
        @assert length(prep_ensemble) == tuple_number "The size of the preparation ensemble does not match the tuple set."
        @assert length(meas_ensemble) == tuple_number "The size of the measurement ensemble does not match the tuple set."
        @assert length(tuple_times) == tuple_number "The number of tuple times does not match the tuple set."
        @assert length(shot_weights) == tuple_number "The number of shot weights does not match the tuple set."
        @assert sum(shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
        @assert all(shot_weights .> 0.0) "The shot weights are not all positive."
        @assert length(experiment_numbers) == tuple_number "The number of experiment numbers does not match the tuple set."
        @assert experiment_number == sum(experiment_numbers) "The experiment number $(experiment_number) does not match the sum of the experiment numbers $(sum(experiment_numbers))."
        @assert size(calculation_times) == (tuple_number, 5) "The calculation times do not match the tuple set."
        @assert ls_type ∈ [:none, :gls, :wls, :ols] "The least squares type $(ls_type) is not supported."
        # Return the design
        return new(
            c,
            full_covariance,
            matrix,
            tuple_set,
            tuple_set_data,
            mapping_ensemble,
            experiment_ensemble,
            covariance_dict_ensemble,
            prep_ensemble,
            meas_ensemble,
            tuple_times,
            shot_weights,
            experiment_numbers,
            experiment_number,
            calculation_times,
            overall_time,
            optimisation_time,
            ls_type,
        )::Design
    end
    # Constructor
    function Design(
        c::T,
        full_covariance::Bool,
        matrix::SparseMatrixCSC{Int32, Int32},
        tuple_set::Vector{Vector{Int}},
        mapping_ensemble::Vector{Vector{Mapping}},
        experiment_ensemble::Vector{Vector{Vector{Int}}},
        covariance_dict_ensemble::Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}},
        prep_ensemble::Vector{Vector{Vector{Layer}}},
        meas_ensemble::Vector{Vector{Layer}},
        tuple_times::Vector{Float64},
        shot_weights::Vector{Float64},
        calculation_times::Matrix{Float64},
        overall_time::Float64;
        ls_type::Symbol = :none,
    ) where {T <: AbstractCircuit}
        # Initialise parameters
        tuple_number = length(tuple_set)
        tuple_set_data = TupleSetData(tuple_set, Vector{Int}[], Int[], Int[])
        experiment_numbers =
            length.([vcat(prep_layer_set...) for prep_layer_set in prep_ensemble])
        experiment_number = length(vcat(vcat(prep_ensemble...)...))
        optimisation_time = 0.0
        # Check parameters
        @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
        @assert length(mapping_ensemble) == tuple_number "The size of the mapping ensemble does not match the tuple set."
        @assert length(experiment_ensemble) == tuple_number "The size of the experiment ensemble does not match the tuple set."
        @assert length(covariance_dict_ensemble) == tuple_number "The size of the covariance dictionary ensemble does not match the tuple set."
        @assert length(prep_ensemble) == tuple_number "The size of the preparation ensemble does not match the tuple set."
        @assert length(meas_ensemble) == tuple_number "The size of the measurement ensemble does not match the tuple set."
        @assert length(tuple_times) == tuple_number "The number of tuple times does not match the tuple set."
        @assert length(shot_weights) == tuple_number "The number of shot weights does not match the tuple set."
        @assert sum(shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
        @assert all(shot_weights .> 0.0) "The shot weights are not all positive."
        @assert length(experiment_numbers) == tuple_number "The number of experiment numbers does not match the tuple set."
        @assert experiment_number == sum(experiment_numbers) "The experiment number $(experiment_number) does not match the sum of the experiment numbers $(sum(experiment_numbers))."
        @assert size(calculation_times) == (tuple_number, 5) "The calculation times do not match the tuple set."
        @assert ls_type ∈ [:none, :gls, :wls, :ols] "The least squares type $(ls_type) is not supported."
        # Return the design
        return new(
            c,
            full_covariance,
            matrix,
            tuple_set,
            tuple_set_data,
            mapping_ensemble,
            experiment_ensemble,
            covariance_dict_ensemble,
            prep_ensemble,
            meas_ensemble,
            tuple_times,
            shot_weights,
            experiment_numbers,
            experiment_number,
            calculation_times,
            overall_time,
            optimisation_time,
            ls_type,
        )::Design
    end
end

function Base.show(io::IO, d::Design)
    return print(
        io,
        "Design for a $(d.c.circuit_param.circuit_name) circuit with $(length(d.tuple_set)) tuples and $(d.experiment_number) experiments.",
    )
end

@struct_hash_equal_isequal Design

"""
    get_pauli_prep_set(gates::Vector{Gate}, n::Int)

Returns all weight-1 Paulis on all `n` qubits, as well as the weight-2 Paulis supported on some gate in `gates`.
"""
function get_pauli_prep_set(gates::Vector{Gate}, n::Int)
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
    pauli_prep_set = Vector{Pauli}(undef, 3n + 9l_2)
    for i in 1:n
        for a in 1:3
            pauli = zeros(Bool, 2n + 1)
            pauli[i] = weight_1_signature[a][1]
            pauli[n + i] = weight_1_signature[a][2]
            pauli_prep_set[3 * (i - 1) + a] = Pauli(pauli, n)
        end
    end
    for i in 1:l_2
        for a in 1:9
            pauli = zeros(Bool, 2n + 1)
            pauli[weight_2_supports[i]] = weight_2_signature[a][1:2]
            pauli[n .+ weight_2_supports[i]] = weight_2_signature[a][3:4]
            pauli_prep_set[3n + 9 * (i - 1) + a] = Pauli(pauli, n)
        end
    end
    return pauli_prep_set::Vector{Pauli}
end

"""
    get_prep_layer(initial::Pauli, initial_support::Vector{Int})

Returns a layer which prepares the Pauli `initial`, supported on the qubits in `initial_support`.
"""
function get_prep_layer(initial::Pauli, initial_support::Vector{Int})
    n = initial.qubit_num
    prep = Vector{Gate}(undef, length(initial_support))
    for (idx, i) in enumerate(initial_support)
        signature = Bool[initial.pauli[i]; initial.pauli[n + i]]
        if signature == Bool[0; 1]
            prep[idx] = Gate("PZ+", 0, [i])
        elseif signature == Bool[1; 0]
            prep[idx] = Gate("PX+", 0, [i])
        elseif signature == Bool[1; 1]
            prep[idx] = Gate("PY+", 0, [i])
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
    get_meas_layer(final::Pauli, final_support::Vector{Int})

Returns a layer which measures the Pauli `final`, supported on the qubits in `final_support`.
"""
function get_meas_layer(final::Pauli, final_support::Vector{Int})
    n = final.qubit_num
    meas = Vector{Gate}(undef, length(final_support))
    for (idx, i) in enumerate(final_support)
        signature = Bool[final.pauli[i]; final.pauli[n + i]]
        if signature == Bool[0; 1]
            meas[idx] = Gate("MZ", 0, [i])
        elseif signature == Bool[1; 0]
            meas[idx] = Gate("MX", 0, [i])
        elseif signature == Bool[1; 1]
            meas[idx] = Gate("MY", 0, [i])
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
    update_design_row!(design_row::Vector{Int}, pauli::Pauli, l::Layer, gate_index::Dict{Gate, Int})

Updates the design matrix row `design_row` according to the Pauli `pauli` for each gate in the layer `l`, using the gate index `gate_index`.
"""
function update_design_row!(
    design_row::Vector{Int},
    pauli::Pauli,
    l::Layer,
    gate_index::Dict{Gate, Int},
)
    n = pauli.qubit_num
    support = get_support(pauli)
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
                @assert ~wrong_Z && ~wrong_X && ~wrong_Y "The Pauli $(pauli) is not prepared or measured by the gate $(gate)."
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
    calc_mapping(initial::Pauli, c::AbstractCircuit)

Returns a [`Mapping`](@ref) object for the Pauli `initial` when mapped by the circuit `c`.
"""
function calc_mapping(initial::Pauli, c::T) where {T <: AbstractCircuit}
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
    initial_support = get_support(initial)
    prep_layer = get_prep_layer(initial, initial_support)
    # Update the design matrix row if appropriate
    if add_prep
        update_design_row!(design_row, initial, prep_layer, gate_index)
    end
    # Prepare the tableau with the initial Pauli
    apply!(t, prep_layer)
    # Store the Pauli in the row corresponding to the first element of the initial support
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    @assert t.tableau[n + initial_support[1], :] == initial.pauli "The initial Pauli has not been appropriately initialised in the tableau."
    # Apply the main circuit
    spread_track = Vector{Vector{Int16}}(undef, 0)
    for l in c.circuit
        # Track the spread of the Pauli
        layer_pauli = Pauli(t.tableau[n + initial_support[1], :], n)
        push!(spread_track, get_support(layer_pauli))
        # Update the design matrix row
        update_design_row!(design_row, layer_pauli, l, gate_index)
        # Apply the circuit layer
        apply!(t, l)
    end
    # Determine the final Pauli to which the initial Pauli is mapped
    final = Pauli(t.tableau[n + initial_support[1], :], n)
    final_support = get_support(final)
    push!(spread_track, final_support)
    # Undo the initial row_sum! operations
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    # Generate the layer which measures the appropriate Pauli
    meas_layer = get_meas_layer(final, final_support)
    if add_meas
        update_design_row!(design_row, final, meas_layer, gate_index)
    end
    # Check that the measurement outcome is appropriate
    measurements = apply!(t, meas_layer; return_measurements = true)
    @assert prod(meas[1] for meas in measurements) == (-1)^final.pauli[2n + 1] "The measurement does not match the final Pauli."
    # Generate a Mapping object
    design_row = convert(SparseVector{Int32, Int32}, dropzeros(sparse(design_row)))
    m = Mapping(initial, final, design_row, spread_track)
    return m::Mapping
end

"""
    calc_mapping_set(c::AbstractCircuit)

Returns a vector of [`Mapping`](@ref) objects for each single-qubit Pauli on the qubits in the circuit `c`, and each two-qubit Pauli supported on some gate in the circuit `c`.
"""
function calc_mapping_set(c::T) where {T <: AbstractCircuit}
    # Determine the Paulis whose mapping set we will calculate
    pauli_prep_set = get_pauli_prep_set(c.gates, c.qubit_num)
    pauli_num = length(pauli_prep_set)
    # Calculate the mapping set
    mapping_set = Vector{Mapping}(undef, pauli_num)
    @threads :static for idx in 1:pauli_num
        mapping_set[idx] = calc_mapping(pauli_prep_set[idx], c)
    end
    # Design matrix row corresponding to the mapping set
    mapping_matrix = convert(SparseMatrixCSC{Int32, Int32}, spzeros(Int32, pauli_num, c.N))
    for (idx, m) in enumerate(mapping_set)
        mapping_matrix[idx, :] = m.design_row
    end
    return (mapping_set::Vector{Mapping}, mapping_matrix::SparseMatrixCSC{Int32, Int32})
end

"""
    calc_consistency_set(mapping_set::Vector{Mapping})

Returns a list for each mapping in `mapping_set` of all other mappings with which it is simultaneously preparable and measurable, and hence tuple-consistent.
"""
function calc_consistency_set(mapping_set::Vector{Mapping})
    L = length(mapping_set)
    n = mapping_set[1].initial.qubit_num
    initial_paulis = [m.initial for m in mapping_set]
    initial_supports = [get_support(initial) for initial in initial_paulis]
    final_paulis = [m.final for m in mapping_set]
    final_supports = [get_support(final) for final in final_paulis]
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
    calc_experiment_set(mapping_set::Vector{Mapping}, consistency_set::Vector{Vector{Int}})

Returns a set of experiments which simultaneously estimate the circuit eigenvalues corresponding to the mappings in `mapping_set`, whose tuple-consistency relations are described by `consistency_set`.
"""
function calc_experiment_set(
    mapping_set::Vector{Mapping},
    consistency_set::Vector{Vector{Int}},
)
    # Final Paulis
    final_paulis = [m.final for m in mapping_set]
    final_supports = [get_support(final) for final in final_paulis]
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
    calc_covariance_dict(c::AbstractCircuit, mapping_set::Vector{Mapping}, experiment_set::Vector{Vector{Int}}, full_covariance::Bool)

Returns a dictionary with the requisite information to calculate the sparse circuit eigenvalue estimator covariance matrix for the arranged circuit `c`, with mappings `mapping_set`, when estimated by the experiments in `experiment_set`.
If `full_covariance` is `true`, then this constructs the full covariance matrix, whereas if `full_covariance` is `false`, then this constructs only the diagonal terms of the covariance matrix.
"""
function calc_covariance_dict(
    c::T,
    mapping_set::Vector{Mapping},
    experiment_set::Vector{Vector{Int}},
    full_covariance::Bool,
) where {T <: AbstractCircuit}
    # Initialise variables
    tuple_experiment_number = length(experiment_set)
    trivial_pauli = Pauli(Bool[0], 0)
    zero_row = SparseVector{Int32, Int32}(0, Int32[], Int32[])
    trivial_row = SparseVector{Int32, Int32}(1, Int32[], Int32[])
    trivial_track = Vector{Int16}[]
    zero_mapping = Mapping(trivial_pauli, trivial_pauli, zero_row, trivial_track)
    trivial_mapping = Mapping(trivial_pauli, trivial_pauli, trivial_row, trivial_track)
    initial_paulis = [m.initial for m in mapping_set]
    # Determine diagonal terms in the covariance matrix
    covariance_dict = Dict{CartesianIndex{2}, Tuple{Mapping, Int}}()
    for j in 1:tuple_experiment_number
        experiment = experiment_set[j]
        # Determine the number of sign configurations as generated by PackCircuit
        experiment_prep_supports =
            [get_support(initial_paulis[pauli]) for pauli in experiment]
        max_prep_support = maximum(length.(experiment_prep_supports))
        @assert max_prep_support ∈ [1; 2] "Currently only supports up to two-qubit Pauli preparations."
        if hasproperty(c, :partition)
            experiment_sign_factor = 2^max_prep_support
        else
            experiment_sign_factor = 1
        end
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
        for j in 1:tuple_experiment_number
            experiment = experiment_set[j]
            S = length(experiment)
            # Determine the number of sign configurations as generated by PackCircuit
            experiment_prep_supports =
                [get_support(initial_paulis[pauli]) for pauli in experiment]
            max_prep_support = maximum(length.(experiment_prep_supports))
            @assert max_prep_support ∈ [1; 2] "Currently only supports up to two-qubit Pauli preparations."
            if hasproperty(c, :partition)
                experiment_sign_factor = 2^max_prep_support
            else
                experiment_sign_factor = 1
            end
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
                                sum_mapping = calc_mapping(sum_pauli, c)
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
    get_experiment_layers(c::AbstractCircuit, mapping_set::Vector{Mapping}, experiment_set::Vector{Vector{Int}})

Returns circuit layers that prepare and measure the initial and final Paulis as given by the mappings in `mapping_set`, for each experiment in `experiment set`.
If `c` has a `partition` field, this is used to generate all necessary sign configurations for the preparation layers.
"""
function get_experiment_layers(
    c::T,
    mapping_set::Vector{Mapping},
    experiment_set::Vector{Vector{Int}},
) where {T <: AbstractCircuit}
    # Generate an appropriate set of circuits based on the experiment set
    n = c.qubit_num
    tuple_experiment_number = length(experiment_set)
    prep_layer_set = Vector{Vector{Layer}}(undef, tuple_experiment_number)
    meas_layer_set = Vector{Layer}(undef, tuple_experiment_number)
    for j in 1:tuple_experiment_number
        # Set up some variables
        experiment = experiment_set[j]
        S = length(experiment)
        initial_set = [mapping_set[pauli].initial for pauli in experiment]
        initial_support_set = [get_support(initial) for initial in initial_set]
        final_set = [mapping_set[pauli].final for pauli in experiment]
        final_support_set = [get_support(final) for final in final_set]
        # Collate the appropriate preparation and measurement circuits
        prep = Vector{Gate}(undef, 0)
        meas = Vector{Gate}(undef, 0)
        for s in 1:S
            append!(prep, get_prep_layer(initial_set[s], initial_support_set[s]).layer)
            append!(meas, get_meas_layer(final_set[s], final_support_set[s]).layer)
        end
        # Remove repeated gates
        unique!(prep)
        unique!(meas)
        # Order the gates by the qubits on which they act
        sort!(prep)
        sort!(meas)
        # Store the measurement layer
        meas_layer_set[j] = Layer(meas, n)
        # Set up the eigenstate sign combinations
        if hasproperty(c, :partition)
            @assert length(unique([c.partition[1]; c.partition[2]])) == c.qubit_num "The partition does not contain the right number of qubits."
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
                prep_layer_set[j] = [Layer(prep, n), Layer(prep_2, n)]
            elseif max_prep_support == 2
                prep_2 = Vector{Gate}(undef, 0)
                prep_3 = Vector{Gate}(undef, 0)
                prep_4 = Vector{Gate}(undef, 0)
                for gate in prep
                    negative_type = replace(gate.type, "+" => "-")
                    negative_gate = Gate(negative_type, gate.index, gate.targets)
                    @assert gate.type != negative_type "The preparation gate $(gate) is unsigned."
                    if gate.targets[1] ∈ c.partition[1]
                        push!(prep_2, gate)
                        push!(prep_3, negative_gate)
                        push!(prep_4, negative_gate)
                    elseif gate.targets[1] ∈ c.partition[2]
                        push!(prep_2, negative_gate)
                        push!(prep_3, gate)
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
                prep_layer_set[j] =
                    [Layer(prep, n), Layer(prep_2, n), Layer(prep_3, n), Layer(prep_4, n)]
            else
                throw(error("Currently only supports up to two-qubit Pauli preparations."))
            end
        else
            # If the circuit cannot partition two qubit gates such that they always act 
            # between two sets of qubits, such as with data (partition[1]) and ancilla
            # (partition[2]) qubits, then we don't bother with the sign configurations.
            prep_layer_set[j] = [Layer(prep, n)]
        end
    end
    return (prep_layer_set::Vector{Vector{Layer}}, meas_layer_set::Vector{Layer})
end

"""
    generate_design(c::AbstractCircuit, tuple_set::Vector{Vector{Int}}; kwargs...)
    generate_design(c::AbstractCircuit, tuple_set_data::TupleSetData; kwargs...)
    generate_design(c::AbstractCircuit; kwargs...)

Returns a [`Design`](@ref) object containing all relevant information describing the experimental design, including the design matrix.

# Arguments

  - `c::AbstractCircuit`: Circuit for which the design matrix is to be generated.
  - `tuple_set::Vector{Vector{Int}}`: Tuple set arranging the circuit layers that is used to generate the experimental design.
  - `tuple_set_data::TupleSetData`: [`TupleSetData`](@ref) object that generates the tuple set.

# Keyword arguments

  - `shot_weights::Union{Vector{Float64}, Nothing} = nothing`: Shot weights for each tuple in the set, which must add to 1. When `nothing`, automatically generates the default shot weights.
  - `full_covariance::Bool = true`: If `true`, generates parameters to construct the full covariance matrix, else if `false`, only generates parameters to construct the terms on the diagonal.
  - `N_warn::Int = 3 * 10^4`: Number of circuit eigenvalues above which to warn the user about certain keyword argument choices.
  - `diagnostics::Bool = false`: Whether to print diagnostic information.
  - `save_data::Bool = false`: Whether to save the design data.
"""
function generate_design(
    c::T,
    tuple_set::Vector{Vector{Int}};
    shot_weights::Union{Vector{Float64}, Nothing} = nothing,
    full_covariance::Bool = true,
    N_warn::Int = 3 * 10^4,
    diagnostics::Bool = false,
    save_data::Bool = false,
    suppress_save_warning::Bool = false,
) where {T <: AbstractCircuit}
    # Set some parameters
    start_time = time()
    tuple_number = length(tuple_set)
    N = c.N
    @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
    if shot_weights !== nothing
        @assert length(shot_weights) == tuple_number "The number of shot weights does not match the tuple set."
        @assert sum(shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
        @assert all(shot_weights .> 0.0) "The shot weights are not all positive."
    end
    # Warn the user if they have unadvisable settings for a large circuit
    if N >= N_warn
        if full_covariance
            @warn "This design is for a very large circuit: generating the full covariance matrix is unadvised."
        end
        if ~diagnostics
            @warn "This design is for a very large circuit: turning on diagnostics is advised."
        end
        if ~save_data && ~suppress_save_warning
            @warn "This design is for a very large circuit: saving the data is advised."
        end
    end
    # Initialise the variables
    design_matrix = convert(SparseMatrixCSC{Int32, Int32}, spzeros(Int32, 0, N))
    mapping_ensemble = Vector{Vector{Mapping}}(undef, tuple_number)
    experiment_ensemble = Vector{Vector{Vector{Int}}}(undef, tuple_number)
    covariance_dict_ensemble =
        Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}(undef, tuple_number)
    prep_ensemble = Vector{Vector{Vector{Layer}}}(undef, tuple_number)
    meas_ensemble = Vector{Vector{Layer}}(undef, tuple_number)
    calculation_times = Matrix{Float64}(undef, tuple_number, 5)
    for idx in 1:tuple_number
        # Determine the design data for the tuple
        circuit_tuple = tuple_set[idx]
        c_tuple = apply_tuple(c, circuit_tuple)
        time_1 = time()
        (mapping_set, mapping_matrix) = calc_mapping_set(c_tuple)
        time_2 = time()
        consistency_set = calc_consistency_set(mapping_set)
        time_3 = time()
        experiment_set = calc_experiment_set(mapping_set, consistency_set)
        time_4 = time()
        covariance_dict =
            calc_covariance_dict(c_tuple, mapping_set, experiment_set, full_covariance)
        time_5 = time()
        (prep_layer_set, meas_layer_set) =
            get_experiment_layers(c_tuple, mapping_set, experiment_set)
        time_6 = time()
        # Track the times taken
        mapping_time = time_2 - time_1
        consistency_time = time_3 - time_2
        pauli_time = time_4 - time_3
        covariance_time = time_5 - time_4
        circuit_time = time_6 - time_5
        if diagnostics
            println(
                "For tuple $(idx), the mappings took $(round(mapping_time, digits = 3)) s, the consistency sets took $(round(consistency_time, digits = 3)) s, packing the Paulis took $(round(pauli_time, digits = 3)) s, the covariance dictionary took $(round(covariance_time, digits = 3)) s, and the circuits took $(round(circuit_time, digits = 3)) s. Since starting, $(round(time_5 - start_time, digits = 3)) s have elapsed.",
            )
        end
        # Store the data
        tuple_set[idx] = circuit_tuple
        design_matrix = vcat(design_matrix, mapping_matrix)
        mapping_ensemble[idx] = mapping_set
        experiment_ensemble[idx] = experiment_set
        covariance_dict_ensemble[idx] = covariance_dict
        prep_ensemble[idx] = prep_layer_set
        meas_ensemble[idx] = meas_layer_set
        calculation_times[idx, :] = [
            mapping_time
            consistency_time
            pauli_time
            covariance_time
            circuit_time
        ]
    end
    # Calculate the tuple times and shot weights
    experiment_numbers =
        length.([vcat(prep_layer_set...) for prep_layer_set in prep_ensemble])
    (tuple_times, tuple_shot_weights) =
        get_tuple_set_params(c, tuple_set, experiment_numbers)
    if shot_weights === nothing
        shot_weights = tuple_shot_weights
    end
    # Save and return the results
    overall_time = time() - start_time
    d = Design(
        c,
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
            "Generated the design for all $(tuple_number) tuples in $(round(overall_time, digits = 3)) s.",
        )
    end
    return d::Design
end

function generate_design(
    c::T,
    tuple_set_data::TupleSetData;
    shot_weights::Union{Vector{Float64}, Nothing} = nothing,
    full_covariance::Bool = true,
    diagnostics::Bool = false,
    save_data::Bool = false,
) where {T <: AbstractCircuit}
    # Generate the design
    tuple_set = get_tuple_set(tuple_set_data)
    d = generate_design(
        c,
        tuple_set;
        shot_weights = shot_weights,
        full_covariance = full_covariance,
        diagnostics = diagnostics,
        save_data = false,
        suppress_save_warning = true,
    )
    @reset d.tuple_set_data = tuple_set_data
    # Save and return the results
    if save_data
        save_design(d)
    end
    return d::Design
end

function generate_design(
    c::T;
    shot_weights::Union{Vector{Float64}, Nothing} = nothing,
    full_covariance::Bool = true,
    diagnostics::Bool = false,
    save_data::Bool = false,
) where {T <: AbstractCircuit}
    # Generate the design
    tuple_set_data = get_tuple_set_data(c)
    d = generate_design(
        c,
        tuple_set_data;
        shot_weights = shot_weights,
        full_covariance = full_covariance,
        diagnostics = diagnostics,
        save_data = save_data,
    )
    # Return the results
    return d::Design
end

"""
    complete_design(d::Design; diagnostics::Bool = false)

Returns a copy of the design `d` where the covariance matrix dictionary generates the full covariance matrix. Prints diagnostics if `diagnostics` is `true`.
"""
function complete_design(d::Design; diagnostics::Bool = false)
    if d.full_covariance
        if diagnostics
            println("The supplied design matrix already has a full covariance matrix.")
        end
        return d::Design
    else
        # Set some parameters
        start_time = time()
        full_covariance = true
        tuple_number = length(d.tuple_set)
        # Regenerate the covariance dictionary set with a full covariance matrix
        covariance_dict_ensemble =
            Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}(undef, tuple_number)
        covariance_times = Vector{Float64}(undef, tuple_number)
        for idx in 1:tuple_number
            circuit_tuple = d.tuple_set[idx]
            c_tuple = apply_tuple(d.c, circuit_tuple)
            mapping_set_tuple = d.mapping_ensemble[idx]
            experiment_set_tuple = d.experiment_ensemble[idx]
            # Calculate the covariance dictionary
            time_4 = time()
            covariance_dict_tuple = calc_covariance_dict(
                c_tuple,
                mapping_set_tuple,
                experiment_set_tuple,
                full_covariance,
            )
            time_5 = time()
            # Track the times taken
            covariance_time = time_5 - time_4
            if diagnostics
                println(
                    "For tuple $(idx), generating the covariance dictionary took $(round(covariance_time, digits = 3)) s. Since starting, $(round(time_5 - start_time, digits = 3)) s have elapsed.",
                )
            end
            # Store the data
            covariance_dict_ensemble[idx] = covariance_dict_tuple
            covariance_times[idx] = covariance_time
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
                "Generated the full covariance dictionaries for all $(tuple_number) tuples in $(round(overall_time, digits = 3)) s.",
            )
        end
        return d_complete::Design
    end
end

"""
    update_noise(d::Design, noise_param::AbstractNoiseParameters)

Returns a copy of `design` where the circuit has been updated with noise generated according to `noise_param`.
"""
function update_noise(d::Design, noise_param::T) where {T <: AbstractNoiseParameters}
    # Generate the noise
    gate_probabilities = get_gate_probabilities(d.c.total_gates, noise_param)
    gate_eigenvalues =
        get_gate_eigenvalues(gate_probabilities, d.c.total_gates, d.c.gate_index, d.c.N)
    # Update the circuit
    d_update = deepcopy(d)
    @reset d_update.c.noise_param = noise_param
    @reset d_update.c.gate_probabilities = gate_probabilities
    @reset d_update.c.gate_eigenvalues = gate_eigenvalues
    return d_update::Design
end
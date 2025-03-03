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

function Base.show(io::IO, p::Pauli)
    print(io, get_pauli_string(p))
    return nothing
end

function Base.:(+)(p_1::Pauli, p_2::Pauli)
    @assert p_1.qubit_num == p_2.qubit_num "The Paulis $(p_1) and $(p_2) do not have the same number of qubits."
    return Pauli(convert(Vector{Bool}, (p_1.pauli .+ p_2.pauli) .% 2), p_1.qubit_num)
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
    print(io, get_pauli_string(m.initial) * " => " * get_pauli_string(m.final))
    return nothing
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
  - `prep_ensemble::Vector{Vector{Layer}}`: Vector of [`Layer`](@ref) objects that prepare qubits in Pauli eigenstates for each experiment in the set, for each tuple in the set.
  - `meas_ensemble::Vector{Vector{Layer}}`: Vector of [`Layer`](@ref) objects that measure qubits in Pauli bases for each experiment in the set, for each tuple in the set.
  - `tuple_times::Vector{Float64}`: Time taken to implement the circuit arranged by each tuple in the set, normalised according to the time factor for the basic tuple set.
  - `shot_weights::Vector{Float64}`: Shot weights for each tuple in the set, which add to 1.
  - `experiment_numbers::Vector{Int}`: Number of experiments for each tuple in the set.
  - `experiment_number::Int`: Total number of experiments.
  - `calculation_times::Matrix{Float64}`: Time taken to generate components of the design for each tuple, which correspond to generating: the mappings, the sets of tuple-consistent Pauli preparations, the experiment sets, and the covariance matrix dictionary.
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
    prep_ensemble::Vector{Vector{Layer}}
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
        prep_ensemble::Vector{Vector{Layer}},
        meas_ensemble::Vector{Vector{Layer}},
        tuple_times::Vector{Float64},
        shot_weights::Vector{Float64},
        experiment_numbers::Vector{Int},
        experiment_number::Integer,
        calculation_times::Matrix{Float64},
        overall_time::Float64,
        optimisation_time::Float64,
        ls_type::Symbol,
    ) where {T <: AbstractCircuit}
        # Check parameters
        tuple_number = length(tuple_set)
        @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
        @assert tuple_set == get_tuple_set(tuple_set_data) "The tuple set does not align with the tuple set data."
        @assert length(mapping_ensemble) == tuple_number "The size of the mapping ensemble does not match the tuple set."
        @assert length(experiment_ensemble) == tuple_number "The size of the experiment ensemble does not match the tuple set."
        @assert length(covariance_dict_ensemble) == tuple_number "The size of the covariance dictionary ensemble does not match the tuple set."
        @assert length(prep_ensemble) == tuple_number "The size of the preparation ensemble does not match the tuple set."
        @assert length(meas_ensemble) == tuple_number "The size of the measurement ensemble does not match the tuple set."
        @assert length.(prep_ensemble) == length.(meas_ensemble) "The preparation and measurement ensembles do not match."
        @assert length(tuple_times) == tuple_number "The number of tuple times does not match the tuple set."
        @assert length(shot_weights) == tuple_number "The number of shot weights does not match the tuple set."
        @assert sum(shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
        @assert all(shot_weights .> 0.0) "The shot weights are not all positive."
        @assert length(experiment_numbers) == tuple_number "The number of experiment numbers does not match the tuple set."
        @assert experiment_number == sum(experiment_numbers) "The experiment number $(experiment_number) does not match the sum of the experiment numbers $(sum(experiment_numbers))."
        @assert size(calculation_times) == (tuple_number, 4) "The calculation times do not match the tuple set."
        @assert ls_type ∈ [:none; :gls; :wls; :ols] "The least squares type $(ls_type) is not supported."
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
        prep_ensemble::Vector{Vector{Layer}},
        meas_ensemble::Vector{Vector{Layer}},
        tuple_times::Vector{Float64},
        shot_weights::Vector{Float64},
        calculation_times::Matrix{Float64},
        overall_time::Float64;
        ls_type::Symbol = :none,
    ) where {T <: AbstractCircuit}
        # Initialise parameters
        tuple_number = length(tuple_set)
        tuple_set_data = TupleSetData(tuple_set, Vector{Int}[], Int[], Tuple{Int, Int}[])
        experiment_numbers = length.(prep_ensemble)
        experiment_number = sum(experiment_numbers)
        optimisation_time = 0.0
        # Check parameters
        @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
        @assert length(mapping_ensemble) == tuple_number "The size of the mapping ensemble does not match the tuple set."
        @assert length(experiment_ensemble) == tuple_number "The size of the experiment ensemble does not match the tuple set."
        @assert length(covariance_dict_ensemble) == tuple_number "The size of the covariance dictionary ensemble does not match the tuple set."
        @assert length(prep_ensemble) == tuple_number "The size of the preparation ensemble does not match the tuple set."
        @assert length(meas_ensemble) == tuple_number "The size of the measurement ensemble does not match the tuple set."
        @assert length.(prep_ensemble) == length.(meas_ensemble) "The preparation and measurement ensembles do not match."
        @assert length(tuple_times) == tuple_number "The number of tuple times does not match the tuple set."
        @assert length(shot_weights) == tuple_number "The number of shot weights does not match the tuple set."
        @assert sum(shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
        @assert all(shot_weights .> 0.0) "The shot weights are not all positive."
        @assert length(experiment_numbers) == tuple_number "The number of experiment numbers does not match the tuple set."
        @assert experiment_number == sum(experiment_numbers) "The experiment number $(experiment_number) does not match the sum of the experiment numbers $(sum(experiment_numbers))."
        @assert size(calculation_times) == (tuple_number, 4) "The calculation times do not match the tuple set."
        @assert ls_type ∈ [:none; :gls; :wls; :ols] "The least squares type $(ls_type) is not supported."
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
    println(
        io,
        "Design for a $(d.c.circuit_param.circuit_name) circuit with $(length(d.tuple_set)) tuples and $(d.experiment_number) experiments.",
    )
    pretty_print(io, d)
    return nothing
end

@struct_hash_equal_isequal Design

"""
    get_eigenvalues(d::Design)
    get_eigenvalues(d::Design, gate_eigenvalues::Vector{Float64})

Returns the circuit eigenvalues of the design `d`, optionally using the gate eigenvalues `gate_eigenvalues`.
"""
function get_eigenvalues(d::Design, gate_eigenvalues::Vector{Float64})
    eigenvalues = exp.(d.matrix * log.(gate_eigenvalues))
    return eigenvalues::Vector{Float64}
end
function get_eigenvalues(d::Design)
    eigenvalues = get_eigenvalues(d, get_gate_eigenvalues(d))
    return eigenvalues::Vector{Float64}
end

"""
    get_eigenvalues_ensemble(d::Design)
    get_eigenvalues_ensemble(d::Design, gate_eigenvalues::Vector{Float64})

Returns the circuit eigenvalue ensemble of the design `d`, optionally using the gate eigenvalues `gate_eigenvalues`.
"""
function get_eigenvalues_ensemble(d::Design, gate_eigenvalues::Vector{Float64})
    # Initialise variables
    tuple_number = length(d.tuple_set)
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([0; mapping_lengths[1:(end - 1)]])
    # Get the eigenvalues and ensemble
    eigenvalues = get_eigenvalues(d, gate_eigenvalues)
    eigenvalues_ensemble = [
        [eigenvalues[mapping_lower[i] + j] for j in 1:mapping_lengths[i]] for
        i in 1:tuple_number
    ]
    @assert eigenvalues == vcat(eigenvalues_ensemble...) "The eigenvalues do not match the ensemble."
    return eigenvalues_ensemble::Vector{Vector{Float64}}
end
function get_eigenvalues_ensemble(d::Design)
    eigenvalues_ensemble = get_eigenvalues_ensemble(d, get_gate_eigenvalues(d))
    return eigenvalues_ensemble::Vector{Vector{Float64}}
end

"""
    get_gate_eigenvalues(d::Design)

Returns the gate eigenvalues of the design `d`.
"""
function get_gate_eigenvalues(d::Design)
    gate_eigenvalues = d.c.gate_eigenvalues
    return gate_eigenvalues::Vector{Float64}
end

"""
    get_marginal_gate_eigenvalues(d::Design)

Returns the gate eigenvalues of the design `d` marginalised over gate orbits.
"""
function get_marginal_gate_eigenvalues(d::Design)
    marginal_gate_eigenvalues =
        get_marginal_gate_eigenvalues(get_gate_eigenvalues(d), d.c.gate_data)
    return marginal_gate_eigenvalues::Vector{Float64}
end

"""
    get_relative_gate_eigenvalues(d::Design)

Returns the gate eigenvalues of the design `d` marginalised over gate orbits which can be estimated to relative precision.
"""
function get_relative_gate_eigenvalues(d::Design)
    relative_gate_eigenvalues =
        get_relative_gate_eigenvalues(get_gate_eigenvalues(d), d.c.gate_data)
    return relative_gate_eigenvalues::Vector{Float64}
end

"""
    get_gate_probabilities(d::Design)

Returns the gate error probabilities of the design `d`.
"""
function get_gate_probabilities(d::Design)
    gate_probabilities = d.c.gate_probabilities
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_marginal_gate_probabilities(d::Design)

Returns the gate error probabilities of the design `d` marginalised over gate orbits.
"""
function get_marginal_gate_probabilities(d::Design)
    marginal_gate_probabilities = get_marginal_gate_probabilities(get_gate_probabilities(d))
    return marginal_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_relative_gate_probabilities(d::Design)

Returns the gate error probabilities of the design `d` marginalised over gate orbits which can be estimated to relative precision.
"""
function get_relative_gate_probabilities(d::Design)
    relative_gate_probabilities =
        get_relative_gate_probabilities(get_gate_probabilities(d), d.c.gate_data)
    return relative_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    get_pauli_prep_set(gates::Vector{Gate}, n::Integer)

Returns all weight-1 Paulis on all `n` qubits, as well as the weight-2 Paulis supported on some gate in `gates`.
"""
function get_pauli_prep_set(gates::Vector{Gate}, n::Integer)
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
    weight_2_supports = Vector{Vector{Int}}()
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
    for (idx, i) in pairs(initial_support)
        signature = Bool[initial.pauli[i]; initial.pauli[n + i]]
        if signature == Bool[0; 1]
            prep[idx] = Gate("PZ", 0, [i])
        elseif signature == Bool[1; 0]
            prep[idx] = Gate("PX", 0, [i])
        elseif signature == Bool[1; 1]
            prep[idx] = Gate("PY", 0, [i])
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
    for (idx, i) in pairs(final_support)
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
    update_design_row!(design_row::Vector{Int}, pauli::Pauli, l::Layer, gate_index_dict::Dict{Gate, Int})

Updates the design matrix row `design_row` according to the Pauli `pauli` for each gate in the layer `l`, using the gate index dictionary `gate_index_dict`.
"""
function update_design_row!(
    design_row::Vector{Int},
    pauli::Pauli,
    l::Layer,
    gate_index_dict::Dict{Gate, Int},
)
    n = pauli.qubit_num
    support = get_support(pauli)
    # Calculate the relevant additions to the design matrix row for each gate in the layer
    for gate in l.layer
        if length(intersect(support, gate.targets)) > 0
            if is_spam(gate)
                # Check that the appropriate preparation or measurement is being performed
                target_1 = gate.targets[1]
                signature = Bool[pauli.pauli[target_1]; pauli.pauli[n + target_1]]
                wrong_Z = (signature == Bool[0; 1] && gate.type ∉ ["PZ"; "MZ"])
                wrong_X = (signature == Bool[1; 0] && gate.type ∉ ["PX"; "MX"])
                wrong_Y = (signature == Bool[1; 1] && gate.type ∉ ["PY"; "MY"])
                @assert ~wrong_Z && ~wrong_X && ~wrong_Y "The Pauli $(pauli) is not prepared or measured by the gate $(gate)."
                @assert signature != Bool[0; 0] "The Pauli $(pauli) is meant to be supported on qubit $(target_1) and yet has signature $(signature) there."
                # Preparations and measurements get unique variables
                pauli_idx = 1
            elseif is_mid_reset(gate)
                # Zero out the design row
                for idx in eachindex(design_row)
                    design_row[idx] = 0
                end
                # Resets get unique variables
                pauli_idx = 1
            elseif is_two_qubit(gate)
                target_1 = gate.targets[1]
                target_2 = gate.targets[2]
                # The index orders 2-qubit Paulis as:
                # (II=0), XI=1,   IX=2,   XX=3
                # ZI=4,   YI=5,   ZX=6,   YX=7
                # IZ=8,   XZ=9,   IY=10,  XY=11
                # ZZ=12,  YZ=13,  ZY=14,  YY=15
                # This is the natural bit string ordering
                pauli_idx =
                    pauli.pauli[target_1] +
                    2 * pauli.pauli[target_2] +
                    4 * pauli.pauli[n + target_1] +
                    8 * pauli.pauli[n + target_2]
            elseif is_one_qubit(gate)
                target_1 = gate.targets[1]
                # The index orders 1-qubit Paulis as:
                # (I=0),  X=1,    Z=2,    Y=3
                # This is the natural bit string ordering
                pauli_idx = pauli.pauli[target_1] + 2 * pauli.pauli[n + target_1]
            else
                throw(error("The gate $(gate) does not operate on either 1 or 2 qubits."))
            end
            design_row[gate_index_dict[gate] + pauli_idx] += 1
        end
    end
    return nothing
end

"""
    calc_mapping(initial::Pauli, gate_index_dict::Dict{Gate, Int}, c::AbstractCircuit)

Returns a [`Mapping`](@ref) object for the Pauli `initial` when mapped by the circuit `c` with gate index dictionary `gate_index_dict`.
"""
function calc_mapping(
    initial::Pauli,
    gate_index_dict::Dict{Gate, Int},
    c::T,
) where {T <: AbstractCircuit}
    # Initialise variables
    n = c.qubit_num
    t = Tableau(n)
    design_row = zeros(Int, c.gate_data.N)
    # Generate the layer which prepares the appropriate Pauli
    @assert initial.pauli[2n + 1] == 0 "This function requires an initial Pauli with positive sign."
    initial_support = get_support(initial)
    @assert length(initial_support) > 0 "The initial Pauli is not supported on any qubits."
    prep_layer = get_prep_layer(initial, initial_support)
    # Update the design matrix row if appropriate
    if c.noisy_prep
        update_design_row!(design_row, initial, prep_layer, gate_index_dict)
    end
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
    spread_track = Vector{Vector{Int16}}()
    tuple_circuit = c.circuit[c.circuit_tuple]
    for l in tuple_circuit
        # Use row_sum! to assemble the Pauli
        for i in eachindex(initial_support[2:end])
            row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
        end
        # Track the spread of the Pauli
        layer_pauli = Pauli(t.tableau[n + initial_support[1], :], n)
        push!(spread_track, get_support(layer_pauli))
        # Update the design matrix row
        update_design_row!(design_row, layer_pauli, l, gate_index_dict)
        # Use row_sum! to disassemble the Pauli
        for i in eachindex(initial_support[2:end])
            row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
        end
        # Apply the circuit layer
        apply!(t, l)
    end
    # Use row_sum! to assemble the Pauli
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    # Determine the final Pauli to which the initial Pauli is mapped
    final = Pauli(t.tableau[n + initial_support[1], :], n)
    final_support = get_support(final)
    push!(spread_track, final_support)
    # Use row_sum! to disassemble the Pauli
    for i in eachindex(initial_support[2:end])
        row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
    end
    # Generate the layer which measures the appropriate Pauli
    meas_layer = get_meas_layer(final, final_support)
    if c.noisy_meas
        update_design_row!(design_row, final, meas_layer, gate_index_dict)
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
    gate_index_dict = get_gate_index_dict(c.gate_data)
    mapping_set = Vector{Mapping}(undef, pauli_num)
    @threads :static for idx in 1:pauli_num
        mapping_set[idx] = calc_mapping(pauli_prep_set[idx], gate_index_dict, c)
    end
    # Design matrix row corresponding to the mapping set
    mapping_matrix =
        convert(SparseMatrixCSC{Int32, Int32}, spzeros(Int32, pauli_num, c.gate_data.N))
    for (idx, m) in pairs(mapping_set)
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
    experiment_set = Vector{Vector{Int}}()
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
            best_pauli_idx = addable_paulis[best_add_idx]
            # Update all relevant sets
            best_unadd_idx = findfirst(x -> x == best_pauli_idx, unadded_paulis)
            push!(experiment, best_pauli_idx)
            union!(set_support, final_supports[best_pauli_idx])
            deleteat!(addable_paulis, best_add_idx)
            intersect!(addable_paulis, consistency_set[best_pauli_idx])
            deleteat!(unadded_paulis, best_unadd_idx)
        end
        # The compatible Paulis are the Paulis compatible with all the Paulis currently added to the set
        compatible_paulis = setdiff(total_paulis, experiment)
        for pauli_idx in experiment
            intersect!(compatible_paulis, consistency_set[pauli_idx])
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
            best_pauli_idx = compatible_paulis[best_comp_idx]
            # Update all relevant sets
            push!(experiment, best_pauli_idx)
            union!(set_support, final_supports[best_pauli_idx])
            deleteat!(compatible_paulis, best_comp_idx)
            intersect!(compatible_paulis, consistency_set[best_pauli_idx])
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

BEWARE: This currently assumes that circuit eigenvalues for Paulis supported on ``n`` qubits require preparing exactly and only ``2^n`` sign configurations.
Errors will occur if this is not the case, and they are not guaranteed to be noisy.
If creating a circuit where this is not the case, you will need to provide new methods for this function.
"""
function calc_covariance_dict(
    c::T,
    mapping_set::Vector{Mapping},
    experiment_set::Vector{Vector{Int}},
    full_covariance::Bool,
) where {T <: AbstractCircuit}
    # Initialise variables
    experiment_number = length(experiment_set)
    trivial_pauli = Pauli(Bool[0], 0)
    zero_row = SparseVector{Int32, Int32}(0, Int32[], Int32[])
    trivial_row = SparseVector{Int32, Int32}(1, Int32[], Int32[])
    trivial_track = Vector{Int16}[]
    zero_mapping = Mapping(trivial_pauli, trivial_pauli, zero_row, trivial_track)
    trivial_mapping = Mapping(trivial_pauli, trivial_pauli, trivial_row, trivial_track)
    initial_paulis = [m.initial for m in mapping_set]
    gate_index_dict = get_gate_index_dict(c.gate_data)
    # Determine diagonal terms in the covariance matrix
    covariance_dict = Dict{CartesianIndex{2}, Tuple{Mapping, Int}}()
    for j in 1:experiment_number
        experiment = experiment_set[j]
        for pauli_idx in experiment
            # Calculate the diagonal term
            diag_pauli_idx = CartesianIndex(pauli_idx, pauli_idx)
            if haskey(covariance_dict, diag_pauli_idx)
                covariance_dict[diag_pauli_idx] = (
                    covariance_dict[diag_pauli_idx][1],
                    covariance_dict[diag_pauli_idx][2] + 1,
                )
            else
                covariance_dict[diag_pauli_idx] = (trivial_mapping, 1)
            end
        end
    end
    # Determine off-diagonal terms in the covariance matrix if required
    # Only store for the upper diagonal as the covariance matrix is symmetric
    if full_covariance
        # Initialise variables
        L = length(c.circuit_tuple)
        reentrant_lock = ReentrantLock()
        design_rows = [m.design_row for m in mapping_set]
        spread_tracks = [m.spread_track for m in mapping_set]
        # Calculate the terms
        for j in 1:experiment_number
            experiment = experiment_set[j]
            S = length(experiment)
            @threads :static for s_1 in 1:S
                # Get the first Pauli
                pauli_idx_1 = experiment[s_1]
                pauli_1 = initial_paulis[pauli_idx_1]
                design_row_1 = design_rows[pauli_idx_1]
                spread_track_1 = spread_tracks[pauli_idx_1]
                for s_2 in (s_1 + 1):S
                    # Get the second Pauli
                    pauli_idx_2 = experiment[s_2]
                    pauli_2 = initial_paulis[pauli_idx_2]
                    design_row_2 = design_rows[pauli_idx_2]
                    spread_track_2 = spread_tracks[pauli_idx_2]
                    # Check if the design matrix row has previously been computed
                    pauli_idx = CartesianIndex(pauli_idx_1, pauli_idx_2)
                    if haskey(covariance_dict, pauli_idx)
                        # If the design matrix row is non-trivial, increment the shots counter
                        lock(reentrant_lock) do
                            if covariance_dict[pauli_idx][1] != zero_mapping
                                covariance_dict[pauli_idx] = (
                                    covariance_dict[pauli_idx][1],
                                    covariance_dict[pauli_idx][2] + 1,
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
                                sum_mapping = calc_mapping(sum_pauli, gate_index_dict, c)
                            else
                                sum_mapping = mapping_set[sum_index]
                            end
                            # Add the mapping if non-trivial
                            lock(reentrant_lock) do
                                if sum_mapping.design_row != design_row_1 + design_row_2
                                    covariance_dict[pauli_idx] = (sum_mapping, 1)
                                else
                                    covariance_dict[pauli_idx] = (zero_mapping, 0)
                                end
                            end
                        else
                            lock(reentrant_lock) do
                                covariance_dict[pauli_idx] = (zero_mapping, 0)
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
    experiment_number = length(experiment_set)
    prep_layer_set = Vector{Layer}(undef, experiment_number)
    meas_layer_set = Vector{Layer}(undef, experiment_number)
    for j in 1:experiment_number
        # Set up some variables
        experiment = experiment_set[j]
        S = length(experiment)
        initial_set = [mapping_set[pauli_idx].initial for pauli_idx in experiment]
        initial_support_set = [get_support(initial) for initial in initial_set]
        final_set = [mapping_set[pauli_idx].final for pauli_idx in experiment]
        final_support_set = [get_support(final) for final in final_set]
        # Collate the appropriate preparation and measurement circuits
        prep = Vector{Gate}()
        meas = Vector{Gate}()
        for s in 1:S
            append!(prep, get_prep_layer(initial_set[s], initial_support_set[s]).layer)
            append!(meas, get_meas_layer(final_set[s], final_support_set[s]).layer)
        end
        # Remove repeated gates
        unique!(prep)
        unique!(meas)
        # Pad the layers with Z preparations and measurements
        non_prep_targets = setdiff(collect(1:n), vcat([gate.targets for gate in prep]...))
        for i in non_prep_targets
            push!(prep, Gate("PZ", 0, [i]))
        end
        non_meas_targets = setdiff(collect(1:n), vcat([gate.targets for gate in meas]...))
        for i in non_meas_targets
            push!(meas, Gate("MZ", 0, [i]))
        end
        # Order the gates by the qubits on which they act
        sort!(prep; by = x -> x.targets)
        sort!(meas; by = x -> x.targets)
        # Store the measurement layer
        meas_layer_set[j] = Layer(meas, n)
        # Set up the eigenstate sign combinations
        prep_layer_set[j] = Layer(prep, n)
    end
    return (prep_layer_set::Vector{Layer}, meas_layer_set::Vector{Layer})
end

"""
    generate_design(c::AbstractCircuit, tuple_set::Vector{Vector{Int}}; kwargs...)
    generate_design(c::AbstractCircuit, tuple_set_data::TupleSetData; kwargs...)
    generate_design(c::AbstractCircuit, d::Design; kwargs...)
    generate_design(c::AbstractCircuit; kwargs...)

Returns a [`Design`](@ref) object containing all relevant information describing the experimental design, including the design matrix.

# Arguments

  - `c::AbstractCircuit`: Circuit for which the design matrix is to be generated.
  - `tuple_set::Vector{Vector{Int}}`: Tuple set arranging the circuit layers that is used to generate the experimental design.
  - `tuple_set_data::TupleSetData`: [`TupleSetData`](@ref) object that generates the tuple set.
  - `d::Design`: Old design object whose parameters are used to generate the design for the circuit `c`.

# Keyword arguments

  - `shot_weights::Union{Vector{Float64}, Nothing} = nothing`: Shot weights for each tuple in the set, which must add to 1. When `nothing`, automatically generates the default shot weights.
  - `N_warn::Integer = 10^5`: Number of circuit eigenvalues above which to warn the user about certain keyword argument choices.
  - `full_covariance::Bool = (c.gate_data.N < N_warn ? true : false)`: If `true`, generates parameters to construct the full covariance matrix, else if `false`, only generates parameters to construct the terms on the diagonal.
  - `add_circuit::Bool = false`: If tuple set data has not been supplied, the circuit itself is added to the repeat tuple set.
  - `repeat_points::Integer = 1`: If tuple set data has not been supplied, the each repeat number is augmented to have `repeat_points` repeat numbers.
  - `weight_experiments::Bool = true`: Whether to weight the shot weights for a tuple by the number of experiments for that tuple.
  - `diagnostics::Bool = false`: Whether to print diagnostic information.
  - `save_data::Bool = false`: Whether to save the design data.
"""
function generate_design(
    c::T,
    tuple_set::Vector{Vector{Int}};
    shot_weights::Union{Vector{Float64}, Nothing} = nothing,
    N_warn::Integer = 10^5,
    full_covariance::Bool = (c.gate_data.N < N_warn ? true : false),
    weight_experiments::Bool = true,
    diagnostics::Bool = false,
    save_data::Bool = false,
    suppress_save_warning::Bool = false,
) where {T <: AbstractCircuit}
    # Set some parameters
    start_time = time()
    tuple_number = length(tuple_set)
    N = c.gate_data.N
    check_tuple!(c, tuple_set)
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
    prep_ensemble = Vector{Vector{Layer}}(undef, tuple_number)
    meas_ensemble = Vector{Vector{Layer}}(undef, tuple_number)
    calculation_times = Matrix{Float64}(undef, tuple_number, 4)
    for i in 1:tuple_number
        # Determine the design data for the tuple
        circuit_tuple = tuple_set[i]
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
        if diagnostics
            println(
                "For tuple $(i), the mappings took $(round(mapping_time, digits = 2)) s, the consistency sets took $(round(consistency_time, digits = 2)) s, packing the Paulis took $(round(pauli_time, digits = 2)) s, and the covariance dictionary took $(round(covariance_time, digits = 2)) s. Since starting, $(round(time_6 - start_time, digits = 2)) s have elapsed.",
            )
        end
        # Store the data
        design_matrix = vcat(design_matrix, mapping_matrix)
        mapping_ensemble[i] = mapping_set
        experiment_ensemble[i] = experiment_set
        covariance_dict_ensemble[i] = covariance_dict
        prep_ensemble[i] = prep_layer_set
        meas_ensemble[i] = meas_layer_set
        calculation_times[i, :] = [
            mapping_time
            consistency_time
            pauli_time
            covariance_time
        ]
        if N >= N_warn
            GC.gc()
        end
    end
    # Calculate the tuple times and shot weights
    experiment_numbers = length.(prep_ensemble)
    (tuple_times, tuple_shot_weights) = get_tuple_set_params(
        c,
        tuple_set,
        experiment_numbers;
        weight_experiments = weight_experiments,
    )
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
            "Generated the design for all $(tuple_number) tuples in $(round(overall_time, digits = 2)) s.",
        )
    end
    return d::Design
end
function generate_design(
    c::T,
    tuple_set_data::TupleSetData;
    shot_weights::Union{Vector{Float64}, Nothing} = nothing,
    N_warn::Integer = 10^5,
    full_covariance::Bool = (c.gate_data.N < N_warn ? true : false),
    weight_experiments::Bool = true,
    diagnostics::Bool = false,
    save_data::Bool = false,
) where {T <: AbstractCircuit}
    # Generate the design
    tuple_set = get_tuple_set(tuple_set_data)
    d = generate_design(
        c,
        tuple_set;
        shot_weights = shot_weights,
        N_warn = N_warn,
        full_covariance = full_covariance,
        weight_experiments = weight_experiments,
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
    N_warn::Integer = 10^5,
    full_covariance::Bool = (c.gate_data.N < N_warn ? true : false),
    add_circuit::Bool = false,
    repeat_points::Integer = 1,
    weight_experiments::Bool = true,
    diagnostics::Bool = false,
    save_data::Bool = false,
) where {T <: AbstractCircuit}
    # Generate the design
    tuple_set_data = get_tuple_set_data(c; add_circuit = add_circuit)
    if repeat_points > 1
        tuple_set_data = get_augmented_tuple_set_data(tuple_set_data, repeat_points)
    end
    d = generate_design(
        c,
        tuple_set_data;
        shot_weights = shot_weights,
        N_warn = N_warn,
        full_covariance = full_covariance,
        weight_experiments = weight_experiments,
        diagnostics = diagnostics,
        save_data = save_data,
    )
    # Return the results
    return d::Design
end
function generate_design(
    c::T,
    d::Design;
    N_warn::Integer = 10^5,
    full_covariance::Bool = (c.gate_data.N < N_warn ? true : false),
    weight_experiments::Bool = true,
    diagnostics::Bool = false,
    save_data::Bool = false,
) where {T <: AbstractCircuit}
    # Generate the design
    d_new = generate_design(
        c,
        d.tuple_set_data;
        shot_weights = d.shot_weights,
        N_warn = N_warn,
        full_covariance = full_covariance,
        weight_experiments = weight_experiments,
        diagnostics = diagnostics,
        save_data = save_data,
    )
    # Return the results
    return d_new::Design
end

"""
    get_full_design(d::Design; diagnostics::Bool = false)

Returns a copy of the design `d` where the covariance matrix dictionary generates the full covariance matrix, printing diagnostics if `diagnostics` is `true`.
"""
function get_full_design(d::Design; diagnostics::Bool = false)
    if d.full_covariance
        if diagnostics
            println("The supplied design already has a full covariance matrix.")
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
        for i in 1:tuple_number
            circuit_tuple = d.tuple_set[i]
            c_tuple = apply_tuple(d.c, circuit_tuple)
            mapping_set_tuple = d.mapping_ensemble[i]
            experiment_set_tuple = d.experiment_ensemble[i]
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
                    "For tuple $(i), generating the covariance dictionary took $(round(covariance_time, digits = 3)) s. Since starting, $(round(time_5 - start_time, digits = 3)) s have elapsed.",
                )
            end
            # Store the data
            covariance_dict_ensemble[i] = covariance_dict_tuple
            covariance_times[i] = covariance_time
        end
        overall_time = time() - start_time
        # Set the new variables 
        d_full = deepcopy(d)
        @reset d_full.full_covariance = full_covariance
        @reset d_full.covariance_dict_ensemble = covariance_dict_ensemble
        @reset d_full.overall_time =
            d_full.overall_time - sum(d_full.calculation_times[:, 4]) + overall_time
        @reset d_full.calculation_times[:, 4] = covariance_times
        if diagnostics
            println(
                "Generated the full covariance dictionaries for all $(tuple_number) tuples in $(round(overall_time, digits = 3)) s.",
            )
        end
        return d_full::Design
    end
end

"""
    get_diag_design(d::Design; diagnostics::Bool = false)

Returns a copy of the design `d` where the covariance matrix dictionary generates a diagonal covariance matrix, printing diagnostics if `diagnostics` is `true`.
"""
function get_diag_design(d::Design; diagnostics::Bool = false)
    if d.full_covariance
        # Set some parameters
        start_time = time()
        full_covariance = false
        tuple_number = length(d.tuple_set)
        # Regenerate the covariance dictionary set with a full covariance matrix
        covariance_dict_ensemble =
            Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}(undef, tuple_number)
        covariance_times = Vector{Float64}(undef, tuple_number)
        for i in 1:tuple_number
            circuit_tuple = d.tuple_set[i]
            c_tuple = apply_tuple(d.c, circuit_tuple)
            mapping_set_tuple = d.mapping_ensemble[i]
            experiment_set_tuple = d.experiment_ensemble[i]
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
                    "For tuple $(i), generating the covariance dictionary took $(round(covariance_time, digits = 3)) s. Since starting, $(round(time_5 - start_time, digits = 3)) s have elapsed.",
                )
            end
            # Store the data
            covariance_dict_ensemble[i] = covariance_dict_tuple
            covariance_times[i] = covariance_time
        end
        overall_time = time() - start_time
        # Set the new variables 
        d_diagonal = deepcopy(d)
        @reset d_diagonal.full_covariance = full_covariance
        @reset d_diagonal.covariance_dict_ensemble = covariance_dict_ensemble
        @reset d_diagonal.overall_time =
            d_diagonal.overall_time - sum(d_diagonal.calculation_times[:, 4]) + overall_time
        @reset d_diagonal.calculation_times[:, 4] = covariance_times
        if diagnostics
            println(
                "Generated the full covariance dictionaries for all $(tuple_number) tuples in $(round(overall_time, digits = 3)) s.",
            )
        end
        return d_diagonal::Design
    else
        if diagnostics
            println("The supplied design already has a diagonal covariance matrix.")
        end
        return d::Design
    end
end

"""
    get_combined_design(d::Design; diagnostics::Bool = false)

Returns a copy of the design `d` where the three parameters describing Pauli X, Y, and Z basis measurements have been combined into a single parameter for each qubit, printing diagnostics if `diagnostics` is `true`.
"""
function get_combined_design(d::Design; diagnostics::Bool = false)
    # Get the measurement aggregated circuit
    c_combined = get_combined_circuit(d.c)
    # Recalculate the mappings, design matrix, and covariance dictionaries
    start_time = time()
    tuple_number = length(d.tuple_set)
    N = c_combined.gate_data.N
    design_matrix = convert(SparseMatrixCSC{Int32, Int32}, spzeros(Int32, 0, N))
    mapping_ensemble = Vector{Vector{Mapping}}(undef, tuple_number)
    covariance_dict_ensemble =
        Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}(undef, tuple_number)
    mapping_times = Vector{Float64}(undef, tuple_number)
    covariance_times = Vector{Float64}(undef, tuple_number)
    for i in 1:tuple_number
        # Calculate the mapping data
        circuit_tuple = d.tuple_set[i]
        c_tuple = apply_tuple(c_combined, circuit_tuple)
        time_1 = time()
        (mapping_set, mapping_matrix) = calc_mapping_set(c_tuple)
        time_2 = time()
        covariance_dict = calc_covariance_dict(
            c_tuple,
            mapping_set,
            d.experiment_ensemble[i],
            d.full_covariance,
        )
        time_3 = time()
        # Track the time taken
        mapping_time = time_2 - time_1
        covariance_time = time_3 - time_2
        if diagnostics
            println(
                "For tuple $(i), generating the mappings took $(round(mapping_time, digits = 3)) s and the covariance dictionary took $(round(covariance_time, digits = 3)) s. Since starting, $(round(time_3 - start_time, digits = 3)) s have elapsed.",
            )
        end
        # Store the data
        design_matrix = vcat(design_matrix, mapping_matrix)
        mapping_ensemble[i] = mapping_set
        covariance_dict_ensemble[i] = covariance_dict
        mapping_times[i] = mapping_time
        covariance_times[i] = covariance_time
    end
    overall_time = time() - start_time
    # Set the new variables
    d_combined = deepcopy(d)
    @reset d_combined.c = c_combined
    @reset d_combined.matrix = design_matrix
    @reset d_combined.mapping_ensemble = mapping_ensemble
    @reset d_combined.covariance_dict_ensemble = covariance_dict_ensemble
    @reset d_combined.overall_time =
        d.overall_time - sum(d.calculation_times[:, 1]) - sum(d.calculation_times[:, 4]) +
        overall_time
    @reset d_combined.calculation_times[:, 1] = mapping_times
    @reset d_combined.calculation_times[:, 4] = covariance_times
    if diagnostics
        println(
            "Generated the combined mappings for all $(tuple_number) tuples in $(round(overall_time, digits = 3)) s.",
        )
    end
    return d_combined::Design
end

"""
    update_noise(d::Design, noise_param::AbstractNoiseParameters)
    update_noise(d::Design, gate_probabilities::Dict{Gate, Vector{Float64}})

Returns a copy of the design `d` where the circuit has been updated with noise generated according to `noise_param`, or the gate probabilities `gate_probabilities`.
"""
function update_noise(d::Design, noise_param::T) where {T <: AbstractNoiseParameters}
    # Update the circuit
    c_update = update_noise(d.c, noise_param)
    if c_update.gate_data.combined && ~d.c.gate_data.combined
        d_update = get_combined_design(d)
    else
        d_update = deepcopy(d)
    end
    @reset d_update.c = c_update
    return d_update::Design
end
function update_noise(d::Design, gate_probabilities::Dict{Gate, Vector{Float64}})
    # Update the circuit
    c_update = update_noise(d.c, gate_probabilities)
    d_update = deepcopy(d)
    @reset d_update.c = c_update
    return d_update::Design
end

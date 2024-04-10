#
function TrivialTupleSet(c::AbstractCircuit)
    trivial_tuple_set = Vector{Int}[[[i] for i in c.unique_layer_indices]; [[]]]
    return trivial_tuple_set
end

#
function TrivialExperimentNumbers(c::AbstractCircuit)
    # Determine the experiment numbers for the trivial tuple set
    trivial_max_targets = Int[
        [
            maximum(length(gate.targets) for gate in c.circuit[i].layer) for
            i in c.unique_layer_indices
        ]
        1
    ]
    # Calculate the experiment numbers for the trivial tuple set
    sign_factor = 2
    pauli_types = 3
    trivial_experiment_numbers = (sign_factor * pauli_types) .^ trivial_max_targets
    return trivial_experiment_numbers::Vector{Int}
end

# 
function TrivialHarmMean(code::Code)
    # Initialise times
    meas_reset_time = code.layer_times[end]
    layer_times = code.layer_times[1:(end - 1)]
    # Initialise the trivial tuple set
    trivial_tuple_set = TrivialTupleSet(code)
    # Determine the times for each of the tuples in the trivial tuple set
    trivial_tuple_times = Vector{Float64}(undef, length(trivial_tuple_set))
    for (idx, tuple) in enumerate(trivial_tuple_set)
        tuple_time = meas_reset_time
        for layer_idx in tuple
            tuple_time += layer_times[layer_idx]
        end
        trivial_tuple_times[idx] = tuple_time
    end
    # Determine the experiment numbers for the trivial tuple set
    trivial_experiment_numbers = TrivialExperimentNumbers(code)
    # Normalise the tuple set times by the haarmonic mean of the trivial tuple set experiment times
    trivial_times_harm_mean =
        sum(trivial_experiment_numbers) /
        sum(trivial_experiment_numbers ./ trivial_tuple_times)
    return trivial_times_harm_mean::Float64
end

# 
function TupleSetParameters(
    code::Code,
    tuple_set::Vector{Vector{Int}},
    experiment_numbers::Vector{Int},
)
    # Initialise times
    meas_reset_time = code.layer_times[end]
    layer_times = code.layer_times[1:(end - 1)]
    # Determine the times for each of the tuples in the tuple set
    unnormalised_tuple_times = Vector{Float64}(undef, length(tuple_set))
    for (idx, tuple) in enumerate(tuple_set)
        tuple_time = meas_reset_time
        for layer_idx in tuple
            tuple_time += layer_times[layer_idx]
        end
        unnormalised_tuple_times[idx] = tuple_time
    end
    # Normalise the tuple set times by the haarmonic mean of the trivial tuple set experiment times
    trivial_times_harm_mean = TrivialHarmMean(code)
    tuple_times = unnormalised_tuple_times / trivial_times_harm_mean
    # Determine the shot weights for the tuple set
    shot_weights =
        (experiment_numbers ./ tuple_times) ./ sum(experiment_numbers ./ tuple_times)
    @assert sum(shot_weights) ≈ 1.0
    return (tuple_times::Vector{Float64}, shot_weights::Vector{Float64})
end

#
function TupleSetData(code::Code)
    # Generate the repeated tuple set data
    if hasproperty(code.code_param, :dynamically_decouple) &&
       code.code_param.dynamically_decouple
        # Repeated tuples for dynamically decoupled circuits
        # Initialise relevant parameters
        types = unique(code.layer_types)
        two_qubit_type = :two_qubit
        dynamical_type = :dynamical
        type_indices = [
            intersect(findall(code.layer_types .== types[idx]), code.unique_layer_indices) for idx in eachindex(types)
        ]
        @assert two_qubit_type ∈ types "The code must have a two-qubit gate layer."
        @assert dynamical_type ∈ types "The code must have a dynamical decoupling gate layer."
        two_qubit_indices = vcat(type_indices[types .== two_qubit_type]...)
        dynamical_indices = vcat(type_indices[types .== dynamical_type]...)
        non_two_qubit_indices = vcat(type_indices[types .!= two_qubit_type]...)
        # This interleaves the two-qubit layers with dynamical decoupling layers but leaves the other layers unchanged
        # The number of repetitions 4 was chosen to ensure that the repeated tuple circuits implement the identity
        # However, this is only necessarily true for the original rotated planar circuit
        non_two_qubit_tuple_set = [[i, i] for i in non_two_qubit_indices]
        two_qubit_tuple_set = reshape(
            [repeat([i, j], 4) for i in two_qubit_indices, j in dynamical_indices],
            (length(two_qubit_indices) * length(dynamical_indices)),
        )
        repeat_tuple_set = vcat(non_two_qubit_tuple_set, two_qubit_tuple_set)
        type_num = length(types)
        repeat_numbers = zeros(Int, type_num)
        non_two_qubit_repeat_indices =
            [findfirst(code.layer_types[i] .== types) for i in non_two_qubit_indices]
        two_qubit_repeat_indices = repeat(
            [findfirst(code.layer_types[i] .== types) for i in two_qubit_indices],
            length(dynamical_indices),
        )
        repeat_indices = vcat(non_two_qubit_repeat_indices, two_qubit_repeat_indices)
    else
        # Repeated tuples for other circuits
        # The repetition of the layers ensures that the repeated tuple circuits implement the identity
        repeat_tuple_set = [[i; i] for i in code.unique_layer_indices]
        types = unique(code.layer_types)
        type_num = length(types)
        repeat_numbers = zeros(Int, type_num)
        repeat_indices =
            [findfirst(code.layer_types[i] .== types) for i in code.unique_layer_indices]
    end
    # Generate the tuple set data
    trivial_tuple_set = TrivialTupleSet(code)
    tuple_set_data =
        TupleSetData(trivial_tuple_set, repeat_tuple_set, repeat_numbers, repeat_indices)
    return tuple_set_data
end

#
function TupleSet(tuple_set_data::TupleSetData)
    # Initialise data
    repeat_tuple_set = tuple_set_data.repeat_tuple_set
    repeat_numbers = tuple_set_data.repeat_numbers
    repeat_indices = tuple_set_data.repeat_indices
    # Generate the tuple set
    repeated_tuple_set = [
        repeat(repeat_tuple_set[idx], repeat_numbers[repeat_indices[idx]]) for
        idx in 1:length(repeat_tuple_set) if repeat_numbers[repeat_indices[idx]] > 0
    ]
    tuple_set = unique([repeated_tuple_set; tuple_set_data.tuple_set])
    return tuple_set::Vector{Vector{Int}}
end

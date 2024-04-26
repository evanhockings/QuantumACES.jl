"""
    TupleSetData

Data parameterising a tuple set.

# Fields

  - `tuple_set::Vector{Vector{Int}}`: The main tuple set, whose tuples are not repeated.
  - `repeat_tuple_set::Vector{Vector{Int}}`: The repeated tuple set, whose tuples are repeated `repeat_numbers` times.
  - `repeat_numbers::Vector{Int}`: The number of repetitions for each tuple in the repeated tuple set `repeat_tuple_set`.
  - `repeat_indices::Vector{Int}`: Indexes the repetition number in `repeat_numbers` that describes the number of repetitions for each tuple in `repeat_tuple_set`.
"""
struct TupleSetData
    tuple_set::Vector{Vector{Int}}
    repeat_tuple_set::Vector{Vector{Int}}
    repeat_numbers::Vector{Int}
    repeat_indices::Vector{Int}
    # Constructor
    function TupleSetData(
        tuple_set::Vector{Vector{Int}},
        repeat_tuple_set::Vector{Vector{Int}},
        repeat_numbers::Vector{Int},
        repeat_indices::Vector{Int},
    )
        @assert all(repeat_numbers .>= 0) "The repeat numbers must be non-negative."
        return new(
            tuple_set,
            repeat_tuple_set,
            repeat_numbers,
            repeat_indices,
        )::TupleSetData
    end
end

@struct_hash_equal_isequal TupleSetData

"""
    get_basic_tuple_set(c::AbstractCircuit)

Returns the basic tuple set for the circuit `c`.
"""
function get_basic_tuple_set(c::T) where {T <: AbstractCircuit}
    basic_tuple_set = Vector{Int}[[[i] for i in c.unique_layer_indices]; [[]]]
    return basic_tuple_set::Vector{Vector{Int}}
end

"""
    get_basic_experiment_numbers(c::AbstractCircuit)

Returns the experiment numbers corresponding to the basic tuple set for the circuit `c`.
"""
function get_basic_experiment_numbers(c::T) where {T <: AbstractCircuit}
    # Determine the experiment numbers for the basic tuple set
    sign_factor = 2
    pauli_types = 3
    basic_max_targets = Int[
        [
            maximum(length(gate.targets) for gate in c.circuit[i].layer) for
            i in c.unique_layer_indices
        ]
        1
    ]
    if hasproperty(c, :partition)
        basic_experiment_numbers = (sign_factor * pauli_types) .^ basic_max_targets
    else
        basic_experiment_numbers = (pauli_types) .^ basic_max_targets
    end
    return basic_experiment_numbers::Vector{Int}
end

"""
    get_basic_times_harm_mean(c::AbstractCircuit)

Returns the harmonic mean of the experiment times corresponding to the basic tuple set for the circuit `c`.
"""
function get_basic_times_harm_mean(c::T) where {T <: AbstractCircuit}
    # Initialise times
    meas_reset_time = c.layer_times[end]
    layer_times = c.layer_times[1:(end - 1)]
    # Initialise the basic tuple set
    basic_tuple_set = get_basic_tuple_set(c)
    # Determine the times for each of the tuples in the basic tuple set
    basic_tuple_times = Vector{Float64}(undef, length(basic_tuple_set))
    for (idx, tuple) in enumerate(basic_tuple_set)
        tuple_time = meas_reset_time
        for layer_idx in tuple
            tuple_time += layer_times[layer_idx]
        end
        basic_tuple_times[idx] = tuple_time
    end
    # Determine the experiment numbers for the basic tuple set
    basic_experiment_numbers = get_basic_experiment_numbers(c)
    # Normalise the tuple set times by the haarmonic mean of the basic tuple set experiment times
    basic_times_harm_mean =
        sum(basic_experiment_numbers) / sum(basic_experiment_numbers ./ basic_tuple_times)
    return basic_times_harm_mean::Float64
end

"""
    get_tuple_set_params(c::AbstractCircuit, tuple_set::Vector{Vector{Int}}, experiment_numbers::Vector{Int})

Returns the time taken to implement each tuple `tuple_times` and the default shot weights for each tuple `shot_weights` corresponding to the circuit `c` with the tuple set `tuple_set` and the experiment numbers `experiment_numbers`.
"""
function get_tuple_set_params(
    c::T,
    tuple_set::Vector{Vector{Int}},
    experiment_numbers::Vector{Int},
) where {T <: AbstractCircuit}
    # Initialise times
    meas_reset_time = c.layer_times[end]
    layer_times = c.layer_times[1:(end - 1)]
    # Determine the times for each of the tuples in the tuple set
    unnormalised_tuple_times = Vector{Float64}(undef, length(tuple_set))
    for (idx, tuple) in enumerate(tuple_set)
        tuple_time = meas_reset_time
        for layer_idx in tuple
            tuple_time += layer_times[layer_idx]
        end
        unnormalised_tuple_times[idx] = tuple_time
    end
    # Normalise the tuple set times by the haarmonic mean of the basic tuple set experiment times
    basic_times_harm_mean = get_basic_times_harm_mean(c)
    tuple_times = unnormalised_tuple_times / basic_times_harm_mean
    # Determine the shot weights for the tuple set
    shot_weights =
        (experiment_numbers ./ tuple_times) ./ sum(experiment_numbers ./ tuple_times)
    @assert sum(shot_weights) ≈ 1.0
    return (tuple_times::Vector{Float64}, shot_weights::Vector{Float64})
end

"""
    get_tuple_set_data(c::AbstractCircuit; init_scaling::Float64 = 0.2)
    get_tuple_set_data(c::AbstractCircuit, tuple_set::Vector{Vector{Int}}; init_scaling::Float64 = 0.2)

Returns the tuple set data corresponding to the circuit `c`, with the non-repeated tuples either being the supplied `tuple_set` or the basic tuple set for `c`.
The repeat numbers are initialised to be inversely proportional to the average noise on the gates in the layers, implicitly assuming depolarising noise, scaled by a factor `init_scaling` which is empirically helpful.
"""
function get_tuple_set_data(
    c::T,
    tuple_set::Vector{Vector{Int}};
    init_scaling::Float64 = 0.2,
) where {T <: AbstractCircuit}
    # Initialise parameters
    r_1 = c.noise_param.r_1
    r_2 = c.noise_param.r_2
    types = unique(c.layer_types)
    type_num = length(types)
    repeat_numbers = zeros(type_num)
    # Generate the repeated tuple set data
    if hasproperty(c.circuit_param, :dynamically_decouple) &&
       c.circuit_param.dynamically_decouple
        # Repeated tuples for dynamically decoupled circuits
        # Initialise relevant parameters
        type_indices = [
            intersect(findall(c.layer_types .== types[idx]), c.unique_layer_indices) for
            idx in eachindex(types)
        ]
        @assert :two_qubit ∈ types "The circuit must have a two-qubit gate layer."
        @assert :dynamical ∈ types "The circuit must have a dynamical decoupling gate layer."
        two_qubit_indices = vcat(type_indices[types .== :two_qubit]...)
        dynamical_indices = vcat(type_indices[types .== :dynamical]...)
        non_two_qubit_indices = vcat(type_indices[types .!= :two_qubit]...)
        # This interleaves the two-qubit layers with dynamical decoupling layers but leaves the other layers unchanged
        non_two_qubit_tuple_set = [[i] for i in non_two_qubit_indices]
        two_qubit_tuple_set = reshape(
            [[i, j, i, j] for i in two_qubit_indices, j in dynamical_indices],
            (length(two_qubit_indices) * length(dynamical_indices)),
        )
        repeat_tuple_set = vcat(non_two_qubit_tuple_set, two_qubit_tuple_set)
        non_two_qubit_repeat_indices =
            [findfirst(c.layer_types[i] .== types) for i in non_two_qubit_indices]
        two_qubit_repeat_indices = repeat(
            [findfirst(c.layer_types[i] .== types) for i in two_qubit_indices],
            length(dynamical_indices),
        )
        repeat_indices = vcat(non_two_qubit_repeat_indices, two_qubit_repeat_indices)
        # Initialise the repeat numbers
        for (idx, type) in enumerate(types)
            if type == :single_qubit || type == :dynamical
                repeat_numbers[idx] = init_scaling * (1 / ((4 / 3) * r_1))
            elseif type == :two_qubit
                repeat_numbers[idx] = init_scaling * (1 / ((4 / 3) * r_1 + (16 / 15) * r_2))
            else
                throw(error("Unsupported gate type $(type)."))
            end
        end
    else
        # Repeated tuples for other circuits
        repeat_tuple_set = [[i] for i in c.unique_layer_indices]
        repeat_indices =
            [findfirst(c.layer_types[i] .== types) for i in c.unique_layer_indices]
        # Initialise the repeat numbers
        for (idx, type) in enumerate(types)
            if type == :single_qubit
                repeat_numbers[idx] = init_scaling * (1 / ((4 / 3) * r_1))
            elseif type == :two_qubit
                repeat_numbers[idx] = init_scaling * (1 / ((16 / 15) * r_2))
            else
                throw(error("Unsupported gate type $(type)."))
            end
        end
    end
    # The repeat number initialisation we use here works well empirically
    # It seems best to have repeat tuples of order 2 repeated an odd number of times
    repeat_numbers = 2 * round.(Int, repeat_numbers / 2) .- 1
    # Generate the tuple set data
    tuple_set_data =
        TupleSetData(tuple_set, repeat_tuple_set, repeat_numbers, repeat_indices)
    return tuple_set_data::TupleSetData
end
function get_tuple_set_data(c::T; init_scaling::Float64 = 0.2) where {T <: AbstractCircuit}
    # Generate the tuple set data
    basic_tuple_set = get_basic_tuple_set(c)
    tuple_set_data = get_tuple_set_data(c, basic_tuple_set; init_scaling = init_scaling)
    return tuple_set_data::TupleSetData
end

"""
    get_tuple_set(tuple_set_data::TupleSetData)

Returns the tuple set corresponding to the data `tuple_set_data`.
"""
function get_tuple_set(tuple_set_data::TupleSetData)
    # Initialise data
    repeat_tuple_set = tuple_set_data.repeat_tuple_set
    repeat_numbers = tuple_set_data.repeat_numbers
    repeat_indices = tuple_set_data.repeat_indices
    @assert length(repeat_tuple_set) == length(repeat_indices) "The number of repeat indices does not match the number of repeated tuples."
    # Generate the tuple set
    repeated_tuple_set = [
        repeat(repeat_tuple_set[idx], repeat_numbers[repeat_indices[idx]]) for
        idx in eachindex(repeat_tuple_set) if repeat_numbers[repeat_indices[idx]] > 0
    ]
    tuple_set = unique([repeated_tuple_set; tuple_set_data.tuple_set])
    return tuple_set::Vector{Vector{Int}}
end

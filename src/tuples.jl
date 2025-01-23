"""
    TupleSetData

Data parameterising a tuple set.

# Fields

  - `tuple_set::Vector{Vector{Int}}`: The main tuple set, whose tuples are not repeated.
  - `repeat_tuple_set::Vector{Vector{Int}}`: The repeat tuple set, whose tuples are repeated according to the `repeat_numbers`.
  - `repeat_numbers::Vector{Int}`: The number of repetitions for tuples in the repeat tuple set `repeat_tuple_set`.
  - `repeat_indices::Vector{Tuple{Int, Int}}`: Indexes pairs (`repeat_tuple`, `repeat_number`) from (`repeat_tuple_set`, `repeat_numbers`) such that the repeated tuple set consists of `repeat_tuple` repeated `repeat_number` times for all indexed pairs.
"""
struct TupleSetData
    tuple_set::Vector{Vector{Int}}
    repeat_tuple_set::Vector{Vector{Int}}
    repeat_numbers::Vector{Int}
    repeat_indices::Vector{Tuple{Int, Int}}
    # Constructor
    function TupleSetData(
        tuple_set::Vector{Vector{Int}},
        repeat_tuple_set::Vector{Vector{Int}},
        repeat_numbers::Vector{Int},
        repeat_indices::Vector{Tuple{Int, Int}},
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
    basic_tuple_set = Vector{Int}[[[]]; [[i] for i in c.unique_layer_indices]]
    return basic_tuple_set::Vector{Vector{Int}}
end

"""
    get_basic_experiment_numbers(c::AbstractCircuit)

Returns the experiment numbers corresponding to the basic tuple set for the circuit `c`.

BEWARE: This currently assumes that circuit eigenvalues for Paulis supported on ``n`` qubits require preparing exactly and only ``2^n`` sign configurations.
Errors will occur if this is not the case, and they are not guaranteed to be noisy.
If creating a circuit where this is not the case, you will need to provide new methods for this function.
"""
function get_basic_experiment_numbers(c::T) where {T <: AbstractCircuit}
    # Determine the experiment numbers for the basic tuple set
    pauli_types = 3
    basic_max_targets = Int[
        1
        [
            maximum(length(gate.targets) for gate in c.circuit[i].layer) for
            i in c.unique_layer_indices
        ]
    ]
    basic_experiment_numbers = pauli_types .^ basic_max_targets
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
    for (idx, tuple) in pairs(basic_tuple_set)
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
    get_tuple_set_params(c::AbstractCircuit, tuple_set::Vector{Vector{Int}}, experiment_numbers::Vector{Int}; weight_experiments::Bool = true)

Returns the time taken to implement each tuple and the default shot weights for each tuple orresponding to the circuit `c` with the tuple set `tuple_set` and the experiment numbers `experiment_numbers`.
If `weight_experiments` is `true`, allocate the same sampling time to each experiment for each tuple; otherwise, allocate the same sampling time to each tuple overall.
The default times allocate the same sampling time to each experiment for each tuple; this seems to be more performant than allocating the same sampling time to each tuple overall, as it allocates more samples to tuples that require more experiments.
"""
function get_tuple_set_params(
    c::T,
    tuple_set::Vector{Vector{Int}},
    experiment_numbers::Vector{Int};
    weight_experiments::Bool = true,
) where {T <: AbstractCircuit}
    # Initialise times
    meas_reset_time = c.layer_times[end]
    layer_times = c.layer_times[1:(end - 1)]
    # Determine the times for each of the tuples in the tuple set
    unnormalised_tuple_times = Vector{Float64}(undef, length(tuple_set))
    for (idx, tuple) in pairs(tuple_set)
        tuple_time = meas_reset_time
        for layer_idx in tuple
            tuple_time += layer_times[layer_idx]
        end
        unnormalised_tuple_times[idx] = tuple_time
    end
    # Normalise the tuple set times by the harmonic mean of the basic tuple set experiment times
    basic_times_harm_mean = get_basic_times_harm_mean(c)
    tuple_times = unnormalised_tuple_times / basic_times_harm_mean
    # Determine the shot weights for the tuple set
    if weight_experiments
        shot_weights =
            (experiment_numbers ./ tuple_times) / sum(experiment_numbers ./ tuple_times)
    else
        shot_weights = (1 ./ tuple_times) / sum(1 ./ tuple_times)
    end
    @assert sum(shot_weights) ≈ 1.0
    return (tuple_times::Vector{Float64}, shot_weights::Vector{Float64})
end

"""
    check_tuple!(c::AbstractCircuit, circuit_tuple::Vector{Int})
    check_tuple!(c::AbstractCircuit, tuple_set::Vector{Vector{Int}})

Checks that the tuple `circuit_tuple`, or each tuple in the tuple set `tuple_set`, is valid for the circuit `c`.
"""
function check_tuple!(c::T, circuit_tuple::Vector{Int}) where {T <: AbstractCircuit}
    # Check the tuple is valid
    tuple_layer_types = c.layer_types[circuit_tuple]
    if :mid_reset ∈ tuple_layer_types && :two_qubit ∈ tuple_layer_types
        throw(
            error(
                "The tuple $(circuit_tuple) contains a mid-circuit reset layer and a two-qubit gate layer, which is not currently supported.",
            ),
        )
    end
    return nothing
end
function check_tuple!(c::T, tuple_set::Vector{Vector{Int}}) where {T <: AbstractCircuit}
    # Check that the tuple set is unique
    @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
    # Check each tuple in the set is valid
    for circuit_tuple in tuple_set
        check_tuple!(c, circuit_tuple)
    end
    return nothing
end

"""
    get_circuit_tuple(c::AbstractCircuit; include_reset::Bool = false)

Returns the circuit tuple corresponding to the circuit `c`, including reset layers if `include_reset` is `true`.
"""
function get_circuit_tuple(c::T; include_reset::Bool = false) where {T <: AbstractCircuit}
    # Generate the circuit tuple
    if include_reset
        ignored_types = []
    else
        ignored_types = [:mid_reset]
    end
    circuit_tuple = [
        c.unique_layer_indices[findfirst(
            c.circuit[idx] == c.circuit[u_idx] for u_idx in c.unique_layer_indices
        )] for idx in c.circuit_tuple if c.layer_types[idx] ∉ ignored_types
    ]
    return circuit_tuple::Vector{Int}
end

"""
    get_tuple_set_data(c::AbstractCircuit; error_target::Real = 0.1, add_circuit::Bool = false)
    get_tuple_set_data(c::AbstractCircuit, tuple_set::Vector{Vector{Int}}; error_target::Real = 0.1, add_circuit::Bool = false)

Returns the tuple set data corresponding to the circuit `c`, with the non-repeated tuples either being the supplied `tuple_set` or the basic tuple set for `c`.
The repeat numbers are initialised to be inversely proportional to the average noise on the gates in the layers, which heuristically roughly target an error rate `error_target`, and the original circuit is added to the repeat tuples if `add_circuit` is `true`.

Note that it appears best to use repeat tuples of order 2 (up to Paulis), with odd repeat numbers where this is efficient for measurement (such as repeating the same circuit layer many times), and even repeat numbers otherwise (such as syndrome extraction circuits).
"""
function get_tuple_set_data(
    c::T,
    tuple_set::Vector{Vector{Int}};
    error_target::Real = 0.1,
    add_circuit::Bool = false,
) where {T <: AbstractCircuit}
    # Initialise parameters
    check_tuple!(c, tuple_set)
    noise_params = c.noise_param.params
    expected_noise = (
        haskey(noise_params, :r_1) &&
        haskey(noise_params, :r_2) &&
        haskey(noise_params, :r_m)
    )
    if expected_noise
        r_1 = noise_params[:r_1]
        r_2 = noise_params[:r_2]
        r_m = noise_params[:r_m]
    else
        println(
            "WARNING: Unexpected noise parameters; the repeat numbers will be initialised without noise information.",
        )
        r_1 = 0.1 / 100
        r_2 = 0.5 / 100
        r_m = 1.0 / 100
    end
    layer_types = c.layer_types
    unique_layer_indices = c.unique_layer_indices
    supported_types = [:single_qubit; :two_qubit; :dynamical; :mid_reset]
    ignored_types = [:mid_reset]
    @assert layer_types ⊆ supported_types "The circuit contains unsupported layer types."
    types = unique(setdiff(layer_types, ignored_types))
    type_num = length(types)
    # Initialise the repeat numbers
    circuit_params = c.circuit_param.params
    is_decoupled = (
        haskey(circuit_params, :dynamically_decouple) &&
        circuit_params[:dynamically_decouple]
    )
    single_qubit_repeat = ceil(Int, (error_target - r_m) / r_1)
    if is_decoupled
        two_qubit_repeat = ceil(Int, (3 * error_target - 2 * r_m) / (r_2 + r_1))
    else
        two_qubit_repeat = ceil(Int, (3 * error_target - 2 * r_m) / r_2)
    end
    dynamical_repeat = ceil(Int, (0.8 * error_target - r_m) / r_1)
    circuit_tuple = get_circuit_tuple(c)
    circuit_error =
        r_1 * sum(layer_types[circuit_tuple] .== :single_qubit) +
        r_2 * sum(layer_types[circuit_tuple] .== :two_qubit) +
        r_1 * sum(layer_types[circuit_tuple] .== :dynamical)
    circuit_repeat = ceil(Int, (3 * error_target - 2 * r_m) / circuit_error)
    repeat_dict = Dict(
        :single_qubit => single_qubit_repeat,
        :two_qubit => two_qubit_repeat,
        :dynamical => dynamical_repeat,
        :circuit => circuit_repeat,
    )
    # Ordinary repeat numbers should be odd
    repeat_numbers = zeros(Int, type_num)
    for idx in eachindex(types)
        repeat_number = repeat_dict[types[idx]]
        if repeat_number >= 3
            if iseven(repeat_number)
                repeat_number -= 1
            end
        else
            repeat_number = 3
        end
        repeat_numbers[idx] = repeat_number
    end
    # Generate the repeat tuple set data
    if is_decoupled
        # Repeated tuples for dynamically decoupled circuits
        # Initialise relevant parameters
        @assert :two_qubit ∈ types "The circuit must have a two-qubit gate layer."
        @assert :dynamical ∈ types "The circuit must have a dynamical decoupling gate layer."
        type_indices = [
            intersect(findall(layer_types .== types[idx]), unique_layer_indices) for
            idx in eachindex(types)
        ]
        non_two_qubit_indices = vcat(type_indices[types .!= :two_qubit]...)
        two_qubit_indices = vcat(type_indices[types .== :two_qubit]...)
        dynamical_indices = vcat(type_indices[types .== :dynamical]...)
        pair_two_qubit_indices = reshape(
            [(i, j) for i in two_qubit_indices, j in dynamical_indices],
            length(two_qubit_indices) * length(dynamical_indices),
        )
        # This interleaves the two-qubit layers with dynamical decoupling layers but leaves the other layers unchanged
        non_two_qubit_tuple_set = [[i] for i in non_two_qubit_indices]
        two_qubit_tuple_set = [[i, j] for (i, j) in pair_two_qubit_indices]
        repeat_tuple_set = vcat(non_two_qubit_tuple_set, two_qubit_tuple_set)
        non_two_qubit_repeat_indices = [
            (idx, findfirst(layer_types[i] .== types)) for
            (idx, i) in pairs(non_two_qubit_indices)
        ]
        two_qubit_repeat_indices = [
            (idx + length(non_two_qubit_indices), findfirst(layer_types[i] .== types))
            for (idx, (i, j)) in pairs(pair_two_qubit_indices)
        ]
        repeat_indices = vcat(non_two_qubit_repeat_indices, two_qubit_repeat_indices)
    else
        # Repeated tuples for other circuits
        type_indices = [i for i in unique_layer_indices if layer_types[i] ∈ types]
        repeat_tuple_set = [[i] for i in type_indices]
        repeat_indices =
            [(idx, findfirst(layer_types[i] .== types)) for (idx, i) in pairs(type_indices)]
    end
    # Append the circuit tuple if appropriate
    if add_circuit
        # Add the circuit type
        push!(types, :circuit)
        type_num = length(types)
        # Add the circuit repeat number, ensuring it is even
        repeat_number = repeat_dict[:circuit]
        if repeat_number >= 2
            if isodd(repeat_number)
                repeat_number -= 1
            end
        else
            repeat_number = 2
        end
        repeat_numbers = vcat(repeat_numbers, repeat_number)
        # Add the circuit tuple
        repeat_tuple_set = vcat(repeat_tuple_set, [circuit_tuple])
        repeat_tuple_num = length(repeat_tuple_set)
        # Add the circuit tuple repeat index
        repeat_indices = vcat(repeat_indices, [(repeat_tuple_num, type_num)])
        @assert repeat_tuple_set[repeat_tuple_num] == circuit_tuple "The circuit tuple must be the last tuple."
        @assert types[type_num] == :circuit "The circuit type must be the last type."
    end
    # Generate the tuple set data
    tuple_set_data = TupleSetData(
        sort(tuple_set; by = x -> length(x)),
        repeat_tuple_set,
        repeat_numbers,
        repeat_indices,
    )
    return tuple_set_data::TupleSetData
end
function get_tuple_set_data(
    c::T;
    error_target::Real = 0.1,
    add_circuit::Bool = false,
) where {T <: AbstractCircuit}
    # Generate the tuple set data
    basic_tuple_set = get_basic_tuple_set(c)
    tuple_set_data = get_tuple_set_data(
        c,
        basic_tuple_set;
        error_target = error_target,
        add_circuit = add_circuit,
    )
    return tuple_set_data::TupleSetData
end

"""
    get_augmented_tuple_set_data(tuple_set_data::TupleSetData, repeat_points::Integer; initial_shrink_factor::Real = 2^(repeat_points - 1))

Returns an augmented version of the tuple set data `tuple_set_data`, where each of the repetition numbers is augmented to have `repeat_points` total repetition numbers, spaced logarithmically between the repetition number and the repetition number shrunk by a factor of `initial_shrink_factor`.
"""
function get_augmented_tuple_set_data(
    tuple_set_data::TupleSetData,
    repeat_points::Integer;
    initial_shrink_factor::Real = 2^(repeat_points - 1),
)
    # Check parameters
    repeat_tuple_set = tuple_set_data.repeat_tuple_set
    repeat_numbers = tuple_set_data.repeat_numbers
    repeat_indices = tuple_set_data.repeat_indices
    @assert all(
        [tuple_idx for (tuple_idx, repeat_idx) in repeat_indices] == collect(1:length(repeat_tuple_set)),
    ) "The supplied tuple set data already appears to be augmented."
    @assert repeat_points > 1 "Must augment each repeat number with at least one additional point."
    # Generate the augmented repeat numbers and indices
    augmented_repeat_numbers_set = Vector{Vector{Int}}(undef, length(repeat_numbers))
    for (repeat_idx, repeat_number) in pairs(repeat_numbers)
        spaced_repeat_numbers =
            ceil.(
                Int,
                exp.(
                    range(
                        log(repeat_number / initial_shrink_factor),
                        log(repeat_number),
                        repeat_points,
                    )
                ),
            )
        # Maintain repeat number parity
        if isodd(repeat_number)
            for idx in eachindex(spaced_repeat_numbers)
                if iseven(spaced_repeat_numbers[idx])
                    spaced_repeat_numbers[idx] -= 1
                end
            end
            spaced_repeat_numbers =
                unique(spaced_repeat_numbers[spaced_repeat_numbers .>= 3])
        else
            for idx in eachindex(spaced_repeat_numbers)
                if isodd(spaced_repeat_numbers[idx])
                    spaced_repeat_numbers[idx] -= 1
                end
            end
            spaced_repeat_numbers =
                unique(spaced_repeat_numbers[spaced_repeat_numbers .>= 2])
        end
        @assert spaced_repeat_numbers[end] == repeat_number
        augmented_repeat_numbers_set[repeat_idx] = spaced_repeat_numbers
    end
    augmented_repeat_numbers = vcat(augmented_repeat_numbers_set...)
    # Generate the augmented repeat indices
    augmented_repeat_numbers_lengths = length.(augmented_repeat_numbers_set)
    augmented_repeat_numbers_lower =
        cumsum([0; augmented_repeat_numbers_lengths[1:(end - 1)]])
    augmented_repeat_indices = Vector{Tuple{Int, Int}}()
    for (tuple_idx, repeat_idx) in repeat_indices
        for point_idx in 1:augmented_repeat_numbers_lengths[repeat_idx]
            push!(
                augmented_repeat_indices,
                (tuple_idx, (augmented_repeat_numbers_lower[repeat_idx] + point_idx)),
            )
        end
    end
    # Return the augmented tuple set data
    augmented_tuple_set_data = TupleSetData(
        sort(tuple_set_data.tuple_set; by = x -> length(x)),
        repeat_tuple_set,
        augmented_repeat_numbers,
        augmented_repeat_indices,
    )
    return augmented_tuple_set_data::TupleSetData
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
    # Generate the tuple set
    repeated_tuple_set = [
        repeat(repeat_tuple_set[tuple_idx], repeat_numbers[repeat_idx]) for
        (tuple_idx, repeat_idx) in repeat_indices if repeat_numbers[repeat_idx] > 0
    ]
    tuple_set = unique([
        tuple_set_data.tuple_set
        repeated_tuple_set
    ])
    return tuple_set::Vector{Vector{Int}}
end

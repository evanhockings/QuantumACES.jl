"""
    optimal_expectation(tuple_set_data::TupleSetData, expectation_dict::Dict{Vector{Int}, Float64}, c::AbstractCircuit; options::OptimOptions = OptimOptions())

Returns the optimised figure of merit, and a dictionary of stored values, for the circuit `c` with tuple set data `tuple_set_data`, with optimised shot weights.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function optimal_expectation(
    tuple_set_data::TupleSetData,
    expectation_dict::Dict{Vector{Int}, Float64},
    c::T;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Get the keyword arguments
    ls_type = options.ls_type
    est_type = options.rep_est_type
    est_weight = options.rep_est_weight
    # Retrieve the NRMSE expectation if it has already been calculated, else calculate it
    repeat_numbers = tuple_set_data.repeat_numbers
    if haskey(expectation_dict, repeat_numbers)
        expectation = expectation_dict[repeat_numbers]
    elseif any(repeat_numbers .< 0)
        expectation = 1e20
        expectation_dict[repeat_numbers] = expectation
    else
        # Generate the tuple set and design
        d = generate_design(c, tuple_set_data; full_covariance = true)
        covariance_log = calc_covariance_log(d)
        # Optimise the shot weights
        # Treat `rep_est_type` and `rep_est_weight` as `est_type` and `est_weight` here
        options_rep = deepcopy(options)
        @reset options_rep.est_type = est_type
        @reset options_rep.est_weight = est_weight
        (d, covariance_log) =
            optimise_weights(d, covariance_log; options = options_rep)[1:2]
        # Calculate the NRMSE expectation
        expectation = calc_ls_moments(
            d,
            covariance_log;
            ls_type = ls_type,
            est_type = est_type,
            est_weight = est_weight,
        )[1]
        expectation_dict[repeat_numbers] = expectation
    end
    return (expectation::Float64, expectation_dict::Dict{Vector{Int}, Float64})
end

"""
    cap_repetition(repeat_number::Int, repeat_tuple_depth::Int, min_depth::Int, max_depth::Int; diagnostics::Bool = false)

Returns a copy of the repeat number `repeat_number` capped according to the minimum and maximum depths `min_depth` and `max_depth`, respectively, given a repeat tuple depth `repeat_tuple_depth`, with the same parity as the original repeat number.
"""
function cap_repetition(
    repeat_number::Int,
    repeat_tuple_depth::Int,
    min_depth::Int,
    max_depth::Int;
    diagnostics::Bool = false,
)
    # Cap the repetition number
    cap_repeat_number = deepcopy(repeat_number)
    min_repeat_number = floor(Int, min_depth / repeat_tuple_depth)
    max_repeat_number = floor(Int, max_depth / repeat_tuple_depth)
    if isodd(repeat_number)
        if cap_repeat_number < 3
            cap_repeat_number = 3
        end
        if cap_repeat_number < min_repeat_number
            cap_repeat_number = min_repeat_number
            if iseven(cap_repeat_number)
                cap_repeat_number += 1
            end
            if diagnostics
                println(
                    "The repetition number $(repeat_number) has been capped according to the minimum depth at $(cap_repeat_number).",
                )
            end
        elseif cap_repeat_number > max_repeat_number
            cap_repeat_number = max_repeat_number
            if iseven(cap_repeat_number)
                cap_repeat_number -= 1
            end
            if diagnostics
                println(
                    "The repetition number $(repeat_number) has been capped according to the maximum depth at $(cap_repeat_number).",
                )
            end
        end
    else
        if cap_repeat_number < 2
            cap_repeat_number = 2
        end
        if cap_repeat_number < min_repeat_number
            cap_repeat_number = min_repeat_number
            if isodd(cap_repeat_number)
                cap_repeat_number += 1
            end
            if diagnostics
                println(
                    "The repetition number $(repeat_number) has been capped according to the minimum depth at $(cap_repeat_number).",
                )
            end
        elseif cap_repeat_number > max_repeat_number
            cap_repeat_number = max_repeat_number
            if isodd(cap_repeat_number)
                cap_repeat_number -= 1
            end
            if diagnostics
                println(
                    "The repetition number $(repeat_number) has been capped according to the maximum depth at $(cap_repeat_number).",
                )
            end
        end
    end
    @assert (iseven(repeat_number) && iseven(cap_repeat_number)) ||
            (isodd(repeat_number) && isodd(cap_repeat_number)) "The capped repeat number must be even if the original repeat number is even, and odd if the original repeat number is odd."
    return cap_repeat_number::Int
end

"""
    step_repetition(tuple_set_data::TupleSetData, expectation_dict::Dict{Vector{Int}, Float64}, step_tracker::Vector{Int}, coordinate_idx::Integer, c::AbstractCircuit; options::OptimOptions = OptimOptions())

Returns the tuple set data, a dictionary of stored figure of merit, and the step tracker after stepping the repetition number for the tuple set data `tuple_set_data` at the index `coordinate_idx` for the circuit `c`.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function step_repetition(
    tuple_set_data::TupleSetData,
    expectation_dict::Dict{Vector{Int}, Float64},
    step_tracker::Vector{Int},
    coordinate_idx::Integer,
    c::T;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Get the keyword arguments
    min_depth = options.min_depth
    max_depth = options.max_depth
    diagnostics = options.rep_diagnostics
    # Calculate the figure of merit for the current repetition numbers
    (expectation, expectation_dict) =
        optimal_expectation(tuple_set_data, expectation_dict, c; options = options)
    # Calculate the figure of merit adding one to the coordinate's repetition number
    tuple_set_data_upper = deepcopy(tuple_set_data)
    tuple_set_data_upper.repeat_numbers[coordinate_idx] += 2
    (upper_expectation, expectation_dict) =
        optimal_expectation(tuple_set_data_upper, expectation_dict, c; options = options)
    # Determine the direction in which to step, if at all
    old_repeat_number = deepcopy(tuple_set_data.repeat_numbers[coordinate_idx])
    if old_repeat_number == 0
        if upper_expectation < expectation
            # Increasing the repetition number decreases the figure of merit
            step_sign = 1
        else
            # The repetition number is a local minimum in the figure of merit
            step_sign = 0
        end
    else
        # Calculate the figure of merit subtracting one from the coordinate's repetition number
        tuple_set_data_lower = deepcopy(tuple_set_data)
        tuple_set_data_lower.repeat_numbers[coordinate_idx] -= 2
        (lower_expectation, expectation_dict) = optimal_expectation(
            tuple_set_data_lower,
            expectation_dict,
            c;
            options = options,
        )
        if lower_expectation > expectation && expectation > upper_expectation
            # Increasing the repetition number decreases the figure of merit
            step_sign = 1
        elseif lower_expectation < expectation && expectation < upper_expectation
            # Decreasing the repetition number decreases the figure of merit
            step_sign = -1
        elseif expectation <= lower_expectation && expectation <= upper_expectation
            # The repetition number is a local minimum in the figure of merit
            step_sign = 0
        else
            # The repetition number is a local maximum in the figure of merit
            if upper_expectation <= lower_expectation
                step_sign = 1
            else
                step_sign = -1
            end
        end
    end
    # Determine the step size
    if step_sign == 0
        step = step_sign
        step_tracker[coordinate_idx] = step_sign
    elseif sign(step_tracker[coordinate_idx]) == step_sign
        # The step sign has not changed
        step = 2 * step_sign * 2^(max(0, abs(step_tracker[coordinate_idx]) - 1))
        step_tracker[coordinate_idx] += step_sign
    else
        # The step sign has changed
        step = 2 * step_sign
        step_tracker[coordinate_idx] = step_sign
    end
    # Step the repetition number and ensure it is within bounds
    @assert iseven(step) "The step size must be even."
    step_repeat_number = old_repeat_number + step
    repeat_tuple_depth = maximum(
        length.(
            tuple_set_data.repeat_tuple_set[[
                tuple_idx for (tuple_idx, idx) in tuple_set_data.repeat_indices if
                idx == coordinate_idx
            ]]
        ),
    )
    cap_repeat_number = cap_repetition(
        step_repeat_number,
        repeat_tuple_depth,
        min_depth,
        max_depth;
        diagnostics = diagnostics,
    )
    tuple_set_data.repeat_numbers[coordinate_idx] = cap_repeat_number
    # Reset the step tracker if the step was capped
    if cap_repeat_number == old_repeat_number
        step_tracker[coordinate_idx] = step_sign
    end
    # Calculate the expectation
    (new_expectation, expectation_dict) =
        optimal_expectation(tuple_set_data, expectation_dict, c; options = options)
    if new_expectation > expectation
        reduced_step_repeat_number = old_repeat_number + 2 * step_sign
        cap_reduced_repeat_number = cap_repetition(
            reduced_step_repeat_number,
            repeat_tuple_depth,
            min_depth,
            max_depth;
            diagnostics = diagnostics,
        )
        tuple_set_data.repeat_numbers[coordinate_idx] = cap_reduced_repeat_number
        step_tracker[coordinate_idx] -= 2 * step_sign
        if diagnostics
            println(
                "Stepping repetition number $(coordinate_idx) worsened the figure of merit; the step size was reduced.",
            )
        end
    end
    return (
        tuple_set_data::TupleSetData,
        expectation_dict::Dict{Vector{Int}, Float64},
        step_tracker::Vector{Int},
    )
end

"""
    optimise_repetitions(c::AbstractCircuit, tuple_set_data::TupleSetData; options::OptimOptions = OptimOptions())

Returns the tuple set data after optimising the repetition numbers in the supplied tuple set data `tuple_set_data` for the circuit `c`.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function optimise_repetitions(
    c::T,
    tuple_set_data::TupleSetData;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Get the keyword arguments
    ls_type = options.ls_type
    est_type = options.rep_est_type
    cycle_convergence_threshold = options.cycle_convergence_threshold
    convergence_cycles = options.convergence_cycles
    min_depth = options.min_depth
    max_depth = options.max_depth
    max_cycles = options.max_cycles
    diagnostics = options.rep_diagnostics
    # Cap the repetition numbers
    for repeat_idx in eachindex(tuple_set_data.repeat_numbers)
        repeat_number = deepcopy(tuple_set_data.repeat_numbers[repeat_idx])
        repeat_tuple_depth = maximum(
            length.(
                tuple_set_data.repeat_tuple_set[[
                    tuple_idx for (tuple_idx, idx) in tuple_set_data.repeat_indices if
                    idx == repeat_idx
                ]]
            ),
        )
        tuple_set_data.repeat_numbers[repeat_idx] = cap_repetition(
            repeat_number,
            repeat_tuple_depth,
            min_depth,
            max_depth;
            diagnostics = diagnostics,
        )
    end
    # Perform cyclic coordinate descent until the repetition numbers converge
    start_time = time()
    if max_cycles > 0
        cycling = true
        cycle = 1
        type_num = length(tuple_set_data.repeat_numbers)
    else
        return tuple_set_data::TupleSetData
    end
    expectation_dict = Dict{Vector{Int}, Float64}()
    expectation_track = Vector{Float64}()
    type_num = length(tuple_set_data.repeat_numbers)
    step_tracker = zeros(Int, type_num)
    old_repeat_numbers = deepcopy(tuple_set_data.repeat_numbers)
    while cycling
        cycle_repeat_numbers = Vector{Vector{Int}}(undef, type_num)
        # Perform a cycle of repetition number steps
        for coordinate_idx in 1:type_num
            (tuple_set_data, expectation_dict, step_tracker) = step_repetition(
                tuple_set_data,
                expectation_dict,
                step_tracker,
                coordinate_idx,
                c;
                options,
            )
            # Update tracking parameters
            repeat_numbers = tuple_set_data.repeat_numbers
            cycle_repeat_numbers[coordinate_idx] = repeat_numbers
        end
        # Check the cycle number and convergence
        expectation = expectation_dict[cycle_repeat_numbers[end]]
        push!(expectation_track, expectation)
        cycling_full_converged = all([
            cycle_repeat_numbers[coordinate_idx] == old_repeat_numbers for
            coordinate_idx in 1:type_num
        ])
        if cycle >= convergence_cycles
            cycling_converged =
                all([
                    abs(1 - expectation_track[idx_2] / expectation_track[idx_1]) <
                    convergence_cycles * cycle_convergence_threshold for
                    idx_1 in (cycle - convergence_cycles + 1):cycle for
                    idx_2 in (cycle - convergence_cycles + 1):cycle
                ]) && (
                    abs(1 - expectation_track[end - 1] / expectation_track[end]) <
                    cycle_convergence_threshold
                )
        else
            cycling_converged = false
        end
        if cycling_full_converged || cycling_converged || cycle >= max_cycles
            cycling = false
            if diagnostics
                if cycling_full_converged
                    println(
                        "Cycling has fully converged after $(cycle) cycles. The $(ls_type) $(est_type) figure of merit is $(round(expectation, sigdigits = 5)), with repetition numbers $(tuple_set_data.repeat_numbers). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                    )
                elseif cycling_converged
                    println(
                        "Cycling has converged after $(cycle) cycles. The $(ls_type) $(est_type) figure of merit is $(round(expectation, sigdigits = 5)), with repetition numbers $(tuple_set_data.repeat_numbers). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                    )
                else
                    println(
                        "The maximum number of cycles $(max_cycles) has been reached without convergence. The $(ls_type) $(est_type) figure of merit is $(round(expectation, sigdigits = 5)), with repetition numbers $(tuple_set_data.repeat_numbers). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                    )
                end
            end
        else
            if diagnostics
                println(
                    "The $(ls_type) $(est_type) figure of merit is $(round(expectation, sigdigits = 5)) after cycle $(cycle), with repetition numbers $(cycle_repeat_numbers[end]). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                )
            end
            old_repeat_numbers = deepcopy(tuple_set_data.repeat_numbers)
            cycle += 1
        end
    end
    return tuple_set_data::TupleSetData
end

"""
    sample_zipf(N::Integer, s::Float64)

Returns a sample from a generalised Zipf distribution supported on 1 to `N` parameterised by the power `s`.
"""
function sample_zipf(N::Integer, s::Float64)
    # Generate the (unnormalised) Zipf PMF
    i_values = collect(1:N)
    zipf_pmf = [1 / i^s for i in i_values]
    # Sample once from the distribution
    zipf_sample = sample(i_values, Weights(zipf_pmf))
    return zipf_sample::Int
end

"""
    tuple_append!(circuit_tuple::Vector{Int}, tuple_length::Integer, s::Float64, unique_indices::Vector{Int})

Appends a random index from `unique_indices` to the tuple `circuit_tuple`, repeated a number of times determined by a Zipf distribution on 1 to `tuple_length` with power `s`.
"""
function tuple_append!(
    circuit_tuple::Vector{Int},
    tuple_length::Integer,
    s::Float64,
    unique_indices::Vector{Int},
)
    # This function should not be used with dynamically decoupled circuits
    # Append a Zipf-distributed number of copies of a random layer index to the tuple
    num_add = sample_zipf(tuple_length, s)
    tuple_add = rand(unique_indices)
    append!(circuit_tuple, repeat([tuple_add], num_add))
    return nothing
end

"""
    tuple_append!(circuit_tuple::Vector{Int}, tuple_length::Integer, s::Float64, unique_indices::Vector{Int}, two_qubit_indices::Vector{Int}, other_indices::Vector{Int})

Appends a random index from `unique_indices` to the tuple `circuit_tuple`, repeated a number of times determined by a Zipf distribution on 1 to `tuple_length` with power `s`.
Ensures that two-qubit layers, whose indices are given by `two_qubit_indices`, are always followed by other layers in the tuple, whose indices are given by `other_indices`.
"""
function tuple_append!(
    circuit_tuple::Vector{Int},
    tuple_length::Integer,
    s::Float64,
    unique_indices::Vector{Int},
    two_qubit_indices::Vector{Int},
    other_indices::Vector{Int},
)
    # This function should only be used with dynamically decoupled circuits
    # Append a Zipf-distributed number of copies, divided by two, of pairs of random layer indices to the tuple
    # Ensure that two-qubit layers are always followed by other layers in the tuple
    if length(circuit_tuple) > 0 && circuit_tuple[end] ∈ two_qubit_indices
        num_add = ceil(Int, sample_zipf(tuple_length, s) / 2)
        tuples_add = [rand(other_indices); rand(unique_indices)]
        append!(circuit_tuple, repeat(tuples_add, num_add))
    else
        num_add = ceil(Int, sample_zipf(tuple_length, s) / 2)
        tuple_add = rand(unique_indices)
        if tuple_add ∈ two_qubit_indices
            tuples_add = [tuple_add; rand(other_indices)]
        else
            tuples_add = [tuple_add; rand(unique_indices)]
        end
        append!(circuit_tuple, repeat(tuples_add, num_add))
    end
    return nothing
end

"""
    random_tuple(c::AbstractCircuit, tuple_length::Integer, s::Float64, mirror::Bool)

Returns a random tuple for the circuit `c` with length `tuple_length`.
The generation is parameterised by the Zipf power `s` and the tuple is mirrored if `mirror` is `true`.
Adds random indices to the tuple, repeated a number of times following a generalised Zipf distribution.
For dynamically decoupled circuits, this function ensures that two-qubit gate layers are always followed by some other layer in the tuple.
"""
function random_tuple(
    c::T,
    tuple_length::Integer,
    s::Float64,
    mirror::Bool,
) where {T <: AbstractCircuit}
    # Set parameters
    two_qubit_type = :two_qubit
    unique_indices = c.unique_layer_indices
    # If the circuit employs dynamical decoupling, ensure the tuples respect that
    tuple_decouple = false
    circuit_params = c.circuit_param.params
    if haskey(circuit_params, :dynamically_decouple) &&
       circuit_params[:dynamically_decouple]
        tuple_decouple = true
        layer_types = c.layer_types
        types = unique(layer_types)
        @assert two_qubit_type ∈ types "The circuit must have a two-qubit gate layer."
        @assert length(types) >= 2 "The circuit must have at least two types of gate layers."
        two_qubit_indices =
            intersect(findall(two_qubit_type .== layer_types), unique_indices)
        other_indices = setdiff(unique_indices, two_qubit_indices)
    end
    circuit_tuple = Vector{Int}()
    if mirror
        # Add the mirror layers
        mirror_length = convert(Int, floor((tuple_length - 1) / 2))
        while length(circuit_tuple) < mirror_length
            if tuple_decouple
                tuple_append!(
                    circuit_tuple,
                    mirror_length,
                    s,
                    unique_indices,
                    two_qubit_indices,
                    other_indices,
                )
            else
                tuple_append!(circuit_tuple, mirror_length, s, unique_indices)
            end
        end
        # Trim extra layers
        circuit_tuple = circuit_tuple[1:mirror_length]
        # Ensure that two-qubit layers aren't repeated in the mirroring
        if tuple_decouple &&
           mirror_length > 0 &&
           circuit_tuple[mirror_length] ∈ two_qubit_indices
            circuit_tuple[mirror_length] = rand(other_indices)
        end
        # Mirror the circuit
        circuit_tuple = [circuit_tuple; reverse(circuit_tuple)]
        # Add the post-mirror layers
        post_mirror_length = tuple_length - 2 * mirror_length
        while length(circuit_tuple) < tuple_length
            if tuple_decouple
                tuple_append!(
                    circuit_tuple,
                    post_mirror_length,
                    s,
                    unique_indices,
                    two_qubit_indices,
                    other_indices,
                )
            else
                tuple_append!(circuit_tuple, post_mirror_length, s, unique_indices)
            end
        end
        # Trim extra layers
        circuit_tuple = circuit_tuple[1:tuple_length]
    else
        # Add the layers
        while length(circuit_tuple) < tuple_length
            if tuple_decouple
                tuple_append!(
                    circuit_tuple,
                    tuple_length,
                    s,
                    unique_indices,
                    two_qubit_indices,
                    other_indices,
                )
            else
                tuple_append!(circuit_tuple, tuple_length, s, unique_indices)
            end
        end
        # Trim extra layers
        circuit_tuple = circuit_tuple[1:tuple_length]
    end
    return circuit_tuple::Vector{Int}
end

"""
    grow_design(d::Design, circuit_tuple::Vector{Int}; weight_experiments::Bool = true)
    grow_design(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}, circuit_tuple::Vector{Int};, weight_experiments::Bool = true)

Returns versions of the design `d` and optionally the circuit log-eigenvalue estimator covariance matrix `covariance_log` after adding the tuple `circuit_tuple` to the tuple set of the design, weighting the shot weights by the experiments if `weight_experiments` is `true` following [`get_tuple_set_params`](@ref).
"""
function grow_design(d::Design, circuit_tuple::Vector{Int}; weight_experiments::Bool = true)
    # Determine the design data for the tuple
    check_tuple!(d.c, circuit_tuple)
    c_tuple = apply_tuple(d.c, circuit_tuple)
    time_1 = time()
    (mapping_set, mapping_matrix) = calc_mapping_set(c_tuple)
    time_2 = time()
    consistency_set = calc_consistency_set(mapping_set)
    time_3 = time()
    experiment_set = calc_experiment_set(mapping_set, consistency_set)
    time_4 = time()
    covariance_dict =
        calc_covariance_dict(c_tuple, mapping_set, experiment_set, d.full_covariance)
    time_5 = time()
    (prep_layer_set, meas_layer_set) =
        get_experiment_layers(c_tuple, mapping_set, experiment_set)
    time_6 = time()
    # Track the times taken
    mapping_time = time_2 - time_1
    consistency_time = time_3 - time_2
    pauli_time = time_4 - time_3
    covariance_time = time_5 - time_4
    calculation_time = [
        mapping_time
        consistency_time
        pauli_time
        covariance_time
    ]
    overall_time = time_6 - time_1
    # Set up design data
    grow_set_data = deepcopy(d.tuple_set_data)
    @reset grow_set_data.tuple_set = [grow_set_data.tuple_set; [circuit_tuple]]
    grow_tuple_set = get_tuple_set(grow_set_data)
    grow_tuple_number = length(grow_tuple_set)
    grow_idx = findfirst(circuit_tuple == grow_tuple for grow_tuple in grow_tuple_set)
    mapping_lengths = length.(d.mapping_ensemble)
    G = sum(mapping_lengths[1:(grow_idx - 1)])
    M = sum(mapping_lengths)
    grow_matrix = vcat(d.matrix[1:G, :], mapping_matrix, d.matrix[(G + 1):M, :])
    grow_mapping_ensemble = Vector{Mapping}[
        d.mapping_ensemble[1:(grow_idx - 1)]
        [mapping_set]
        d.mapping_ensemble[grow_idx:end]
    ]
    grow_experiment_ensemble = Vector{Vector{Int}}[
        d.experiment_ensemble[1:(grow_idx - 1)]
        [experiment_set]
        d.experiment_ensemble[grow_idx:end]
    ]
    grow_covariance_dict_ensemble = [
        d.covariance_dict_ensemble[1:(grow_idx - 1)]
        covariance_dict
        d.covariance_dict_ensemble[grow_idx:end]
    ]
    grow_prep_ensemble = Vector{Layer}[
        d.prep_ensemble[1:(grow_idx - 1)]
        [prep_layer_set]
        d.prep_ensemble[grow_idx:end]
    ]
    grow_meas_ensemble = Vector{Layer}[
        d.meas_ensemble[1:(grow_idx - 1)]
        [meas_layer_set]
        d.meas_ensemble[grow_idx:end]
    ]
    grow_experiment = length(prep_layer_set)
    grow_experiment_numbers = [
        d.experiment_numbers[1:(grow_idx - 1)]
        grow_experiment
        d.experiment_numbers[grow_idx:end]
    ]
    grow_calculation_times = vcat(
        d.calculation_times[1:(grow_idx - 1), :],
        calculation_time',
        d.calculation_times[grow_idx:end, :],
    )
    # Create the shot weights and tuple times
    (grow_tuple_times, default_shot_weights) = get_tuple_set_params(
        d.c,
        grow_tuple_set,
        grow_experiment_numbers;
        weight_experiments = weight_experiments,
    )
    growless_indices = setdiff(1:grow_tuple_number, grow_idx)
    grow_shot_weights = deepcopy(default_shot_weights)
    grow_shot_weights[growless_indices] =
        d.shot_weights * sum(default_shot_weights[growless_indices])
    # Construct the new design
    d_grow = Design(
        d.c,
        d.full_covariance,
        grow_matrix,
        grow_tuple_set,
        grow_set_data,
        grow_mapping_ensemble,
        grow_experiment_ensemble,
        grow_covariance_dict_ensemble,
        grow_prep_ensemble,
        grow_meas_ensemble,
        grow_tuple_times,
        grow_shot_weights,
        grow_experiment_numbers,
        d.experiment_number + grow_experiment,
        grow_calculation_times,
        d.overall_time + overall_time,
        d.optimisation_time,
        d.ls_type,
    )
    return d_grow::Design
end
function grow_design(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int},
    circuit_tuple::Vector{Int};
    weight_experiments::Bool = true,
)
    # Grow the design
    d_grow = grow_design(d, circuit_tuple; weight_experiments = weight_experiments)
    # Grow the covariance matrix
    mapping_lengths = length.(d.mapping_ensemble)
    grow_tuple_number = length(d_grow.tuple_set)
    grow_lengths = length.(d_grow.mapping_ensemble)
    grow_idx = findfirst(circuit_tuple == grow_tuple for grow_tuple in d_grow.tuple_set)
    growless_indices = setdiff(1:grow_tuple_number, grow_idx)
    @assert grow_lengths[growless_indices] == mapping_lengths
    @assert d_grow.tuple_set[growless_indices] == d.tuple_set
    d_extra = Design(
        d.c,
        d.full_covariance,
        d_grow.matrix[
            (sum(grow_lengths[1:(grow_idx - 1)]) + 1):sum(grow_lengths[1:grow_idx]),
            :,
        ],
        [circuit_tuple],
        TupleSetData([circuit_tuple], Vector{Int}[], Int[], Tuple{Int, Int}[]),
        [d_grow.mapping_ensemble[grow_idx]],
        [d_grow.experiment_ensemble[grow_idx]],
        [d_grow.covariance_dict_ensemble[grow_idx]],
        [d_grow.prep_ensemble[grow_idx]],
        [d_grow.meas_ensemble[grow_idx]],
        [d_grow.tuple_times[grow_idx]],
        [1.0],
        [d_grow.experiment_numbers[grow_idx]],
        d_grow.experiment_numbers[grow_idx],
        convert(Matrix{Float64}, d_grow.calculation_times[grow_idx, :]'),
        sum(d_grow.calculation_times[grow_idx, :]),
        0.0,
        :none,
    )
    # Construct the covariance matrix of the extra tuple
    covariance_log_extra = calc_covariance_log(d_extra)
    # Unweight the extra covariance matrix
    @assert d_grow.tuple_times[grow_idx] == d_extra.tuple_times[1]
    covariance_log_extra_unweighted = covariance_log_extra / d_extra.tuple_times[1]
    # Unweight the original covariance matrix
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
    covariance_log_unweighted = covariance_log * shot_weights_factor_inv
    # Construct the new unweighted covariance matrix
    M = sum(mapping_lengths)
    G = sum(mapping_lengths[1:(grow_idx - 1)])
    @assert nnz(covariance_log_unweighted[1:G, (G + 1):M]) == 0
    @assert nnz(covariance_log_unweighted[(G + 1):M, 1:G]) == 0
    covariance_log_grow_unweighted = blockdiag(
        covariance_log_unweighted[1:G, 1:G],
        covariance_log_extra_unweighted,
        covariance_log_unweighted[(G + 1):M, (G + 1):M],
    )
    # Reweight the new covariance matrix
    shot_weights_grow_factor =
        get_shot_weights_factor(d_grow.shot_weights, d_grow.tuple_times, grow_lengths)
    covariance_log_grow = covariance_log_grow_unweighted * shot_weights_grow_factor
    return (d_grow::Design, covariance_log_grow::SparseMatrixCSC{Float64, Int})
end

"""
    prune_design(d::Design, prune_idx::Integer)
    prune_design(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}, prune_idx::Integer; update_weights::Bool = true)

Returns versions of the design `d` and optionally the circuit log-eigenvalue estimator covariance matrix `covariance_log` after removing the tuple at index `prune_idx` from the tuple set of the design.
If `update_weights` is `false`, do not update the shot weight factor for the covariance matrix.
"""
function prune_design(d::Design, prune_idx::Integer)
    # Determine the indices of the tuples and eigenvalues to keep
    pruned_indices = setdiff(1:length(d.tuple_set), prune_idx)
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    prune_eigenvalue_indices =
        vcat([collect(mapping_lower[i]:mapping_upper[i]) for i in pruned_indices]...)
    # Prune the tuple from the tuple set data
    prune_tuple = d.tuple_set[prune_idx]
    @assert prune_tuple ∈ d.tuple_set_data.tuple_set "The pruned tuple is a repeated tuple, which cannot be pruned, or not in the tuple set data."
    prune_set_data = deepcopy(d.tuple_set_data)
    @reset prune_set_data.tuple_set = setdiff(d.tuple_set_data.tuple_set, [prune_tuple])
    @assert d.tuple_set[pruned_indices] == get_tuple_set(prune_set_data) "The pruned tuple set is inconsistent with the pruned tuple set data."
    # Update the shot weights
    prune_shot_weights =
        d.shot_weights[pruned_indices] / sum(d.shot_weights[pruned_indices])
    # Construct the new design
    d_prune = Design(
        d.c,
        d.full_covariance,
        d.matrix[prune_eigenvalue_indices, :],
        d.tuple_set[pruned_indices],
        prune_set_data,
        d.mapping_ensemble[pruned_indices],
        d.experiment_ensemble[pruned_indices],
        d.covariance_dict_ensemble[pruned_indices],
        d.prep_ensemble[pruned_indices],
        d.meas_ensemble[pruned_indices],
        d.tuple_times[pruned_indices],
        prune_shot_weights,
        d.experiment_numbers[pruned_indices],
        sum(d.experiment_numbers[pruned_indices]),
        d.calculation_times[pruned_indices, :],
        d.overall_time - sum(d.calculation_times[prune_idx, :]),
        d.optimisation_time,
        d.ls_type,
    )
    return d_prune::Design
end
function prune_design(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int},
    prune_idx::Integer;
    update_weights::Bool = true,
)
    # Prune the design
    d_prune = prune_design(d, prune_idx)
    # Update the covariance matrix
    pruned_indices = setdiff(1:length(d.tuple_set), prune_idx)
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    prune_eigenvalue_indices =
        vcat([collect(mapping_lower[i]:mapping_upper[i]) for i in pruned_indices]...)
    if update_weights
        # Unweight the covariance matrix
        shot_weights_factor_inv =
            get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
        covariance_log_unweighted = covariance_log * shot_weights_factor_inv
        # Prune the unweighted covariance matrix
        covariance_log_prune_unweighted =
            covariance_log_unweighted[prune_eigenvalue_indices, prune_eigenvalue_indices]
        # Reweight the covariance matrix
        prune_lengths = mapping_lengths[pruned_indices]
        shot_weights_prune_factor = get_shot_weights_factor(
            d_prune.shot_weights,
            d_prune.tuple_times,
            prune_lengths,
        )
        covariance_log_prune = covariance_log_prune_unweighted * shot_weights_prune_factor
    else
        covariance_log_prune =
            covariance_log[prune_eigenvalue_indices, prune_eigenvalue_indices]
    end
    return (d_prune::Design, covariance_log_prune::SparseMatrixCSC{Float64, Int})
end

"""
    grow_design_excursion(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; options::OptimOptions = OptimOptions())

Returns versions of the design `d` and circuit log-eigenvalue estimator covariance matrix `covariance_log` after optimising the tuple set of the design with an excursion that grows the tuple set.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function grow_design_excursion(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    basic_tuple_number = length(get_tuple_set(get_tuple_set_data(d.c)))
    ls_type = options.ls_type
    est_type = options.tuple_est_type
    est_weight = options.tuple_est_weight
    max_tuple_number =
        basic_tuple_number + options.extra_tuple_number + options.excursion_length
    max_tuple_length = options.max_tuple_length
    tuple_length_zipf_power = options.tuple_length_zipf_power
    repeat_zipf_powers = options.repeat_zipf_powers
    mirror_values = options.mirror_values
    trial_factor = options.trial_factor
    grow_greedy = options.grow_greedy
    weight_experiments = options.weight_experiments
    seed = options.seed
    diagnostics = options.tuple_diagnostics
    # Initialise the trial tuples
    if length(d.tuple_set) >= max_tuple_number
        growing = false
        trial_tuple_set = Vector{Vector{Int}}()
    else
        growing = true
        # Generate the set of trial circuit tuples
        if seed !== nothing
            Random.seed!(seed)
        end
        trial_num = trial_factor * (max_tuple_number - length(d.tuple_set))
        trial_tuple_set = Vector{Vector{Int}}()
        while length(trial_tuple_set) < trial_num
            try
                # Check the random tuple to make sure it is supported before adding it to the set
                tuple_length = sample_zipf(max_tuple_length, tuple_length_zipf_power)
                s = rand(repeat_zipf_powers)
                mirror = rand(mirror_values)
                trial_tuple = random_tuple(d.c, tuple_length, s, mirror)
                check_tuple!(d.c, trial_tuple)
                if trial_tuple ∉ d.tuple_set && trial_tuple ∉ trial_tuple_set
                    push!(trial_tuple_set, trial_tuple)
                end
            catch
            end
        end
        @assert length(trial_tuple_set) == trial_num
        if seed !== nothing
            Random.seed!()
        end
    end
    # Calculate the figure of merit if growing greedily
    if grow_greedy
        expectation = calc_ls_moments(
            d,
            covariance_log;
            ls_type = ls_type,
            est_type = est_type,
            est_weight = est_weight,
        )[1]
    end
    # Add each of the tuples from the trial set that improve the figure of merit
    while growing
        # Try adding a tuple to the design
        circuit_tuple = pop!(trial_tuple_set)
        (d_trial, covariance_log_trial) = grow_design(
            d,
            covariance_log,
            circuit_tuple;
            weight_experiments = weight_experiments,
        )
        if grow_greedy
            # Calculate the figure of merit
            expectation_trial = calc_ls_moments(
                d_trial,
                covariance_log_trial;
                ls_type = ls_type,
                est_type = est_type,
                est_weight = est_weight,
            )[1]
            # Add the tuple if it improves the figure of merit
            if expectation_trial < expectation
                d = d_trial
                covariance_log = covariance_log_trial
                expectation = expectation_trial
                if diagnostics
                    println(
                        "The $(ls_type) $(est_type) figure of merit is $(round(expectation, sigdigits = 6)) with $(length(d.tuple_set)) tuples in the set.",
                    )
                end
            end
        else
            d = d_trial
            covariance_log = covariance_log_trial
            if diagnostics
                println("There are now $(length(d.tuple_set)) tuples in the set.")
            end
        end
        if length(d.tuple_set) >= max_tuple_number || length(trial_tuple_set) == 0
            growing = false
        end
    end
    @assert length(d.tuple_set) <= max_tuple_number
    if diagnostics
        if length(d.tuple_set) == max_tuple_number
            println("Grew the tuple set to $(max_tuple_number) tuples.")
        else
            println(
                "Unable to find enough tuples in the $(trial_num) tested to grow the tuple set of $(length(d.tuple_set)) tuples to $(max_tuple_number).",
            )
        end
    end
    return (d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
end

"""
    prune_design_excursion(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; options::OptimOptions = OptimOptions())

Returns versions of the design `d` and circuit log-eigenvalue estimator covariance matrix `covariance_log` after optimising the tuple set of the design with an excursion that prunes the tuple set.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function prune_design_excursion(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    basic_tuple_number = length(get_tuple_set(get_tuple_set_data(d.c)))
    ls_type = options.ls_type
    est_type = options.tuple_est_type
    est_weight = options.tuple_est_weight
    max_tuple_number = basic_tuple_number + options.extra_tuple_number
    diagnostics = options.tuple_diagnostics
    # Set some parameters
    if length(d.tuple_set) >= max_tuple_number
        pruning = true
    else
        pruning = false
    end
    # Greedily prune tuples from the trial set according to the figure of merit
    expectation = calc_ls_moments(
        d,
        covariance_log;
        ls_type = ls_type,
        est_type = est_type,
        est_weight = est_weight,
    )[1]
    while pruning
        # Try removing each of the tuples from the design
        tuple_number = length(d.tuple_set)
        expectation_trials = Array{Float64}(undef, tuple_number)
        for i in 1:tuple_number
            if d.tuple_set[i] ∈ d.tuple_set_data.tuple_set
                (d_trial, covariance_log_trial) = prune_design(d, covariance_log, i)
                # Calculate the figure of merit
                # Sometimes removing a tuple can cause the design to be less than full rank
                # If we get an error, we set the figure of merit to a large number
                try
                    expectation_trial = calc_ls_moments(
                        d_trial,
                        covariance_log_trial;
                        ls_type = ls_type,
                        est_type = est_type,
                        est_weight = est_weight,
                    )[1]
                    expectation_trials[i] = expectation_trial
                catch
                    @debug "Error in calculating the figure of merit after removing the $(t)th tuple; the design matrix is probably no longer full-rank."
                    expectation_trials[i] = 1e20
                end
            else
                expectation_trials[i] = 1e20
            end
        end
        # Prune the tuple it improves the figure of merit or if the desired tuple number has not been reached
        (expectation_min, t_min) = findmin(expectation_trials)
        if expectation_min < 1e20 &&
           (expectation_min < expectation || length(d.tuple_set) > max_tuple_number)
            (d, covariance_log) = prune_design(d, covariance_log, t_min)
            expectation = expectation_min
            if diagnostics
                println(
                    "The $(ls_type) $(est_type) figure of merit is $(round(expectation, sigdigits = 6)) with $(length(d.tuple_set)) tuples in the set.",
                )
            end
        else
            pruning = false
        end
    end
    @assert length(d.tuple_set) <= max_tuple_number
    if diagnostics
        println("Pruned the tuple set until only $(length(d.tuple_set)) remain.")
    end
    return (d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
end

"""
    optimise_tuple_set(d::Design, covariance_log::SparseMatrixCSC{Float64, Int}; options::OptimOptions = OptimOptions())

Returns versions of the design `d` and circuit log-eigenvalue estimator covariance matrix `covariance_log` after optimising the tuple set of the design with repeated excursions that grow and prune the tuple set.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function optimise_tuple_set(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    ls_type = options.ls_type
    est_type = options.tuple_est_type
    est_weight = options.tuple_est_weight
    max_excursions = options.max_excursions
    excursions_unchanged = options.excursions_unchanged
    seed = options.seed
    diagnostics = options.tuple_diagnostics
    # Generate the requisite random seeds using the fixed seed
    start_time = time()
    if seed !== nothing
        Random.seed!(seed)
    end
    seeds = rand(UInt64, max_excursions)
    if seed !== nothing
        Random.seed!()
    end
    # Repeatedly perform excursions
    if max_excursions > 0
        excursioning = true
        unchanged_excursions = 0
        excursion = 1
        expectation_initial = 10^10
    else
        excursioning = false
    end
    while excursioning
        # Grow the design
        d_initial = deepcopy(d)
        covariance_log_initial = deepcopy(covariance_log)
        seeded_options = deepcopy(options)
        @reset seeded_options.seed = seeds[excursion]
        (d, covariance_log) =
            grow_design_excursion(d, covariance_log; options = seeded_options)
        # Prune the design
        (d, covariance_log) = prune_design_excursion(d, covariance_log; options = options)
        # Calculate the figure of merit
        expectation = calc_ls_moments(
            d,
            covariance_log;
            ls_type = ls_type,
            est_type = est_type,
            est_weight = est_weight,
        )[1]
        # Determine whether to stop excursioning
        if expectation > expectation_initial
            d = d_initial
            covariance_log = covariance_log_initial
            if diagnostics
                println(
                    "Excursion $(excursion) worsened the design and was reset, resulting in the same figure of merit $(round(expectation_initial, sigdigits = 6)). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                )
            end
        else
            expectation_initial = expectation
            if diagnostics
                println(
                    "Excursion $(excursion) complete, obtaining a design with figure of merit $(round(expectation, sigdigits = 6)). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                )
            end
        end
        # Determine whether to stop excursioning
        if d == d_initial
            unchanged_excursions += 1
        else
            unchanged_excursions = 0
        end
        if excursion == max_excursions || unchanged_excursions >= excursions_unchanged
            excursioning = false
        else
            excursion += 1
        end
    end
    optimisation_time = time() - start_time
    @reset d.optimisation_time = d.optimisation_time + optimisation_time
    if diagnostics
        println(
            "Design optimisation complete. The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
        )
    end
    return (d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
end

"""
    optimise_design(c::AbstractCircuit; options::OptimOptions = OptimOptions())
    optimise_design(c::AbstractCircuit, tuple_set_data::TupleSetData; options::OptimOptions = OptimOptions())

Returns an optimised experimental design for the circuit `c` initialised with the tuple set data `tuple_set_data`.
The optimisation is parameterised by the [`OptimOptions`](@ref) object `options`.
"""
function optimise_design(
    c::T,
    tuple_set_data::TupleSetData;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Get the keyword arguments
    save_data = options.save_data
    repeat_points = options.repeat_points
    initial_shrink_factor = options.initial_shrink_factor
    rep_diagnostics = options.rep_diagnostics
    weight_experiments = options.weight_experiments
    tuple_diagnostics = options.tuple_diagnostics
    # Optimise the repetitions
    time_1 = time()
    tuple_set_data = optimise_repetitions(c, tuple_set_data; options = options)
    time_2 = time()
    if rep_diagnostics
        println(
            "The time taken to optimise the repetitions is $(round(time_2 - time_1, digits = 3)) s.",
        )
    end
    # Generate the design
    d = generate_design(
        c,
        tuple_set_data;
        full_covariance = true,
        weight_experiments = weight_experiments,
    )
    covariance_log = calc_covariance_log(d)
    # Optimise the tuples in the design
    time_3 = time()
    (d, covariance_log) = optimise_tuple_set(d, covariance_log; options = options)
    time_4 = time()
    if tuple_diagnostics
        println(
            "The time taken to optimise the tuple set is $(round(time_4 - time_3, digits = 3)) s, and the overall time elapsed is $(round(time_4 - time_1, digits = 3)) s.",
        )
    end
    # Augment the tuple set
    if repeat_points > 1
        augmented_tuple_set_data = get_augmented_tuple_set_data(
            d.tuple_set_data,
            repeat_points;
            initial_shrink_factor = initial_shrink_factor,
        )
        d = generate_design(
            c,
            augmented_tuple_set_data;
            full_covariance = true,
            weight_experiments = weight_experiments,
        )
        covariance_log = calc_covariance_log(d)
    end
    # Optimise the shot weights
    # Treat `tuple_diagnostics`, rather than `grad_diagnostics`, as the diagnostics here
    options_grad = deepcopy(options)
    @reset options_grad.grad_diagnostics = tuple_diagnostics
    d = optimise_weights(d, covariance_log; options = options_grad)[1]
    time_5 = time()
    if tuple_diagnostics
        println(
            "The time taken to optimise the shot weights is $(round(time_5 - time_4, digits = 3)) s, and the overall time elapsed is $(round(time_5 - time_1, digits = 3)) s.",
        )
    end
    # Sort the tuples in the tuple set
    optimisation_time = time_5 - time_1
    @reset d.optimisation_time = optimisation_time
    if save_data
        save_design(d)
    end
    return d::Design
end
function optimise_design(
    c::T;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Generate the tuple set data
    error_target = options.error_target
    add_circuit = options.add_circuit
    tuple_set_data =
        get_tuple_set_data(c; error_target = error_target, add_circuit = add_circuit)
    # Optimise the design
    d = optimise_design(c, tuple_set_data; options = options)
    return d::Design
end

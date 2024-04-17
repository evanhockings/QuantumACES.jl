#
function optimal_expectation(
    tuple_set_data::TupleSetData,
    expectation_dict::Dict{Vector{Int}, Float64},
    c::T;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Get the keyword arguments
    ls_type = options.ls_type
    # Retrieve the NRMSE expectation if it has already been calculated, else calculate it
    repeat_numbers = tuple_set_data.repeat_numbers
    if haskey(expectation_dict, repeat_numbers)
        expectation = expectation_dict[repeat_numbers]
    elseif any(repeat_numbers .< 0)
        expectation = 1e20
        expectation_dict[repeat_numbers] = expectation
    else
        # Generate the tuple set and design
        d = generate_design(c, tuple_set_data)
        covariance_log = calc_covariance_log(d)
        # Optimise the design
        (d, covariance_log) = optimise_weights(d, covariance_log; options = options)[1:2]
        # Calculate the NRMSE expectation
        expectation = calc_ls_moments(d, covariance_log, ls_type)[1]
        expectation_dict[repeat_numbers] = expectation
    end
    return (expectation::Float64, expectation_dict::Dict{Vector{Int}, Float64})
end

function step_repetitions(
    tuple_set_data::TupleSetData,
    expectation_dict::Dict{Vector{Int}, Float64},
    step_tracker::Vector{Int},
    coordinate_idx::Int,
    c::T;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Calculate the figure of merit for the current repetition numbers
    (expectation, expectation_dict) =
        optimal_expectation(tuple_set_data, expectation_dict, c; options = options)
    # Calculate the figure of merit adding one to the coordinate's repetition number
    tuple_set_data_upper = deepcopy(tuple_set_data)
    tuple_set_data_upper.repeat_numbers[coordinate_idx] += 2
    (upper_expectation, expectation_dict) =
        optimal_expectation(tuple_set_data_upper, expectation_dict, c; options = options)
    # Determine the direction in which to step, if at all
    if tuple_set_data.repeat_numbers[coordinate_idx] == 0
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
    # Step the repetition number
    tuple_set_data.repeat_numbers[coordinate_idx] =
        max(1, tuple_set_data.repeat_numbers[coordinate_idx] + step)
    (expectation, expectation_dict) =
        optimal_expectation(tuple_set_data, expectation_dict, c; options = options)
    return (
        tuple_set_data::TupleSetData,
        expectation_dict::Dict{Vector{Int}, Float64},
        step_tracker::Vector{Int},
    )
end

function optimise_repetitions(
    c::T,
    tuple_set_data::TupleSetData;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Get the keyword arguments
    max_cycles = options.max_cycles
    diagnostics = options.rep_diagnostics
    # Perform cyclic coordinate descent until the repetition numbers converge
    start_time = time()
    if max_cycles > 0
        cycling = true
        cycle = 1
        type_num = length(tuple_set_data.repeat_numbers)
        if tuple_set_data.repeat_numbers == zeros(Int, type_num)
            tuple_set_data.repeat_numbers = ones(Int, type_num)
        end
    else
        cycling = false
    end
    expectation_dict = Dict{Vector{Int}, Float64}()
    type_num = length(tuple_set_data.repeat_numbers)
    step_tracker = zeros(Int, type_num)
    old_repeat_numbers = deepcopy(tuple_set_data.repeat_numbers)
    while cycling
        cycle_repeat_numbers = Vector{Vector{Int}}(undef, type_num)
        # Perform a cycle of repetition number steps
        for coordinate_idx in 1:type_num
            (tuple_set_data, expectation_dict, step_tracker) = step_repetitions(
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
        cycling_converged = all([
            cycle_repeat_numbers[coordinate_idx] == old_repeat_numbers for
            coordinate_idx in 1:type_num
        ])
        if cycling_converged || cycle >= max_cycles
            cycling = false
            if diagnostics
                if cycling_converged
                    println(
                        "Cycling has converged after $(cycle) cycles. The figure of merit is $(round(expectation, sigdigits = 5)) with repetition numbers $(cycle_repeat_numbers[end]). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                    )
                else
                    println(
                        "The maximum number of cycles $(max_cycles) has been reached without convergence. The figure of merit is $(round(expectation, sigdigits = 5)) with repetition numbers $(cycle_repeat_numbers[end]). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                    )
                end
            end
        else
            if diagnostics
                println(
                    "The figure of merit is $(round(expectation, sigdigits = 5)) after cycle $(cycle), with repetition numbers $(cycle_repeat_numbers[end]). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
                )
            end
            old_repeat_numbers = deepcopy(tuple_set_data.repeat_numbers)
            cycle += 1
        end
    end
    return tuple_set_data::TupleSetData
end

"""
    sample_zipf(N::Int, s::Float64)

Sample from a generalised Zipf distribution supported on 1 to N with parameter s.
"""
function sample_zipf(N::Int, s::Float64)
    # Generate the (unnormalised) Zipf PMF
    i_values = collect(1:N)
    zipf_pmf = [1 / i^s for i in i_values]
    # Sample once from the distribution
    zipf_sample = sample(i_values, Weights(zipf_pmf))
    return zipf_sample::Int
end

#
function tuple_append!(
    circuit_tuple::Vector{Int},
    tuple_length::Int,
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

#
function tuple_append!(
    circuit_tuple::Vector{Int},
    tuple_length::Int,
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
    random_tuple(c::T, tuple_length::Int, s::Float64, mirror::Bool)

Generates a random tuple, or arrangement with repetition, for the `unique_layer_indices` of a circuit whose length is `tuple_length`. Adds random layers to the tuple, with the number of copies following a generalised Zipf distribution; when the parameter `s` is `Inf`, this only adds one copy, and 2 is another common choice. If `mirror`, mirrors the first `floor((tuple_length - 1) / 2)` layers of the circuit.
"""
function random_tuple(
    c::T,
    tuple_length::Int,
    s::Float64,
    mirror::Bool,
) where {T <: AbstractCircuit}
    # Set parameters
    two_qubit_type = :two_qubit
    unique_indices = c.unique_layer_indices
    # If the circuit employs dynamical decoupling, ensure the tuples respect that
    tuple_decouple = false
    if hasproperty(c.circuit_param, :dynamically_decouple) &&
       c.circuit_param.dynamically_decouple
        tuple_decouple = true
        layer_types = c.layer_types
        types = unique(layer_types)
        @assert two_qubit_type ∈ types "The circuit must have a two-qubit gate layer."
        @assert length(types) >= 2 "The circuit must have at least two types of gate layers."
        two_qubit_indices =
            intersect(findall(two_qubit_type .== layer_types), unique_indices)
        other_indices = setdiff(unique_indices, two_qubit_indices)
    end
    circuit_tuple = Vector{Int}(undef, 0)
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

#
function grow_design(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int},
    circuit_tuple::Vector{Int},
)
    # Determine the design data for the tuple
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
    circuit_time = time_6 - time_5
    calculation_time = [
        mapping_time
        consistency_time
        pauli_time
        covariance_time
        circuit_time
    ]
    overall_time = time_6 - time_1
    # Construct the new design
    grow_tuple_set = [d.tuple_set; [circuit_tuple]]
    grow_set_data = deepcopy(d.tuple_set_data)
    grow_set_data.tuple_set = [grow_set_data.tuple_set; [circuit_tuple]]
    grow_experiment = length(vcat(prep_layer_set...))
    grow_experiment_numbers = [d.experiment_numbers; grow_experiment]
    (grow_tuple_times, default_shot_weights) =
        get_tuple_set_params(d.c, grow_tuple_set, grow_experiment_numbers)
    grow_shot_weights = deepcopy(default_shot_weights)
    grow_shot_weights[1:(end - 1)] = d.shot_weights * sum(default_shot_weights[1:(end - 1)])
    d_grow = Design(
        d.c,
        d.full_covariance,
        vcat(d.matrix, mapping_matrix),
        grow_tuple_set,
        grow_set_data,
        [d.mapping_ensemble; [mapping_set]],
        [d.experiment_ensemble; [experiment_set]],
        [d.covariance_dict_ensemble; [covariance_dict]],
        [d.prep_ensemble; [prep_layer_set]],
        [d.meas_ensemble; [meas_layer_set]],
        grow_tuple_times,
        grow_shot_weights,
        grow_experiment_numbers,
        d.experiment_number + grow_experiment,
        vcat(d.calculation_times, calculation_time'),
        d.overall_time + overall_time,
        d.optimisation_time,
        d.ls_type,
    )
    # Grow the covariance matrix
    d_extra = Design(
        d.c,
        d.full_covariance,
        mapping_matrix,
        [circuit_tuple],
        TupleSetData([circuit_tuple], Vector{Int}[], Int[], Int[]),
        [mapping_set],
        [experiment_set],
        [covariance_dict],
        [prep_layer_set],
        [meas_layer_set],
        [grow_tuple_times[end]],
        [1.0],
        [grow_experiment],
        grow_experiment,
        convert(Matrix{Float64}, calculation_time'),
        overall_time,
        0.0,
        :none,
    )
    # Construct the covariance matrix of the extra tuple
    covariance_log_extra = calc_covariance_log(d_extra)
    # Unweight the extra covariance matrix
    covariance_log_extra_unweighted = covariance_log_extra / d_extra.tuple_times[end]
    # Unweight the original covariance matrix
    mapping_lengths = length.(d.mapping_ensemble)
    shot_weights_factor_inv =
        get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
    covariance_log_unweighted = covariance_log * shot_weights_factor_inv
    # Construct the new unweighted covariance matrix
    covariance_log_grow_unweighted =
        blockdiag(covariance_log_unweighted, covariance_log_extra_unweighted)
    # Reweight the new covariance matrix
    grow_lengths = length.(d_grow.mapping_ensemble)
    shot_weights_grow_factor =
        get_shot_weights_factor(d_grow.shot_weights, d_grow.tuple_times, grow_lengths)
    covariance_log_grow = covariance_log_grow_unweighted * shot_weights_grow_factor
    return (d_grow::Design, covariance_log_grow::SparseMatrixCSC{Float64, Int})
end

#
function prune_design(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int},
    prune_idx::Int;
    update_weights::Bool = true,
)
    # Determine the indices of the tuples and eigenvalues to keep
    prune_indices = setdiff(1:length(d.tuple_set), prune_idx)
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    prune_eigenvalue_indices =
        vcat([collect(mapping_lower[i]:mapping_upper[i]) for i in prune_indices]...)
    # Prune the tuple from the tuple set data
    prune_tuple = d.tuple_set[prune_idx]
    prune_set_data = deepcopy(d.tuple_set_data)
    if prune_tuple ∈ d.tuple_set_data.tuple_set
        prune_set_data.tuple_set = setdiff(d.tuple_set_data.tuple_set, [prune_tuple])
    elseif prune_tuple ∈ d.tuple_set
        repeat_prune_idx = findfirst(
            prune_tuple == tuple for
            tuple in setdiff(d.tuple_set, d.tuple_set_data.tuple_set)
        )
        repeat_prune_indices =
            setdiff(1:length(d.tuple_set_data.repeat_tuple_set), repeat_prune_idx)
        # Note that if this rendered some repeat number unused, it is not pruned
        prune_set_data.repeat_tuple_set =
            d.tuple_set_data.repeat_tuple_set[repeat_prune_indices]
        prune_set_data.repeat_indices =
            d.tuple_set_data.repeat_indices[repeat_prune_indices]
    else
        throw(error("The pruned tuple does not appear in the tuple set data."))
    end
    @assert d.tuple_set[prune_indices] == get_tuple_set(prune_set_data)
    # Update the shot weights
    prune_weights = d.shot_weights[prune_indices] / sum(d.shot_weights[prune_indices])
    # Construct the new design
    d_prune = Design(
        d.c,
        d.full_covariance,
        d.matrix[prune_eigenvalue_indices, :],
        d.tuple_set[prune_indices],
        prune_set_data,
        d.mapping_ensemble[prune_indices],
        d.experiment_ensemble[prune_indices],
        d.covariance_dict_ensemble[prune_indices],
        d.prep_ensemble[prune_indices],
        d.meas_ensemble[prune_indices],
        d.tuple_times[prune_indices],
        prune_weights,
        d.experiment_numbers[prune_indices],
        sum(d.experiment_numbers[prune_indices]),
        d.calculation_times[prune_indices, :],
        d.overall_time - sum(d.calculation_times[prune_idx, :]),
        d.optimisation_time,
        d.ls_type,
    )
    # Update the covariance matrix
    if update_weights
        # Unweight the covariance matrix
        shot_weights_factor_inv =
            get_shot_weights_factor_inv(d.shot_weights, d.tuple_times, mapping_lengths)
        covariance_log_unweighted = covariance_log * shot_weights_factor_inv
        # Prune the unweighted covariance matrix
        covariance_log_prune_unweighted =
            covariance_log_unweighted[prune_eigenvalue_indices, prune_eigenvalue_indices]
        # Reweight the covariance matrix
        prune_lengths = mapping_lengths[prune_indices]
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
"""
function grow_design_excursion(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    ls_type = options.ls_type
    max_tuple_number = options.max_tuple_number + options.excursion_length
    max_tuple_length = options.max_tuple_length
    tuple_length_zipf_power = options.tuple_length_zipf_power
    repeat_zipf_powers = options.repeat_zipf_powers
    mirror_values = options.mirror_values
    trial_factor = options.trial_factor
    grow_greedy = options.grow_greedy
    seed = options.seed
    diagnostics = options.tuple_diagnostics
    # Initialise the trial tuples
    if length(d.tuple_set) >= max_tuple_number
        growing = false
        trial_tuple_set = Vector{Vector{Int}}(undef, 0)
    else
        growing = true
        # Generate the set of trial circuit tuples
        if seed !== nothing
            Random.seed!(seed)
        end
        trial_num = trial_factor * (max_tuple_number - length(d.tuple_set))
        trial_tuple_set = Vector{Vector{Int}}(undef, 0)
        while length(trial_tuple_set) < trial_num
            tuple_length = sample_zipf(max_tuple_length, tuple_length_zipf_power)
            s = rand(repeat_zipf_powers)
            mirror = rand(mirror_values)
            trial_tuple = random_tuple(d.c, tuple_length, s, mirror)
            if trial_tuple ∉ d.tuple_set && trial_tuple ∉ trial_tuple_set
                push!(trial_tuple_set, trial_tuple)
            end
        end
        @assert length(trial_tuple_set) == trial_num
        if seed !== nothing
            Random.seed!()
        end
    end
    # Calculate the figure of merit if growing greedily
    if grow_greedy
        expectation = calc_ls_moments(d, covariance_log, ls_type)[1]
    end
    # Add each of the tuples from the trial set that improve the figure of merit
    while growing
        # Try adding a tuple to the design
        circuit_tuple = pop!(trial_tuple_set)
        (d_trial, covariance_log_trial) = grow_design(d, covariance_log, circuit_tuple)
        if grow_greedy
            # Calculate the figure of merit
            expectation_trial = calc_ls_moments(d_trial, covariance_log_trial, ls_type)[1]
            # Add the tuple if it improves the figure of merit
            if expectation_trial < expectation
                d = d_trial
                covariance_log = covariance_log_trial
                expectation = expectation_trial
                if diagnostics
                    println(
                        "The figure of merit is $(round(expectation, sigdigits = 6)) with $(length(d.tuple_set)) tuples in the set.",
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

#
function prune_design_excursion(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    ls_type = options.ls_type
    max_tuple_number = options.max_tuple_number
    diagnostics = options.tuple_diagnostics
    # Set some parameters
    if length(d.tuple_set) >= max_tuple_number
        pruning = true
    else
        pruning = false
    end
    # Greedily prune tuples from the trial set according to the figure of merit
    expectation = calc_ls_moments(d, covariance_log, ls_type)[1]
    while pruning
        # Try removing each of the tuples from the design
        tuple_number = length(d.tuple_set)
        expectation_trial = Array{Float64}(undef, tuple_number)
        for idx in 1:tuple_number
            (d_trial, covariance_log_trial) = prune_design(d, covariance_log, idx)
            # Calculate the figure of merit
            # Sometimes removing a tuple can cause the design to be less than full rank
            # If we get an error, we set the figure of merit to a large number
            try
                expectation_trial[idx] =
                    calc_ls_moments(d_trial, covariance_log_trial, ls_type)[1]
            catch
                @debug "Error in calculating the figure of merit after removing the $(t)th tuple; the design matrix is probably no longer full-rank."
                expectation_trial[idx] = 1e20
            end
        end
        # Prune the tuple it improves the figure of merit or if the desired tuple number has not been reached
        (expectation_min, t_min) = findmin(expectation_trial)
        if expectation_min < 1e20 &&
           (expectation_min < expectation || length(d.tuple_set) > max_tuple_number)
            (d, covariance_log) = prune_design(d, covariance_log, t_min)
            expectation = expectation_min
            if diagnostics
                println(
                    "The figure of merit is $(round(expectation, sigdigits = 6)) with $(length(d.tuple_set)) tuples in the set.",
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

#
function optimise_tuple_set(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int};
    options::OptimOptions = OptimOptions(),
)
    # Get the keyword arguments
    ls_type = options.ls_type
    excursion_number = options.excursion_number
    seed = options.seed
    diagnostics = options.tuple_diagnostics
    # Generate the requisite random seeds using the fixed seed
    start_time = time()
    if seed !== nothing
        Random.seed!(seed)
    end
    seeds = rand(UInt64, excursion_number)
    if seed !== nothing
        Random.seed!()
    end
    # Repeatedly perform excursions
    if excursion_number > 0
        excursioning = true
        excursion = 1
    else
        excursioning = false
    end
    while excursioning
        # Grow the design
        seeded_options = deepcopy(options)
        @reset seeded_options.seed = seeds[excursion]
        (d, covariance_log) =
            grow_design_excursion(d, covariance_log; options = seeded_options)
        # Prune the design
        (d, covariance_log) = prune_design_excursion(d, covariance_log; options = options)
        if diagnostics
            expectation = calc_ls_moments(d, covariance_log, ls_type)[1]
            println(
                "Excursion $(excursion) complete, obtaining a design with figure of merit $(round(expectation, sigdigits = 6)). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
            )
        end
        if excursion == excursion_number
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

function optimise_design(
    c::T,
    tuple_set_data::TupleSetData;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Get the keyword arguments
    save_data = options.save_data
    rep_diagnostics = options.rep_diagnostics
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
    d = generate_design(c, tuple_set_data)
    covariance_log = calc_covariance_log(d)
    # Optimise the tuples in the design
    (d, covariance_log) = optimise_tuple_set(d, covariance_log; options = options)
    time_3 = time()
    if tuple_diagnostics
        println(
            "The time taken to optimise the tuple set is $(round(time_3 - time_2, digits = 3)) s, and the overall time elapsed is $(round(time_3 - time_1, digits = 3)) s.",
        )
    end
    # Optimise the shot weights
    # Treat `tuple_diagnostics` as the diagnostics here, rather than `grad_diagnostics`
    options_copy = deepcopy(options)
    @reset options_copy.grad_diagnostics = tuple_diagnostics
    d = optimise_weights(d, covariance_log; options = options_copy)[1]
    time_4 = time()
    if tuple_diagnostics
        println(
            "The time taken to optimise the shot weights is $(round(time_4 - time_3, digits = 3)) s, and the overall time elapsed is $(round(time_4 - time_1, digits = 3)) s.",
        )
    end
    optimisation_time = time_4 - time_1
    @reset d.optimisation_time = optimisation_time
    if save_data
        save_design(d)
    end
    return d::Design
end

#
function optimise_design(
    c::T;
    options::OptimOptions = OptimOptions(),
) where {T <: AbstractCircuit}
    # Generate the tuple set data
    tuple_set_data = get_tuple_set_data(c)
    # Optimise the design
    d = optimise_design(c, tuple_set_data; options = options)
    return d::Design
end

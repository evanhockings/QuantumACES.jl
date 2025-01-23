"""
    OptimOptions

Keyword arguments for [`optimise_design`](@ref), and specifically the optimisation functions within it, including for the gradient descent function [`optimise_weights`](@ref), the repetition number optimisation function [`optimise_repetitions`](@ref), and the tuple set optimisation function [`optimise_tuple_set`](@ref).

# General options

  - `ls_type::Symbol = :gls`: Type of least squares estimator for which we optimise the design, which can be `:gls`, `:wls`, or `:ols`.
  - `est_type::Symbol = :prod`: Type of estimator for which we optimise the design, which can be `:ordinary`, `:marginal`, `:relative`, or `:sum` or `:prod`, which optimise for the arithmetic or geometric mean, respectively, of the `:ordinary` and `:relative` estimator figures of merit weighted by `est_weight`.
  - `est_weight::Real = 0.5`: Weighting of the `:ordinary` estimator figure of merit in the arithmetic or geometric mean when `est_type` is `:sum` or `:prod`, respectively, such that `:relative` is weighted by `1 - est_weight`.
  - `save_data::Bool = false`: Whether to automatically save the optimised design.

# Gradient descent options

  - `learning_rate::Real = (ls_type == :ols ? 1 : 10.0^(3/4))`: Learning rate for the gradient descent algorithm.
  - `momentum::Real = 0.99`: Momentum for the gradient descent algorithm.
  - `learning_rate_scale_factor::Real = 10.0^(1/4)`: Factor by which to reduce the learning rate if the gradient descent algorithm repeatedly steps in directions that reduce the figure of merit.
  - `max_steps::Integer = 1000`: Maximum number of gradient descent steps to take.
  - `convergence_threshold::Real = 1e-5`: Convergence threshold for the gradient descent algorithm.
  - `convergence_steps::Integer = 5`: Number of steps over which to check convergence.
  - `grad_diagnostics::Bool = false`: Whether to display gradient descent diagnostics.

# Reptition number optimisation options

  - `rep_est_type::Symbol = :relative`: Type of estimator for which we optimise the repetition number, which can be `:ordinary`, `:marginal`, `:relative`, or `:sum` or `:prod`, which optimise for the arithmetic or geometric mean, respectively, of the `:ordinary` and `:relative` estimator figures of merit weighted by `rep_est_weight`.
  - `rep_est_weight::Real = 0.5`: Weighting of the `:ordinary` estimator figure of merit in the arithmetic or geometric mean when `rep_est_type` is `:sum` or `:prod`, respectively, such that `:relative` is weighted by `1 - rep_est_weight`.
  - `error_target::Real = 0.1`: Controls initialisation of the repetition numbers, which heuristically roughly target this error rate; note that setting this to `0` initialises repetition numbers as 0, and all such repetition numbers are then set to 2 before optimisation.
  - `add_circuit::Bool = true`: Whether to add the circuit itself to the repeat tuple set.
  - `cycle_convergence_threshold::Real = 1e-5`: Convergence threshold for the cyclic coordinate descent algorithm for optimising repetition numbers.
  - `convergence_cycles::Integer = 5`: Number of steps over which to check convergence for the cyclic coordinate descent algorithm for optimising repetition numbers.
  - `min_depth::Integer = 0`: Minimum depth of repeated tuples.
  - `max_depth::Integer = 1000`: Maximum depth of repeated tuples.
  - `max_cycles::Integer = 100`: Maximum number of cycles to use in the cyclic coordinate descent algorithm for optimising repetition numbers.
  - `repeat_points::Integer = 1`: Each repeat number is augmented to have `repeat_points` repeat numbers, logarithmically spaced with the smallest repeat number which is shrunk by a factor `initial_shrink_factor`.
  - `initial_shrink_factor::Real = 2^(repeat_points - 1)`: Factor by which the smallest repeat number in an augmented tuple set is shrunk.
  - `rep_diagnostics::Bool = true`: Whether to display repetition number optimisation diagnostics.

# Tuple set optimisation options

  - `tuple_est_type::Symbol = :ordinary`: Type of estimator for which we optimise the tuple set, which can be `:ordinary`, `:marginal`, `:relative`, `:sum`, or `:prod`, which optimise for the arithmetic or geometric mean, respectively, of the `:ordinary` and `:relative` estimator figures of merit weighted by `tuple_est_weight`.
  - `tuple_est_weight::Real = 0.5`: Weighting of the `:ordinary` estimator figure of merit in the arithmetic or geometric mean when `tuple_est_type` is `:sum` or `:prod`, respectively, such that `:relative` is weighted by `1 - tuple_est_weight`.
  - `max_excursions::Integer = 100`: Number of excurisons used to optimise the tuple set.
  - `excursions_unchanged::Integer = 3`: Number of excursions that must not change the tuple set before the optimisation routine terminates.
  - `excursion_length::Integer = 5`: Number of tuples added by each excursion.
  - `extra_tuple_number::Integer = 5`: Number of tuples beyond the number in the basic tuple set added by the tuple set optimisation procedure.
  - `max_tuple_length::Integer = 20`: Maximum length of random tuples.
  - `tuple_length_zipf_power::Real = 1`: Zipf power to use when Zipf-randomly choosing the length of random tuples.
  - `repeat_zipf_powers::Vector{Float64} = [Inf; 2.0]`: Zipf power to use, chosen uniformly at random from the vector, when Zipf-randomly choosing how many times to repeat entries that will be appended to the end of a random tuple during its generation.
  - `mirror_values::Vector{Bool} = [false; true]`: Whether to mirror the tuple, chosen uniformly at random from the vector, when generating random tuples.
  - `trial_factor::Integer = 20`: Number of random tuples trialled for each tuple the excursion needs to add to the tuple set to grow it.
  - `grow_greedy::Bool = true`: Whether the excursions add tuples to the set greedily according to the figure of merit, or to add them even if this reduces the figure of merit.
  - `weight_experiments::Bool = false`: Whether to weight the shot weights for a tuple by the number of experiments for that tuple; while usually `true`, `false` performs better for tuple set optimisation.
  - `seed::Union{UInt64, Nothing} = nothing`: Seed used to randomly generate tuples.
  - `tuple_diagnostics::Bool = true`: Whether to display tuple set optimisation diagnostics.
"""
struct OptimOptions
    # General options
    ls_type::Symbol
    est_type::Symbol
    est_weight::Float64
    save_data::Bool
    # Gradient descent options
    learning_rate::Float64
    momentum::Float64
    learning_rate_scale_factor::Float64
    max_steps::Int
    convergence_threshold::Float64
    convergence_steps::Int
    grad_diagnostics::Bool
    # Reptition options
    rep_est_type::Symbol
    rep_est_weight::Float64
    error_target::Float64
    add_circuit::Bool
    cycle_convergence_threshold::Float64
    convergence_cycles::Int
    min_depth::Int
    max_depth::Int
    max_cycles::Int
    repeat_points::Int
    initial_shrink_factor::Real
    rep_diagnostics::Bool
    # Tuple set options
    tuple_est_type::Symbol
    tuple_est_weight::Float64
    max_excursions::Int
    excursions_unchanged::Int
    excursion_length::Int
    extra_tuple_number::Int
    max_tuple_length::Int
    tuple_length_zipf_power::Float64
    repeat_zipf_powers::Vector{Float64}
    mirror_values::Vector{Bool}
    trial_factor::Int
    grow_greedy::Bool
    weight_experiments::Bool
    seed::Union{UInt64, Nothing}
    tuple_diagnostics::Bool
    # Default constructor
    function OptimOptions(
        ls_type::Symbol,
        est_type::Symbol,
        est_weight::Float64,
        save_data::Bool,
        learning_rate::Float64,
        momentum::Float64,
        learning_rate_scale_factor::Float64,
        max_steps::Integer,
        convergence_threshold::Float64,
        convergence_steps::Integer,
        grad_diagnostics::Bool,
        rep_est_type::Symbol,
        rep_est_weight::Float64,
        error_target::Float64,
        add_circuit::Bool,
        cycle_convergence_threshold::Float64,
        convergence_cycles::Integer,
        min_depth::Integer,
        max_depth::Integer,
        max_cycles::Integer,
        repeat_points::Integer,
        initial_shrink_factor::Real,
        rep_diagnostics::Bool,
        tuple_est_type::Symbol,
        tuple_est_weight::Float64,
        max_excursions::Integer,
        excursions_unchanged::Integer,
        excursion_length::Integer,
        extra_tuple_number::Integer,
        max_tuple_length::Integer,
        tuple_length_zipf_power::Float64,
        repeat_zipf_powers::Vector{Float64},
        mirror_values::Vector{Bool},
        trial_factor::Integer,
        grow_greedy::Bool,
        weight_experiments::Bool,
        seed::Union{UInt64, Nothing},
        tuple_diagnostics::Bool,
    )
        @assert ls_type ∈ [:gls; :wls; :ols] "Must use a valid least squares type."
        @assert est_type ∈ [:ordinary; :marginal; :relative; :sum; :prod] "Must use a valid estimator type."
        @assert est_weight >= 0.0 && est_weight <= 1.0 "The estimator weight must be non-negative and less than or equal to 1."
        @assert learning_rate > 0.0 "The learning rate must be positive."
        @assert momentum >= 0.0 && momentum <= 1.0 "The momentum must be non-negative and less than or equal to 1."
        @assert learning_rate_scale_factor > 1.0 "The learning rate scale factor must be greater than 1."
        @assert max_steps >= 0 "The number of gradient descent steps must be a non-negative integer."
        @assert convergence_threshold > 0.0 "The convergence threshold must be positive."
        @assert convergence_steps >= 1 "The number of steps over which to check convergence must be a positive integer."
        @assert rep_est_type ∈ [:ordinary; :marginal; :relative; :sum; :prod] "Must use a valid estimator type."
        @assert rep_est_weight >= 0.0 && rep_est_weight <= 1.0 "The estimator weight must be non-negative and less than or equal to 1."
        @assert error_target >= 0.0 "The repeat number initialisation scaling must be non-negative."
        @assert cycle_convergence_threshold > 0.0 "The convergence threshold for the cyclic coordinate descent algorithm must be positive."
        @assert convergence_cycles >= 1 "The number of steps over which to check convergence for the cyclic coordinate descent algorithm must be a positive integer."
        @assert min_depth >= 0 "The minimum repeated tuple depth must be a non-negative integer."
        @assert max_depth >= 0 "The maximum repeated tuple depth must be a non-negative integer."
        @assert max_depth >= min_depth "The maximum repeated tuple depth must be at least the minimum repeated tuple depth."
        @assert max_cycles >= 0 "The number of cycles must be a non-negative integer."
        @assert repeat_points >= 1 "The number of repeat points must be a positive integer."
        @assert initial_shrink_factor >= 1 "The initial shrink factor must be at least one."
        @assert tuple_est_type ∈ [:ordinary; :marginal; :relative; :sum; :prod] "Must use a valid estimator type."
        @assert tuple_est_weight >= 0.0 && tuple_est_weight <= 1.0 "The estimator weight must be non-negative and less than or equal to 1."
        @assert max_excursions >= 0 "The maximum excursion number must be a non-negative integer."
        @assert excursions_unchanged >= 1 "The number of excursions that must not change the tuple set before the optimisation routine terminates must be a positive integer."
        @assert excursion_length >= 1 "The excursion length must be a positive integer."
        @assert extra_tuple_number > 0 "The maximum number of tuples must be positive, and if the number is too small the optimisation routine will not produce a full-rank design."
        @assert max_tuple_length >= 2 "The maximum tuple length must be an integer that is at least 2."
        @assert tuple_length_zipf_power >= 0.0 "The tuple length Zipf power must be a non-negative float."
        @assert repeat_zipf_powers == unique(repeat_zipf_powers) "The repeat Zipf powers must be a vector of unique floats."
        @assert all(repeat_zipf_powers .>= 0.0) "The repeat Zipf powers must be non-negative."
        @assert mirror_values == unique(mirror_values) "The mirror values flags must be a vector of unique booleans."
        @assert trial_factor >= 1 "The trial factor must be a positive integer."
        new(
            ls_type,
            est_type,
            est_weight,
            save_data,
            learning_rate,
            momentum,
            learning_rate_scale_factor,
            max_steps,
            convergence_threshold,
            convergence_steps,
            grad_diagnostics,
            rep_est_type,
            rep_est_weight,
            error_target,
            add_circuit,
            cycle_convergence_threshold,
            convergence_cycles,
            min_depth,
            max_depth,
            max_cycles,
            repeat_points,
            initial_shrink_factor,
            rep_diagnostics,
            tuple_est_type,
            tuple_est_weight,
            max_excursions,
            excursions_unchanged,
            excursion_length,
            extra_tuple_number,
            max_tuple_length,
            tuple_length_zipf_power,
            repeat_zipf_powers,
            mirror_values,
            trial_factor,
            grow_greedy,
            weight_experiments,
            seed,
            tuple_diagnostics,
        )
    end
    # Keyword constructor
    function OptimOptions(;
        ls_type::Symbol = :gls,
        est_type::Symbol = :prod,
        est_weight::Real = 0.5,
        save_data::Bool = false,
        learning_rate::Real = (ls_type == :ols ? 1 : (10.0)^(3 / 4)),
        momentum::Real = 0.99,
        learning_rate_scale_factor::Real = 10.0^(1 / 4),
        max_steps::Integer = 1000,
        convergence_threshold::Real = 1e-5,
        convergence_steps::Integer = 5,
        grad_diagnostics::Bool = false,
        rep_est_type::Symbol = :relative,
        rep_est_weight::Real = 0.5,
        error_target::Real = 0.1,
        add_circuit::Bool = true,
        cycle_convergence_threshold::Real = 1e-5,
        convergence_cycles::Integer = 8,
        min_depth::Integer = 0,
        max_depth::Integer = 1000,
        max_cycles::Integer = 100,
        repeat_points::Integer = 1,
        initial_shrink_factor::Real = 2^(repeat_points - 1),
        rep_diagnostics::Bool = true,
        tuple_est_type::Symbol = :ordinary,
        tuple_est_weight::Real = 0.5,
        max_excursions::Integer = 100,
        excursions_unchanged::Integer = 3,
        excursion_length::Integer = 5,
        extra_tuple_number::Integer = 5,
        max_tuple_length::Integer = 20,
        tuple_length_zipf_power::Real = 1,
        repeat_zipf_powers::Vector{Float64} = [Inf; 2.0],
        mirror_values::Vector{Bool} = [false; true],
        trial_factor::Integer = 20,
        grow_greedy::Bool = true,
        weight_experiments::Bool = false,
        seed::Union{UInt64, Nothing} = nothing,
        tuple_diagnostics::Bool = true,
    )
        @assert ls_type ∈ [:gls; :wls; :ols] "Must use a valid least squares type."
        @assert est_type ∈ [:ordinary; :marginal; :relative; :sum; :prod] "Must use a valid estimator type."
        @assert est_weight >= 0.0 && est_weight <= 1.0 "The estimator weight must be non-negative and less than or equal to 1."
        @assert learning_rate > 0.0 "The learning rate must be positive."
        @assert momentum >= 0.0 && momentum <= 1.0 "The momentum must be non-negative and less than or equal to 1."
        @assert learning_rate_scale_factor > 1.0 "The learning rate scale factor must be greater than 1."
        @assert max_steps >= 0 "The number of gradient descent steps must be a non-negative integer."
        @assert convergence_threshold > 0.0 "The convergence threshold must be positive."
        @assert convergence_steps >= 1 "The number of steps over which to check convergence must be a positive integer."
        @assert rep_est_type ∈ [:ordinary; :marginal; :relative; :sum; :prod] "Must use a valid estimator type."
        @assert rep_est_weight >= 0.0 && rep_est_weight <= 1.0 "The estimator weight must be non-negative and less than or equal to 1."
        @assert error_target >= 0.0 "The repeat number initialisation scaling must be non-negative."
        @assert cycle_convergence_threshold > 0.0 "The convergence threshold for the cyclic coordinate descent algorithm must be positive."
        @assert convergence_cycles >= 1 "The number of steps over which to check convergence for the cyclic coordinate descent algorithm must be a positive integer."
        @assert min_depth >= 0 "The minimum repeated tuple depth must be a non-negative integer."
        @assert max_depth >= 0 "The maximum repeated tuple depth must be a non-negative integer."
        @assert max_depth >= min_depth "The maximum repeated tuple depth must be at least the minimum repeated tuple depth."
        @assert max_cycles >= 0 "The number of cycles must be a non-negative integer."
        @assert repeat_points >= 1 "The number of repeat points must be a positive integer."
        @assert initial_shrink_factor >= 1 "The initial shrink factor must be at least one."
        @assert tuple_est_type ∈ [:ordinary; :marginal; :relative; :sum; :prod] "Must use a valid estimator type."
        @assert tuple_est_weight >= 0.0 && tuple_est_weight <= 1.0 "The estimator weight must be non-negative and less than or equal to 1."
        @assert max_excursions >= 0 "The maximum excursion number must be a non-negative integer."
        @assert excursions_unchanged >= 1 "The number of excursions that must not change the tuple set before the optimisation routine terminates must be a positive integer."
        @assert excursion_length >= 1 "The excursion length must be a positive integer."
        @assert extra_tuple_number > 0 "The maximum number of tuples must be positive, and if the number is too small the optimisation routine will not produce a full-rank design."
        @assert max_tuple_length >= 2 "The maximum tuple length must be an integer that is at least 2."
        @assert tuple_length_zipf_power >= 0.0 "The tuple length Zipf power must be a non-negative float."
        @assert repeat_zipf_powers == unique(repeat_zipf_powers) "The repeat Zipf powers must be a vector of unique floats."
        @assert all(repeat_zipf_powers .>= 0.0) "The repeat Zipf powers must be non-negative."
        @assert mirror_values == unique(mirror_values) "The mirror values flags must be a vector of unique booleans."
        @assert trial_factor >= 1 "The trial factor must be a positive integer."
        return new(
            ls_type,
            est_type,
            est_weight,
            save_data,
            learning_rate,
            momentum,
            learning_rate_scale_factor,
            max_steps,
            convergence_threshold,
            convergence_steps,
            grad_diagnostics,
            rep_est_type,
            rep_est_weight,
            error_target,
            add_circuit,
            cycle_convergence_threshold,
            convergence_cycles,
            min_depth,
            max_depth,
            max_cycles,
            repeat_points,
            initial_shrink_factor,
            rep_diagnostics,
            tuple_est_type,
            tuple_est_weight,
            max_excursions,
            excursions_unchanged,
            excursion_length,
            extra_tuple_number,
            max_tuple_length,
            tuple_length_zipf_power,
            repeat_zipf_powers,
            mirror_values,
            trial_factor,
            grow_greedy,
            weight_experiments,
            seed,
            tuple_diagnostics,
        )
    end
end

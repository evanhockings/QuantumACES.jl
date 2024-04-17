# This provides keywords for design optimsation functions in weights.jl and optimise.jl
struct OptimOptions
    # General options
    ls_type::Symbol
    save_data::Bool
    # Gradient descent options
    learning_rate::Float64
    momentum::Float64
    learning_rate_scale_factor::Float64
    shot_weights_clip::Float64
    max_steps::Int
    convergence_threshold::Float64
    convergence_steps::Int
    grad_diagnostics::Bool
    # Reptition options
    max_cycles::Int
    rep_diagnostics::Bool
    # Tuple set options
    excursion_number::Int
    excursion_length::Int
    max_tuple_number::Int
    max_tuple_length::Int
    tuple_length_zipf_power::Float64
    repeat_zipf_powers::Vector{Float64}
    mirror_values::Vector{Bool}
    trial_factor::Int
    grow_greedy::Bool
    seed::Union{UInt64, Nothing}
    tuple_diagnostics::Bool
    # Default constructor
    function OptimOptions(
        ls_type::Symbol,
        save_data::Bool,
        learning_rate::Float64,
        momentum::Float64,
        learning_rate_scale_factor::Float64,
        shot_weights_clip::Float64,
        max_steps::Int,
        convergence_threshold::Float64,
        convergence_steps::Int,
        grad_diagnostics::Bool,
        max_cycles::Int,
        rep_diagnostics::Bool,
        excursion_number::Int,
        excursion_length::Int,
        max_tuple_number::Int,
        max_tuple_length::Int,
        tuple_length_zipf_power::Float64,
        repeat_zipf_powers::Vector{Float64},
        mirror_values::Vector{Bool},
        trial_factor::Int,
        grow_greedy::Bool,
        seed::Union{UInt64, Nothing},
        tuple_diagnostics::Bool,
    )
        @assert ls_type ∈ [:gls; :wls; :ols] "Must use a valid least squares type."
        @assert learning_rate > 0.0 "The learning rate must be positive."
        @assert momentum >= 0.0 && momentum <= 1.0 "The momentum must be non-negative and less than or equal to 1."
        @assert learning_rate_scale_factor > 1.0 "The learning rate scale factor must be greater than 1."
        @assert shot_weights_clip >= 0.0 "The minimum shot weight must be non-negative."
        @assert max_steps >= 0 "The number of gradient descent steps must be a non-negative integer."
        @assert convergence_threshold > 0.0 "The convergence threshold must be positive."
        @assert convergence_steps >= 1 "The number of steps over which to check convergence must be a positive integer."
        @assert max_cycles >= 0 "The number of cycles must be a non-negative integer."
        @assert excursion_number >= 0 "The excursion number must be a non-negative integer."
        @assert excursion_length >= 1 "The excursion length must be a positive integer."
        @assert max_tuple_number > 0 "The maximum number of tuples must be positive, and if the number is too small the optimisation routine will not produce a full-rank design."
        @assert max_tuple_length >= 2 "The maximum tuple length must be an integer that is at least 2."
        @assert tuple_length_zipf_power >= 0.0 "The tuple length Zipf power must be a non-negative float."
        @assert repeat_zipf_powers == unique(repeat_zipf_powers) "The repeat Zipf powers must be a vector of unique floats."
        @assert all(repeat_zipf_powers .>= 0.0) "The repeat Zipf powers must be non-negative."
        @assert mirror_values == unique(mirror_values) "The mirror values flags must be a vector of unique booleans."
        @assert trial_factor >= 1 "The trial factor must be a positive integer."
        new(
            ls_type,
            save_data,
            learning_rate,
            momentum,
            learning_rate_scale_factor,
            shot_weights_clip,
            max_steps,
            convergence_threshold,
            convergence_steps,
            grad_diagnostics,
            max_cycles,
            rep_diagnostics,
            excursion_number,
            excursion_length,
            max_tuple_number,
            max_tuple_length,
            tuple_length_zipf_power,
            repeat_zipf_powers,
            mirror_values,
            trial_factor,
            grow_greedy,
            seed,
            tuple_diagnostics,
        )
    end
    # Keyword constructor
    function OptimOptions(;
        ls_type::Symbol = :wls,
        save_data::Bool = false,
        learning_rate::Float64 = (ls_type == :ols ? 1.0 : (10.0)^(3 / 4)),
        momentum::Float64 = 0.99,
        learning_rate_scale_factor::Float64 = 10.0^(1 / 4),
        shot_weights_clip::Float64 = 1e-5,
        max_steps::Int = 200,
        convergence_threshold::Float64 = 1e-5,
        convergence_steps::Int = 5,
        grad_diagnostics::Bool = false,
        max_cycles::Int = 50,
        rep_diagnostics::Bool = true,
        excursion_number::Int = 5,
        excursion_length::Int = 5,
        max_tuple_number::Int = 35,
        max_tuple_length::Int = 20,
        tuple_length_zipf_power::Float64 = 1.0,
        repeat_zipf_powers::Vector{Float64} = [Inf, 2.0],
        mirror_values::Vector{Bool} = [false; true],
        trial_factor::Int = 20,
        grow_greedy::Bool = true,
        seed::Union{UInt64, Nothing} = nothing,
        tuple_diagnostics::Bool = true,
    )
        @assert ls_type ∈ [:gls; :wls; :ols] "Must use a valid least squares type."
        @assert learning_rate > 0.0 "The learning rate must be positive."
        @assert momentum >= 0.0 && momentum <= 1.0 "The momentum must be non-negative and less than or equal to 1."
        @assert learning_rate_scale_factor > 1.0 "The learning rate scale factor must be greater than 1."
        @assert shot_weights_clip >= 0.0 "The minimum shot weight must be non-negative."
        @assert max_steps >= 0 "The number of gradient descent steps must be a non-negative integer."
        @assert convergence_threshold > 0.0 "The convergence threshold must be positive."
        @assert convergence_steps >= 1 "The number of steps over which to check convergence must be a positive integer."
        @assert max_cycles >= 0 "The number of cycles must be a non-negative integer."
        @assert excursion_number >= 0 "The excursion number must be a non-negative integer."
        @assert excursion_length >= 1 "The excursion length must be a positive integer."
        @assert max_tuple_number > 0 "The maximum number of tuples must be positive, and if the number is too small the optimisation routine will not produce a full-rank design."
        @assert max_tuple_length >= 2 "The maximum tuple length must be an integer that is at least 2."
        @assert tuple_length_zipf_power >= 0.0 "The tuple length Zipf power must be a non-negative float."
        @assert repeat_zipf_powers == unique(repeat_zipf_powers) "The repeat Zipf powers must be a vector of unique floats."
        @assert all(repeat_zipf_powers .>= 0.0) "The repeat Zipf powers must be non-negative."
        @assert mirror_values == unique(mirror_values) "The mirror values flags must be a vector of unique booleans."
        @assert trial_factor >= 1 "The trial factor must be a positive integer."
        return new(
            ls_type,
            save_data,
            learning_rate,
            momentum,
            learning_rate_scale_factor,
            shot_weights_clip,
            max_steps,
            convergence_threshold,
            convergence_steps,
            grad_diagnostics,
            max_cycles,
            rep_diagnostics,
            excursion_number,
            excursion_length,
            max_tuple_number,
            max_tuple_length,
            tuple_length_zipf_power,
            repeat_zipf_powers,
            mirror_values,
            trial_factor,
            grow_greedy,
            seed,
            tuple_diagnostics,
        )
    end
end

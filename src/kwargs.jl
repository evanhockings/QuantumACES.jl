# This provides a number of functions that parse keyword arguments
# These simplify the otherwise messy keyword argument passing in gradient.jl and optimise.jl

#
function LeastSquaresType(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :ls_type)
        ls_type = kwarg_dict[:ls_type]
        @assert ls_type ∈ [:gls; :wls; :ols] "Must use a valid least squares type."
    else
        ls_type = :wls
    end
    return ls_type::Symbol
end

#
function LearningRate(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :learning_rate)
        learning_rate = kwarg_dict[:learning_rate]
        @assert learning_rate > 0.0 "The learning rate must be positive."
    else
        learning_rate = 1e-1
    end
    return learning_rate
end

#
function Momentum(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :momentum)
        momentum = kwarg_dict[:momentum]
        @assert momentum >= 0.0 && momentum <= 1.0 "The momentum must be non-negative and less than or equal to 1."
    else
        momentum = 0.9
    end
    return momentum
end

#
function ConvergenceThreshold(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :convergence_threshold)
        convergence_threshold = kwarg_dict[:convergence_threshold]
        @assert convergence_threshold > 0.0 "The convergence threshold must be positive."
    else
        convergence_threshold = 1e-5
    end
    return convergence_threshold
end

#
function ConvergenceSteps(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :convergence_steps)
        convergence_steps = kwarg_dict[:convergence_steps]
        @assert typeof(convergence_steps) <: Integer && convergence_steps >= 1 "The number of steps over which to check convergence must be a positive integer."
    else
        convergence_steps = 5
    end
    return convergence_steps::Int
end

#
function MaxSteps(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :max_steps)
        max_steps = kwarg_dict[:max_steps]
        @assert typeof(max_steps) <: Integer && max_steps >= 0 "The number of gradient descent steps must be a non-negative integer."
    else
        max_steps = 200
    end
    return max_steps::Int
end

#
function LearningRateScaleFactor(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :learning_rate_scale_factor)
        learning_rate_scale_factor = kwarg_dict[:learning_rate_scale_factor]
        @assert learning_rate_scale_factor > 1.0 "The learning rate scale factor must be greater than 1."
    else
        learning_rate_scale_factor = 2.0
    end
    return learning_rate_scale_factor
end

#
function ShotWeightsNoise(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :shot_weights_noise)
        shot_weights_noise = kwarg_dict[:shot_weights_noise]
        @assert shot_weights_noise >= 0.0 "The shot weight noise scaling factor must be non-negative."
    else
        shot_weights_noise = 0.0
    end
    return shot_weights_noise
end

#
function ShotWeightsClip(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :shot_weights_clip)
        shot_weights_clip = kwarg_dict[:shot_weights_clip]
        @assert shot_weights_clip >= 0.0 "The minimum shot weight must be non-negative."
    else
        shot_weights_clip = 1e-5
    end
    return shot_weights_clip
end

#
function MaxCycles(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :max_cycles)
        max_cycles = kwarg_dict[:max_cycles]
        @assert typeof(max_cycles) <: Integer && max_cycles >= 0 "The number of cycles must be a non-negative integer."
    else
        max_cycles = 100
    end
    return max_cycles::Int
end

#
function Diagnostics(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :diagnostics)
        diagnostics = kwarg_dict[:diagnostics]
        @assert typeof(diagnostics) <: Bool "The diagnostics flag must be a boolean."
    else
        diagnostics = nothing
    end
    return diagnostics
end

#
function GradDiagnostics(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :grad_diagnostics)
        grad_diagnostics = kwarg_dict[:grad_diagnostics]
        @assert typeof(grad_diagnostics) <: Bool "The gradient diagnostics flag must be a boolean."
    else
        grad_diagnostics = nothing
    end
    return grad_diagnostics
end

#
function SaveData(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :save_data)
        save_data = kwarg_dict[:save_data]
        @assert typeof(save_data) <: Bool "The save data flag must be a boolean."
    else
        save_data = nothing
    end
    return save_data
end

#
function TupleNumber(kwarg_dict::Dict{Symbol, Any}, code::Code)
    if haskey(kwarg_dict, :tuple_num)
        tuple_num = kwarg_dict[:tuple_num]
        @assert typeof(tuple_num) <: Integer && tuple_num >= length(code.circuit) + 3 "The number of tuples must be a positive integer that is at least two larger than the number required for a trivial full-rank design."
    else
        tuple_num = 5 * length(code.unique_layer_indices)
    end
    return tuple_num::Int
end

#
function ExcursionLength(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :excursion_len)
        excursion_len = kwarg_dict[:excursion_len]
        @assert typeof(excursion_len) <: Integer && excursion_len >= 1 "The excursion length must be a positive integer."
    else
        excursion_len = 10
    end
    return excursion_len::Int
end

#
function ExcursionNumber(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :excursion_num)
        excursion_num = kwarg_dict[:excursion_num]
        @assert typeof(excursion_num) <: Integer && excursion_num >= 0 "The excursion number must be a non-negative integer."
    else
        excursion_num = 3
    end
    return excursion_num::Int
end

#
function MaxTupleLength(kwarg_dict::Dict{Symbol, Any}, code::Code)
    if haskey(kwarg_dict, :max_tuple_len)
        max_tuple_len = kwarg_dict[:max_tuple_len]
        @assert typeof(max_tuple_len) <: Integer && max_tuple_len >= 2 "The maximum tuple length must be an integer greater than or equal to 2."
    else
        max_tuple_len = 2 * length(code.circuit)
    end
    return max_tuple_len::Int
end

#
function TupleLengthZipfPower(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :tuple_length_zipf_power)
        tuple_length_zipf_power = kwarg_dict[:tuple_length_zipf_power]
        @assert typeof(tuple_length_zipf_power) <: Float64 && tuple_length_zipf_power >= 0.0 "The tuple length Zipf power must be a non-negative float."
    else
        tuple_length_zipf_power = 1.0
    end
    return tuple_length_zipf_power::Float64
end

#
function RepeatZipfPowers(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :repeat_zipf_powers)
        repeat_zipf_powers = kwarg_dict[:repeat_zipf_powers]
        @assert typeof(repeat_zipf_powers) <: Vector{Float64} &&
                length(repeat_zipf_powers) >= 1 &&
                repeat_zipf_powers == unique(repeat_zipf_powers) "The repeat Zipf powers must be a vector of unique floats."
    else
        repeat_zipf_powers = [Inf, 2.0]
    end
    return repeat_zipf_powers::Vector{Float64}
end

#
function MirrorValues(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :mirror_values)
        mirror_values = kwarg_dict[:mirror_values]
        @assert typeof(mirror_values) <: Vector{Bool} &&
                length(mirror_values) >= 1 &&
                mirror_values == unique(mirror_values) "The mirror values flags must be a vector of unique booleans."
    else
        mirror_values = [false; true]
    end
    return mirror_values::Vector{Bool}
end

#
function TrialFactor(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :trial_factor)
        trial_factor = kwarg_dict[:trial_factor]
        @assert typeof(trial_factor) <: Integer && trial_factor >= 1 "The trial factor must be a positive integer."
    else
        trial_factor = 20
    end
    return trial_factor::Int
end

#
function GrowGreedy(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :grow_greedy)
        grow_greedy = kwarg_dict[:grow_greedy]
        @assert typeof(grow_greedy) <: Bool "The greedy growth flag must be a boolean."
    else
        grow_greedy = true
    end
    return grow_greedy::Bool
end

#
function Seed(kwarg_dict::Dict{Symbol, Any})
    if haskey(kwarg_dict, :seed)
        seed = kwarg_dict[:seed]
        @assert typeof(seed) <: UInt64 "The seed must be a UInt64."
    else
        seed = nothing
    end
    return seed
end

"""
    LognormalParameters

Parameterises a log-normally random Pauli noise model.

# Fields

  - `params::Dict{Symbol, Any}`: Dictionary of the noise parameters described below.
  - `noise_name::String`: Noise parameter name for saving data.

# Parameters

  - `r_1::Float64`: Average single-qubit gate entanglement infidelity, the sum of all 3 non-identity Pauli error probabilities.
  - `r_2::Float64`: Average two-qubit gate entanglement infidelity, the sum of all 15 non-identity Pauli error probabilities.
  - `r_m::Float64`: Average measurement error probability.
  - `r_im::Float64`: Average mid-circuit measurement idle entanglement infidelity.
  - `r_r::Float64`: Average mid-circuit reset error probability.
  - `r_1_std_log::Float64`: Approximate standard deviation of the logarithm of the single-qubit gate entanglement infidelity.
  - `r_2_std_log::Float64`: Approximate standard deviation of the logarithm of the two-qubit gate entanglement infidelity.
  - `r_m_std_log::Float64`: Approximate standard deviation of the logarithm of the measurement error probability.
  - `r_im_std_log::Float64`: Approximate standard deviation of the logarithm of the mid-circuit measurement idle entanglement infidelity.
  - `r_r_std_log::Float64`: Approximate standard deviation of the logarithm of the mid-circuit reset error probability.
  - `seed::UInt64`: Random seed used to generate the noise.
  - `combined::Bool`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
"""
struct LognormalParameters <: AbstractNoiseParameters
    params::Dict{Symbol, Any}
    noise_name::String
    # Default constructor
    function LognormalParameters(params::Dict{Symbol, Any}, noise_name::String)
        # Check noise parameters are present
        @assert haskey(params, :r_1) "The single-qubit gate entanglement infidelity is missing."
        @assert haskey(params, :r_2) "The two-qubit gate entanglement infidelity is missing."
        @assert haskey(params, :r_m) "The measurement entanglement infidelity is missing."
        @assert haskey(params, :r_1_std_log) "The standard deviation of the logarithm of the single-qubit gate entanglement infidelity is missing."
        @assert haskey(params, :r_2_std_log) "The standard deviation of the logarithm of the two-qubit gate entanglement infidelity is missing."
        @assert haskey(params, :r_m_std_log) "The standard deviation of the logarithm of the measurement entanglement infidelity is missing."
        @assert haskey(params, :seed) "The random seed is missing."
        @assert haskey(params, :combined) "The combined flag is missing."
        r_1 = params[:r_1]
        r_2 = params[:r_2]
        r_m = params[:r_m]
        r_im = params[:r_im]
        r_r = params[:r_r]
        r_1_std_log = params[:r_1_std_log]
        r_2_std_log = params[:r_2_std_log]
        r_m_std_log = params[:r_m_std_log]
        r_im_std_log = params[:r_im_std_log]
        r_r_std_log = params[:r_r_std_log]
        seed = params[:seed]
        combined = params[:combined]
        # Check some conditions
        @assert (r_1 >= 0) && (r_1 <= 3 / 4) "The single-qubit gate entanglement infidelity $(r_1) is out of bounds."
        @assert (r_2 >= 0) && (r_2 <= 15 / 16) "The two-qubit gate entanglement infidelity $(r_2) is out of bounds."
        @assert (r_m >= 0) && (r_m <= 1 / 2) "The measurement entanglement infidelity $(r_m) is out of bounds."
        @assert (r_im >= 0) && (r_im <= 3 / 4) "The mid-circuit measurement idle entanglement infidelity $(r_im) is out of bounds."
        @assert (r_r >= 0) && (r_r <= 1 / 2) "The mid-circuit reset entanglement infidelity $(r_r) is out of bounds."
        @assert r_1_std_log >= 0 "The standard deviation of the logarithm of the single-qubit gate entanglement infidelity $(r_1_std_log) is out of bounds."
        @assert r_2_std_log >= 0 "The standard deviation of the logarithm of the two-qubit gate entanglement infidelity $(r_2_std_log) is out of bounds."
        @assert r_m_std_log >= 0 "The standard deviation of the logarithm of the measurement entanglement infidelity $(r_m_std_log) is out of bounds."
        @assert r_im_std_log >= 0 "The standard deviation of the logarithm of the mid-circuit measurement idle entanglement infidelity $(r_im_std_log) is out of bounds."
        @assert r_r_std_log >= 0 "The standard deviation of the logarithm of the mid-circuit reset entanglement infidelity $(r_r_std_log) is out of bounds."
        @assert typeof(seed) == UInt64 "The random seed $(seed) is not a UInt64."
        @assert typeof(combined) == Bool "The combined flag $(combined) is not a Bool."
        # Return parameters with the appropriate name
        sigdigits = 3
        new_noise_name = "lognormal_$(round(r_1; sigdigits = sigdigits))_$(round(r_2; sigdigits = sigdigits))_$(round(r_m; sigdigits = sigdigits))_$(round(r_im; sigdigits = sigdigits))_$(round(r_r; sigdigits = sigdigits))_$(round(r_1_std_log; sigdigits = sigdigits))_$(round(r_2_std_log; sigdigits = sigdigits))_$(round(r_m_std_log; sigdigits = sigdigits))_$(round(r_im_std_log; sigdigits = sigdigits))_$(round(r_r_std_log; sigdigits = sigdigits))_$(seed)_$(combined)"
        return new(params, new_noise_name)::LognormalParameters
    end
end

@struct_hash_equal_isequal LognormalParameters

"""
    get_log_param(r_1::Real, r_2::Real, r_m::Real, std_log::Real; seed::Union{UInt64, Nothing} = nothing, combined::Bool = false)
    get_log_param(; r_1::Real, r_2::Real, r_m::Real, r_im::Real = r_m, r_r::Real = r_m, r_1_std_log::Real, r_2_std_log::Real, r_m_std_log::Real, r_im_std_log::Real = r_m_std_log, r_r_std_log::Real = r_m_std_log, seed::Union{UInt64, Nothing} = nothing, combined::Bool = false)

Return a [`LognormalParameters`](@ref) object that parameterises a log-normally distributed random Pauli noise model.

Optionally specify all standard deviations as `std_log`.

# Arguments

  - `r_1::Float64`: Average single-qubit gate entanglement infidelity, the sum of all 3 non-identity Pauli error probabilities.
  - `r_2::Float64`: Average two-qubit gate entanglement infidelity, the sum of all 15 non-identity Pauli error probabilities.
  - `r_m::Float64`: Average measurement error probability.
  - `r_im::Float64`: Average mid-circuit measurement idle entanglement infidelity.
  - `r_r::Float64`: Average mid-circuit reset error probability.
  - `r_1_std_log::Float64`: Approximate standard deviation of the logarithm of the single-qubit gate entanglement infidelity.
  - `r_2_std_log::Float64`: Approximate standard deviation of the logarithm of the two-qubit gate entanglement infidelity.
  - `r_m_std_log::Float64`: Approximate standard deviation of the logarithm of the measurement error probability.
  - `r_im_std_log::Float64`: Approximate standard deviation of the logarithm of the mid-circuit measurement idle entanglement infidelity.
  - `r_r_std_log::Float64`: Approximate standard deviation of the logarithm of the mid-circuit reset error probability.
  - `seed::UInt64`: Random seed used to generate the noise.
  - `combined::Bool`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
"""
function get_log_param(;
    r_1::Real,
    r_2::Real,
    r_m::Real,
    r_im::Real = r_m,
    r_r::Real = r_m,
    r_1_std_log::Real,
    r_2_std_log::Real,
    r_m_std_log::Real,
    r_im_std_log::Real = r_m_std_log,
    r_r_std_log::Real = r_m_std_log,
    seed::Union{UInt64, Nothing} = nothing,
    combined::Bool = false,
)
    if seed === nothing
        seed = rand(UInt64)
    end
    params = Dict{Symbol, Any}(
        :r_1 => r_1,
        :r_2 => r_2,
        :r_m => r_m,
        :r_im => r_im,
        :r_r => r_r,
        :r_1_std_log => r_1_std_log,
        :r_2_std_log => r_2_std_log,
        :r_m_std_log => r_m_std_log,
        :r_im_std_log => r_im_std_log,
        :r_r_std_log => r_r_std_log,
        :seed => seed,
        :combined => combined,
    )
    log_param = LognormalParameters(params, "lognormal")
    return log_param::LognormalParameters
end
function get_log_param(
    r_1::Real,
    r_2::Real,
    r_m::Real,
    std_log::Real;
    seed::Union{UInt64, Nothing} = nothing,
    combined::Bool = false,
)
    log_param = get_log_param(;
        r_1 = r_1,
        r_2 = r_2,
        r_m = r_m,
        r_1_std_log = std_log,
        r_2_std_log = std_log,
        r_m_std_log = std_log,
        seed = seed,
        combined = combined,
    )
    return log_param::LognormalParameters
end

"""
    init_gate_probabilities(total_gates::Vector{Gate}, log_param::LognormalParameters)

Returns a dictionary of the Pauli error probabilities for each gate in `total_gates` generated according to the noise parameters `log_param`, ordered lexicographically following [`get_paulis`](@ref).
"""
function init_gate_probabilities(total_gates::Vector{Gate}, log_param::LognormalParameters)
    # Extract the parameters for generating the noise
    r_1 = log_param.params[:r_1]
    r_2 = log_param.params[:r_2]
    r_m = log_param.params[:r_m]
    r_im = log_param.params[:r_im]
    r_r = log_param.params[:r_r]
    r_1_std_log = log_param.params[:r_1_std_log]
    r_2_std_log = log_param.params[:r_2_std_log]
    r_m_std_log = log_param.params[:r_m_std_log]
    r_im_std_log = log_param.params[:r_im_std_log]
    r_r_std_log = log_param.params[:r_r_std_log]
    seed = log_param.params[:seed]
    combined = log_param.params[:combined]
    # We approximate the sum of log-normal random variables as a log-normal random variable
    # We match the mean and variance of the sum
    # We ensure that the overall gate infidelity has the specified standard deviation
    p_1_std_log = sqrt(log(1 + 3 * (exp(r_1_std_log^2) - 1)))
    p_2_std_log = sqrt(log(1 + 15 * (exp(r_2_std_log^2) - 1)))
    p_m_std_log = r_m_std_log
    p_im_std_log = sqrt(log(1 + 3 * (exp(r_im_std_log^2) - 1)))
    p_r_std_log = r_r_std_log
    p_1_mean_log = log(r_1 / 3) - p_1_std_log^2 / 2
    p_2_mean_log = log(r_2 / 15) - p_2_std_log^2 / 2
    p_m_mean_log = log(r_m) - p_m_std_log^2 / 2
    p_im_mean_log = log(r_im / 3) - p_im_std_log^2 / 2
    p_r_mean_log = log(r_r) - p_r_std_log^2 / 2
    # Fix the random seed
    Random.seed!(seed)
    # Generate the noise
    spam_gates = Gate[]
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for (idx, gate) in pairs(total_gates)
        if is_spam(gate)
            gate_probs = [exp(p_m_std_log * randn() + p_m_mean_log)]
        elseif is_mid_reset(gate)
            gate_probs = [exp(p_r_std_log * randn() + p_r_mean_log)]
        elseif is_meas_idle(gate)
            gate_probs = exp.(p_im_std_log * randn(3) .+ p_im_mean_log)
        elseif is_mid_meas(gate)
            throw(error("Mid-circuit measurement gates, such as $(gate), are unsupported."))
        else
            gate_support_size = length(gate.targets)
            if gate_support_size == 1
                gate_probs = exp.(p_1_std_log * randn(3) .+ p_1_mean_log)
            elseif gate_support_size == 2
                gate_probs = exp.(p_2_std_log * randn(15) .+ p_2_mean_log)
            else
                throw(error("The gate $(gate) is unsupported."))
            end
        end
        @assert sum(gate_probs) < 1 "The unnormalised gate probabilities $(gate_probs) sum to more than 1; change the input parameters."
        if combined && is_spam(gate)
            if gate ∉ spam_gates
                # Generate the SPAM gates
                spam_z = Gate(gate.type[1] * "Z", 0, gate.targets)
                spam_x = Gate(gate.type[1] * "X", 0, gate.targets)
                spam_y = Gate(gate.type[1] * "Y", 0, gate.targets)
                # Check the ordering is as expected
                @assert gate == spam_z &&
                        total_gates[idx + 1] == spam_x &&
                        total_gates[idx + 2] == spam_y "The SPAM gates $(spam_z), $(spam_x), and $(spam_y) are not in the expected order."
                # Set the SPAM error probabilities
                gate_probs = [exp(p_m_std_log * randn() + p_m_mean_log)]
                spam_probability = [1 - sum(gate_probs); gate_probs]
                gate_probabilities[spam_z] = spam_probability
                gate_probabilities[spam_x] = spam_probability
                gate_probabilities[spam_y] = spam_probability
                # Add the SPAM gates to the list
                append!(spam_gates, [spam_z; spam_x; spam_y])
            end
        else
            gate_probabilities[gate] = [1 - sum(gate_probs); gate_probs]
        end
    end
    # Reset the random seed
    Random.seed!()
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

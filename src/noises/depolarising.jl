"""
    DepolarisingParameters

Parameterises a depolarising Pauli noise model.

# Fields

  - `params::Dict{Symbol, Any}`: Dictionary of the noise parameters described below.
  - `noise_name::String`: Noise parameter name for saving data.

# Parameters

  - `r_1::Real`: Single-qubit gate entanglement infidelity, the sum of all 3 non-identity Pauli error probabilities.
  - `r_2::Real`: Two-qubit gate entanglement infidelity, the sum of all 15 non-identity Pauli error probabilities.
  - `r_m::Real`: Measurement error probability.
  - `r_im::Real`: Mid-circuit measurement idle entanglement infidelity.
  - `r_r::Real`: Mid-circuit reset error probability.
  - `combined::Bool`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
"""
struct DepolarisingParameters <: AbstractNoiseParameters
    params::Dict{Symbol, Any}
    noise_name::String
    # Default constructor
    function DepolarisingParameters(params::Dict{Symbol, Any}, noise_name::String)
        # Check noise parameters are present
        @assert haskey(params, :r_1) "The single-qubit gate entanglement infidelity is missing."
        @assert haskey(params, :r_2) "The two-qubit gate entanglement infidelity is missing."
        @assert haskey(params, :r_m) "The measurement entanglement infidelity is missing."
        @assert haskey(params, :r_im) "The mid-circuit measurement idle entanglement infidelity is missing."
        @assert haskey(params, :r_r) "The mid-circuit reset entanglement infidelity is missing."
        @assert haskey(params, :combined) "The combined flag is missing."
        r_1 = params[:r_1]
        r_2 = params[:r_2]
        r_m = params[:r_m]
        r_im = params[:r_im]
        r_r = params[:r_r]
        combined = params[:combined]
        # Check some conditions
        @assert (r_1 >= 0) && (r_1 <= 3 / 4) "The single-qubit gate entanglement infidelity $(r_1) is out of bounds."
        @assert (r_2 >= 0) && (r_2 <= 15 / 16) "The two-qubit gate entanglement infidelity $(r_2) is out of bounds."
        @assert (r_m >= 0) && (r_m <= 1 / 2) "The measurement entanglement infidelity $(r_m) is out of bounds."
        @assert (r_im >= 0) && (r_im <= 3 / 4) "The mid-circuit measurement idle entanglement infidelity $(r_im) is out of bounds."
        @assert (r_r >= 0) && (r_r <= 1 / 2) "The mid-circuit reset entanglement infidelity $(r_r) is out of bounds."
        @assert typeof(combined) == Bool "The combined flag $(combined) is not a Bool."
        # Return parameters with the appropriate name
        sigdigits = 3
        new_noise_name = "depolarising_$(round(r_1; sigdigits = sigdigits))_$(round(r_2; sigdigits = sigdigits))_$(round(r_m; sigdigits = sigdigits))_$(round(r_im; sigdigits = sigdigits))_$(round(r_r; sigdigits = sigdigits))_$(combined)"
        return new(params, new_noise_name)::DepolarisingParameters
    end
end

@struct_hash_equal_isequal DepolarisingParameters

"""
    get_dep_param(r_1::Real, r_2::Real, r_m::Real)
    get_dep_param(; r_1::Real, r_2::Real, r_m::Real, r_im::Real = r_m, r_r::Real = r_m)

Return a [`DepolarisingParameters`](@ref) object that parameterises a depolarising Pauli noise model.

# Arguments

  - `r_1::Real`: Single-qubit gate entanglement infidelity, the sum of all 3 non-identity Pauli error probabilities.
  - `r_2::Real`: Two-qubit gate entanglement infidelity, the sum of all 15 non-identity Pauli error probabilities.
  - `r_m::Real`: Measurement error probability.
  - `r_im::Real`: Mid-circuit measurement idle entanglement infidelity.
  - `r_r::Real`: Mid-circuit reset error probability.
  - `combined::Bool`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
"""
function get_dep_param(;
    r_1::Real,
    r_2::Real,
    r_m::Real,
    r_im::Real = r_m,
    r_r::Real = r_m,
    combined::Bool = false,
)
    params = Dict{Symbol, Any}(
        :r_1 => r_1,
        :r_2 => r_2,
        :r_m => r_m,
        :r_im => r_im,
        :r_r => r_r,
        :combined => combined,
    )
    dep_param = DepolarisingParameters(params, "depolarising")
    return dep_param::DepolarisingParameters
end
function get_dep_param(r_1::Real, r_2::Real, r_m::Real; combined::Bool = false)
    dep_param = get_dep_param(; r_1 = r_1, r_2 = r_2, r_m = r_m, combined = combined)
    return dep_param::DepolarisingParameters
end

"""
    init_gate_probabilities(total_gates::Vector{Gate}, dep_param::DepolarisingParameters)

Returns a dictionary of the Pauli error probabilities for each gate in `total_gates` generated according to the noise parameters `dep_param`, ordered lexicographically following [`get_paulis`](@ref).
"""
function init_gate_probabilities(
    total_gates::Vector{Gate},
    dep_param::DepolarisingParameters,
)
    # Extract the parameters for generating the noise
    p_1 = dep_param.params[:r_1] / 3
    p_2 = dep_param.params[:r_2] / 15
    p_m = dep_param.params[:r_m]
    p_r = dep_param.params[:r_r]
    p_im = dep_param.params[:r_im] / 3
    # Generate the noise
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate in total_gates
        if is_spam(gate)
            gate_probs = [p_m]
        elseif is_mid_reset(gate)
            gate_probs = [p_r]
        elseif is_meas_idle(gate)
            gate_probs = p_im * ones(3)
        elseif is_mid_meas(gate)
            throw(error("Mid-circuit measurement gates, such as $(gate), are unsupported."))
        else
            gate_support_size = length(gate.targets)
            if gate_support_size == 1
                gate_probs = p_1 * ones(3)
            elseif gate_support_size == 2
                gate_probs = p_2 * ones(15)
            else
                throw(error("The gate $(gate) is unsupported."))
            end
        end
        @assert sum(gate_probs) < 1 "The unnormalised gate probabilities $(gate_probs) sum to more than 1; change the input parameters."
        gate_probabilities[gate] = [1 - sum(gate_probs); gate_probs]
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

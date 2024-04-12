"""
Depolarising Pauli noise. Do all the documentation here.
"""
struct DepolarisingParameters <: AbstractNoiseParameters
    # Single-qubit gate entanglement infidelity
    # This is the sum of all 3 non-identity Pauli error probabilities
    r_1::Float64
    # Two-qubit gate entanglement infidelity
    # This is the sum of all 15 non-identity Pauli error probabilities
    r_2::Float64
    # Measurement entanglement infidelity
    # This is the measurement error probability
    r_m::Float64
    # Name of the noise for saving data
    noise_name::String
    # Default constructor
    function DepolarisingParameters(
        r_1::Float64,
        r_2::Float64,
        r_m::Float64;
        noise_name::String = "depolarising",
    )
        # Check the nosie parameters
        @assert (r_1 >= 0) && (r_1 <= 3 / 4) "The single-qubit gate entanglement infidelity $(r_1) is out of bounds."
        @assert (r_2 >= 0) && (r_2 <= 15 / 16) "The two-qubit gate entanglement infidelity $(r_2) is out of bounds."
        @assert (r_m >= 0) && (r_m <= 1 / 2) "The measurement entanglement infidelity $(r_m) is out of bounds."
        # Return parameters
        return new(r_1, r_2, r_m, noise_name)::DepolarisingParameters
    end
end

@struct_hash_equal_isequal DepolarisingParameters

function get_gate_probabilities(
    total_gates::Vector{Gate},
    noise_param::DepolarisingParameters,
)
    # Extract the parameters for generating the noise
    p_1 = noise_param.r_1 / 3
    p_2 = noise_param.r_2 / 15
    p_m = noise_param.r_m
    # Generate the noise
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate in total_gates
        if gate.type ∈ ["MZ", "MX", "MY"]
            probability = [p_m]
        elseif length(gate.targets) == 1
            probability = p_1 * ones(3)
        elseif length(gate.targets) == 2
            probability = p_2 * ones(15)
        else
            throw(error("The gate $(gate) is unsupported."))
        end
        @assert sum(probability) < 1 "The probabilities $(probability) sum to more than 1; change the input parameters."
        gate_probabilities[gate] = [1 - sum(probability); probability]
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
Log-normal Pauli noise. Do all the documentation here.
"""
struct LognormalParameters <: AbstractNoiseParameters
    # Mean of the single-qubit gate entanglement infidelity
    r_1::Float64
    # Mean of the two-qubit gate entanglement infidelity
    r_2::Float64
    # Mean of the measurement entanglement infidelity
    r_m::Float64
    # Approximate standard deviation of the logarithm of the entanglement infidelity
    total_std_log::Float64
    # Random seed
    seed::UInt64
    # Name of the noise for saving data
    noise_name::String
    # Constructor
    function LognormalParameters(
        r_1::Float64,
        r_2::Float64,
        r_m::Float64,
        total_std_log::Float64;
        seed::Union{UInt64, Nothing} = nothing,
        noise_name::String = "lognormal",
    )
        # Check the nosie parameters
        @assert (r_1 >= 0) && (r_1 <= 3 / 4) "The single-qubit gate entanglement infidelity $(r_1) is out of bounds."
        @assert (r_2 >= 0) && (r_2 <= 15 / 16) "The two-qubit gate entanglement infidelity $(r_2) is out of bounds."
        @assert (r_m >= 0) && (r_m <= 1 / 2) "The measurement entanglement infidelity $(r_m) is out of bounds."
        # Randomly set the seed if one isn't supplied
        if seed === nothing
            seed = rand(UInt64)
        end
        # Return parameters
        return new(r_1, r_2, r_m, total_std_log, seed, noise_name)::LognormalParameters
    end
end

@struct_hash_equal_isequal LognormalParameters

function get_gate_probabilities(total_gates::Vector{Gate}, noise_param::LognormalParameters)
    # Extract the parameters for generating the noise
    r_1 = noise_param.r_1
    r_2 = noise_param.r_2
    r_m = noise_param.r_m
    total_std_log = noise_param.total_std_log
    seed = noise_param.seed
    # We approximate the sum of log-normal random variables as a log-normal random variable # with the same mean and standard deviation, in order to ensure that the sum has a mean
    # that is approximately independent of the standard deviation, total_std_log
    p_1_std_log = sqrt(log(1 + 3 * (exp(total_std_log^2) - 1)))
    p_2_std_log = sqrt(log(1 + 15 * (exp(total_std_log^2) - 1)))
    p_m_std_log = total_std_log
    p_1_mean_log = log(r_1 / 3) - p_1_std_log^2 / 2
    p_2_mean_log = log(r_2 / 15) - p_2_std_log^2 / 2
    p_m_mean_log = log(r_m) - p_m_std_log^2 / 2
    # Fix the random seed
    Random.seed!(seed)
    # Generate the noise
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate in total_gates
        if gate.type ∈ ["MZ", "MX", "MY"]
            probability = [exp(p_m_std_log * randn() + p_m_mean_log)]
        elseif length(gate.targets) == 1
            probability = exp.(p_1_std_log * randn(3) .+ p_1_mean_log)
        elseif length(gate.targets) == 2
            probability = exp.(p_2_std_log * randn(15) .+ p_2_mean_log)
        else
            throw(error("The gate $(gate) is unsupported."))
        end
        @assert sum(probability) < 1 "The probabilities $(probability) sum to more than 1; change the input parameters."
        gate_probabilities[gate] = [1 - sum(probability); probability]
    end
    # Reset the random seed
    Random.seed!()
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

#
function get_gate_eigenvalues(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    total_gates::Vector{Gate},
    gate_index::Dict{Gate, Int},
    N::Int,
)
    # Generate the Walsh-Hadamard transform matrices
    W_1 = wht_matrix(1)
    W_2 = wht_matrix(2)
    # Generate the gate eigenvalues
    gate_eigenvalues = Vector{Float64}(undef, N)
    for gate in total_gates
        if gate.type ∈ ["MZ", "MX", "MY"]
            gate_eig = 1 - 2 * gate_probabilities[gate][2]
            gate_eigenvalues[gate_index[gate] + 1] = gate_eig
        else
            if length(gate.targets) == 1
                gate_eig = W_1 * gate_probabilities[gate]
            elseif length(gate.targets) == 2
                gate_eig = W_2 * gate_probabilities[gate]
            else
                throw(error("Unsupported gate $(gate)."))
            end
            gate_eigenvalues[(gate_index[gate] + 1):(gate_index[gate] + (4^length(
                gate.targets,
            ) - 1))] = gate_eig[2:end]
        end
    end
    return gate_eigenvalues::Vector{Float64}
end

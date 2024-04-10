"""
    DepolarisingNoise(d::Design, p_1::Float64, p_2::Float64, m::Float64)

Generate depolarising noise of per-Pauli strength p_1 for 1-qubit gates and per-Pauli strength p_2 for 2-qubit gates, and of strength p_m for measurements.
"""
function DepolarisingNoise(
    total_gates::Vector{Gate},
    p_1::Float64,
    p_2::Float64,
    p_m::Float64,
)
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
    LogNormalNoise(code::Code, log_p_1_mean::Float64, log_p_1_std::Float64, log_p_2_mean::Float64, log_p_2_std::Float64, log_m_mean::Float64, log_m_std::Float64; seed::Union{UInt64, Nothing} = nothing)

Generate log-normally distributed Pauli noise where the supplied means and standard deviations describe the normal distributions for individual Pauli errors for one- and two-qubit gates as well as measurements.
"""
function LogNormalNoise(
    total_gates::Vector{Gate},
    p_1_mean_log::Float64,
    p_1_std_log::Float64,
    p_2_mean_log::Float64,
    p_2_std_log::Float64,
    p_m_mean_log::Float64,
    p_m_std_log::Float64,
    seed::UInt64,
)
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
function GenerateNoise(total_gates::Vector{Gate}, noise_param::AbstractNoiseParameters)
    if typeof(noise_param) == DepolarisingParameters
        # Extract the parameters
        r_1 = noise_param.r_1
        r_2 = noise_param.r_2
        r_m = noise_param.r_m
        # Determine the appropriate noise parameters
        p_1 = r_1 / 3
        p_2 = r_2 / 15
        p_m = r_m
        # Generate the noise
        gate_probabilities = DepolarisingNoise(total_gates, p_1, p_2, p_m)
    elseif typeof(noise_param) == LogNormalParameters
        # Extract the parameters
        r_1 = noise_param.r_1
        r_2 = noise_param.r_2
        r_m = noise_param.r_m
        total_std_log = noise_param.total_std_log
        seed = noise_param.seed
        # Determine the appropriate noise parameters
        # Ensure that the entanglement infidelity of all of the gate types is approximately total_std_log
        p_1_std_log = sqrt(log(1 + 3 * (exp(total_std_log^2) - 1)))
        p_2_std_log = sqrt(log(1 + 15 * (exp(total_std_log^2) - 1)))
        p_m_std_log = total_std_log
        p_1_mean_log = log(r_1 / 3) - p_1_std_log^2 / 2
        p_2_mean_log = log(r_2 / 15) - p_2_std_log^2 / 2
        p_m_mean_log = log(r_m) - p_m_std_log^2 / 2
        # Generate the noise
        gate_probabilities = LogNormalNoise(
            total_gates,
            p_1_mean_log,
            p_1_std_log,
            p_2_mean_log,
            p_2_std_log,
            p_m_mean_log,
            p_m_std_log,
            seed,
        )
    else
        throw(error("Unsupported noise parameter type $(typeof(noise_param))."))
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

#
function GateEigenvalues(
    gate_probabilities::Dict{Gate, Vector{Float64}},
    total_gates::Vector{Gate},
    gate_index::Dict{Gate, Int},
    N::Int,
)
    # Generate the Walsh-Hadamard transform matrices
    W_1 = WHTMatrix(1)
    W_2 = WHTMatrix(2)
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

#
function Update(c::Code, noise_param::AbstractNoiseParameters)
    # Generate the noise
    gate_probabilities = GenerateNoise(c.total_gates, noise_param)
    gate_eigenvalues = GateEigenvalues(gate_probabilities, c.total_gates, c.gate_index, c.N)
    # Update the circuit
    c_update = deepcopy(c)
    @reset c_update.noise_param = noise_param
    @reset c_update.gate_probabilities = gate_probabilities
    @reset c_update.gate_eigenvalues = gate_eigenvalues
    return c_update::Code
end

#
function Update(d::Design, noise_param::AbstractNoiseParameters)
    # Generate the noise
    gate_probabilities = GenerateNoise(d.code.total_gates, noise_param)
    gate_eigenvalues =
        GateEigenvalues(gate_probabilities, d.code.total_gates, d.code.gate_index, d.code.N)
    # Update the circuit
    d_update = deepcopy(d)
    @reset d_update.code.noise_param = noise_param
    @reset d_update.code.gate_probabilities = gate_probabilities
    @reset d_update.code.gate_eigenvalues = gate_eigenvalues
    return d_update::Design
end

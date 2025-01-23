using QuantumACES, Test
# Set up parameters
p = 0.025 / 100
m = 2.0 / 100
r_1 = 3 * p
r_2 = 6 * p + 9 * p^2
r_m = m
budget_set = [10^7; 2 * 10^7]
repetitions = 3
seed = UInt(0)
z_score_cutoff_upper = 3.5
z_score_cutoff_lower = -3.5
# Example circuit parameters
struct ExampleParameters <: AbstractCircuitParameters
    params::Dict{Symbol, Any}
    circuit_name::String
    # Default constructor
    function ExampleParameters(params::Dict{Symbol, Any}, circuit_name::String)
        # Check circuit parameters are present
        @assert haskey(params, :pad_identity) "The pad identity flag is missing."
        @assert haskey(params, :layer_time_dict) "The layer time dictionary is missing."
        pad_identity = params[:pad_identity]
        layer_time_dict = params[:layer_time_dict]
        # Check some conditions
        @assert haskey(layer_time_dict, :single_qubit) "The layer time dictionary must contain the key :single_qubit."
        @assert haskey(layer_time_dict, :two_qubit) "The layer time dictionary must contain the key :two_qubit."
        @assert haskey(layer_time_dict, :meas_reset) "The layer time dictionary must contain the key :meas_reset."
        @assert layer_time_dict[:single_qubit] > 0.0 "The single-qubit layer time must be positive."
        @assert layer_time_dict[:two_qubit] > 0.0 "The two-qubit layer time must be positive."
        @assert layer_time_dict[:meas_reset] > 0.0 "The measurement and reset time must be positive."
        # Set the circuit name
        new_circuit_name = "example_circuit"
        if pad_identity != true
            circuit_name *= "_no_pad_identity"
        end
        # Return parameters
        return new(params, new_circuit_name)::ExampleParameters
    end
end
# Construct the example parameters
function get_example_param(;
    pad_identity = true,
    single_qubit_time::Real = 29,
    two_qubit_time::Real = 29,
    meas_reset_time::Real = 660,
)
    # Construct the layer time dictionary
    layer_time_dict = Dict{Symbol, Float64}(
        :single_qubit => single_qubit_time,
        :two_qubit => two_qubit_time,
        :meas_reset => meas_reset_time,
    )
    # Construct the circuit parameters
    params = Dict{Symbol, Any}(
        :pad_identity => pad_identity,
        :layer_time_dict => layer_time_dict,
    )
    # Return parameters
    example_param = ExampleParameters(params, "example_circuit")
    return example_param::ExampleParameters
end
# Construct the example circuit
function example_circuit(example_param::ExampleParameters)
    # Set up variables
    pad_identity = example_param.params[:pad_identity]
    layer_time_dict = example_param.params[:layer_time_dict]
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    # Generate the circuit
    qubit_num = 3
    circuit = [
        Layer([Gate("CZ", 0, [2; 3])], qubit_num),
        Layer([Gate("CZ", 0, [1; 2]), Gate("H", 0, [3])], qubit_num),
        Layer([Gate("H", 0, [1]), Gate("S", 0, [2]), Gate("H", 0, [3])], qubit_num),
    ]
    layer_types = [two_qubit_type, two_qubit_type, single_qubit_type]
    layer_times = get_layer_times(layer_types, layer_time_dict)
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [pad_layer(l) for l in circuit]
    end
    return (
        circuit::Vector{Layer},
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
    )
end
# Construct the circuit struct
function QuantumACES.get_circuit(
    example_param::ExampleParameters,
    noise_param::T;
    noisy_prep::Bool = false,
    noisy_meas::Bool = true,
    combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] :
                     false,
    strict::Bool = false,
) where {T <: AbstractNoiseParameters}
    # Construct the circuit
    (circuit, layer_types, layer_times) = example_circuit(example_param)
    c = get_circuit(
        circuit,
        layer_types,
        layer_times,
        noise_param;
        circuit_param = example_param,
        noisy_prep = noisy_prep,
        noisy_meas = noisy_meas,
        combined = combined,
        strict = strict,
    )
    return c::Circuit
end
# Phenomenological noise model parameters
struct PhenomenologicalParameters <: AbstractNoiseParameters
    params::Dict{Symbol, Any}
    noise_name::String
    # Default constructor
    function PhenomenologicalParameters(params::Dict{Symbol, Any}, noise_name::String)
        # Check noise parameters are present
        @assert haskey(params, :p) "The phenomenological gate error probability is missing."
        @assert haskey(params, :m) "The measurement error probability is missing."
        @assert haskey(params, :combined) "The combined flag is missing."
        p = params[:p]
        m = params[:m]
        combined = params[:combined]
        # Check some conditions
        @assert (p >= 0) && (p <= 1 / 10) "The phenomenological gate error probability $(p) is out of bounds."
        @assert (m >= 0) && (m <= 1 / 2) "The phenomenological measurement error probability $(m) is out of bounds."
        @assert typeof(combined) == Bool "The combined flag $(combined) is not a Bool."
        # Return parameters with the appropriate name
        sigdigits = 3
        new_noise_name = "phenomenological_$(round(p; sigdigits = sigdigits))_$(round(m; sigdigits = sigdigits))_$(combined)"
        return new(params, new_noise_name)::PhenomenologicalParameters
    end
end
# Construct the phenomenological noise model parameters
function get_phen_param(p::Float64, m::Float64; combined::Bool = false)
    params = Dict{Symbol, Any}(:p => p, :m => m, :combined => combined)
    phen_param = PhenomenologicalParameters(params, "phenomenological")
    return phen_param::PhenomenologicalParameters
end
# Construct the gate probabilities
function QuantumACES.init_gate_probabilities(
    total_gates::Vector{Gate},
    phen_param::PhenomenologicalParameters,
)
    # Extract the parameters for generating the noise
    p = phen_param.params[:p]
    m = phen_param.params[:m]
    im = phen_param.params[:m] / 3
    # Determine the weight of the error corresponding to each gate error probability
    one_qubit_support_size = ones(3)
    n = 2
    two_qubit_support_size = Vector{Int}()
    bit_array = BitArray(undef, 2n + 1)
    for bit_array.chunks[1] in 1:(4^n - 1)
        two_qubit_pauli = Pauli(convert(Vector{Bool}, bit_array), n)
        push!(two_qubit_support_size, length(get_support(two_qubit_pauli)))
    end
    @assert sum(two_qubit_support_size .== 1) == 6
    @assert sum(two_qubit_support_size .== 2) == 9
    @assert length(two_qubit_support_size) == 15
    # Generate the noise
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate in total_gates
        if is_spam(gate) || is_mid_meas_reset(gate)
            gate_probs = [m]
        elseif is_meas_idle(gate)
            gate_probs = im * one_qubit_support_size
        else
            gate_support_size = length(gate.targets)
            if gate_support_size == 1
                gate_probs = p .^ one_qubit_support_size
            elseif gate_support_size == 2
                gate_probs = p .^ two_qubit_support_size
            else
                throw(error("The gate $(gate) is unsupported."))
            end
        end
        @assert sum(gate_probs) < 1 "The probabilities $(gate_probs) sum to more than 1; change the input parameters."
        gate_probabilities[gate] = [1 - sum(gate_probs); gate_probs]
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end
# Test creating an example circuit and phenomenological noise model
@testset "Example circuit and phenomenological noise" begin
    # Set up noise models
    phen_param = get_phen_param(p, m)
    dep_param = get_dep_param(r_1, r_2, r_m)
    # Generate the circuit
    example_param = get_example_param()
    circuit_example = get_circuit(example_param, dep_param)
    # Generate the experimental design
    d = generate_design(circuit_example)
    display(d)
    # Update the noise to the phenomenological noise model
    d_phen = update_noise(d, phen_param)
    merit_phen = calc_merit(d_phen)
    display(merit_phen)
    # Simulate ACES experiments
    aces_data = simulate_aces(d_phen, budget_set; repetitions = repetitions, seed = seed)
    pretty_print(aces_data, merit_phen; projected = true)
    # Test that the simulations agree sufficiently with the predicted distributions
    noise_score_coll = get_noise_score(aces_data, merit_phen)
    for noise_score in noise_score_coll
        @test is_score_expected(noise_score, z_score_cutoff_lower, z_score_cutoff_upper)
    end
end

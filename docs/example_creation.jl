using QuantumACES

# Example circuit parameters
struct ExampleParameters <: AbstractCircuitParameters
    pad_identity::Bool
    layer_time_dict::Dict{Symbol, Float64}
    circuit_name::String
end

# Construct the example parameters
function get_example_param(;
    pad_identity = true,
    single_qubit_time::Float64 = 29.0,
    two_qubit_time::Float64 = 29.0,
    meas_reset_time::Float64 = 660.0,
)
    # Create the example parameters
    @assert single_qubit_time > 0.0 "The single-qubit layer time must be positive."
    @assert two_qubit_time > 0.0 "The two-qubit layer time must be positive."
    @assert meas_reset_time > 0.0 "The measurement and reset layer time must be positive."
    layer_time_dict = Dict(
        :single_qubit => single_qubit_time,
        :two_qubit => two_qubit_time,
        :meas_reset => meas_reset_time,
    )
    circuit_name = "example_circuit"
    if pad_identity != true
        circuit_name *= "_pad_identity_$(pad_identity)"
    end
    example_param = ExampleParameters(pad_identity, layer_time_dict, circuit_name)
    return example_param::ExampleParameters
end

# Construct the example circuit
function example_circuit(example_param::ExampleParameters)
    # Set up variables
    pad_identity = example_param.pad_identity
    layer_time_dict = example_param.layer_time_dict
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
        qubit_num::Int,
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
    )
end

# Construct the circuit struct
function QuantumACES.get_circuit(
    example_param::ExampleParameters,
    noise_param::T;
    add_prep::Bool = false,
    add_meas::Bool = true,
) where {T <: AbstractNoiseParameters}
    # Construct the circuit
    (circuit, qubit_num, layer_types, layer_times) = example_circuit(example_param)
    circuit_tuple = collect(1:length(circuit))
    # Prepare the circuit and generate additional parameters
    (
        labelled_circuit,
        unique_layer_indices,
        gates,
        total_gates,
        gate_index,
        N,
        gate_probabilities,
        gate_eigenvalues,
    ) = prepare_circuit(
        circuit,
        qubit_num,
        layer_types,
        layer_times,
        noise_param;
        add_prep = add_prep,
        add_meas = add_meas,
    )
    # Return the circuit
    c = Circuit(
        example_param,
        labelled_circuit,
        circuit_tuple,
        qubit_num,
        unique_layer_indices,
        layer_types,
        layer_times,
        gates,
        total_gates,
        gate_index,
        N,
        noise_param,
        gate_probabilities,
        gate_eigenvalues,
        add_prep,
        add_meas,
    )
    return c::Circuit
end

# Phenomenological noise model parameters
struct PhenomenologicalParameters <: AbstractNoiseParameters
    p::Float64
    m::Float64
    noise_name::String
end

# Construct the phenomenological noise model parameters
function get_phen_param(p::Float64, m::Float64)
    @assert (p >= 0) && (p <= 1 / 10) "The phenomenological gate error probability $(p) is out of bounds."
    @assert (m >= 0) && (m <= 1 / 2) "The phenomenological measurement error probability $(m) is out of bounds."
    noise_name = "phenomenological_$(round(p; sigdigits = 4))_$(round(m; sigdigits = 4))"
    return PhenomenologicalParameters(p, m, noise_name)::PhenomenologicalParameters
end

# Construct the gate probabilities
function QuantumACES.get_gate_probabilities(
    total_gates::Vector{Gate},
    noise_param::PhenomenologicalParameters,
)
    # Extract the parameters for generating the noise
    p = noise_param.p
    m = noise_param.m
    # Determine the weight of the error corresponding to each gate error probability
    one_qubit_support_size = ones(3)
    n = 2
    two_qubit_support_size = Vector{Int}(undef, 0)
    bit_array = BitArray(undef, 2n + 1)
    for bit_array.chunks[1] in 1:(2^(2n) - 1)
        two_qubit_pauli = Pauli(convert(Vector{Bool}, bit_array), n)
        push!(two_qubit_support_size, length(get_support(two_qubit_pauli)))
    end
    @assert sum(two_qubit_support_size .== 1) == 6
    @assert sum(two_qubit_support_size .== 2) == 9
    # Generate the noise
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate in total_gates
        if gate.type ∈ ["MZ", "MX", "MY"]
            probability = [m]
        elseif length(gate.targets) == 1
            probability = p .^ one_qubit_support_size
        elseif length(gate.targets) == 2
            probability = p .^ two_qubit_support_size
        else
            throw(error("The gate $(gate) is unsupported."))
        end
        @assert sum(probability) < 1 "The probabilities $(probability) sum to more than 1; change the input parameters."
        gate_probabilities[gate] = [1 - sum(probability); probability]
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end

# Set up noise models
p = 0.025 / 100
m = 2.0 / 100
r_1 = 3 * p
r_2 = 6 * p + 9 * p^2
r_m = m
phen_param = get_phen_param(p, m)
dep_param = get_dep_param(r_1, r_2, r_m)
# Generate the circuit
example_param = get_example_param()
circuit_example = get_circuit(example_param, dep_param)
# Optimise the experimental design
seed = UInt(0)
d = optimise_design(circuit_example; options = OptimOptions(; ls_type = :gls, seed = seed))
pretty_print(d)
# Updat the noise to the phenomenological noise model
d_phen = update_noise(d, phen_param)
merit_set_dep = calc_merit_set(d)
merit_set_phen = calc_merit_set(d_phen)
# Simulate ACES experiments
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_phen, budget_set; repetitions = repetitions, seed = seed)
pretty_print(aces_data, merit_set_phen)

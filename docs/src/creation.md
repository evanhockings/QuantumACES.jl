# Creating Circuits and Noise Models

QuantumACES makes it easy to create new circuits and noise models.
At a high level, we create new parameter types for the circuit or noise model and then create new methods for the functions [`get_circuit`](@ref) and [`init_gate_probabilities`](@ref), respectively, that take these parameter types as arguments.
Then [`get_circuit`](@ref) uses the circuit and noise model parameters to create a [`Circuit`](@ref) object.

The [`Circuit`](@ref) object enables the functionality we saw in [Example usage](@ref), with two exceptions.
First, calculating performance predictions requires the circuit to be parameterised by a `dist` parameter, typically the distance of the code underlying the syndrome extraction circuit.
Second, simulating memory experiments in Stim requires the circuit to be a syndrome extraction annotated with the appropriate information.
The [`Circuit`](@ref) object contains an `extra_fields` field which is a dictionary that can store additional parameters to enable functionality such as this.
In particular, syndrome extraction circuits must store a [`CodeParameters`](@ref)object in this dictionary, which is then used to generate the detectors in Stim for memory experiment circuits, enabling decoding.

Let us begin by creating a new circuit, following the example circuit shown in Figure 2 of [arXiv:2404.06545](https://arxiv.org/abs/2404.06545).
The first step is to create a parameter object for the circuit, which must be a subtype of [`AbstractCircuitParameters`](@ref) and contain the necessary fields `params` and `circuit_name`, using the `StructEquality` package to automatically generate hash and equality relationships for the object to enable comparisons for the resulting [`Circuit`](@ref) objects.

```julia
using QuantumACES, StructEquality
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

@struct_hash_equal_isequal ExampleParameters
```

We need a function to construct the parameter object.

```julia
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
```

And we need a function to create the circuit from the parameter object.

```julia
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
    extra_fields = Dict{Symbol, Any}()
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [pad_layer(l) for l in circuit]
    end
    return (
        circuit::Vector{Layer},
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
        extra_fields::Dict{Symbol, Any},
    )
end
```

Finally, we create a function that generates the circuit in the form of a [`Circuit`](@ref) object.
We do this by adding a method to [`get_circuit`](@ref) which uses the new parameter object as an argument, leveraging this circuit creation function.

```julia
function QuantumACES.get_circuit(
    example_param::ExampleParameters,
    noise_param::T;
    noisy_prep::Bool = false,
    noisy_meas::Bool = true,
    combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] :
                     false,
    strict::Bool = false,
) where {T <: AbstractNoiseParameters}
    (circuit, layer_types, layer_times, extra_fields) = example_circuit(example_param)
    c = get_circuit(
        circuit,
        layer_types,
        layer_times,
        noise_param;
        circuit_param = example_param,
        extra_fields = extra_fields,
        noisy_prep = noisy_prep,
        noisy_meas = noisy_meas,
        combined = combined,
        strict = strict,
    )
    return c::Circuit
end
```

Circuit creation works similarly for more complicated circuits, such as the syndrome extraction circuits generated for [`RotatedPlanarParameters`](@ref), [`UnrotatedPlanarParameters`](@ref), and [`HeavyHexParameters`](@ref).
These each have specialised methods added to [`get_circuit`](@ref) which uses these parameter structs as arguments.

If you wish to add additional information to a [`Circuit`](@ref) object, store it in the `extra_fields` field, which for example can contain a [`CodeParameters`](@ref) object that enables the construction and simulation of a memory circuit through [`get_stim_memory_circuit`](@ref) and [`simulate_memory`](@ref).

Next, we will create a phenomenological noise model where each Pauli error probability has some probability `p` of occurring, so that two-qubit errors have a probability `p^2`.
As with the circuit, we begin by creating a parameter object for the noise model, which must be a subtype of [`AbstractNoiseParameters`](@ref) and contains the noise parameters in the necessary field `params`, as well as a name in the necessary field `noise_name`.
Including a `combined` field in the noise model specifies to the circuit and experimental design whether to treat Pauli X, Y, and Z basis SPAM noise as the same.

```julia
struct PhenomenologicalParameters <: AbstractNoiseParameters
    params::Dict{Symbol, Any}
    noise_name::String
    # Default constructor
    function PhenomenologicalParameters(params::Dict{Symbol, Any}, noise_name::String)
        # Check noise parameters are present
        @assert haskey(params, :p) "The phenomenological gate error probability is missing."
        @assert haskey(params, :m) "The measurement error probability is missing."
        @assert haskey(params, :m_r) "The measurement reset error probability is missing."
        @assert haskey(params, :m_i) "The measurement idle error probability is missing."
        @assert haskey(params, :combined) "The combined flag is missing."
        p = params[:p]
        m = params[:m]
        m_r = params[:m_r]
        m_i = params[:m_i]
        combined = params[:combined]
        # Check some conditions
        @assert (p >= 0) && (p <= 1 / 10) "The phenomenological gate error probability $(p) is out of bounds."
        @assert (m >= 0) && (m <= 1 / 2) "The phenomenological measurement error probability $(m) is out of bounds."
        @assert (m_r >= 0) && (m_r <= 1 / 2) "The phenomenological measurement reset error probability $(m_r) is out of bounds."
        @assert (m_i >= 0) && (m_i <= 1 / 4) "The phenomenological measurement idle error probability $(m_i) is out of bounds."
        @assert typeof(combined) == Bool "The combined flag $(combined) is not a Bool."
        # Return parameters with the appropriate name
        sigdigits = 3
        new_noise_name = "phenomenological_$(round(p; sigdigits = sigdigits))_$(round(m; sigdigits = sigdigits))_$(round(m_r; sigdigits = sigdigits))_$(round(m_i; sigdigits = sigdigits))_$(combined)"
        return new(params, new_noise_name)::PhenomenologicalParameters
    end
end
```

We need a function to construct the parameter object.

```julia
function get_phen_param(
    p::Float64,
    m::Float64;
    m_r::Real = m,
    m_i::Real = m / 3,
    combined::Bool = false,
)
    params =
        Dict{Symbol, Any}(:p => p, :m => m, :m_r => m_r, :m_i => m_i, :combined => combined)
    phen_param = PhenomenologicalParameters(params, "phenomenological")
    return phen_param::PhenomenologicalParameters
end
```

And we need a function to create the noise model for a set of gates from the parameter object.
As with the circuit, we add a method to [`init_gate_probabilities`](@ref) which uses the new parameter object as an argument.

```julia
function QuantumACES.init_gate_probabilities(
    total_gates::Vector{Gate},
    phen_param::PhenomenologicalParameters,
)
    # Set up variables
    p = phen_param.params[:p]
    m = phen_param.params[:m]
    m_r = phen_param.params[:m_r]
    m_i = phen_param.params[:m_i]
    # Determine the weights of the Pauli errors
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
        if is_spam(gate)
            gate_probs = [m]
        elseif is_mid_meas_reset(gate)
            gate_probs = [m_r]
        elseif is_meas_idle(gate)
            gate_probs = m_i .^ one_qubit_support_size
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
        @assert sum(gate_probs) < 1
        gate_probabilities[gate] = [1 - sum(gate_probs); gate_probs]
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end
```

Noise model creation should follow the form of this example, much like the noise models [`DepolarisingParameters`](@ref) and [`LognormalParameters`](@ref).

Now we are ready to reprise the basic usage example for this example circuit.
First, parameterise depolarising and phenomenological noise models.

```julia
p = 0.025 / 100
m = 2.0 / 100
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
phen_param = get_phen_param(p, m)
dep_param = get_dep_param(r_1, r_2, r_m)
```

Then construct the circuit.

```julia
example_param = get_example_param()
circuit_example = get_circuit(example_param, dep_param)
```

The default optimisation target is the generalised least squares (GLS) estimator, as it performs best, even if it cannot be scaled up to very large numbers of qubits.
This is not an issue here because because the circuit acts on only three qubits.
QuantumACES is geared towards syndrome extraction circuits, which typically are performant when repeated an even number of times.
This is not guaranteed for the example circuit, so we set `add_circuit` to be `false`, though note this is the default behaviour.

```julia
add_circuit = false
repeat_points = 3
seed = UInt(0)
d = optimise_design(
    circuit_example;
    options = OptimOptions(;
        add_circuit = add_circuit,
        repeat_points = repeat_points,
        seed = seed,
    ),
)
merit = calc_merit(d)
```

Now create a randomised experimental design for the circuit.
This is not necessary for simulations, but is when implementing the experimental design on an actual quantum device.

```julia
min_randomisations = 50
target_shot_budget = 10^7
experiment_shots = 512
d_rand = generate_rand_design(
    d,
    min_randomisations,
    target_shot_budget,
    experiment_shots;
    seed = seed,
)
d_shot = get_design(d_rand)
merit_shot = calc_merit(d_shot)
```

Create a copy of the optimised design that associates phenomenological noise with the circuit to compare the predicted performance of the experimental design with depolarising and phenomenological noise.
In particular, we can predict the expectation and mean of the normalised root-mean-square (RMS) error between the estimated and true gate eigenvalues: either for all of the gate eigenvalues; or after marginalising eigenvalues across the Pauli orbits for each gate; or for those marginalised eigenvalues which can be estimated to relative precision, namely those not associated with state preparation and measurement (SPAM) noise.

```julia
d_phen = update_noise(d_shot, phen_param)
merit_phen = calc_merit(d_phen)
```

We can also simulate noise characterisation experiments with this experimental design and phenomenological noise, and compare the performance to predictions by computing z-scores for the normalised RMS error with respect to the predicted expectation and variance.
Note that the generalised least squares (GLS) estimator is the most performant and the focus here.
It is implemented as an iterative feasible generalised least squares (FGLS) method.

```julia
budget_set = [10^6; 10^7; 10^8]
repetitions = 10
aces_data = simulate_aces(d_phen, budget_set; repetitions = repetitions, seed = seed)
pretty_print(aces_data, merit_phen)
```

As before, note that the distribution of the normalised RMS error between the estimated and true gate eigenvalues is not quite normally distributed.
Hence the z-scores shown here, which are normalised by the predicted performance of the experimental design, will not quite be normally distributed.

Now we can examine one of the tuples in the design and the corresponding experiments used to estimate its circuit eigenvalues.

```julia
example_tuple = [2; 3; 2; 1; 1]
idx = 9
@assert d_shot.tuple_set[idx] == example_tuple
experiment_set = d_shot.experiment_ensemble[idx]
example_mappings =
    [d_shot.mapping_ensemble[idx][experiment] for experiment in experiment_set]
```

We can also examine slices from the gate eigenvalue estimator covariance matrix corresponding to a particular gate.

```julia
gls_covariance = calc_gls_covariance(d_shot)
gls_marginal_covariance = get_marginal_gate_covariance(d_shot, gls_covariance)
h_22_gate_index = d_shot.c.gate_data.gate_indices[4]
h_22_indices = h_22_gate_index.indices
h_22_marg_indices = h_22_gate_index.marg_indices
display(gls_covariance[h_22_indices, h_22_indices])
display(gls_marginal_covariance[h_22_marg_indices, h_22_marg_indices])
```

For the Hadamard gate considered here, Pauli X and Z form one orbit, and Pauli Y forms the other.
We see that the X and Z eigenvalues have large variance, but also large negative covariance, such that upon marginalisation over gate orbits, we see that the marginal orbit eigenvalue has small variance comparable to that of the Y eigenvalue.
This demonstrates that ACES implicitly performs relative precision estimation of these marginal orbit eigenvalues.

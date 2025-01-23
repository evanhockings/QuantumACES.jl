# Package Guide

## Introduction

`QuantumACES.jl` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits with averaged circuit eigenvalue sampling (ACES).
It is particularly interested in characterising the noise associated with fault-tolerant gadgets in the context of topological quantum error correcting codes, such as surface code syndrome extraction circuits.
It interfaces with [Stim](https://github.com/quantumlib/Stim) and [PyMatching](https://github.com/oscarhiggott/PyMatching) for stabiliser circuit simulations and decoding of syndrome extraction circuits, respectively, and with [Qiskit](https://github.com/Qiskit/qiskit) for implementation on quantum devices.

The methods used in this package are based on those detailed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545), and the code generating the data for that paper can be found in the `scalable_aces` folder on the [scalable_aces](https://github.com/evanhockings/QuantumACES.jl/tree/scalable_aces) branch, though the code uses an older version of the package.
These methods build on the original ACES protocol presented in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).

Typical usage of this package involves a few steps:

  - Construct the circuit and noise model for which you aim to perform an ACES noise characterisation experiment, using either provided functions or your own.
  - Optimise an experimental design, in general for a depolarising noise model with roughly the same average error rates as the noise you aim to characterise.
  - Simulate noise characterisation experiments with the optimised experimental design across a range of specified measurement budgets.

Some other things you can do with this package include:

  - Predict the performance of an experimental design of a syndrome extraction circuit as a function of the code distance.
  - Simulate memory experiments using syndrome extraction circuits of quantum error correcting codes.
  - Export an experimental design to Qiskit circuits, enabling characterisation of physical quantum devices.

## Installation and setup

To install this package, run the following command in the Julia REPL.

```
] add QuantumACES
```

CAUTION: This package uses [PythonCall](https://github.com/JuliaPy/PythonCall.jl) to call a number of Python packages.
If PythonCall and these packages are not configured correctly, associated functions will not work.
The packages attempts to load the following Python packages:

  - [Stim](https://github.com/quantumlib/Stim), installed with `pip install stim`.
  - [PyMatching](https://github.com/oscarhiggott/PyMatching), installed with `pip install pymatching`.
  - [BeliefMatching](https://github.com/oscarhiggott/BeliefMatching), installed with `pip install beliefmatching`.
  - [Qiskit](https://github.com/Qiskit/qiskit), installed with `pip install qiskit`.
  - [Aer](https://github.com/Qiskit/qiskit-aer), installed with `pip install qiskit-aer`.

By default, PythonCall creates its own Python environment, but you may wish to use an existing Python installation.

One helpful method for managing Python versions is [pyenv](https://github.com/pyenv/pyenv), or for Windows, [pyenv-win](https://github.com/pyenv-win/pyenv-win); these are analogous to [Juliaup](https://github.com/JuliaLang/juliaup) for Julia.
The following assumes you are using pyenv or pyenv-win.

On Windows, to instruct PythonCall to use the Python version set by pyenv, configure PythonCall's environment variables by adding the following to your `~/.julia/config/startup.jl` file

```julia
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
python_exe = readchomp(`cmd /C pyenv which python`)
ENV["JULIA_PYTHONCALL_EXE"] = python_exe
```

On Unix systems, shell commands are parsed directly by Julia and appear to be unaware of your PATH variable, and I am not sure how to work around this.
Therefore, you may need to manually supply `python_exe` for the Python version `<version>` as

```julia
python_exe = homedir() * "/.pyenv/versions/<version>/bin/python"
```

This package uses [SCS.jl](https://github.com/jump-dev/SCS.jl) to project probability distributions into the simplex in the Mahalanobis distance.
This projection performs best, at least at small scales, when single-threaded by setting `ENV["OMP_NUM_THREADS"]="1"` in `~/.julia/config/startup.jl`, but this causes large performance regressions in other linear algebra routines in the package, so the setting is only advised when benchmarking the time taken to perform ACES noise estimation.

## Example usage

Beware that the examples shown below can take over an hour to run.
Ensure that Julia is set up to use as many threads as your CPU can handle.

First parameterise a depolarising noise model with single-qubit gate infidelity `r_1`, two-qubit gate infidelity `r_2`, and measurement infidelity `r_m`, and a log-normal random Pauli noise model with the same gate infidelities and a standard deviation of the underlying normal distributions `total_std_log`, alongside a random seed used when generating the noise model.

```julia
using QuantumACES
r_1 = 0.05 / 100
r_2 = 0.4 / 100
r_m = 0.8 / 100
total_std_log = 0.5
seed = UInt(0)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
```

Then generate the syndrome extraction circuit for a distance `dist` (rotated) surface code.

```julia
dist = 3
rotated_param = get_rotated_param(dist)
circuit_dep = get_circuit(rotated_param, dep_param)
circuit_log = get_circuit(rotated_param, log_param)
```

Next, generate an experimental design for this circuit.

```julia
d = generate_design(circuit_dep)
display(d)
```

Alternatively, optimise an experimental design to improve its sample efficiency, configuring the optimisation with the parameters associated with [`OptimOptions`](@ref).

```julia
d = optimise_design(circuit_dep; options = OptimOptions(; seed = seed))
display(d)
```

Create a copy of the optimised design that associates log-normal random Pauli noise with the circuit, and simulate `repetitions` rounds of an ACES noise characterisation experiment across all of the supplied measurement budgets in `budget_set`, which are measurement shots normalised by the time taken to perform the experiment.

```julia
d_log = update_noise(d, log_param)
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_log, budget_set; repetitions = repetitions, seed = seed)
```

We can compare the performance to performance predictions at the largest measurement budget, although we note that the z-scores will not quite be normally distributed as the underlying distribution is not quite normal.

```julia
merit_log = calc_merit(d_log)
display(merit_log)
pretty_print(aces_data, merit_log)
noise_score_coll = get_noise_score(aces_data, merit_log)
gls_z_scores = [noise_score.gls_z_score for noise_score in noise_score_coll]
display(gls_z_scores)
```

Next, calculate the performance scaling of this design as a fuction of the code distance, for both depolarising noise and over random instances of log-normal random Pauli noise.

```julia
dist_max = 7
merit_scaling = calc_merit_scaling(d, dist_max)
ensemble_scaling = calc_ensemble_scaling(d_log, dist_max; seed = seed)
```

Next, transfer the optimised experimental design to the syndrome extraction circuit for a much larger distance surface code with log-normal random Pauli noise.
Disable full covariance matrix generation, as inverting such large covariance matrices is very computationally expensive.

```julia
dist_big = 13
rotated_param_big = get_rotated_param(dist_big)
circuit_big = get_circuit(rotated_param_big, log_param)
d_big = generate_design(
    circuit_big,
    d.tuple_set_data;
    shot_weights = d.shot_weights,
    full_covariance = false,
    diagnostics = true,
)
```

Now we can simulate a large-scale ACES noise characterisation experiment across the supplied measurement budgets.
Make sure to split projection of the gate error probabilities into the simplex across each of the gates, rather than doing all gates collectively, as the latter is very computationally expensive.

```julia
aces_data_big = simulate_aces(d_big, budget_set; seed = seed, split = true)
```

Now we can compare the performance to predictions across the measurement budgets for weighted least squares, though here the noise estimates for generalised least squares will fall back on weighted least squares as the full covariance matrix information has not been generated.
Note that we would not expect the z-scores here to actually correspond to a normal distribution as the underlying distribution is not quite normal, and there is a substantive amount of uncertainty associated with the fit and the random noise model.

```julia
scaling_fit = get_scaling_fit(merit_scaling)
ensemble_fit = get_ensemble_fit(ensemble_scaling; precision = 1e-1)
wls_pred_expectation = ensemble_fit.wls_expectation_model(dist_big)
wls_pred_variance = ensemble_fit.wls_variance_model(dist_big)
wls_z_scores_big = [
    (noise_error.wls_nrmse .- wls_pred_expectation) / sqrt(wls_pred_variance) for
    noise_error in aces_data_big.noise_error_coll[1, :]
]
display(wls_z_scores_big)
wls_pred_relative_expectation = ensemble_fit.wls_relative_expectation_model(dist_big)
wls_pred_relative_variance = ensemble_fit.wls_relative_variance_model(dist_big)
wls_relative_z_scores_big = [
    (noise_error.wls_relative_nrmse .- wls_pred_relative_expectation) /
    sqrt(wls_pred_relative_variance) for
    noise_error in aces_data_big.noise_error_coll[1, :]
]
display(wls_relative_z_scores_big)
```

To implement this experimental design on an actual quantum device, we need to first randomly compile the experiments by creating a randomised design.

```julia
min_randomisations = 40
target_shot_budget = 10^7
experiment_shots = 512
d_rand = generate_rand_design(
    d_log,
    min_randomisations,
    target_shot_budget,
    experiment_shots;
    seed = seed,
)
```

We can calculate the figure of merit of the randomised design.

```julia
d_shot = get_design(d_rand)
merit_shot = calc_merit(d_shot)
display(d_shot)
display(merit_shot)
```

Then we can simultaneously generate Stim and Qiskit circuits corresponding to the randomised design, while mapping the qubits onto a particular part of the device.
The Qiskit circuit will act on `qiskit_qubit_num` qubits and the qubits will be mapped onto Qiskit qubits by `qiskit_qubit_map`, noting that Qiskit indexes qubits from 0.

```julia
qiskit_qubit_num = 17
qiskit_qubit_map = collect(0:(qiskit_qubit_num - 1))
(stim_ensemble, qiskit_ensemble) =
    get_stim_qiskit_ensemble(d_rand, qiskit_qubit_num, qiskit_qubit_map)
```

Now we can simulate the randomised experimental design with Stim.

```julia
simulate_stim_ensemble(d_rand, stim_ensemble, experiment_shots; seed = seed)
rand_noise_est = estimate_stim_ensemble(d_rand, experiment_shots; simulation_seed = seed)
rand_noise_error = get_noise_error(d_rand, rand_noise_est)
rand_noise_score = get_noise_score(rand_noise_error, merit_shot)
```

Finally, we can examine decoding a memory experiment using this syndrome extraction circuit, comparing decoding with the true gate probabilities, the estimated gate probabilities, and depolarising noise.

```julia
rounds = dist
shots = 10^6
decoder_gate_probabilities = [
    circuit_log.gate_probabilities,
    rand_noise_est.gls_gate_probabilities,
    circuit_dep.gate_probabilities,
]
memory_data = simulate_memory(
    circuit_log,
    rounds,
    shots;
    seed = seed,
    decoder_gate_probabilities = decoder_gate_probabilities,
)
display(memory_data)
```

We can also check the Z and X distances of the code, which conventionally in this package are the vertical and horizontal distances, respectively.

```julia
(z_dist, x_dist) = calc_memory_distances(circuit_dep)
```

## Circuit and noise model creation

This package also supports creating your own circuits and noise models.

Let us begin by creating a new circuit, following the example circuit shown in Figure 2 of [arXiv:2404.06545](https://arxiv.org/abs/2404.06545).
The first step is to create a parameter struct for the circuit, which must be a subtype of [`AbstractCircuitParameters`](@ref) and contain the necessary fields `params` and `circuit_name`, using the `StructEquality` package to automatically generate hash and equality relationships for the struct to enable comparisons.

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

We need a function to construct the parameter struct.

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

And we need a function to create the circuit from the parameter struct.

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
```

Finally, we create a function that generates the circuit in the form of a [`Circuit`](@ref) object.
We do this by adding a method to [`get_circuit`](@ref) which uses the new parameter struct as an argument.

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
```

Circuit creation works similarly for more complicated circuits, such as the syndrome extraction circuits generated for [`RotatedPlanarParameters`](@ref), [`UnrotatedPlanarParameters`](@ref), and [`HeavyHexParameters`](@ref).
These each have specialised methods added to [`get_circuit`](@ref) which uses these parameter structs as arguments.
If you wish to add additional information to a [`Circuit`](@ref) object, store it in the `extra_fields` field, which for example can contain a [`CodeParameters`](@ref) object that enables the construction and simulation of a memory circuit through [`get_stim_memory_circuit`](@ref) and [`simulate_memory`](@ref).

Next, we will create a phenomenological noise model where each Pauli error probability has some probability `p` of occurring, so that two-qubit errors have a probability `p^2`.
As with the circuit, we begin by creating a parameter struct for the noise model, which must be a subtype of [`AbstractNoiseParameters`](@ref) and contains the noise parameters in the necessary field `params`, as well as a name in the necessary field `noise_name`.
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
```

We need a function to construct the parameter struct.

```julia
function get_phen_param(p::Float64, m::Float64; combined::Bool = false)
    params = Dict{Symbol, Any}(:p => p, :m => m, :combined => combined)
    phen_param = PhenomenologicalParameters(params, "phenomenological")
    return phen_param::PhenomenologicalParameters
end
```

And we need a function to create the noise model for a set of gates from the parameter struct.
As with the circuit, we add a method to [`init_gate_probabilities`](@ref) which uses the new parameter struct as an argument.

```julia
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
This package is geared towards syndrome extraction circuits, which typically are performant when repeated an even number of times; this is not guaranteed for the example circuit, so we set `add_circuit = false`.

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
display(d)
display(merit)
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
display(d_shot)
display(merit_shot)
```

Create a copy of the optimised design that associates phenomenological noise with the circuit to compare the predicted performance of the experimental design with depolarising and phenomenological noise.
In particular, we can predict the expectation and mean of the normalised root-mean-square (RMS) error between the estimated and true gate eigenvalues: either for all of the gate eigenvalues; or after marginalising eigenvalues across the Pauli orbits for each gate; or for those marginalised eigenvalues which can be estimated to relative precision, namely those not associated with state preparation and measurement (SPAM) noise.

```julia
d_phen = update_noise(d_shot, phen_param)
merit_phen = calc_merit(d_phen)
display(merit_phen)
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
display(example_mappings)
```

We can also examine slices from the gate eigenvalue estimator covariance matrix corresponding to a particular gate.

```julia
gls_covariance = calc_gls_covariance(d_shot)
gls_marginal_covariance = get_marginal_gate_covariance(d_shot, gls_covariance)
h_22_gate_index = d_shot.c.gate_data.gate_indices[4]
h_22_indices = h_22_gate_index.indices
display(gls_covariance[h_22_indices, h_22_indices])
h_22_marg_indices = h_22_gate_index.marg_indices
display(gls_marginal_covariance[h_22_marg_indices, h_22_marg_indices])
```

For the Hadamard gate considered here, Pauli X and Z form one orbit, and Pauli Y forms the other.
We see that the X and Z eigenvalues have large variance, but also large negative covariance, such that upon marginalisation over gate orbits, we see that the marginal orbit eigenvalue has small variance comparable to that of the Y eigenvalue.
This demonstrates that ACES implicitly performs relative precision estimation of these marginal orbit eigenvalues.

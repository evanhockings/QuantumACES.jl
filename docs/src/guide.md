# Package Guide

## Introduction

`QuantumACES` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits with averaged circuit eigenvalue sampling (ACES).
It focuses on the context of quantum error correction and fault-tolerant circuits and, in particular, on the syndrome extraction circuits of topological quantum error correcting codes.
It interfaces with [Stim](https://github.com/quantumlib/Stim) for stabiliser circuit simulation, [PyMatching](https://github.com/oscarhiggott/PyMatching) and [BeliefMatching](https://github.com/oscarhiggott/BeliefMatching) for decoding, and [Qiskit](https://github.com/Qiskit/qiskit) for implementation on quantum devices.

Typical usage of QuantumACES involves first doing the following:

  - Construct the circuit and the noise model that you aim to characterise, either using existing functions or your own.
  - Optimise an ACES experimental design for noise characterisation of a small-scale instance of the circuit, typically according to a deterministic noise model, such as depolarising noise, with roughly the same average error rates as the noise you aim to characterise.

This subsequently enables:

  - Transferring the optimised experimental design to larger-scale instances of the circuit, including with different noise models.
  - Simulate noise characterisation experiments with ACES experimental designs, including at large scales, using Stim.
  - Calculating performance predictions for experimental designs at small scales and fitting the performance predictions, in particular for syndrome extraction circuits as a function of the distance of the underlying code, to predict performance at large scales.
  - Simulating memory experiments for syndrome extraction circuits using Stim, and then decoding with PyMatching or BeliefMatching with decoder priors informed by a range of noise models, including ACES noise estimates.
  - Creating Pauli frame randomised ACES experimental designs, exporting them to Qiskit circuits, and processing the results, enabling implementation on quantum devices.

The methods used in this package are based on [arXiv:2404.06545](https://arxiv.org/abs/2404.06545) and [arXiv:2502.21044](https://arxiv.org/abs/2502.21044), and they build on the original ACES protocol introduced in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).
The code for [arXiv:2404.06545](https://arxiv.org/abs/2404.06545) can be found in the `scalable_aces` folder on the [scalable_aces](https://github.com/evanhockings/QuantumACES.jl/tree/scalable_aces) branch.
The code for [arXiv:2502.21044](https://arxiv.org/abs/2502.21044) can be found in the `aces_decoding` folder on the [aces_decoding](https://github.com/evanhockings/QuantumACES.jl/tree/aces_decoding) branch.

If you find this package helpful for your research, please cite it using the supplied `CITATION.cff` file, and consider citing the associated papers if appropriate.
If you wish to contribute to this package, please refer to the `CONTRIBUTING.md` file.

## Installation and setup

To install this package, run the following command in the Julia REPL.

```
] add QuantumACES
```

BEWARE: This package uses [PythonCall](https://github.com/JuliaPy/PythonCall.jl) to call a number of Python packages.
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

## Package usage

Beware that the examples shown below can take a long time to run.
Ensure that Julia is set up to use as many threads as your CPU can handle.

First parameterise a depolarising noise model with single-qubit gate infidelity `r_1`, two-qubit gate infidelity `r_2`, and measurement infidelity `r_m`, and a log-normal random Pauli noise model with the same gate infidelities and a standard deviation of the underlying normal distributions `total_std_log`, specifying the seed `seed` for reproducibility.

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

Similarly, we create circuit parameters for the syndrome extraction circuit of a distance `dist` (rotated) surface code.

```julia
dist = 3
rotated_param = get_rotated_param(dist)
```

Next, we create versions of this circuit with both noise models.
```julia
circuit_dep = get_circuit(rotated_param, dep_param)
circuit_log = get_circuit(rotated_param, log_param)
```

Now we can generate an experimental design for this circuit.

```julia
d = generate_design(circuit_dep)
```

Alternatively, we can optimise an experimental design to improve its sample efficiency, configuring the optimisation with the parameters associated with [`OptimOptions`](@ref).

```julia
d = optimise_design(circuit_dep; options = OptimOptions(; seed = seed))
```

This experimental design can be transferred to the circuit with the log-normal Pauli noise model.

```julia
d_log = generate_design(circuit_log, d)
```

If we only wish to update the noise model, however, we can do this more efficiently.

```julia
d_log = update_noise(d, circuit_log)
```

Now we can simulate `repetitions` instances of ACES noise characterisation across a range of measurement budgets `budget_set`, which are measurement shots normalised by the time taken to perform the experiment.

```julia
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_log, budget_set; repetitions = repetitions, seed = seed)
```

We can compare the performance to predictions, although we note that the z-scores will not quite be normally distributed as the underlying distribution is not quite normal.

```julia
merit_log = calc_merit(d_log)
pretty_print(aces_data, merit_log)
```

We can also simulate ACES noise characterisation at scale.
First create a new design at a large code distance `dist_big`.
Setting `full_covariance` to be false means only the diagonal circuit eigenvalue estimator covariance matrix is generated, which saves a substantial amount of time.
It also prevents the design from attempting to perform generalised least squares (GLS) with the full covariance matrix, which can consume large amounts of memory at large scales, restricting the design to weighted least squares (WLS).

```julia
dist_big = 13
rotated_param_big = get_rotated_param(dist_big)
circuit_big = get_circuit(rotated_param_big, dep_param)
circuit_big_log = get_circuit(rotated_param_big, log_param)
d_big = generate_design(circuit_big_log, d; full_covariance = false, diagnostics = true)
```

Now simulate this new design, setting `split` to be `true` to avoid memory issues by splitting projection of the gate error probabilities into the simplex across each of the gates, rather than doing all gates collectively.

```julia
aces_data_big = simulate_aces(d_big, budget_set; seed = seed, split = true)
```

It is expensive to directly calculate the performance of the experimental design at this scale.
Instead, we calculate the performance scaling of the experimental design at small code distances and then extrapolate.
We can do this for depolarising noise, and for an average over instances of log-normal Pauli noise, calculating up to `dist_max`, and then extracting fits.

```julia
dist_max = 7
merit_scaling = calc_merit_scaling(d, dist_max)
scaling_fit = get_scaling_fit(merit_scaling)
ensemble_scaling = calc_ensemble_scaling(d_log, dist_max; seed = seed)
ensemble_fit = get_ensemble_fit(ensemble_scaling)
```

This allows us to predict expectations and variances, and compare them to the true values, for both the ordinary figure of merit and the relative precision figure of merit.
These are not exactly z-scores in particular because the simulation was for a single instance of log-normal Pauli noise, whereas the predictions are averaged over instances.
In practice, performance appears to be self-averaging so prediction works well at scale.

```julia
wls_pred_expectation = ensemble_fit.wls_expectation_model(dist_big)
wls_pred_variance = ensemble_fit.wls_variance_model(dist_big)
wls_scores_big = [
    (noise_error.wls_nrmse .- wls_pred_expectation) / sqrt(wls_pred_variance) for
    noise_error in aces_data_big.noise_error_coll[1, :]
]
wls_pred_relative_expectation = ensemble_fit.wls_relative_expectation_model(dist_big)
wls_pred_relative_variance = ensemble_fit.wls_relative_variance_model(dist_big)
wls_relative_scores_big = [
    (noise_error.wls_relative_nrmse .- wls_pred_relative_expectation) /
    sqrt(wls_pred_relative_variance) for
    noise_error in aces_data_big.noise_error_coll[1, :]
]
```

We can now use Stim to simulate a memory experiment with `big_rounds` rounds, sampling `big_shots` shots.
We inform the decoder, PyMatching by default, with a range of noise models including our noise estimates.

```julia
big_rounds = dist_big
big_shots = 5 * 10^6
decoder_gate_probabilities = [
    circuit_big_log.gate_probabilities
    circuit_big.gate_probabilities
    [noise_est.wls_gate_probabilities 
    for noise_est in aces_data_big.noise_est_coll[1, :]
    ]
]
decoder_labels = [
    "True"
    "Depolarising"
    ["ACES S=$(budget)" for budget in budget_set]
]
big_memory_data = simulate_memory(circuit_big_log, big_rounds, big_shots;
    seed = seed,
    decoder_gate_probabilities = decoder_gate_probabilities,
    decoder_labels = decoder_labels,
    diagnostics = true,
)
big_memory_summary = get_memory_summary(big_memory_data)
```

To implement this experimental design on an actual quantum device, we need to first construct a Pauli frame randomised version of the experimental design and generate corresponding Qiskit circuits.
We specify a minimum number of randomisations `min_randomisations` and a target shot budget `target_shot_budget`, and `experiment_shots` shots per randomised experiment.

```julia
min_randomisations = 64
target_shot_budget = 5 * 10^6
experiment_shots = 64
d_rand = generate_rand_design(
    d_log,
    min_randomisations,
    target_shot_budget,
    experiment_shots;
    seed = seed,
)
```

This modifies the shot weights, so we can calculate the merit of this new design.

```julia
d_shot = get_design(d_rand)
merit_shot = calc_merit(d_shot)
```

Now we simultaneously generate ensembles of Stim and Qiskit circuits that implement this experimental design.
The Qiskit circuits act on `qiskit_qubit_num` qubits and `qiskit_qubit_map` maps QuantumACES qubit indices to Qiskit qubit indices, noting that Julia indexes from 1 whereas Python indexes from 0.

```julia
qiskit_qubit_num = 17
qiskit_qubit_map = collect(0:(qiskit_qubit_num - 1))
(stim_ensemble, qiskit_ensemble) =
    get_stim_qiskit_ensemble(d_rand, qiskit_qubit_num, qiskit_qubit_map)
```

We only simulate in Stim as the Qiskit stabiliser circuit simulator is much slower.

```julia
simulate_stim_ensemble(d_rand, stim_ensemble, experiment_shots; seed = seed)
rand_noise_est = estimate_stim_ensemble(d_rand, experiment_shots; simulation_seed = seed)
rand_noise_error = get_noise_error(d_rand, rand_noise_est)
rand_noise_score = get_noise_score(rand_noise_error, merit_shot)
```

Suppose we then run the Qiskit circuits on a quantum device.
The results must be stored in an appropriate folder to be processed.
Given a prefix `backend`, which typically describes the device on which the circuits are run, the results must be stored relative to the current directory in a folder whose name is given by `qiskit_results_folder`.

```julia
backend = "backend"
d_rand_filename = rand_design_filename(d_rand)
@assert d_rand_filename[(end - 4):end] == ".jld2"
qiskit_results_folder = "data/$(backend)_$(d_rand_filename[1:(end - 5)])"
```

The ensemble `qiskit_ensemble` is a vector containing vectors of Qiskit circuits, each of which comprise a job.
Each job should be stored as a pickle in `qiskit_results_folder` with prefix `prefix`, followed by an underscore and the job index, starting from 1.

```julia
prefix = "job"
example_job_1_filename = "$(qiskit_results_folder)/$(prefix)_1.pickle"
```

Then it is simple to process the data and estimate the noise.

```julia
process_qiskit_ensemble(
    d_rand,
    qiskit_qubit_num,
    qiskit_qubit_map,
    experiment_shots;
    backend = backend,
    prefix = prefix,
)
noise_est =
    estimate_qiskit_ensemble(d_rand, qiskit_qubit_map, experiment_shots; backend = backend)
```

Finally, we can analyse the consistency of this noise estimate with our noise model.

```julia
model_violation = get_model_violation(d_shot, noise_est)
```

This quantity is a z-score which is approximately normally distributed if the circuit-level Pauli noise model is upheld, as it is in simulation.
Note the model violation score is substantially larger when calculated for the projected noise estimates, that is, the noise estimates after projecting the Pauli error probabilities into the probability simplex even in simulation.
By default, then, this function calculates the model violation for the unprojected noise estimates.

We can also create a version of the experimental design corresponding to a combined noise model for which Pauli ``X``, ``Z``, and ``Y`` basis SPAM noise are combined into a single parameter for each qubit, so that for ``n`` qubits we have ``n`` SPAM noise parameters.
Previously, we considered ``3n`` SPAM noise parameters for each Pauli basis, which we will call the ordinary noise model.
Then we can estimate the noise with the combined noise model and calculate its model violation.

```julia
d_comb = get_combined_design(d_shot)
comb_noise_est = estimate_gate_noise(d_comb, noise_est)
comb_model_violation = get_model_violation(d_comb, comb_noise_est)
```

We might want to perform model selection with the Akaike information criterion (AIC) or the Bayesian information criterion (BIC), which are straightforward to calculate.

```julia
aic = get_aic(d_shot, noise_est)
bic = get_bic(d_shot, noise_est)
comb_aic = get_aic(d_comb, comb_noise_est)
comb_bic = get_bic(d_comb, comb_noise_est)
```

The preferred model is the one that minimises the AIC or BIC, depending on the metric of choice.
In practice, the combined noise model tends not to be formally preferred but is nevertheless more parsimonious as the SPAM noise estimates in the ordinary noise model differ across Pauli bases more than can reasonably be expected.
Therefore we tend not to use this model selection procedure.

## Circuit and noise model creation

`QuantumACES` makes it easy to create new circuits and noise models.
At a high level, we create new parameter types for the circuit or noise model and then create new methods for the functions [`get_circuit`](@ref) and [`init_gate_probabilities`](@ref), respectively, that take these parameter types as arguments.
Then [`get_circuit`](@ref) uses the circuit and noise model parameters to create a [`Circuit`](@ref) object.

The [`Circuit`](@ref) object enables the functionality we saw in [Package usage](@ref), with two exceptions.
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
`QuantumACES` is geared towards syndrome extraction circuits, which typically are performant when repeated an even number of times.
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

# Package Guide

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

## Example usage

Beware that the examples shown below can take a long time to run.
Ensure that Julia is set up to use as many threads as your CPU can handle.
For example usage that is faster to run, see [Creating Circuits and Noise Models](@ref).

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

There are a number of options that allow you to reduce the optimisation time.
For example, we can disable cyclic coordinate descent optimisation of the circuit depth of repeated tuples in the design by setting `max_cycles = 0`.
We can also allow the greedy search over ordinary tuples to terminate once they are left unchanged by single excursion in the search by setting `excursions_unchanged = 1`.

```julia
d = optimise_design(
    circuit_dep;
    options = OptimOptions(; max_cycles = 0, excursions_unchanged = 1, seed = seed),
)
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

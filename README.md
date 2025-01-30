# QuantumACES.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://evanhockings.github.io/QuantumACES.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://evanhockings.github.io/QuantumACES.jl/dev/)
[![Build Status](https://github.com/evanhockings/QuantumACES.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/evanhockings/QuantumACES.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/evanhockings/QuantumACES.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/evanhockings/QuantumACES.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Unitary Foundation](https://img.shields.io/badge/Supported%20By-Unitary%20Foundation-FFFF00.svg)](https://unitary.foundation)

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

## Example usage

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

By default, PythonCall creates its own Python environment, but you may wish to [configure](https://juliapy.github.io/PythonCall.jl/stable/pythoncall/#pythoncall-config) it to use an existing Python installation.

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

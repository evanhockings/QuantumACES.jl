# QuantumACES.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://evanhockings.github.io/QuantumACES.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://evanhockings.github.io/QuantumACES.jl/dev/)
[![Build Status](https://github.com/evanhockings/QuantumACES.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/evanhockings/QuantumACES.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/evanhockings/QuantumACES.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/evanhockings/QuantumACES.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Unitary Foundation](https://img.shields.io/badge/Supported%20By-Unitary%20Foundation-FFFF00.svg)](https://unitary.foundation)

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

## Example usage

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

Alternatively, we can optimise an experimental design to improve its sample efficiency, configuring the optimisation with the parameters associated with `OptimOptions`.

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

## Attribution

The methods used in this package are based on [arXiv:2404.06545](https://arxiv.org/abs/2404.06545) and [arXiv:2502.21044](https://arxiv.org/abs/2502.21044), and they build on the original ACES protocol introduced in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).

The code for [arXiv:2404.06545](https://arxiv.org/abs/2404.06545) can be found in the `scalable_aces` folder on the [scalable_aces](https://github.com/evanhockings/QuantumACES.jl/tree/scalable_aces) branch.

The code for [arXiv:2502.21044](https://arxiv.org/abs/2502.21044) can be found in the `aces_decoding` folder on the [aces_decoding](https://github.com/evanhockings/QuantumACES.jl/tree/aces_decoding) branch.

If you find this package helpful for your research, please cite it using the supplied `CITATION.cff` file, and consider citing the associated papers if appropriate.
If you wish to contribute to this package, please refer to the `CONTRIBUTING.md` file.

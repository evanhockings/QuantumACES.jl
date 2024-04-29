# QuantumACES.jl

{::comment}[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://evanhockings.github.io/QuantumACES.jl/stable/){:/comment}
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://evanhockings.github.io/QuantumACES.jl/dev/)
[![Build Status](https://github.com/evanhockings/QuantumACES.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/evanhockings/QuantumACES.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/evanhockings/QuantumACES.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/evanhockings/QuantumACES.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

`QuantumACES.jl` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits with averaged circuit eigenvalue sampling (ACES).
It is particularly interested in characterising the noise associated with fault-tolerant gadgets in the context of topological quantum error correcting codes, such as surface code syndrome extraction circuits.

The methods used in this package are detailed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545).{::comment}, and the code generating the data for this paper can be found in the `scalable_aces` folder on the [scalable_aces](https://github.com/evanhockings/QuantumACES.jl/tree/scalable_aces) branch.{:/comment}
These methods build on the original ACES protocol presented in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).
This package relies on [Stim](https://github.com/quantumlib/Stim) for stabiliser circuit simulations.

## Example usage

First parameterise a depolarising noise model with single-qubit gate infidelity `r_1`, two-qubit gate infidelity `r_2`, and measurement infidelity `r_m`, and a log-normal random Pauli noise model with the same gate infidelities and a standard deviation of the underlying normal distributions `total_std_log`, alongside a random seed used when generating the noise model.

```julia
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
```

Then generate the syndrome extraction circuit for a distance `dist` (rotated) surface code.

```julia
dist = 3
rotated_param = get_rotated_param(dist)
rotated_planar = get_circuit(rotated_param, dep_param)
```

Optimise an experimental design for these parameters.

```julia
d = optimise_design(rotated_planar, options = OptimOptions(; ls_type = :wls, seed = seed))
```

Create a copy of the optimised design that associates log-normal random Pauli noise with the circuit, and simulate `repetitions` rounds of an ACES noise characterisation experiment across all of the supplied measurement budgets in `budget_set`, which are measurement shots normalised by the time taken to perform the experiment.

```julia
d_log = update_noise(d, log_param)
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_log, budget_set; repetitions = repetitions)
```

We can compare the performance to performance predictions at the largest measurement budget, although we note that the z-scores will not be normally distributed as the underlying distribution is not quite normal.

```julia
wls_merit_log = calc_wls_merit(d_log)
fgls_z_scores =
    (aces_data.fgls_gate_norm_coll[:, 3] .- wls_merit_log.expectation) /
    sqrt(wls_merit_log.variance)
```

Next, calculate the performance scaling of this design as a fuction of the code distance up to some large distance `dist_max` with the weighted least squares (WLS) estimator, for both depolarising and log-normal random Pauli noise.

```julia
dist_max = 9
dep_planar_scaling = calc_depolarising_planar_scaling(d, dist_max; ls_type = :wls)
log_planar_scaling = calc_lognormal_planar_scaling(d_log, dist_max; ls_type = :wls, seed = seed)
```

Next, transfer the optimised experimental design to the syndrome extraction circuit for a distance `dist_big` surface code with log-normal random Pauli noise.

```julia
dist_big = 13
rotated_param_big = get_rotated_param(dist_big)
rotated_planar_big = get_circuit(rotated_param_big, log_param)
d_big = generate_design(rotated_planar_big, d.tuple_set_data)
```

We can simulate a large-scale ACES noise characterisation experiment across the supplied measurement budgets.

```julia
budget_set_big = [10^6; 10^7; 10^8; 10^9]
aces_data_big = simulate_aces(d_big, budget_set_big)
```

Finally, we compare the performance to predictions across the measurement budgets, although note that we would not expect the z-scores here to actually correspond to a normal distribution as the underlying distribution is not quite normal, and there is a substantive amount of uncertainty associated with the fit.

```julia
pred_expectation = log_planar_scaling.expectation_fit(dist_big)
pred_variance = log_planar_scaling.variance_fit(dist_big)
fgls_z_scores_big =
    (aces_data_big.fgls_gate_norm_coll[1, :] .- pred_expectation) / sqrt(pred_variance)
```

## Installation and setup

This is not currently a registered package, so to add it you can run

```
julia> # press ] to enter the Pkg REPL

pkg> add https://github.com/evanhockings/QuantumACES.jl
```

This package relies on the Python package [Stim](https://github.com/quantumlib/Stim) to perform stabiliser simulations.
It calls stim with [PythonCall](https://github.com/JuliaPy/PythonCall.jl), which can be a little tricky to set up.
One helpful method for managing Python versions is [pyenv](https://github.com/pyenv/pyenv), or for Windows, [pyenv-win](https://github.com/pyenv-win/pyenv-win), which is analogous to [Juliaup](https://github.com/JuliaLang/juliaup) for Julia.

On Windows, to instruct PythonCall to use the Python version set by pyenv, configure PythonCall's environment variables by adding the following to your `~/.julia/config/startup.jl` file

```
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
python_exe = readchomp(`cmd /C pyenv which python`)
ENV["JULIA_PYTHONCALL_EXE"] = python_exe
```

On Unix systems, shell commands are parsed directly by Julia and appear to be unaware of your PATH variable.
I am not sure how to fix this, so you may need to manually supply `python_exe` for the Python version `<version>` as

```
python_exe = homedir() * "/.pyenv/versions/<version>/bin/python"
```

Then ensure Stim is installed by running `pip install stim` in your terminal.

# ACES.jl

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

`ACES.jl` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits with averaged circuit eigenvalue sampling (ACES).
It is particularly interested in characterising the noise associated with fault-tolerant gadgets in the context of topological quantum error-correcting codes, such as syndrome extraction circuits.

The methods used in this package are detailed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545), building on the original ACES method laid out in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).
This package relies on [Stim](https://github.com/quantumlib/Stim) for stabiliser circuit simulations.

## Example usage

Parameterise a depolarising noise model with single-qubit gate infidelity `r_1`, two-qubit gate infidelity `r_2`, and measurement infidelity `r_m`.

```
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
dep_param = get_dep_param(r_1, r_2, r_m)
```

Then generate the syndrome extraction circuit for a distance `dist` (rotated) surface code.

```
dist = 3
rotated_param = get_rotated_param(dist)
rotated_planar = get_circuit(rotated_param, dep_param)
```

Optimise an experimental design for these parameters.

```
d = optimise_design(rotated_planar)
```

Calculate the performance scaling of this design as a fuction of the code distance up to some large distance `dist_max` with the weighted least squares (WLS) estimator.

```
dist_max = 9
dep_planar_scaling = calc_depolarising_planar_scaling(d, dist_max; ls_type = :wls)
```

Next, transfer the optimised experimental design to the syndrome extraction circuit for a distance `dist_max` surface code.

```
rotated_param_big = get_rotated_param(dist_max)
rotated_planar_big = get_circuit(rotated_param_big, dep_param)
d_big = generate_design(rotated_planar_big, d.tuple_set_data)
```

Finally, simulate an ACES noise characterisation experiment across all of the supplied measurement budgets in `budget_set`, which are measurement shots normalised by the time taken to perform the experiment.

```
budget_set = [10^6; 10^7; 10^8]
aces_data_big = simulate_aces(d_big, budget_set)
```

## Installation and setup

This is not currently a registered package, so to add it you can run

```
julia> # press ] to enter the Pkg REPL

pkg> add https://github.com/evanhockings/ACES.jl
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

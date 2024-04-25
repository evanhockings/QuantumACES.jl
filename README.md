# ACES.jl

`ACES.jl` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits.
It is particularly interested in characterising the noise associated with fault-tolerant gadgets in the context of topological quantum error-correcting codes, such as syndrome extraction circuits.

The methods used in this package are detailed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545), building on the original method laid out in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).
This package relies on [Stim](https://github.com/quantumlib/Stim) for stabiliser circuit simulations.

## Example usage

Parameterise a depolarising noise model with

```
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
dep_param = get_dep_param(r_1, r_2, r_m)
```

Generate the syndrome extraction circuit of a rotated planar code with

```
dist = 3
rotated_param = get_rotated_param(dist)
rotated_planar = get_circuit(rotated_param, dep_param)
```

Now optimise an experimental design for these parameters using the supplied least squares estimator

```
d = optimise_design(rotated_planar)
```

Use the optimised design to create a design for a larger code circuit

```
dist_big = 9
rotated_param_big = get_rotated_param(dist_big)
rotated_planar_big = get_circuit(rotated_param_big, dep_param)
d_big = generate_design(rotated_planar_big, d.tuple_set_data)
```

Finally, simulate ACES noise characterisation using the specified number of measurement shots with 

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
It uses [PythonCall](https://github.com/JuliaPy/PythonCall.jl) to interface with Stim, and this can be a little tricky to set up.
One helpful method for managing Python versions is [pyenv](https://github.com/pyenv/pyenv), or for Windows, [pyenv-win](https://github.com/pyenv-win/pyenv-win), which is analogous to [Juliaup](https://github.com/JuliaLang/juliaup) in Julia.

On Windows, to instruct PythonCall to use the Python version set by pyenv on Windows, configure PythonCall's environment variables by adding the following to your `~/.julia/config/startup.jl` file

```
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
python_exe = readchomp(`cmd /C pyenv which python`)
ENV["JULIA_PYTHONCALL_EXE"] = python_exe
```

On Unix systems, however, shell commands are parsed directly by Julia and are unaware of your PATH variable.
I am not sure how to fix this, so you may need to manually supply `python_exe` for the Python version `<version>` as

```
python_exe = homedir() * "/.pyenv/versions/<version>/bin/python"
```

Then ensure Stim is installed by running `pip install stim`.

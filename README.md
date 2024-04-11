# AveragedCircuitEigenvalueSampling.jl

**WARNING**: This package is currently in a prerelease state.
The API is NOT STABLE or appropriately documented.
The next release will change the entire API.
Please wait for it before building on this package.

`AveragedCircuitEigenvalueSampling.jl` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits.
It is particularly interested in characterising the noise associated with fault-tolerant gadgets in the context of topological quantum error-correcting codes, such as syndrome extraction circuits.

The methods used in this package are detailed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545), building on the original method laid out in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).
This package relies on [Stim](https://github.com/quantumlib/Stim) for stabiliser circuit simulations.

## Example usage

Parameterise a depolarising noise model with

```
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
dep_param = DepolarisingParameters(r_1, r_2, r_m)
```

Generate the syndrome extraction circuit of a rotated planar code with

```
dist = 3
rotated_param = RotatedPlanarParameters(dist)
code = Code(rotated_param, dep_param)
```

Now optimise an experimental design for these parameters using the supplied least squares estimator

```
d = OptimiseDesign(code; ls_type = :wls)
```

Finally, simulate ACES noise characterisation using the specified number of measurement shots with 

```
shots_set = [10^6; 10^7; 10^8]
aces_data = SimulateACES(d, shots_set)
```

## Installation

This is not yet a registered package, and so must be downloaded manually.

If you encounter issues installing [Stim](https://github.com/quantumlib/Stim), a Python package used by PythonCall, try typing

```
julia> using CondaPkg

julia> # press ] to enter the Pkg REPL

pkg> conda resolve
```

Alternatively, install the package manually with

```
pkg> conda pip_add stim
```

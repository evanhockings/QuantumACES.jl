# ACES.jl

`ACES.jl` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits with averaged circuit eigenvalue sampling (ACES).
It is particularly interested in characterising the noise associated with fault-tolerant gadgets in the context of topological quantum error-correcting codes, such as syndrome extraction circuits.

The methods used in this package are detailed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545), building on the original ACES method laid out in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).
This package relies on [Stim](https://github.com/quantumlib/Stim) for stabiliser circuit simulations.

The [Package Guide](@ref) shows how to start using `ACES.jl`.

The [Public API](@ref) section documents the public functions and types of `ACES.jl`.

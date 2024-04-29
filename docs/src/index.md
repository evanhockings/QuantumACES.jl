# QuantumACES.jl

`QuantumACES.jl` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits with averaged circuit eigenvalue sampling (ACES).
It is particularly interested in characterising the noise associated with fault-tolerant gadgets in the context of topological quantum error correcting codes, such as surface code syndrome extraction circuits.

The methods used in this package are detailed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545), and the code generating the data for this paper can be found in the `scalable_aces` folder on the [scalable_aces](https://github.com/evanhockings/QuantumACES.jl/tree/scalable_aces) branch.
These methods build on the original ACES protocol presented in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).
This package relies on [Stim](https://github.com/quantumlib/Stim) for stabiliser circuit simulations.

The [Package Guide](@ref) describes how to start using `QuantumACES.jl`.

The [Public API](@ref) section documents the public functions and types of `QuantumACES.jl`.

# QuantumACES.jl

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

## Index

The [Package Guide](@ref) describes how you can start using `QuantumACES`, and [Creating Circuits and Noise Models](@ref) described the creation of new circuits and noise models to which the methods of `QuantumACES` can be applied.

The [Public API](@ref) documents the public functions and types, and the internal API is also documented.

The [Package Performance](@ref) section describes some important performance optimisations in `QuantumACES`.

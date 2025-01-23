---
title: "QuantumACES.jl: design noise characterisation experiments for quantum computers"
tags:
  - Julia
  - Quantum computing
  - Quantum error correction
authors:
  - name: Evan T. Hockings
    orcid: 0000-0002-5890-5705
    affiliation: 1, 2
affiliations:
  - name: School of Physics, The University of Sydney, Sydney, NSW 2006, Australia
    index: 1
  - name: ARC Centre of Excellence for Engineered Quantum Systems
    index: 2
date: 23 January 2025
bibliography: paper.bib
---

# Summary

QuantumACES.jl is a Julia&nbsp;[@bezanson_julia_2017] package for designing, simulating, and performing scalable Pauli noise characterisation experiments for quantum computers.
Noise in quantum devices is the key obstacle to large-scale quantum computation.
Consequently, quantum computers will require fault-tolerant architectures that replace physical qubits and operations with redundantly-encoded logical equivalents, which entails regularly measuring the parity checks of quantum error correcting codes&nbsp;[@shor_faulttolerant_1996; @gottesman_introduction_2010; @aliferis_introduction_2013].
Detailed characterisation of the noise in prototype quantum computers can demonstrate correct device performance and identify poorly functioning qubits and operations.
It can also enable noise-aware decoding of quantum error correcting codes by calibrating the decoder prior on which physical configurations of errors are more and less likely&nbsp;[@tuckett_faulttolerant_2020; @chen_calibrated_2022; @sundaresan_demonstrating_2023; @tiurev_correcting_2023; @higgott_improved_2023].

QuantumACES.jl uses the framework of averaged circuit eigenvalue sampling (ACES)&nbsp;[@flammia_averaged_2022] to design experiments that characterise Pauli noise in stabiliser circuits, following the theory and protocol outlined in&nbsp;[@hockings_scalable_2024].
Stabiliser circuits are a restricted class of quantum circuits that admit efficient classical simulation&nbsp;[@aaronson_improved_2004; @gottesman_stabilizer_1997], including with Pauli noise.
Quantum noise is tailored into Pauli noise by techniques such as Pauli frame randomisation&nbsp;[@knill_quantum_2005], randomised compiling&nbsp;[@wallman_noise_2016], or quantum error correction itself&nbsp;[@beale_coherence_2018], and the theory of quantum error correction and fault tolerance generally relies on modelling noise as Pauli noise&nbsp;[@terhal_quantum_2015].

QuantumACES.jl contains routines for designing optimised designs for noise characterisation experiments, given an arbitrary stabiliser circuit and Pauli noise model, using functions that precisely predict the performance of these experimental designs.
It has built-in circuits and noise models and also allows users to define their own.
These noise characterisation experiments are simulated with the open-source Python package Stim&nbsp;[@gidney_stim_2021], a fast simulator for stabiliser circuits with Pauli noise.

In a typical fault-tolerant quantum computing architecture, the bulk of the physical qubits and gate operations are dedicated to performing the syndrome extraction circuits that measure the parity checks of quantum error correcting codes.
These syndrome extraction circuits, which are stabiliser circuits, are therefore the key target for noise characterisation experiments.
QuantumACES.jl is tailored to characterising Pauli noise in syndrome extraction circuits, particularly for topological quantum error correcting codes such as the surface code&nbsp;[@bravyi_quantum_1998; @dennis_topological_2002; @kitaev_faulttolerant_2003; @fowler_surface_2012].
It leverages the fact that the simple structures of the syndrome extraction circuits of topological quantum codes remain similar across code sizes, enabling the optimised experimental design for the syndrome extraction circuit of a small-scale code to be transferred to syndrome extraction circuits of larger-scale versions of the same code.
QuantumACES.jl is capable of calculating and precisely fitting the performance scaling of these experimental designs as a function of the code size, enabling performance predictions at scales where explicit calculation becomes intractable.

Moreover, QuantumACES.jl makes it easy to simulate memory experiments for its syndrome extraction circuits with Stim, and supports decoding these experiments with the open-source Python packages PyMatching&nbsp;[@higgott_pymatching_2022; @higgott_sparse_2025] and BeliefMatching&nbsp;[@higgott_improved_2023].
It also contains an interface with the open-source Python package Qiskit&nbsp;[@javadi-abhari_quantum_2024], enabling the export of experimental designs and circuits to Qiskit circuits which can then be implemented on quantum devices to characterise noise in real quantum hardware.

# Statement of need

The utility of detailed and scalable Pauli noise characterisation methods grows as experimental progress pushes quantum devices towards scales of hundreds of qubits and initial demonstrations of fault tolerance.
QuantumACES.jl enables noise characterisation in this context, as demonstrated in&nbsp;[@hockings_scalable_2024].
While there are several software packages for benchmarking and noise characterisation, there are no open-source packages capable of detailed and scalable Pauli noise characterisation of quantum devices.
Forest-Benchmarking&nbsp;[@combes_forest_2019], is an open-source Python package containing many routines for quantum characterisation, verification, and validation (QCVV), but its detailed noise characterisation techniques are not scalable.
Gate set tomography (GST)&nbsp;[@nielsen_gate_2021], is a principled and extremely detailed noise characterisation protocol implemented by the open-source Python package pyGSTi&nbsp;[@nielsen_pygsti_2022], but it is limited to characterising extremely small numbers of qubits.
Cycle error reconstruction (CER)&nbsp;[@carignan-dugas_error_2023] is the noise characterisation protocol whose capabilities are most similar to ACES, but it is implemented by the commercial software True-Q&nbsp;[@beale_trueq_2020].

# Acknowledgements

This work was supported by the Australian Research Council Centre of Excellence for Engineered Quantum Systems (CE170100009), the U.S. Army Research Office (W911NF-21-1-0001, W911NF-23-S-0004), and the Unitary Foundation.

# References

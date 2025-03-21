---
title: "QuantumACES.jl: design noise characterisation experiments for quantum computers"
tags:
  - Julia
  - Quantum computing
  - Quantum error correction
authors:
  - name: Evan T. Hockings
    orcid: 0000-0002-5890-5705
    affiliation: Australian Research Council (ARC) Centre of Excellence for Engineered Quantum Systems, School of Physics, The University of Sydney, Sydney, New South Wales 2006, Australia
date: 22 March 2025
bibliography: paper.bib
---

# Summary

QuantumACES is a Julia [@bezanson_julia_2017] package for designing, simulating, and performing scalable Pauli noise characterisation experiments for quantum computers.
Noise in quantum devices is the key obstacle to large-scale quantum computation.
Consequently, quantum computers will require fault-tolerant architectures that replace physical qubits and operations with redundantly encoded logical equivalents, which entails regularly measuring the parity checks of quantum error correcting codes [@shor_faulttolerant_1996; @gottesman_introduction_2010; @aliferis_introduction_2013].
Decoders process the resulting error syndrome data and attempt to infer the most likely underlying physical errors in the quantum device.
Subsequently, they determine correction operations that attempt to preserve logical information.
The noise estimates produced by QuantumACES can inform decoders of the likelihood of physical error configurations in a quantum device [@tuckett_faulttolerant_2020; @chen_calibrated_2022; @sundaresan_demonstrating_2023; @tiurev_correcting_2023; @higgott_improved_2023], enabling noise-aware decoding.
They can also be used to generate simulated data for training more accurate decoders, such as in [@bausch_learning_2024], to verify appropriate device calibration, and to inform the co-design of quantum error correcting codes, decoders, fault-tolerant circuits, and quantum devices.

QuantumACES designs experiments to characterise Pauli noise in stabiliser circuits within the framework of averaged circuit eigenvalue sampling (ACES) [@flammia_averaged_2022], following the theory and protocol outlined in [@hockings_scalable_2025].
Stabiliser circuits are a restricted class of quantum circuits that admit efficient classical simulation [@aaronson_improved_2004; @gottesman_stabilizer_1997], including with Pauli noise.
Quantum noise is tailored into Pauli noise by techniques such as Pauli frame randomisation [@knill_quantum_2005], randomised compiling [@wallman_noise_2016], or quantum error correction itself [@beale_coherence_2018].
Additionally, the theory of quantum error correction and fault tolerance generally relies on modelling noise as Pauli noise [@terhal_quantum_2015].

QuantumACES contains routines for optimising designs for noise characterisation experiments, given an arbitrary stabiliser circuit and Pauli noise model, using functions that precisely predict the performance of these experimental designs.
It has built-in circuits and noise models and also allows users to define their own.
These noise characterisation experiments are simulated with the open-source Python package Stim [@gidney_stim_2021], a fast simulator for stabiliser circuits with Pauli noise.

In a typical fault-tolerant quantum computing architecture, the bulk of the physical qubits and gate operations are dedicated to performing the syndrome extraction circuits that measure the parity checks of quantum error correcting codes.
These syndrome extraction circuits, which are stabiliser circuits, are therefore the key target for noise characterisation experiments.
QuantumACES is tailored to characterising Pauli noise in syndrome extraction circuits, particularly for topological quantum error correcting codes such as the surface code [@bravyi_quantum_1998; @dennis_topological_2002; @kitaev_faulttolerant_2003; @fowler_surface_2012].
It leverages the fact that the simple structures of the syndrome extraction circuits of topological quantum codes remain similar across code sizes.
This enables optimised large-scale noise characterisation experiments that use experimental designs optimised at small scales.
QuantumACES is capable of calculating and precisely fitting the performance scaling of these experimental designs as a function of the code size, enabling performance predictions at scales where explicit calculation becomes intractable.

Moreover, QuantumACES supports the simulation and decoding of memory experiments for syndrome extraction circuits with Stim and the open-source Python packages PyMatching [@higgott_pymatching_2022; @higgott_sparse_2025] and BeliefMatching [@higgott_improved_2023], respectively.
It also provides an interface with the open-source Python package Qiskit [@javadi-abhari_quantum_2024], enabling the export of experimental designs to Qiskit circuits that can then be implemented to characterise noise in real quantum devices.

# Statement of need

The utility of detailed and scalable Pauli noise characterisation methods grows as experimental progress pushes quantum devices towards scales of hundreds of qubits and initial demonstrations of fault tolerance.
QuantumACES enables noise characterisation and noise-aware decoding in this context, as demonstrated in [@hockings_scalable_2025] and [@hockings_improving_2025], respectively.
While there are several software packages for benchmarking and noise characterisation, there are no open-source packages capable of detailed and scalable Pauli noise characterisation of quantum devices.
Forest-Benchmarking [@combes_forest_2019] is an open-source Python package containing many routines for quantum characterisation, verification, and validation, but its detailed noise characterisation techniques are not scalable.
Gate set tomography [@nielsen_gate_2021] is a principled and extremely detailed noise characterisation protocol implemented by the open-source Python package pyGSTi [@nielsen_pygsti_2022], but it is limited to characterising very small numbers of qubits.
Cycle error reconstruction [@carignan-dugas_error_2023] is the noise characterisation protocol whose capabilities are most similar to ACES, but it is implemented by the commercial software True-Q [@beale_trueq_2020].

# Acknowledgements

This work was supported by the Australian Research Council Centre of Excellence for Engineered Quantum Systems (CE170100009), the U.S. Army Research Office (W911NF-21-1-0001, W911NF-23-S-0004), and the Unitary Foundation.

# References

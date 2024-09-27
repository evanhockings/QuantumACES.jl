---
title: "QuantumACES.jl: Design noise characterisation experiments for quantum computers"
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
date: 16 September 2024
bibliography: paper.bib
---

# Summary

`QuantumACES.jl` is a Julia [@bezanson_julia_2017] package for designing, optimising, and simulating scalable and performant noise characterisation experiments for quantum computers.

Noise in quantum devices is the key obstacle to large-scale quantum computation.
Consequently, quantum computers will require fault-tolerant architectures that encode quantum information in quantum error correcting codes [@shor_faulttolerant_1996; @gottesman_introduction_2010].
Detailed characterisation of the noise in prototype quantum computers can demonstrate that the device is performing correctly, and improve the performance of error correction by identifying which physical configurations of errors are more and less likely [@sundaresan_demonstrating_2023].

`QuantumACES.jl` designs experiments with the framework of averaged circuit eigenvalue sampling (ACES) [@flammia_averaged_2022] that characterise Pauli noise in stabiliser circuits, following theory outlined in [@hockings_scalable_2024].
Stabiliser circuits are a restricted class of quantum circuits that admit efficient classical simulation [@aaronson_improved_2004; @gottesman_stabilizer_1997].

In a typical fault-tolerant quantum computing architecture, syndrome extraction circuits, which are stabiliser circuits that measure the parity checks of quantum error correcting codes, comprise the bulk of the physical qubits and gate operations.
Moreover, the theory of quantum error correction and fault-tolerance generally relies on modelling noise as Pauli noise [@terhal_quantum_2015].
Quantum noise is tailored into Pauli noise by techniques such as Pauli frame randomisation or randomised compiling [@knill_quantum_2005; @wallman_noise_2016], or quantum error correction itself [@beale_coherence_2018].

Accordingly, `QuantumACES.jl` focuses on characterising Pauli noise in syndrome extraction circuits and fault-tolerant logical circuits.
It contains routines for creating and optimising experimental designs according to their predicted performance.
The optimised experimental design for the syndrome extraction circuit of a small-scale quantum error correcting code can be applied to that circuit for larger-scale versions of the same code, with predictable performance scaling.
End-to-end simulation of a noise characterisation experiment leverages the open-source Python package `Stim` [@gidney_stim_2021] to rapidly simulate stabiliser circuits with Pauli noise, allowing verification of these performance predictions.

# Statement of need

As experimental progress pushes quantum devices towards scales of hundreds of qubits and initial demonstrations of fault-tolerance, the utility of detailed and scalable noise characterisation methods also grows.
`QuantumACES.jl` enables noise characterisation in this context, as demonstrated in [@hockings_scalable_2024].
While there are a number of software packages for benchmarking and noise characterisation, there are no open-source packages capable of detailed and scalable noise characterisation of quantum devices.
`Forest-Benchmarking` [@combes_forest_2019] is an open-source Python package containing a large number of routines for quantum characterisation, verification, and validation (QCVV), but its detailed noise characterisation techniques are not scalable.
Gate set tomography (GST) [@nielsen_gate_2021] is a principled and detailed noise characterisation technique implemented by the open-source Python package `pyGSTi` [@nielsen_pygstio_2022], but it is only capable of characterising small numbers of qubits.
Cycle error reconstruction (CER) [@carignan-dugas_error_2023] is a noise characterisation technique whose capabilities are most similar to ACES, but it is implemented by the commercial software True-Q [@beale_trueq_2020].

# Acknowledgements

ETH is supported by an Australian Government Research Training Program (RTP) Scholarship.
This work was supported by the Australian Research Council Centre of Excellence for Engineered Quantum Systems (EQUS, CE170100009).

# References
# Scalable noise characterisation of syndrome extraction circuits with averaged circuit eigenvalue sampling

This repository contains the code used to generate the results for [arXiv:2404.06545](https://arxiv.org/abs/2404.06545).

For these files, the prefix `rot` refers to the rotated surface code, whereas the prefix `unrot` refers to the unrotated surface code, and the Jupyter notebooks corresponding to each file plots and displays the results.
- `optimise` files optimise designs at a range of depolarising noise strengths.
- `scaling` files calculate the scaling behaviour of the optimised design for depolarising and log-normal Pauli noise.
- `simulate` files perform many small-scale simulations for optimised and basic designs, both for depolarising noise and the seed-0 instance of log-normal Pauli noise.
- `simulate_big` files perform large-scale simulation for optimised and basic designs, both for depolarising noise and the seed-0 instance of log-normal Pauli noise.
- `runfiles` files run all of the above files.

Moreover, the `toy_design.ipynb` notebook examines relative precision estimation in the context of a toy experimental design, and the `google_data.ipynb` notebook examines the error probabilities of the Google quantum device used in their surface code scaling experiment.

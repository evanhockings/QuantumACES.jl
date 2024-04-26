# Package Guide

## Index

```@index
Pages = ["manual.md"]
```

## Introduction


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

## Creating new circuits


## Creating new noise models


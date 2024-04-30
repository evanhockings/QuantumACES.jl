# Package Guide

## Introduction

`QuantumACES.jl` is a package for designing and simulating scalable and performant Pauli noise characterisation experiments for stabiliser circuits with averaged circuit eigenvalue sampling (ACES).
It is particularly interested in characterising the noise associated with fault-tolerant gadgets in the context of topological quantum error correcting codes, such as surface code syndrome extraction circuits.

The methods used in this package are detailed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545), and the code generating the data for this paper can be found in the `scalable_aces` folder on the [scalable_aces](https://github.com/evanhockings/QuantumACES.jl/tree/scalable_aces) branch.
These methods build on the original ACES protocol presented in [arXiv:2108.05803](https://arxiv.org/abs/2108.05803).
This package relies on [Stim](https://github.com/quantumlib/Stim) for stabiliser circuit simulations.

Typical usage of this package involves a few steps:

  - Constructing the circuit and noise model for which you aim to perform an ACES noise characterisation experiment, using either provided functions or your own.
  - Optimise an experimental design for a depolarising noise model with the same average error rates as the noise model you aim to characterise.
      - Optional: Predict the performance of the experimental design, including its scaling as a function of circuit parameters, for example as a function of the code distance in the case of surface code syndrome extraction circuits.
  - Simulate noise characterisation experiments across a range of specified measurement budgets that use the optimised experimental design. 

## Installation and setup

To install this package, run the following command in the Julia REPL.

```
] add QuantumACES
```

This package relies on the Python package [Stim](https://github.com/quantumlib/Stim) to perform stabiliser circuit simulations.
It calls Stim with [PythonCall](https://github.com/JuliaPy/PythonCall.jl).
By default, PythonCall creates its own Python environment, but you may wish to use an existing Python installation.

One helpful method for managing Python versions is [pyenv](https://github.com/pyenv/pyenv), or for Windows, [pyenv-win](https://github.com/pyenv-win/pyenv-win); these are analogous to [Juliaup](https://github.com/JuliaLang/juliaup) for Julia.
The following assumes you are using pyenv or pyenv-win.

On Windows, to instruct PythonCall to use the Python version set by pyenv, configure PythonCall's environment variables by adding the following to your `~/.julia/config/startup.jl` file

```julia
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
python_exe = readchomp(`cmd /C pyenv which python`)
ENV["JULIA_PYTHONCALL_EXE"] = python_exe
```

On Unix systems, shell commands are parsed directly by Julia and appear to be unaware of your PATH variable, and I am not sure how to work around this.
Therefore, you may need to manually supply `python_exe` for the Python version `<version>` as

```julia
python_exe = homedir() * "/.pyenv/versions/<version>/bin/python"
```

Then ensure Stim is installed by running `pip install stim` in your terminal.

## Example usage

To see a full example showing how this package can be used, see the code that generates the data for [arXiv:2404.06545](https://arxiv.org/abs/2404.06545) in the `scalable_aces` folder on the [scalable_aces](https://github.com/evanhockings/QuantumACES.jl/tree/scalable_aces) branch.

Beware that the examples shown below may take a number of hours to run.
Ensure that Julia is set up to use as many threads as your CPU can handle.

First parameterise a depolarising noise model with single-qubit gate infidelity `r_1`, two-qubit gate infidelity `r_2`, and measurement infidelity `r_m`, and a log-normal random Pauli noise model with the same gate infidelities and a standard deviation of the underlying normal distributions `total_std_log`, alongside a random seed used when generating the noise model.

```julia
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
```

Then generate the syndrome extraction circuit for a distance `dist` (rotated) surface code.

```julia
dist = 3
rotated_param = get_rotated_param(dist)
rotated_planar = get_circuit(rotated_param, dep_param)
```

Next, generate an experimental design for this circuit.

```julia
d = generate_design(rotated_planar)
```

Alternatively, optimise an experimental design to improve its sample efficiency, configuring the optimisation with the parameters associated with [`OptimOptions`](@ref).

```julia
d = optimise_design(rotated_planar, options = OptimOptions(; ls_type = :wls, seed = seed))
```

Create a copy of the optimised design that associates log-normal random Pauli noise with the circuit, and simulate `repetitions` rounds of an ACES noise characterisation experiment across all of the supplied measurement budgets in `budget_set`, which are measurement shots normalised by the time taken to perform the experiment.

```julia
d_log = update_noise(d, log_param)
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_log, budget_set; repetitions = repetitions)
```

We can compare the performance to performance predictions at the largest measurement budget, although we note that the z-scores will not be normally distributed as the underlying distribution is not quite normal.

```julia
wls_merit_log = calc_wls_merit(d_log)
fgls_z_scores =
    (aces_data.fgls_gate_norm_coll[:, 3] .- wls_merit_log.expectation) /
    sqrt(wls_merit_log.variance)
```

Next, calculate the performance scaling of this design as a fuction of the code distance up to some large distance `dist_max` with the weighted least squares (WLS) estimator, for both depolarising and log-normal random Pauli noise.

```julia
dist_max = 9
dep_planar_scaling = calc_depolarising_planar_scaling(d, dist_max; ls_type = :wls)
log_planar_scaling = calc_lognormal_planar_scaling(d_log, dist_max; ls_type = :wls, seed = seed)
```

Next, transfer the optimised experimental design to the syndrome extraction circuit for a distance `dist_big` surface code with log-normal random Pauli noise.

```julia
dist_big = 13
rotated_param_big = get_rotated_param(dist_big)
rotated_planar_big = get_circuit(rotated_param_big, log_param)
d_big = generate_design(rotated_planar_big, d.tuple_set_data)
```

We can simulate a large-scale ACES noise characterisation experiment across the supplied measurement budgets.

```julia
budget_set_big = [10^6; 10^7; 10^8; 10^9]
aces_data_big = simulate_aces(d_big, budget_set_big)
```

Finally, we compare the performance to predictions across the measurement budgets, although note that we would not expect the z-scores here to actually correspond to a normal distribution as the underlying distribution is not quite normal, and there is a substantive amount of uncertainty associated with the fit.

```julia
pred_expectation = log_planar_scaling.expectation_fit(dist_big)
pred_variance = log_planar_scaling.variance_fit(dist_big)
fgls_z_scores_big =
    (aces_data_big.fgls_gate_norm_coll[1, :] .- pred_expectation) / sqrt(pred_variance)
```

## More advanced usage

This package also supports creating your own circuits and noise models.

Let us begin by creating a new circuit, following the example circuit shown in Figure 2 of [arXiv:2404.06545](https://arxiv.org/abs/2404.06545).
The first step is to create a parameter struct for the circuit, which must be a subtype of [`AbstractCircuitParameters`](@ref) and contain the necessary fields `layer_time_dict` and `circuit_name`.

```julia
struct ExampleParameters <: AbstractCircuitParameters
    pad_identity::Bool
    layer_time_dict::Dict{Symbol, Float64}
    circuit_name::String
end
```

We need a function to construct the parameter struct.

```julia
function get_example_parameters(;
    pad_identity = true,
    single_qubit_time::Float64 = 29.0,
    two_qubit_time::Float64 = 29.0,
    meas_reset_time::Float64 = 660.0,
)
    # Create the example parameters
    @assert single_qubit_time > 0.0 "The single-qubit layer time must be positive."
    @assert two_qubit_time > 0.0 "The two-qubit layer time must be positive."
    @assert meas_reset_time > 0.0 "The measurement and reset layer time must be positive."
    layer_time_dict = Dict(
        :single_qubit => single_qubit_time,
        :two_qubit => two_qubit_time,
        :meas_reset => meas_reset_time,
    )
    circuit_name = "example_circuit"
    if pad_identity != true
        circuit_name *= "_pad_identity_$(pad_identity)"
    end
    example_param = ExampleParameters(pad_identity, layer_time_dict, "example_circuit")
    return example_param::ExampleParameters
end
```

And we need a function to create the circuit from the parameter struct.

```julia
function example_circuit(example_param::ExampleParameters)
    # Set up variables
    pad_identity = example_param.pad_identity
    layer_time_dict = example_param.layer_time_dict
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    # Generate the circuit
    qubit_num = 3
    circuit = [
        Layer([Gate("CZ", 0, [2; 3])], qubit_num),
        Layer([Gate("CX", 0, [1; 2]), Gate("H", 0, [3])], qubit_num),
        Layer([Gate("H", 0, [1]), Gate("S", 0, [2]), Gate("H", 0, [3])], qubit_num),
    ]
    layer_types = [two_qubit_type, two_qubit_type, single_qubit_type]
    layer_times = get_layer_times(layer_types, layer_time_dict)
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [pad_layer(l) for l in circuit]
    end
    return (
        circuit::Vector{Layer},
        qubit_num::Int,
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
    )
end
```

Finally, we create a function that generates the circuit in the form of a [`Circuit`](@ref) struct.
We do this by adding a method to [`get_circuit`](@ref) which uses the new parameter struct as an argument.
A helpful function for this is [`prepare_circuit`](@ref), which deals with much of the busywork associated with generating the circuit.

```julia
function QuantumACES.get_circuit(
    example_param::ExampleParameters,
    noise_param::T;
    add_prep::Bool = false,
    add_meas::Bool = true,
) where {T <: AbstractNoiseParameters}
    # Construct the circuit
    (circuit, qubit_num, layer_types, layer_times) = example_circuit(example_param)
    circuit_tuple = collect(1:length(circuit))
    # Prepare the circuit and generate additional parameters
    (
        labelled_circuit,
        unique_layer_indices,
        gates,
        total_gates,
        gate_index,
        N,
        gate_probabilities,
        gate_eigenvalues,
    ) = prepare_circuit(
        circuit,
        qubit_num,
        layer_types,
        layer_times,
        noise_param;
        add_prep = add_prep,
        add_meas = add_meas,
    )
    # Return the circuit
    c = Circuit(
        example_param,
        labelled_circuit,
        circuit_tuple,
        qubit_num,
        unique_layer_indices,
        layer_types,
        layer_times,
        gates,
        total_gates,
        gate_index,
        N,
        noise_param,
        gate_probabilities,
        gate_eigenvalues,
        add_prep,
        add_meas,
    )
    return c::Circuit
end
```

Make sure to follow the form of this example function when creating your own circuit.
If you want to include more data with your circuit, create a new struct that is a subtype of [`AbstractCircuit`](@ref), and add the necessary fields, namely all of those present in [`Circuit`](@ref).

Next, we will create a phenomenological noise model where each Pauli error probability has some probability `p` of occurring, so that two-qubit errors have a probability `p^2`.
As with the circuit, we begin by creating a parameter struct for the noise model, which must be a subtype of [`AbstractNoiseParameters`](@ref) and contain the necessary field `noise_name`.

```julia
struct PhenomenologicalParameters <: AbstractNoiseParameters
    p::Float64
    m::Float64
    noise_name::String
end
```

We need a function to construct the parameter struct.

```julia
function get_phen_param(p::Float64, m::Float64)
    @assert (p >= 0) && (p <= 1 / 10) "The phenomenological gate error probability $(p) is out of bounds."
    @assert (m >= 0) && (m <= 1 / 2) "The phenomenological measurement error probability $(m) is out of bounds."
    noise_name = "phenomenological_$(round(p; sigdigits = 4))_$(round(m; sigdigits = 4))"
    return PhenomenologicalParameters(p, m, noise_name)::PhenomenologicalParameters
end
```

And we need a function to create the noise model for a set of gates from the parameter struct.
As with the circuit, we add a method to [`get_gate_probabilities`](@ref) which uses the new parameter struct as an argument.

```julia
function QuantumACES.get_gate_probabilities(
    total_gates::Vector{Gate},
    noise_param::PhenomenologicalParameters,
)
    # Extract the parameters for generating the noise
    p = noise_param.p
    m = noise_param.m
    # Determine the weight of the error corresponding to each gate error probability
    one_qubit_support_size = ones(3)
    n = 2
    two_qubit_support_size = Vector{Int}(undef, 0)
    bit_array = BitArray(undef, 2n + 1)
    for bit_array.chunks[1] in 1:(2^(2n) - 1)
        two_qubit_pauli = Pauli(convert(Vector{Bool}, bit_array), n)
        push!(two_qubit_support_size, length(get_support(two_qubit_pauli)))
    end
    @assert sum(two_qubit_support_size .== 1) == 6
    @assert sum(two_qubit_support_size .== 2) == 9
    # Generate the noise
    gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate in total_gates
        if gate.type ∈ ["MZ", "MX", "MY"]
            probability = [m]
        elseif length(gate.targets) == 1
            probability = p .^ one_qubit_support_size
        elseif length(gate.targets) == 2
            probability = p .^ two_qubit_support_size
        else
            throw(error("The gate $(gate) is unsupported."))
        end
        @assert sum(probability) < 1 "The probabilities $(probability) sum to more than 1; change the input parameters."
        gate_probabilities[gate] = [1 - sum(probability); probability]
    end
    return gate_probabilities::Dict{Gate, Vector{Float64}}
end
```

Make sure to follow the form of this example function when creating your own noise model.

Now we are ready to reprise the basic usage example for this example circuit.
First, parameterise depolarising and phenomenological noise models, making sure to create a depolarising noise model with the same gate infidelity as the phenomenological noise model.

```julia
p = 0.025 / 100
m = 2.0 / 100
r_1 = 3 * p
r_2 = 6 * p + 9 * p^2
r_m = m
phen_param = get_phen_param(p, m)
dep_param = get_dep_param(r_1, r_2, r_m)
```

Then construct the circuit

```
example_param = get_example_parameters()
circuit_example = get_circuit(example_param, dep_param)
```

Optimise the experimental design for the generalised least squares (GLS) estimator as we are not interested in scaling this experimental design to large numbers of qubits, configuring the optimisation with the parameters associated with [`OptimOptions`](@ref).
This is because the circuit acts on only three qubits and, unlike the surface code syndrome extraction circuits, does not form a family of circuits.

```julia
seed = UInt(0)
d = optimise_design(circuit_example; options = OptimOptions(; ls_type = :gls, seed = seed))
```

Create a copy of the optimised design that associates phenomenological noise with the circuit, and compare the predicted performance of the experimental design with depolarising and phenomenological noise.

```julia
d_phen = update_noise(d, phen_param)
merit_dep = calc_gls_merit(d)
merit_phen = calc_gls_merit(d_phen)
```

We can also simulate the performance of the experimental design with phenomenological noise.

```julia
budget_set = [10^6; 10^7; 10^8]
repetitions = 20
aces_data = simulate_aces(d_phen, budget_set; repetitions = repetitions, seed = seed)
```

Finally, compare the performance to predictions at the largest measurement budget.

```julia
fgls_z_scores_phen =
    (aces_data.fgls_gate_norm_coll[:, 3] .- merit_phen.expectation) /
    sqrt(merit_phen.variance)
```

As before, note that the distribution of the normalised RMS error between the estimated and true gate eigenvalues is not quite normally distributed.
Hence the z-scores shown here, which are normalised by the predicted performance of the experimental design, will not quite be normally distributed.

For a fuller understanding of the methods used in this package, refer to [arXiv:2404.06545](https://arxiv.org/abs/2404.06545).

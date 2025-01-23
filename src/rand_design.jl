"""
    RandDesign

Pauli frame randomised experimental ensemble for an experimental design.

# Fields

  - `d::Design`: Experimental design.
  - `pauli_randomisation_ensemble::Vector{Vector{Vector{Vector{String}}}}`: For each tuple, experiment, and randomisation, a vector of Pauli frame randomisations between each layer in the randomised circuit.
  - `pauli_sign_ensemble::Vector{Vector{Vector{Vector{Bool}}}}`: For each tuple, experiment, and randomisation, a vector of signs for the Paulis estimated by the randomised experiment.
  - `covariance_pauli_sign_ensemble::Vector{Vector{Vector{Vector{Bool}}}}`: For each tuple, experiment, and randomisation, a vector of signs for the covariance Paulis estimated by the randomised experiment.
  - `job_indices::Vector{Vector{NTuple{3, Int}}}`: For each job, a vector of indices (tuple, experiment, randomisation) for each circuit in the job.
  - `randomisations::Vector{Int}`: Number of randomisations for the experiment set corresponding to each tuple in the design.
  - `shot_budget::Int`: Shot budget for the randomised experimental ensemble.
  - `experiment_shots::Int`: Number of shots for each experiment in the randomised experimental ensemble.
  - `seed::UInt64`: Random seed used to generate the randomised experimental ensemble.
"""
struct RandDesign
    d::Design
    pauli_randomisation_ensemble::Vector{Vector{Vector{Vector{String}}}}
    pauli_sign_ensemble::Vector{Vector{Vector{Vector{Bool}}}}
    covariance_pauli_sign_ensemble::Vector{Vector{Vector{Vector{Bool}}}}
    job_indices::Vector{Vector{NTuple{3, Int}}}
    randomisations::Vector{Int}
    shot_budget::Int
    experiment_shots::Int
    seed::UInt64
end

function Base.show(io::IO, d_rand::RandDesign)
    println(
        io,
        "Randomised experimental design for a $(d_rand.d.c.circuit_param.circuit_name) circuit with $(length(d_rand.d.tuple_set)) tuples and $(d_rand.d.experiment_number) experiments measured across $(sum(d_rand.randomisations .* d_rand.d.experiment_numbers)) randomised experiments.",
    )
    pretty_print(io, d_rand.d)
    return nothing
end

@struct_hash_equal_isequal RandDesign

"""
    get_pauli_layer(pauli_string::String, n::Int16)

Returns a layer of Pauli gates corresponding to the supplied Pauli string `pauli_string` on `n` qubits.`
"""
function get_pauli_layer(pauli_string::String, n::Int16)
    @assert length(pauli_string) == n "The Pauli layer has the wrong number of qubits."
    pauli_layer = Layer([Gate(string(pauli_string[i]), -1, [i]) for i in 1:n], n)
    return pauli_layer::Layer
end

"""
    get_randomised_circuit(prep_layer::Layer, tuple_circuit::Vector{Layer}, meas_layer::Layer, pauli_randomisation::Vector{String})

Returns the randomised circuit corresponding to the preparation layer `prep_layer`, tuple circuit `tuple_circuit`, measurement layer `meas_layer`, interspersed with the Pauli frame randomisation layers specified by `pauli_randomisation`.
"""
function get_randomised_circuit(
    prep_layer::Layer,
    tuple_circuit::Vector{Layer},
    meas_layer::Layer,
    pauli_randomisation::Vector{String},
)
    # Initialise variables
    n = prep_layer.qubit_num
    @assert all(l.qubit_num == n for l in tuple_circuit) "The circuit layers do not have the same qubit number."
    @assert meas_layer.qubit_num == n "The circuit layers do not have the same qubit number."
    L = length(tuple_circuit)
    @assert length(pauli_randomisation) == L + 2 "The Pauli randomisation does not have the correct number of layers."
    randomised_circuit = Vector{Layer}(undef, 2 * L + 4)
    # Set up the randomisation and preparation layer
    randomised_circuit[1] = get_pauli_layer(pauli_randomisation[1], n)
    randomised_circuit[2] = prep_layer
    for idx in 1:L
        # Combine the randomisation layer with the circuit layer
        randomised_circuit[2 * idx + 1] = get_pauli_layer(pauli_randomisation[idx + 1], n)
        randomised_circuit[2 * idx + 2] = tuple_circuit[idx]
    end
    # Set up the randomisation and measurement layer
    randomised_circuit[2 * L + 3] = get_pauli_layer(pauli_randomisation[L + 2], n)
    randomised_circuit[2 * L + 4] = meas_layer
    return randomised_circuit::Vector{Layer}
end

"""
    calc_mapping_signs(experiment_mappings::Vector{Mapping}, randomised_circuit::Vector{Layer})

Returns the signs of the Paulis corresponding to the mappings in `experiment_mappings` for the randomised circuit `randomised_circuit`.
"""
function calc_mapping_signs(
    experiment_mappings::Vector{Mapping},
    randomised_circuit::Vector{Layer},
)
    # Initialise variables
    n = randomised_circuit[1].qubit_num
    R = length(randomised_circuit)
    @assert R >= 4 "The randomised circuit does not have enough layers."
    @assert all(l.qubit_num == n for l in randomised_circuit) "The randomised circuit layers do not have the same qubit number."
    E = length(experiment_mappings)
    initial_support_set = [get_support(m.initial) for m in experiment_mappings]
    final_support_set = [get_support(m.final) for m in experiment_mappings]
    # Set up the tableau
    t = Tableau(n)
    # Apply the randomisation and preparation layers
    apply!(t, randomised_circuit[1])
    apply!(t, randomised_circuit[2])
    # Check the initial Paulis
    for idx in 1:E
        m = experiment_mappings[idx]
        initial_support = initial_support_set[idx]
        # Use row_sum! to assemble the Pauli
        for i in eachindex(initial_support[2:end])
            row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
        end
        @assert t.tableau[n + initial_support[1], 1:(2n)] == m.initial.pauli[1:(2n)] "The initial Pauli has not been appropriately initialised in the tableau."
        # Use row_sum! to disassemble the Pauli
        for i in eachindex(initial_support[2:end])
            row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
        end
    end
    # Apply the circuit and randomisation layers
    for l in randomised_circuit[3:(R - 1)]
        apply!(t, l)
    end
    # Check the final Paulis and record their signs
    pauli_signs = Vector{Bool}(undef, E)
    for idx in 1:E
        m = experiment_mappings[idx]
        initial_support = initial_support_set[idx]
        # Use row_sum! to assemble the Pauli
        for i in eachindex(initial_support[2:end])
            row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
        end
        @assert t.tableau[n + initial_support[1], 1:(2n)] == m.final.pauli[1:(2n)] "The final Pauli has not been appropriately initialised in the tableau."
        pauli_signs[idx] = t.tableau[n + initial_support[1], 2n + 1]
        # Use row_sum! to disassemble the Pauli
        for i in eachindex(initial_support[2:end])
            row_sum!(t, n + initial_support[1], n + initial_support[i + 1])
        end
    end
    # Apply the measurement layer
    measurements = apply!(t, randomised_circuit[R]; return_measurements = true)
    for idx in 1:E
        final_support = final_support_set[idx]
        @assert prod(meas[1] for meas in measurements[final_support]) ==
                (-1)^pauli_signs[idx] "The measurement signs do not match the expected Pauli signs."
    end
    return pauli_signs::Vector{Bool}
end

"""
    get_randomisations(d::Design, min_randomisations::Integer, target_shot_budget::Integer, experiment_shots::Integer; p_norm::Real = 2)

Returns the number of randomisations for each tuple in the experimental design `d`, given a target shot budget `target_shot_budget`, number of shots for each experiment `experiment_shots`, and minimum number of randomisations `min_randomisations`.
The number of randomisations is greedily optimised to minimise the `p_norm`-norm to the shot weights associated with the experimental design.

Typically, the 2-norm works better for smaller shot budgets, and the 1-norm works better for larger shot budgets.
The number of randomisations corresponds to unnormalised shot weights `randomisations .* d.experiment_numbers`.
"""
function get_randomisations(
    d::Design,
    min_randomisations::Integer,
    target_shot_budget::Integer,
    experiment_shots::Integer;
    p_norm::Real = 2,
)
    # Check that the minimum number of randomisations is appropriate
    tuple_number = length(d.tuple_set)
    experiment_numbers = d.experiment_numbers
    # Set up the number of randomisations
    randomisations = min_randomisations * ones(Int, tuple_number)
    shot_budget = sum(experiment_shots * randomisations .* experiment_numbers)
    if shot_budget > target_shot_budget
        @warn "The supplied experiment shots $(experiment_shots) with the minimum number of randomisations $(min_randomisations) produces a minimum shot budget $(shot_budget) that exceeds the target $(target_shot_budget)."
        adding = false
    else
        adding = true
    end
    # Iteratively add randomisations in order to minimise the 2-norm to the design's shot weights
    while adding
        trial_norms = Vector{Float64}(undef, tuple_number)
        for i in 1:tuple_number
            trial_randomisations = deepcopy(randomisations)
            trial_randomisations[i] += 1
            trial_unnormalised_shot_weights = trial_randomisations .* experiment_numbers
            trial_shot_weights =
                trial_unnormalised_shot_weights / sum(trial_unnormalised_shot_weights)
            trial_norms[i] = norm(trial_shot_weights - d.shot_weights, p_norm)
        end
        add_randomisations = deepcopy(randomisations)
        add_randomisations[findmin(trial_norms)[2]] += 1
        add_shot_budget = sum(experiment_shots * add_randomisations .* experiment_numbers)
        if add_shot_budget > target_shot_budget
            adding = false
        else
            randomisations = add_randomisations
        end
    end
    shot_budget = sum(experiment_shots * randomisations .* experiment_numbers)
    return (randomisations::Vector{Int}, shot_budget::Integer)
end

"""
    generate_rand_design(d::Design, randomisations::Vector{Int}, shot_budget::Integer, experiment_shots::Integer; kwargs...)
    generate_rand_design(d::Design, min_randomisations::Integer, target_shot_budget::Integer, experiment_shots::Integer; kwargs...)

Returns a [`RandDesign`](@ref) object describing a Pauli frame randomised experimental ensemble for the design `d`.
Supply either the number of randomisations `randomisations` and shot budget `shot_budget` output by [`get_randomisations`](@ref), as well as the number of shots for each experiment `experiment_shots`, or the arguments of that function, namely the minimum number of randomisations `min_randomisations`, target shot budget `target_shot_budget`, and number of shots for each experiment `experiment_shots`.

# Arguments

  - `d::Design`: Experimental design.
  - `randomisations::Vector{Int}`: Number of randomisations for each tuple in the design.
  - `shot_budget::Integer`: Shot budget for the randomised experimental ensemble.
  - `min_randomisations::Integer`: Minimum number of randomisations for each tuple in the design.
  - `target_shot_budget::Integer`: Target shot budget for the randomised experimental ensemble.
  - `experiment_shots::Integer`: Number of shots for each experiment in the randomised experimental ensemble.

# Keyword arguments

  - `job_circuit_number::Integer = 200`: Number of circuits in each job.
  - `seed::Union{UInt64, Nothing} = nothing`: Random seed for the randomisation.
  - `diagnostics::Bool = false`: Whether to print diagnostics.
  - `save_data::Bool = false`: Whether to save the randomised design.
"""
function generate_rand_design(
    d::Design,
    randomisations::Vector{Int},
    shot_budget::Integer,
    experiment_shots::Integer;
    job_circuit_number::Integer = 200,
    seed::Union{UInt64, Nothing} = nothing,
    diagnostics::Bool = false,
    save_data::Bool = false,
)
    # Initialise variables and generate random seeds
    start_time = time()
    n = d.c.qubit_num
    tuple_number = length(d.tuple_set)
    experiment_numbers = d.experiment_numbers
    @assert all(randomisations .> 1) "Each tuple must have more than 1 randomisation."
    if seed === nothing
        seed = rand(UInt64)
    end
    Random.seed!(seed)
    seeds = [
        [rand(UInt64, randomisations[i]) for j in 1:experiment_numbers[i]] for
        i in 1:tuple_number
    ]
    perm_seed = rand(UInt64)
    Random.seed!(perm_seed)
    tuple_experiments = randomisations .* experiment_numbers
    total_experiments = sum(tuple_experiments)
    total_experiments_perm = randperm(total_experiments)
    Random.seed!()
    # Generate the Pauli frame randomisation layers
    pauli_randomisation_ensemble =
        Vector{Vector{Vector{Vector{String}}}}(undef, tuple_number)
    ensemble_indices = Vector{NTuple{3, Int}}()
    for (i, circuit_tuple) in pairs(d.tuple_set)
        tuple_circuit = d.c.circuit[circuit_tuple]
        L = length(tuple_circuit)
        pauli_randomisation_ensemble[i] =
            Vector{Vector{Vector{String}}}(undef, experiment_numbers[i])
        for j in 1:experiment_numbers[i]
            pauli_randomisation_ensemble[i][j] =
                Vector{Vector{String}}(undef, randomisations[i])
            for r in 1:randomisations[i]
                # Get the random Pauli frame for each layer
                Random.seed!(seeds[i][j][r])
                pauli_randomisation_ensemble[i][j][r] =
                    [get_random_pauli(n) for idx in 1:(L + 2)]
                Random.seed!()
                push!(ensemble_indices, (i, j, r))
            end
        end
    end
    # Calculate the signs of the Paulis and covariance Paulis in each randomised experiment
    experiment_ensemble = get_experiment_data(d)[1]
    covariance_keys_ensemble = get_covariance_mapping_data(d)[1]
    covariance_experiment_ensemble = get_covariance_experiment_data(d)[1]
    covariance_mapping_lengths = length.(covariance_keys_ensemble)
    pauli_sign_ensemble = Vector{Vector{Vector{Vector{Bool}}}}(undef, tuple_number)
    covariance_pauli_sign_ensemble =
        Vector{Vector{Vector{Vector{Bool}}}}(undef, tuple_number)
    for i in 1:tuple_number
        tuple_start = time()
        circuit_tuple = d.tuple_set[i]
        tuple_circuit = d.c.circuit[circuit_tuple]
        mapping_set = d.mapping_ensemble[i]
        covariance_dict = d.covariance_dict_ensemble[i]
        covariance_keys = covariance_keys_ensemble[i]
        pauli_sign_ensemble[i] = Vector{Vector{Vector{Bool}}}(undef, experiment_numbers[i])
        has_covariance = (covariance_mapping_lengths[i] > 0)
        if has_covariance
            covariance_pauli_sign_ensemble[i] =
                Vector{Vector{Vector{Bool}}}(undef, experiment_numbers[i])
        else
            covariance_pauli_sign_ensemble[i] = Vector{Vector{Vector{Bool}}}()
        end
        for j in 1:experiment_numbers[i]
            prep_layer = d.prep_ensemble[i][j]
            meas_layer = d.meas_ensemble[i][j]
            experiment = experiment_ensemble[i][j]
            pauli_sign_ensemble[i][j] = Vector{Vector{Bool}}(undef, randomisations[i])
            if has_covariance
                covariance_experiment = covariance_experiment_ensemble[i][j]
                covariance_pauli_sign_ensemble[i][j] =
                    Vector{Vector{Bool}}(undef, randomisations[i])
            end
            for r in 1:randomisations[i]
                # Get the randomised circuit
                pauli_randomisation = pauli_randomisation_ensemble[i][j][r]
                randomised_circuit = get_randomised_circuit(
                    prep_layer,
                    tuple_circuit,
                    meas_layer,
                    pauli_randomisation,
                )
                # Calculate the Pauli signs
                experiment_mappings = mapping_set[experiment]
                E = length(experiment_mappings)
                if has_covariance
                    covariance_experiment_mappings = [
                        covariance_dict[covariance_key][1] for
                        covariance_key in covariance_keys[covariance_experiment]
                    ]
                    pair_mappings = [experiment_mappings; covariance_experiment_mappings]
                else
                    pair_mappings = experiment_mappings
                end
                pair_signs = calc_mapping_signs(pair_mappings, randomised_circuit)
                pauli_sign_ensemble[i][j][r] = pair_signs[1:E]
                if has_covariance
                    covariance_pauli_sign_ensemble[i][j][r] = pair_signs[(E + 1):end]
                end
            end
        end
        if diagnostics
            println(
                "Calculated the Pauli signs for tuple $(i) in $(round(time() - tuple_start, digits = 2)) s. The overall time elapsed is $(round(time() - start_time, digits = 2)) s.",
            )
        end
    end
    # Organise the randomised ensemble into jobs
    @assert length(ensemble_indices) == total_experiments "The number of experiments $(length(ensemble_indices)) does not match the expected number $(total_experiments)."
    job_number = cld(total_experiments, job_circuit_number)
    job_indices = Vector{Vector{NTuple{3, Int}}}(undef, job_number)
    for job_idx in 1:job_number
        if job_idx == job_number
            job_circuits = total_experiments - job_circuit_number * (job_number - 1)
        else
            job_circuits = job_circuit_number
        end
        job_indices[job_idx] = Vector{NTuple{3, Int}}(undef, job_circuits)
        for circuit_idx in 1:job_circuits
            total_experiments_idx =
                total_experiments_perm[(job_idx - 1) * job_circuit_number + circuit_idx]
            (i, j, r) = ensemble_indices[total_experiments_idx]
            job_indices[job_idx][circuit_idx] = (i, j, r)
        end
    end
    # Create the randomised design
    overall_time = time() - start_time
    d_rand = RandDesign(
        d,
        pauli_randomisation_ensemble,
        pauli_sign_ensemble,
        covariance_pauli_sign_ensemble,
        job_indices,
        randomisations,
        shot_budget,
        experiment_shots,
        seed,
    )
    if save_data
        save_rand_design(d_rand)
    end
    if diagnostics
        println(
            "Generated the randomised design for all $(tuple_number) tuples in $(round(overall_time, digits = 2)) s.",
        )
    end
    return d_rand::RandDesign
end
function generate_rand_design(
    d::Design,
    min_randomisations::Integer,
    target_shot_budget::Integer,
    experiment_shots::Integer;
    p_norm::Real = 2,
    job_circuit_number::Integer = 200,
    seed::Union{UInt64, Nothing} = nothing,
    diagnostics::Bool = false,
    save_data::Bool = false,
)
    # Set up the number of randomisations
    (randomisations, shot_budget) = get_randomisations(
        d,
        min_randomisations,
        target_shot_budget,
        experiment_shots;
        p_norm = p_norm,
    )
    # Generate the randomised design
    d_rand = generate_rand_design(
        d,
        randomisations,
        shot_budget,
        experiment_shots;
        job_circuit_number = job_circuit_number,
        seed = seed,
        diagnostics = diagnostics,
        save_data = save_data,
    )
    return d_rand::RandDesign
end

"""
    get_design(d_rand::RandDesign)

Returns the design corresponding to the randomised design `d_rand`.
"""
function get_design(d_rand::RandDesign)
    # Set the shot weights according to the randomisations
    d = deepcopy(d_rand.d)
    unnormalised_shot_weights = d_rand.randomisations .* d.experiment_numbers
    @reset d.shot_weights = unnormalised_shot_weights / sum(unnormalised_shot_weights)
    return d::Design
end

"""
    get_meas_budget(d_rand::RandDesign)

Returns the measurement budget corresponding to the randomised design `d_rand`.
"""
function get_meas_budget(d_rand::RandDesign)
    # Get the shot budget
    experiment_shots = d_rand.experiment_shots
    randomisations = d_rand.randomisations
    experiment_numbers = d_rand.d.experiment_numbers
    shot_budget = sum(experiment_shots * randomisations .* experiment_numbers)
    @assert shot_budget == d_rand.shot_budget
    # Get the tuple times factor
    unnormalised_shot_weights = randomisations .* experiment_numbers
    shot_weights = unnormalised_shot_weights / sum(unnormalised_shot_weights)
    tuple_times = d_rand.d.tuple_times
    times_factor = sum(shot_weights .* tuple_times)
    # Return the measurement budget
    meas_budget = shot_budget * times_factor
    return meas_budget::Float64
end

"""
    append_qiskit_gate!(qiskit_circuit, gate_type::String, gate_targets::Vector{Vector{Int}})

Appends the gates of the type `gate_type` to the Qiskit circuit `qiskit_circuit` with the targets specified by `gate_targets`.
"""
function append_qiskit_gate!(
    qiskit_circuit,
    gate_type::String,
    gate_targets::Vector{Vector{Int}},
)
    # Append the gates to the Qiskit circuit
    if gate_type == "CX" || gate_type == "CNOT"
        @assert length(gate_targets) == 2
        qiskit_circuit.cx(gate_targets[1], gate_targets[2])
    elseif gate_type == "CZ"
        @assert length(gate_targets) == 2
        qiskit_circuit.cz(gate_targets[1], gate_targets[2])
    elseif gate_type == "H"
        @assert length(gate_targets) == 1
        qiskit_circuit.h(gate_targets[1])
    elseif gate_type == "S"
        @assert length(gate_targets) == 1
        qiskit_circuit.s(gate_targets[1])
    elseif gate_type == "S_DAG"
        @assert length(gate_targets) == 1
        qiskit_circuit.sdg(gate_targets[1])
    elseif gate_type == "I" || gate_type == "IM"
        @assert length(gate_targets) == 1
        qiskit_circuit.id(gate_targets[1])
    elseif gate_type == "X"
        @assert length(gate_targets) == 1
        qiskit_circuit.x(gate_targets[1])
    elseif gate_type == "Z"
        @assert length(gate_targets) == 1
        qiskit_circuit.z(gate_targets[1])
    elseif gate_type == "Y"
        @assert length(gate_targets) == 1
        qiskit_circuit.y(gate_targets[1])
    elseif gate_type ∈ ["PZ"; "PX"; "PY"]
        @assert length(gate_targets) == 1
        if gate_type == "PZ"
        elseif gate_type == "PX"
            qiskit_circuit.h(gate_targets[1])
        elseif gate_type == "PY"
            qiskit_circuit.h(gate_targets[1])
            qiskit_circuit.s(gate_targets[1])
        else
            throw(error("Invalid preparation gate type $(gate_type)."))
        end
    elseif gate_type ∈ ["MZ"; "MX"; "MY"]
        @assert length(gate_targets) == 1
        if gate_type == "MZ"
        elseif gate_type == "MX"
            qiskit_circuit.h(gate_targets[1])
        elseif gate_type == "MY"
            qiskit_circuit.sdg(gate_targets[1])
            qiskit_circuit.h(gate_targets[1])
        else
            throw(error("Invalid measurement gate type $(gate_type)."))
        end
    elseif gate_type == "M"
        @assert length(gate_targets) == 1
        throw(error("Mid-circuit measurement gates are not currently supported."))
        qiskit_circuit.measure(gate_targets[1], gate_targets[1])
    elseif gate_type == "R"
        @assert length(gate_targets) == 1
        qiskit_circuit.reset(gate_targets[1])
    elseif gate_type == "II"
        @assert length(gate_targets) == 2
        qiskit_circuit.id(gate_targets[1])
        qiskit_circuit.id(gate_targets[2])
    elseif gate_type ∈ ["XX"; "SQRT_XX"; "SQRT_XX_DAG"]
        @assert length(gate_targets) == 2
        if gate_type == "XX"
            rot_angle = π
        elseif gate_type == "SQRT_XX"
            rot_angle = π / 2
        elseif gate_type == "SQRT_XX_DAG"
            rot_angle = -π / 2
        end
        qiskit_circuit.rxx(rot_angle, gate_targets[1], gate_targets[2])
    elseif gate_type ∈ ["ZZ", "SQRT_ZZ", "SQRT_ZZ_DAG"]
        @assert length(gate_targets) == 2
        if gate_type == "ZZ"
            rot_angle = π
        elseif gate_type == "SQRT_ZZ"
            rot_angle = π / 2
        elseif gate_type == "SQRT_ZZ_DAG"
            rot_angle = -π / 2
        end
        qiskit_circuit.rzz(rot_angle, gate_targets[1], gate_targets[2])
    elseif gate_type ∈ ["YY"; "SQRT_YY"; "SQRT_YY_DAG"]
        @assert length(gate_targets) == 2
        if gate_type == "YY"
            rot_angle = π
        elseif gate_type == "SQRT_YY"
            rot_angle = π / 2
        elseif gate_type == "SQRT_YY_DAG"
            rot_angle = -π / 2
        end
        qiskit_circuit.ryy(rot_angle, gate_targets[1], gate_targets[2])
    else
        throw(error("Invalid gate type $(gate_type)."))
    end
    return nothing
end

"""
    append_qiskit!(qiskit_circuit::Py, l::Layer, qiskit_qubit_map::Vector{Int})

Appends the gates in the layer `l` to the Qiskit circuit `qiskit_circuit`, using the mapping to Qiskit qubits given by `qiskit_qubit_map`.
"""
function append_qiskit!(qiskit_circuit::Py, l::Layer, qiskit_qubit_map::Vector{Int})
    # Get the unique gate types in the layer
    gate_types = unique(gate.type for gate in l.layer)
    for gate_type in gate_types
        # Get the targets for the gates of this type
        type_gates = [gate for gate in l.layer if gate.type == gate_type]
        target_num = length(type_gates[1].targets)
        @assert all(length(gate.targets) == target_num for gate in type_gates) "The gates of type $(gate_type) do not have the same number of targets."
        gate_targets = Vector{Vector{Int}}(undef, target_num)
        for target_idx in 1:target_num
            gate_targets[target_idx] = Vector{Int}(undef, length(type_gates))
        end
        for (gate_idx, gate) in pairs(type_gates)
            for target_idx in 1:target_num
                gate_targets[target_idx][gate_idx] =
                    qiskit_qubit_map[gate.targets[target_idx]]
            end
        end
        # Append all gates of this type to the circuit
        append_qiskit_gate!(qiskit_circuit, gate_type, gate_targets)
    end
    return nothing
end

"""
    get_stim_qiskit_circuit(randomised_circuit::Vector{Layer}, gate_probabilities::Dict{Gate, Vector{Float64}}, qiskit_qubit_num::Integer, qiskit_qubit_map::Vector{Int}, noisy_prep::Bool, noisy_meas::Bool; extra_fields::Dict{Symbol, Any} = Dict{Symbol, Any}())

Generate the Stim circuit string and Qiskit circuit for a Pauli frame randomised circuit `randomised_circuit`.

The Stim circuit features noise on each of the non-randomisation gates specified by `gate_probabilities`, as well as noisy preparations if `noisy_prep` is `true`, and noisy measurements if `noisy_meas` is `true`.
Qubit coordinates are specified by `extra_fields` if it contains a [`CodeParameters`](@ref) object.

The Qiskit circuit acts on `qiskit_qubit_num` qubits, which in general must be at least `maximum(qiskit_qubit_map) + 1`, and the mapping of the qubits upon which the circuit acts to the qubits in the Qiskit circuit is controlled by `qiskit_qubit_map`, which is a vector whose length is the number of qubits upon which the circuit acts.
"""
function get_stim_qiskit_circuit(
    randomised_circuit::Vector{Layer},
    gate_probabilities::Dict{Gate, Vector{Float64}},
    qiskit_qubit_num::Integer,
    qiskit_qubit_map::Vector{Int},
    noisy_prep::Bool,
    noisy_meas::Bool;
    extra_fields::Dict{Symbol, Any} = Dict{Symbol, Any}(),
    reset_type::Symbol = :reset,
)
    # Check variables
    @assert reset_type ∈ [:reset; :meas; :meas_reset] "The reset type must be either `:reset`, `:meas`, or `:meas_reset`."
    # Stim orders its one-qubit Paulis as
    # X, Y, Z
    # We order our one-qubit Paulis as
    # X, Z, Y
    # To go from the Stim ordering to our ordering, and vice versa, we have the indexing
    # 1, 3, 2.
    # We add 1 to ignore the normalising probability of I.
    order_1 = [1; 3; 2] .+ 1
    # Stim orders its two-qubit Paulis as
    # IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ.
    # We order our two-qubit Paulis as
    # XI, IX, XX, ZI, YI, ZX, YX, IZ, XZ, IY, XY, ZZ, YZ, ZY, YY.
    # To go from the Stim ordering to our ordering, we have the indexing
    # 4, 1, 5, 12, 8, 13, 9, 3, 7, 2, 6, 15, 11, 14, 10.
    # To go from our ordering to the Stim ordering, we have the indexing
    # 2, 10, 8, 1, 3, 11, 9, 5, 7, 15, 13, 4, 6, 14, 12.
    # Here we need the latter ordering.
    # We add 1 to ignore the normalising probability of II.
    order_2 = [2; 10; 8; 1; 3; 11; 9; 5; 7; 15; 13; 4; 6; 14; 12] .+ 1
    # Measurement and preparation are both only associated with a single error probability
    # This is the second index, as the first index is the probability of no error.
    # Randomised circuits must have 2L+4 layers, where the tuple circuit has L layers
    # The odd numbered layers correspond to the Pauli frame randomisation layers
    # The even numbered layer correspond to the circuit layers, with the firdst being preparation, the next L being the tuple circuit, and the last being measurement
    R = length(randomised_circuit)
    L = convert(Int, (R - 4) // 2)
    @assert R == 2 * L + 4 "The randomised circuit does not have the correct number of layers."
    # Generate the gate list for Stim and the circuit for Qiskit
    @assert all(qiskit_qubit_num >= layer.qubit_num for layer in randomised_circuit) "The Qiskit qubit number $(qiskit_qubit_num) is insufficient for the randomised circuit."
    @assert qiskit_qubit_num >= maximum(qiskit_qubit_map) + 1 "The Qiskit qubit number $(qiskit_qubit_num) is insufficient for the maximum qubit index $(maximum(qiskit_qubit_map)) in the qubit map, noting that Qiskit indexes from 0."
    qiskit_circuit = qiskit.QuantumCircuit(qiskit_qubit_num, qiskit_qubit_num)
    string_vector = Vector{String}(undef, R)
    # Generate the circuit gate layers
    for (i, l) in pairs(randomised_circuit)
        # Append the layer to the Stim circuit
        # Identify the layer type
        is_prep_layer = (i == 2)
        is_gate_layer = (i != 2 && i != R && iseven(i))
        is_meas_layer = (i == R)
        is_randomisation_layer = (isodd(i))
        @assert is_prep_layer || is_gate_layer || is_meas_layer || is_randomisation_layer "Invalid layer index $(i)."
        # Noise channels are only added to the non-random layers and preparation if noisy_prep is true
        k = length(l.layer)
        if (is_prep_layer && noisy_prep) || is_gate_layer
            noise_string_vector = Vector{String}(undef, k)
        end
        if is_gate_layer
            meas_noise_string_vector = Vector{String}(undef, k)
        end
        # Preparation layers are separated into Hadamard and phase preparation gates
        # All other layers are not separated in the same manner
        if is_prep_layer
            hadamard_prep_string_vector = Vector{String}(undef, k)
            phase_prep_string_vector = Vector{String}(undef, k)
        else
            gate_string_vector = Vector{String}(undef, k)
        end
        # Implement the layer of gates
        if is_prep_layer
            for (j, gate) in pairs(l.layer)
                @assert is_state_prep(gate) "Invalid preparation gate type $(gate)."
                # Generate the Stim preparation gate string
                if noisy_prep
                    noise_string_vector[j] = "X_ERROR($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                end
                if gate.type == "PZ"
                    hadamard_prep_string_vector[j] = ""
                    phase_prep_string_vector[j] = ""
                elseif gate.type == "PX"
                    hadamard_prep_string_vector[j] = "H $(gate.targets[1])\n"
                    phase_prep_string_vector[j] = ""
                elseif gate.type == "PY"
                    hadamard_prep_string_vector[j] = "H $(gate.targets[1])\n"
                    phase_prep_string_vector[j] = "S $(gate.targets[1])\n"
                else
                    throw(error("Invalid preparation gate type $(gate.type)."))
                end
            end
        elseif is_meas_layer
            for (j, gate) in pairs(l.layer)
                @assert is_state_meas(gate) "Invalid measurement gate type $(gate)."
                # Generate the Stim measurement gate string
                if noisy_meas
                    gate_string_vector[j] = "$(gate.type)($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                else
                    gate_string_vector[j] = "$(gate.type) $(gate.targets[1])\n"
                end
            end
        else
            for (j, gate) in pairs(l.layer)
                # Generate the Stim gate string
                if is_one_qubit(gate)
                    if is_gate_layer
                        noise_string_vector[j] = ""
                        meas_noise_string_vector[j] = ""
                        if is_mid_meas(gate)
                            gate_string_vector[j] = "$(gate.type)($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                            meas_noise_string_vector[j] = "X_ERROR($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                        elseif is_mid_reset(gate)
                            if reset_type == :reset
                                gate_string_vector[j] = "$(gate.type) $(gate.targets[1])\n"
                            elseif reset_type == :meas || reset_type == :meas_reset
                                # Get the measurement error from the relevant SPAM error if possible
                                spam_meas_gate = Gate("MZ", 0, [gate.targets[1]])
                                spam_prep_gate = Gate("PZ", 0, [gate.targets[1]])
                                if haskey(gate_probabilities, spam_meas_gate)
                                    spam_gate = spam_meas_gate
                                elseif haskey(gate_probabilities, spam_prep_gate)
                                    spam_gate = spam_prep_gate
                                else
                                    spam_gate = gate
                                end
                                gate_string_vector[j] = "$(reset_type == :meas ? "M" : "MR")($(gate_probabilities[spam_gate][2])) $(gate.targets[1])\n"
                            else
                                throw(error("Unsupported reset type $(reset_type)."))
                            end
                            meas_noise_string_vector[j] = "X_ERROR($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                        else
                            noise_string_vector[j] = "PAULI_CHANNEL_1($(join(gate_probabilities[gate][order_1], ", "))) $(gate.targets[1])\n"
                            gate_string_vector[j] = "$((is_meas_idle(gate) ? "I" : gate.type)) $(gate.targets[1])\n"
                        end
                    else
                        @assert is_pauli(gate) "A Pauli frame randomisation layer contains a non-Pauli gate $(gate)."
                        gate_string_vector[j] = "$((is_meas_idle(gate) ? "I" : gate.type)) $(gate.targets[1])\n"
                    end
                elseif is_two_qubit(gate; stim_supported = true)
                    @assert is_gate_layer "A Pauli frame randomisation layer contains a two-qubit gate $(gate)."
                    noise_string_vector[j] = "PAULI_CHANNEL_2($(join(gate_probabilities[gate][order_2], ", "))) $(join(gate.targets, " "))\n"
                    gate_string_vector[j] = "$(gate.type) $(join(gate.targets, " "))\n"
                    meas_noise_string_vector[j] = ""
                else
                    throw(error("Invalid gate $(gate)."))
                end
            end
        end
        # Join and store the Stim gate strings
        if is_prep_layer
            if noisy_prep
                # Noisy preparation layers, separaing Hadamard and phase preparation gates
                string_vector[i] =
                    "TICK\n" *
                    join(noise_string_vector) *
                    join(hadamard_prep_string_vector) *
                    join(phase_prep_string_vector)
            else
                # Preparation layers, separaing Hadamard and phase preparation gates
                string_vector[i] =
                    "TICK\n" *
                    join(hadamard_prep_string_vector) *
                    join(phase_prep_string_vector)
            end
        elseif is_gate_layer
            # Noise is only added to the non-randomised layers
            string_vector[i] =
                "TICK\n" *
                join(noise_string_vector) *
                join(gate_string_vector) *
                join(meas_noise_string_vector)

        else
            string_vector[i] = "TICK\n" * join(gate_string_vector)
        end
        # Append the layer to the Qiskit circuit
        append_qiskit!(qiskit_circuit, l, qiskit_qubit_map)
        # Apply barriers between the Pauli frame randomised layer triples, and measurement at the end of the circuit
        if is_prep_layer || is_gate_layer
            qiskit_circuit.barrier(qiskit_qubit_map)
        elseif is_meas_layer
            qiskit_circuit.barrier(qiskit_qubit_map)
            qiskit_circuit.measure(qiskit_qubit_map, qiskit_qubit_map)
        end
    end
    # Generate the qubit coordinates
    if haskey(extra_fields, :code_param)
        code_param = extra_fields[:code_param]
        @assert length(code_param.qubits) == length(qiskit_qubit_map) "The qubit map does not have the correct number of qubits."
        insert!(string_vector, 1, get_stim_qubits(code_param))
    end
    # Generate the Stim circuit string
    stim_circuit_string = join(string_vector)
    return (stim_circuit_string::String, qiskit_circuit::Py)
end

"""
    get_stim_qiskit_ensemble(d_rand::RandDesign, qiskit_qubit_num::Integer, qiskit_qubit_map::Vector{Int})

Returns the Stim and Qiskit circuits for the randomised ensemble `d_rand`.

The Qiskit circuits act on `qiskit_qubit_num` qubits, which in general must be at least `maximum(qiskit_qubit_map) + 1`, and the mapping of the qubits upon which the circuits act to the qubits in the Qiskit circuits is controlled by `qiskit_qubit_map`, which is a vector whose length is the number of qubits upon which the circuits act.
"""
function get_stim_qiskit_ensemble(
    d_rand::RandDesign,
    qiskit_qubit_num::Integer,
    qiskit_qubit_map::Vector{Int},
)
    # Set up variables
    d = d_rand.d
    gate_probabilities = d.c.gate_probabilities
    noisy_prep = d.c.noisy_prep
    noisy_meas = d.c.noisy_meas
    extra_fields = d.c.extra_fields
    # Generate the Stim and Qiskit circuits for the randomised ensemble
    job_number = length(d_rand.job_indices)
    job_circuit_numbers = length.(d_rand.job_indices)
    stim_ensemble = Vector{Vector{String}}(undef, job_number)
    qiskit_ensemble = pylist([])
    for job_idx in 1:job_number
        stim_ensemble[job_idx] = Vector{String}(undef, job_circuit_numbers[job_idx])
        qiskit_experiment_set = pylist([])
        for circuit_idx in 1:job_circuit_numbers[job_idx]
            # Get the randomised circuit
            (i, j, r) = d_rand.job_indices[job_idx][circuit_idx]
            circuit_tuple = d.tuple_set[i]
            tuple_circuit = d.c.circuit[circuit_tuple]
            prep_layer = d.prep_ensemble[i][j]
            meas_layer = d.meas_ensemble[i][j]
            pauli_randomisation = d_rand.pauli_randomisation_ensemble[i][j][r]
            randomised_circuit = get_randomised_circuit(
                prep_layer,
                tuple_circuit,
                meas_layer,
                pauli_randomisation,
            )
            # Get the Stim and Qiskit circuits
            (stim_string, qiskit_circuit) = get_stim_qiskit_circuit(
                randomised_circuit,
                gate_probabilities,
                qiskit_qubit_num,
                qiskit_qubit_map,
                noisy_prep,
                noisy_meas;
                extra_fields = extra_fields,
            )
            stim_ensemble[job_idx][circuit_idx] = stim_string
            qiskit_experiment_set.append(qiskit_circuit)
        end
        qiskit_ensemble.append(qiskit_experiment_set)
    end
    return (stim_ensemble::Vector{Vector{String}}, qiskit_ensemble::Py)
end

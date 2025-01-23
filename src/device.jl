"""
    simulate_stim_ensemble(d_rand::RandDesign, stim_ensemble::Vector{Vector{String}}, experiment_shots::Integer; seed::Union{UInt64, Nothing} = nothing, diagnostics::Bool = true)

Saves data simulated with Stim for the randomised experimental design `d_rand` with an associated ensemble of Stim circuits `stim_ensemble` produced by [`get_stim_qiskit_ensemble`](@ref), with `experiment_shots` shots for each circuit.
The simulation uses the random seed `seed` and displays diagnostics if `diagnostics` is `true`.
"""
function simulate_stim_ensemble(
    d_rand::RandDesign,
    stim_ensemble::Vector{Vector{String}},
    experiment_shots::Integer;
    seed::Union{UInt64, Nothing} = nothing,
    diagnostics::Bool = true,
)
    # Initialise variables
    start_time = time()
    job_number = length(d_rand.job_indices)
    job_circuit_numbers = length.(d_rand.job_indices)
    @assert (job_number == length(stim_ensemble)) &&
            (job_circuit_numbers == length.(stim_ensemble)) "The number of jobs and circuits in each job in the randomised design are not consistent with the supplied Stim ensemble."
    # Generate the requisite random seeds for Stim
    if seed !== nothing
        Random.seed!(seed)
        backend = "stim_$(seed)"
    else
        backend = "stim"
    end
    stim_seeds = Vector{Vector{UInt64}}(undef, job_number)
    for job_idx in 1:job_number
        stim_seeds[job_idx] = rand(UInt64, job_circuit_numbers[job_idx])
    end
    Random.seed!()
    # Simulate the circuits with Stim
    for job_idx in 1:job_number
        job_counts = Vector{Matrix{UInt8}}(undef, job_circuit_numbers[job_idx])
        for circuit_idx in 1:job_circuit_numbers[job_idx]
            job_counts[circuit_idx] = stim_sample(
                stim_ensemble[job_idx][circuit_idx],
                experiment_shots;
                stim_seed = stim_seeds[job_idx][circuit_idx],
            )
        end
        save_rand_design_job(d_rand, backend, job_counts, job_idx)
        if diagnostics
            println(
                "Simulated sampling job $(job_idx) of $(job_number). The time elapsed since simulation started is $(round(time() - start_time, digits = 3)) s.",
            )
        end
    end
    return nothing
end

"""
    simulate_qiskit_ensemble(d_rand::RandDesign, qiskit_ensemble::Py, experiment_shots::Integer; backend::String = "qiskit", prefix::String = "result", simulator::Py = aer.AerSimulator(; method = "stabilizer"), diagnostics::Bool = true)

Saves data simulated with Qiskit for the randomised experimental design `d_rand` with an associated ensemble of Qiskit circuits `qiskit_ensemble` produced by [`get_stim_qiskit_ensemble`](@ref), with `experiment_shots` shots for each circuit.
The saved data is labelled according to the supplied `backend` and `prefix`, and simulated using the supplied `simulator`, displaying diagnostics if `diagnostics` is `true.
"""
function simulate_qiskit_ensemble(
    d_rand::RandDesign,
    qiskit_ensemble::Py,
    experiment_shots::Integer;
    backend::String = "qiskit",
    prefix::String = "result",
    simulator::Py = aer.AerSimulator(; method = "stabilizer"),
    diagnostics::Bool = true,
)
    # Check that the relevant folders are present
    foldername = backend * "_" * rand_design_filename(d_rand)
    @assert foldername[(end - 4):end] == ".jld2"
    foldername = foldername[1:(end - 5)]
    if ~isdir("data")
        mkdir("data")
    end
    if ~isdir("data/$(foldername)")
        mkdir("data/$(foldername)")
    end
    # Simulate with the Qiskit stabiliser simulator
    start_time = time()
    job_number = length(d_rand.job_indices)
    for job_idx in 1:job_number
        # Simulate the job with Qiskit
        py_job_result =
            simulator.run(qiskit_ensemble[job_idx - 1]; shots = experiment_shots).result()
        # Save the Qiskit job data
        filename = "$(foldername)/$(prefix)_$(job_idx).pickle"
        open(pwd() * "/data/" * filename, "w") do f
            pickle.dump(py_job_result, f)
        end
        if diagnostics
            println(
                "Simulated sampling job $(job_idx) of $(job_number). The time elapsed since simulation started is $(round(time() - start_time, digits = 3)) s.",
            )
        end
    end
    return nothing
end

"""
    process_qiskit_ensemble(d_rand::RandDesign, qiskit_qubit_num::Integer, qiskit_qubit_map::Vector{Int}, experiment_shots::Integer; backend::String = "qiskit", prefix::String = "result")

Saves a processed version of saved Qiskit data for the randomised experimental design `d_rand`, produced for example by [`simulate_qiskit_ensemble`](@ref), whose associated ensemble of Qiskit circuits act on `qiskit_qubit_num` qubits and are indexed by the mapping `qiskit_qubit_map`, and which were run with `experiment_shots` shots for each circuit.
The saved data is loaded according to the labelling implied by the supplied `backend` and `prefix`.
"""
function process_qiskit_ensemble(
    d_rand::RandDesign,
    qiskit_qubit_num::Integer,
    qiskit_qubit_map::Vector{Int},
    experiment_shots::Integer;
    backend::String = "qiskit",
    prefix::String = "result",
)
    # Initialise filenaming
    foldername = backend * "_" * rand_design_filename(d_rand)
    @assert foldername[(end - 4):end] == ".jld2"
    foldername = foldername[1:(end - 5)]
    # Convert the Qiskit data into the Stim counts format
    job_number = length(d_rand.job_indices)
    for job_idx in 1:job_number
        # Load the Qiskit job data
        filename = "$(foldername)/$(prefix)_$(job_idx).pickle"
        py_job_result = open(pwd() * "/data/" * filename, "r") do f
            pickle.load(f)
        end
        # Format the data in a Julia-compatible format
        qiskit_job_counts = Vector{Dict{String, Int}}()
        try
            qiskit_job_counts =
                pyconvert(Vector{Dict{String, Int}}, py_job_result.get_counts())
        catch
            qiskit_job_counts = [
                pyconvert(Dict{String, Int}, py_circuit_result.get_counts()) for
                py_circuit_result in py_job_result
            ]
        end
        job_counts = [
            parse_qiskit_dict(
                counts,
                experiment_shots,
                qiskit_qubit_num;
                qiskit_qubit_map = qiskit_qubit_map,
            ) for counts in qiskit_job_counts
        ]
        # Save the data
        save_rand_design_job(d_rand, backend, job_counts, job_idx)
    end
    return nothing
end

"""
    estimate_stim_ensemble(d_rand::RandDesign, experiment_shots::Integer; fail_jobs::Vector{Int} = Int[], simulation_seed::Union{UInt64, Nothing} = nothing, ls_type::Symbol = :fgls, split::Bool = false, min_eigenvalue::Real = 0.1)

Returns estimated circuit eigenvalues and gate noise estimates from simulated Stim data for the randomised experimental design `d_rand` with `experiment_shots` shots per circuit produced by [`simulate_stim_ensemble`](@ref).
Avoids attempting to load failed jobs specified by `fail_jobs`.
The simulated data is loaded according to the label implied by the supplied `simulation_seed`, and the estimation procedure uses the least squares estimator `ls_type`, splitting the gate eigenvalue estimator projection across the gates if `split` is `true`, and includes only circuit eigenvalues which are at least `min_eigenvalue`.
"""
function estimate_stim_ensemble(
    d_rand::RandDesign,
    experiment_shots::Integer;
    fail_jobs::Vector{Int} = Int[],
    simulation_seed::Union{UInt64, Nothing} = nothing,
    min_eigenvalue::Real = 0.1,
    split::Bool = false,
    precision::Real = 1e-8,
)
    # Initialise variables
    d = d_rand.d
    tuple_number = length(d.tuple_set)
    if simulation_seed !== nothing
        backend = "stim_$(simulation_seed)"
    else
        backend = "stim"
    end
    (experiment_ensemble, meas_indices_ensemble) = get_experiment_data(d)[1:2]
    pauli_sign_ensemble = d_rand.pauli_sign_ensemble
    mapping_lengths = length.(meas_indices_ensemble)
    @assert mapping_lengths == length.(d.mapping_ensemble)
    (covariance_experiment_ensemble, covariance_meas_indices_ensemble) =
        get_covariance_experiment_data(d)[1:2]
    covariance_pauli_sign_ensemble = d_rand.covariance_pauli_sign_ensemble
    covariance_mapping_lengths = length.(covariance_meas_indices_ensemble)
    # Estimate the eigenvalues
    job_indices = d_rand.job_indices
    job_number = length(job_indices)
    job_circuit_numbers = length.(job_indices)
    @assert all(fail_jobs .> 1) && all(fail_jobs .<= job_number) "The failed job indices are not consistent with the number of jobs."
    success_jobs = setdiff(1:job_number, fail_jobs)
    est_eigenvalues_experiment_ensemble =
        [[Float64[] for j in 1:mapping_lengths[i]] for i in 1:tuple_number]
    est_covariance_experiment_ensemble =
        [[Float64[] for j in 1:covariance_mapping_lengths[i]] for i in 1:tuple_number]
    for job_idx in success_jobs
        job_counts = load_rand_design_job(d_rand, backend, job_idx)
        @assert length(job_counts) == job_circuit_numbers[job_idx] "The number of circuits in the job does not match the number of counts."
        for circuit_idx in 1:job_circuit_numbers[job_idx]
            # Obtain the job indices
            (i, j, r) = job_indices[job_idx][circuit_idx]
            # Use the data to contribute to the eigenvalue estimates
            experiment = experiment_ensemble[i][j]
            experiment_meas_indices = meas_indices_ensemble[i][experiment]
            experiment_pauli_signs = pauli_sign_ensemble[i][j][r]
            est_eigenvalues_experiment = estimate_experiment(
                job_counts[circuit_idx],
                experiment_meas_indices,
                experiment_pauli_signs,
                experiment_shots,
            )
            for (idx, pauli_idx) in pairs(experiment)
                push!(
                    est_eigenvalues_experiment_ensemble[i][pauli_idx],
                    est_eigenvalues_experiment[idx],
                )
            end
            has_covariance = (covariance_mapping_lengths[i] > 0)
            if has_covariance
                covariance_experiment = covariance_experiment_ensemble[i][j]
                covariance_experiment_meas_indices =
                    covariance_meas_indices_ensemble[i][covariance_experiment]
                covariance_experiment_pauli_signs = covariance_pauli_sign_ensemble[i][j][r]
                est_covariance_eigenvalues_experiment = estimate_experiment(
                    job_counts[circuit_idx],
                    covariance_experiment_meas_indices,
                    covariance_experiment_pauli_signs,
                    experiment_shots,
                )
                for (idx, pauli_idx) in pairs(covariance_experiment)
                    push!(
                        est_covariance_experiment_ensemble[i][pauli_idx],
                        est_covariance_eigenvalues_experiment[idx],
                    )
                end
            end
        end
    end
    # Get the design and estimate the noise
    noise_est = estimate_gate_noise(
        d_rand,
        est_eigenvalues_experiment_ensemble,
        est_covariance_experiment_ensemble;
        min_eigenvalue = min_eigenvalue,
        split = split,
        precision = precision,
    )
    return noise_est::NoiseEstimate
end

"""
    estimate_qiskit_ensemble(d_rand::RandDesign, qiskit_qubit_map::Vector{Int}, experiment_shots::Integer; fail_jobs::Vector{Int} = Int[], backend::String = "qiskit", ls_type::Symbol = :fgls, split::Bool = false, min_eigenvalue::Real = 0.1)

Returns estimated circuit eigenvalues, and gate noise estimates from Qiskit data processed with [`process_qiskit_ensemble`] for the randomised experimental design `d_rand` with `experiment_shots` shots per circuit.
Avoids attempting to load failed jobs specified by `fail_jobs`.
The simulated data is loaded according to the label implied by the supplied `backend`, and the estimation procedure uses the least squares estimator `ls_type`, splitting the gate eigenvalue estimator projection across the gates if `split` is `true`, and includes only circuit eigenvalues which are at least `min_eigenvalue`.
"""
function estimate_qiskit_ensemble(
    d_rand::RandDesign,
    qiskit_qubit_map::Vector{Int},
    experiment_shots::Integer;
    fail_jobs::Vector{Int} = Int[],
    backend::String = "qiskit",
    min_eigenvalue::Real = 0.1,
    split::Bool = false,
    precision::Real = 1e-8,
)
    # Initialise variables
    d = d_rand.d
    tuple_number = length(d.tuple_set)
    (experiment_ensemble, meas_indices_ensemble) =
        get_experiment_data(d; qiskit_qubit_map = qiskit_qubit_map)[1:2]
    pauli_sign_ensemble = d_rand.pauli_sign_ensemble
    mapping_lengths = length.(meas_indices_ensemble)
    @assert mapping_lengths == length.(d.mapping_ensemble)
    (covariance_experiment_ensemble, covariance_meas_indices_ensemble) =
        get_covariance_experiment_data(d; qiskit_qubit_map = qiskit_qubit_map)[1:2]
    covariance_pauli_sign_ensemble = d_rand.covariance_pauli_sign_ensemble
    covariance_mapping_lengths = length.(covariance_meas_indices_ensemble)
    # Estimate the eigenvalues
    job_indices = d_rand.job_indices
    job_number = length(job_indices)
    job_circuit_numbers = length.(job_indices)
    @assert all(fail_jobs .> 1) && all(fail_jobs .<= job_number) "The failed job indices are not consistent with the number of jobs."
    success_jobs = setdiff(1:job_number, fail_jobs)
    est_eigenvalues_experiment_ensemble =
        [[Float64[] for j in 1:mapping_lengths[i]] for i in 1:tuple_number]
    est_covariance_experiment_ensemble =
        [[Float64[] for j in 1:covariance_mapping_lengths[i]] for i in 1:tuple_number]
    for job_idx in success_jobs
        job_counts = load_rand_design_job(d_rand, backend, job_idx)
        @assert length(job_counts) == job_circuit_numbers[job_idx] "The number of circuits in the job does not match the number of counts."
        for circuit_idx in 1:job_circuit_numbers[job_idx]
            # Obtain the job indices
            (i, j, r) = job_indices[job_idx][circuit_idx]
            # Use the data to contribute to the eigenvalue estimates
            experiment = experiment_ensemble[i][j]
            experiment_meas_indices = meas_indices_ensemble[i][experiment]
            experiment_pauli_signs = pauli_sign_ensemble[i][j][r]
            est_eigenvalues_experiment = estimate_experiment(
                job_counts[circuit_idx],
                experiment_meas_indices,
                experiment_pauli_signs,
                experiment_shots,
            )
            for (idx, pauli_idx) in pairs(experiment)
                push!(
                    est_eigenvalues_experiment_ensemble[i][pauli_idx],
                    est_eigenvalues_experiment[idx],
                )
            end
            has_covariance = (covariance_mapping_lengths[i] > 0)
            if has_covariance
                covariance_experiment = covariance_experiment_ensemble[i][j]
                covariance_experiment_meas_indices =
                    covariance_meas_indices_ensemble[i][covariance_experiment]
                covariance_experiment_pauli_signs = covariance_pauli_sign_ensemble[i][j][r]
                est_covariance_eigenvalues_experiment = estimate_experiment(
                    job_counts[circuit_idx],
                    covariance_experiment_meas_indices,
                    covariance_experiment_pauli_signs,
                    experiment_shots,
                )
                for (idx, pauli_idx) in pairs(covariance_experiment)
                    push!(
                        est_covariance_experiment_ensemble[i][pauli_idx],
                        est_covariance_eigenvalues_experiment[idx],
                    )
                end
            end
        end
    end
    # Get the design and estimate the noise
    noise_est = estimate_gate_noise(
        d_rand,
        est_eigenvalues_experiment_ensemble,
        est_covariance_experiment_ensemble;
        min_eigenvalue = min_eigenvalue,
        split = split,
        precision = precision,
    )
    return noise_est::NoiseEstimate
end

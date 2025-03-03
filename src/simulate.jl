"""
    ACESData

ACES noise characterisation experiment simulation results.

# Fields

  - `d::Design`: Experimental design.
  - `noise_est_coll::Matrix{NoiseEstimate}`: Noise estimates for each repetition and measurement budget.
  - `noise_error_coll::Matrix{NoiseError}`: Noise errors for each repetition and measurement budget.
  - `budget_set::Vector{Int}`: Measurement budgets.
  - `shots_set::Vector{Int}`: Measurement shots corresponding to the measurement budgets in `budget_set`.
  - `repetitions::Int`: Repetitions of simulating an ACES noise characterisation experiment.
  - `seed::UInt64`: Seed used to generate the random seeds for each repetition.
  - `seeds::Vector{UInt64}`: Seeds for each of the repetitions.
  - `calculation_times::Matrix{Float64}`: Time taken to simulate sampling, and to estimate the gate eigenvalues and project them into the probability simplex, for each repetition.
  - `overall_time::Float64`: Overall time taken to simulate ACES across all repetitions.
"""
struct ACESData
    d::Design
    noise_est_coll::Matrix{NoiseEstimate}
    noise_error_coll::Matrix{NoiseError}
    budget_set::Vector{Int}
    shots_set::Vector{Int}
    repetitions::Int
    seed::UInt64
    seeds::Vector{UInt64}
    calculation_times::Matrix{Float64}
    overall_time::Float64
end

function Base.show(io::IO, a::ACESData)
    print(
        io,
        "ACES data from a design for a $(a.d.c.circuit_param.circuit_name) circuit with $(length(a.d.tuple_set)) tuples and $(a.d.experiment_number) experiments.",
    )
    return nothing
end

@struct_hash_equal_isequal ACESData

"""
    stim_sample(stim_circuit_string::String, shots::Integer; stim_seed::Union{UInt64, Nothing} = nothing)

Returns bit-packaged measurement outcomes from simulating the circuit `stim_circuit_string` over `shots` measurement shots using the Python package Stim.
While the seed `stim_seed` can be fixed, Stim guarantees inconsistency when using the same seed on different versions.
"""
function stim_sample(
    stim_circuit_string::String,
    shots::Integer;
    stim_seed::Union{UInt64, Nothing} = nothing,
)
    # As of v0.9.22, PythonCall is capable of running Python on thread 1 without segfaulting due to garbage collection issues
    @assert Threads.threadid() == 1 "This function cannot be multithreaded."
    # Multithreaded Python is also possible but requires messing with the Python GIL
    # Generate the Stim sampler
    stim_circuit_sampler =
        stim.Circuit(stim_circuit_string).compile_sampler(; seed = stim_seed)
    # Sample the shots
    stim_samples = pyconvert(
        Matrix{UInt8},
        stim_circuit_sampler.sample(; shots = shots, bit_packed = true),
    )
    return stim_samples::Matrix{UInt8}
end

"""
    simulate_stim_estimate(d::Design, shots_set::Vector{Int}; kwargs...)

Returns simulated estimated circuit eigenvalues, covariance circuit eigenvalues, and corresponding circuit eigenvalue pairs, for the experimental design `d` across each of the measurement shots in `shots_set`.

# Keyword arguments

  - `seed::Union{UInt64, Nothing} = nothing`: Seed controlling the random seeds for Stim.
  - `max_samples::Integer = 10^9`: Maximum number of samples to take in a single Stim simulation.
  - `detailed_diagnostics::Bool = false`: Whether to print diagnostic information about the simulation.
"""
function simulate_stim_estimate(
    d::Design,
    shots_set::Vector{Int};
    seed::Union{UInt64, Nothing} = nothing,
    max_samples::Integer = 10^9,
    detailed_diagnostics::Bool = false,
)
    # Initialise variables
    start_time = time()
    tuple_number = length(d.tuple_set)
    experiment_numbers = d.experiment_numbers
    n = d.c.qubit_num
    gate_probabilities = d.c.gate_probabilities
    noisy_prep = d.c.noisy_prep
    noisy_meas = d.c.noisy_meas
    extra_fields = d.c.extra_fields
    budget_count = length(shots_set)
    # Divide the shots between the experiments
    @assert sum(d.shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
    @assert all(d.shot_weights .> 0.0) "The shot weights are not all positive."
    shots_divided_float = [
        shots_set[s] * d.shot_weights[t] / experiment_numbers[t] for
        t in 1:tuple_number, s in 1:budget_count
    ]
    @assert vec(sum(shots_divided_float .* experiment_numbers; dims = 1)) ≈ shots_set "The shots have not been divided correctly."
    shots_divided = ceil.(Int, shots_divided_float)
    shots_maximum = vec(maximum(shots_divided; dims = 2))
    # Generate the requisite random seeds for Stim
    if seed !== nothing
        Random.seed!(seed)
    end
    batched_shots = any(shots_maximum * n .> max_samples)
    if batched_shots
        stim_seeds = [
            [
                rand(UInt64, length(batch_shots(shots_maximum[i], n, max_samples))) for
                j in 1:experiment_numbers[i]
            ] for i in 1:tuple_number
        ]
    else
        stim_seeds = [rand(UInt64, experiment_numbers[i]) for i in 1:tuple_number]
    end
    if seed !== nothing
        Random.seed!()
    end
    # Initialise data
    (experiment_ensemble, meas_indices_ensemble, pauli_sign_ensemble) =
        get_experiment_data(d)
    mapping_lengths = length.(meas_indices_ensemble)
    @assert mapping_lengths == length.(pauli_sign_ensemble) &&
            mapping_lengths == length.(d.mapping_ensemble)
    (
        covariance_experiment_ensemble,
        covariance_meas_indices_ensemble,
        covariance_pauli_sign_ensemble,
    ) = get_covariance_experiment_data(d)
    covariance_mapping_lengths = length.(covariance_meas_indices_ensemble)
    @assert covariance_mapping_lengths == length.(covariance_pauli_sign_ensemble)
    # Estimate the eigenvalues
    est_eigenvalues_experiment_ensemble_set = [
        [[Float64[] for j in 1:mapping_lengths[i]] for i in 1:tuple_number] for
        s in 1:budget_count
    ]
    est_covariance_experiment_ensemble_set = [
        [[Float64[] for j in 1:covariance_mapping_lengths[i]] for i in 1:tuple_number]
        for s in 1:budget_count
    ]
    for i in 1:tuple_number
        # Initialise tuple circuit
        circuit_tuple = d.tuple_set[i]
        tuple_circuit_string = get_stim_circuit(
            d.c.circuit[circuit_tuple],
            gate_probabilities,
            noisy_prep,
            noisy_meas,
        )
        has_covariance = (covariance_mapping_lengths[i] > 0)
        for j in 1:experiment_numbers[i]
            # Initialise preparation and measurement layers
            prep_layer_string = get_stim_circuit(
                [d.prep_ensemble[i][j]],
                gate_probabilities,
                noisy_prep,
                noisy_meas;
                extra_fields = extra_fields,
            )
            meas_layer_string = get_stim_circuit(
                [d.meas_ensemble[i][j]],
                gate_probabilities,
                noisy_prep,
                noisy_meas,
            )
            # Generate the circuit
            stim_circuit_string =
                prep_layer_string * tuple_circuit_string * meas_layer_string
            # Sample the circuit with Stim
            if batched_shots
                # If there are too many shots to sample at once, divide them into batches
                shot_batches = batch_shots(shots_maximum[i], n, max_samples)
                batches = length(shot_batches)
                stim_counts_coll = Vector{Matrix{UInt8}}(undef, batches)
                for b in 1:batches
                    stim_counts_coll[b] = stim_sample(
                        stim_circuit_string,
                        shot_batches[b];
                        stim_seed = stim_seeds[i][j][b],
                    )
                end
                stim_counts = vcat(stim_counts_coll...)
            else
                # Sample all of the shots at once
                stim_counts = stim_sample(
                    stim_circuit_string,
                    shots_maximum[i];
                    stim_seed = stim_seeds[i][j],
                )
            end
            # Estimate the eigenvalues and covariance eigenvalues from the counts data
            experiment = experiment_ensemble[i][j]
            experiment_meas_indices = meas_indices_ensemble[i][experiment]
            experiment_pauli_signs = pauli_sign_ensemble[i][experiment]
            if has_covariance
                covariance_experiment = covariance_experiment_ensemble[i][j]
                covariance_experiment_meas_indices =
                    covariance_meas_indices_ensemble[i][covariance_experiment]
                covariance_experiment_pauli_signs =
                    covariance_pauli_sign_ensemble[i][covariance_experiment]
            end
            for s in 1:budget_count
                shots_value = shots_divided[i, s]
                est_eigenvalues_experiment = estimate_experiment(
                    stim_counts,
                    experiment_meas_indices,
                    experiment_pauli_signs,
                    shots_value,
                )
                for (idx, pauli_idx) in pairs(experiment)
                    push!(
                        est_eigenvalues_experiment_ensemble_set[s][i][pauli_idx],
                        est_eigenvalues_experiment[idx],
                    )
                end
                if has_covariance
                    est_covariance_eigenvalues_experiment = estimate_experiment(
                        stim_counts,
                        covariance_experiment_meas_indices,
                        covariance_experiment_pauli_signs,
                        shots_value,
                    )
                    for (idx, pauli_idx) in pairs(covariance_experiment)
                        push!(
                            est_covariance_experiment_ensemble_set[s][i][pauli_idx],
                            est_covariance_eigenvalues_experiment[idx],
                        )
                    end
                end
            end
            if detailed_diagnostics
                println(
                    "Simulated sampling experiment $(j) of $(experiment_numbers[i]) for tuple $(i) of $(tuple_number). The time elapsed since simulation started is $(round(time() - start_time, digits = 3)) s.",
                )
            end
        end
    end
    return (
        est_eigenvalues_experiment_ensemble_set::Vector{Vector{Vector{Vector{Float64}}}},
        est_covariance_experiment_ensemble_set::Vector{Vector{Vector{Vector{Float64}}}},
    )
end

"""
    simulate_aces(d::Design, budget_set::Vector{Int}; kwargs...)

Returns simulated ACES noise characterisation experiment results as an [`ACESData`](@ref) object for the experimental design `d` across each of the measurement budgets in `budget_set`.

WARNING: Seeding has the same features as in Stim.
The behaviour of the same random seed will differ across different versions of Stim.
Also, when measurement shots are sampled in batches, which occurs when `max_samples` is exceeded, the results will differ from when all shots are sampled at once.

# Arguments

  - `d::Design`: Experimental design.
  - `budget_set::Vector{Int}`: Measurement budgets for which to simulate ACES.

# Keyword arguments

  - `repetitions::Integer = 1`: Number of simulation repetitions.
  - `seed::Union{UInt64, Nothing} = nothing`: Seed to use for random number generation.
  - `N_warn::Integer = 10^5`: Number of circuit eigenvalues above which to warn the user about certain keyword argument choices.
  - `max_samples::Integer = 10^9`: Maximum number of Stim samples collected in a single simulation.
  - `min_eigenvalue::Real = 0.1`: Circuit eigenvalues below this threshold are omitted from the estimation procedure.
  - `clip_warning::Bool = false`: Whether to warn the user about clipped eigenvalues.
  - `N_warn_split::Integer = 5 * 10^3`: Number of circuit eigenvalues above which to warn the user about if `split` is false.
  - `split::Bool = (d.c.gate_data.N < N_warn_split ? false : true)`: Whether to split the gate eigenvalue projection across gates or simultaneously project all gate eigenvalues.
  - `precision::Real = 1e-8`: Precision of the solver for projecting noise estimates into the probability simplex.
  - `diagnostics::Bool = true`: Whether to print diagnostics.
  - `detailed_diagnostics::Bool = false`: Whether to print detailed diagnostics for Stim simulations.
  - `save_data::Bool = false`: Whether to save the data.
  - `save_interval::Integer = 50`: Repetition interval at which to save the data.
  - `clear_design::Bool = false`: Whether to clear the saved design data after saving the full simulation data.
"""
function simulate_aces(
    d::Design,
    budget_set::Vector{Int};
    repetitions::Integer = 1,
    seed::Union{UInt64, Nothing} = nothing,
    N_warn::Integer = 10^5,
    max_samples::Integer = 10^9,
    min_eigenvalue::Real = 0.1,
    clip_warning::Bool = false,
    N_warn_split::Integer = 5 * 10^3,
    split::Bool = (d.c.gate_data.N < N_warn_split ? false : true),
    precision::Real = 1e-8,
    diagnostics::Bool = true,
    detailed_diagnostics::Bool = false,
    save_data::Bool = false,
    save_interval::Integer = 50,
    clear_design::Bool = false,
)
    # Warn the user if they have inadvisable settings for a large circuit
    N = d.c.gate_data.N
    @assert N == size(d.matrix, 2)
    if N >= N_warn_split
        if ~split
            @warn "This ACES simulation is for a large circuit: splitting the gate eigenvalue projection across gates is advised."
        end
    end
    if N >= N_warn
        if ~diagnostics
            @warn "This ACES simulation is for a very large circuit: turning on diagnostics is advised."
        end
        if ~save_data
            @warn "This ACES simulation is for a very large circuit: saving the data is advised."
        end
    end
    # Normalise the sampled shot count by the amount of time taken to perform the circuits
    tuple_number = length(d.tuple_set)
    budget_count = length(budget_set)
    times_factor = sum(d.shot_weights .* d.tuple_times)
    shots_set = round.(Int, budget_set / times_factor)
    # Generate the requisite random seeds using the fixed seed
    if seed === nothing
        seed = rand(UInt64)
    end
    Random.seed!(seed)
    seeds = rand(UInt64, repetitions)
    Random.seed!()
    # Initialise data storage
    noise_est_coll = Matrix{NoiseEstimate}(undef, repetitions, budget_count)
    noise_error_coll = Matrix{NoiseError}(undef, repetitions, budget_count)
    # Time data
    calculation_times = Matrix{Float64}(undef, repetitions, 2)
    overall_time = 0.0
    # Load saved data
    saved_repetitions = 0
    filename = aces_data_filename(d, budget_set, seed)
    if isfile(pwd() * "/data/" * filename)
        saved_aces_data = load_aces(d, budget_set, seed)
        if saved_aces_data.repetitions >= repetitions
            if saved_aces_data.repetitions > repetitions
                @warn "The saved data has more repetitions $(saved_aces_data.repetitions) than the requested number of repetitions $(repetitions)."
            elseif diagnostics
                println(
                    "The saved data already contains the requested number $(repetitions) of simulated repetitions of ACES.",
                )
            end
            return saved_aces_data::ACESData
        end
        d_saved = saved_aces_data.d
        same_data =
            (d_saved.matrix == d.matrix) &&
            (d_saved.tuple_set == d.tuple_set) &&
            (d_saved.tuple_set_data == d.tuple_set_data) &&
            (d_saved.shot_weights ≈ d.shot_weights) &&
            (saved_aces_data.seed == seed) &&
            (saved_aces_data.seeds == seeds[1:(saved_aces_data.repetitions)])
        # Only used the saved data if it has the same parameters
        if same_data
            saved_repetitions = saved_aces_data.repetitions
            noise_est_coll[1:saved_repetitions, :] = saved_aces_data.noise_est_coll
            noise_error_coll[1:saved_repetitions, :] = saved_aces_data.noise_error_coll
            calculation_times[1:saved_repetitions, :] = saved_aces_data.calculation_times
            overall_time = saved_aces_data.overall_time
            if diagnostics
                println(
                    "Loading $(saved_repetitions) of $(repetitions) repetitions of simulated ACES data.",
                )
            end
        else
            @warn "Attempted and failed to load and reuse the saved data; its parameters differ from the supplied data."
        end
    end
    repetitions_start = time()
    for idx in (saved_repetitions + 1):repetitions
        # Simulate ACES circuits and process the data to estimate the circuit eigenvalues
        simulate_start = time()
        (est_eigenvalues_experiment_ensemble_set, est_covariance_experiment_ensemble_set) =
            simulate_stim_estimate(
                d,
                shots_set;
                seed = seeds[idx],
                max_samples = max_samples,
                detailed_diagnostics = detailed_diagnostics,
            )
        if N >= N_warn
            GC.gc()
        end
        simulate_time = time() - simulate_start
        if diagnostics
            println(
                "Sampling $(round(maximum(shots_set), sigdigits = 4)) shots, corresponding to the measurement budget $(round(maximum(budget_set), sigdigits = 4)), divided between $(d.experiment_number) experiments for the $(tuple_number) tuples, and then estimating the circuit eigenvalues for each of the $(budget_count) measurement budgets, took $(round(simulate_time, digits = 3)) s.",
            )
        end
        # Estimate the gate eigenvalues and probabilities for each of the shots
        estimate_start = time()
        for s in 1:budget_count
            noise_est = estimate_gate_noise(
                d,
                est_eigenvalues_experiment_ensemble_set[s],
                est_covariance_experiment_ensemble_set[s],
                shots_set[s];
                min_eigenvalue = min_eigenvalue,
                clip_warning = clip_warning,
                split = split,
                split_warning = false,
                precision = precision,
            )
            noise_error = get_noise_error(d, noise_est)
            noise_est_coll[idx, s] = noise_est
            noise_error_coll[idx, s] = noise_error
            if N >= N_warn
                GC.gc()
            end
        end
        estimate_time = time() - estimate_start
        calculation_times[idx, :] = [simulate_time; estimate_time]
        overall_time += simulate_time + estimate_time
        if diagnostics
            if diagnostics
                println(
                    "Estimating the gate eigenvalues and projecting the estimates into the probability simplex for each of the $(budget_count) measurement budgets took $(round(estimate_time, digits = 3)) s.",
                )
            end
        end
        if (idx % save_interval == 0) && idx != repetitions
            if save_data
                aces_data = ACESData(
                    d,
                    noise_est_coll[1:idx, :],
                    noise_error_coll[1:idx, :],
                    budget_set,
                    shots_set,
                    idx,
                    seed,
                    seeds[1:idx],
                    calculation_times[1:idx, :],
                    overall_time,
                )
                save_aces(aces_data)
            end
        end
        if diagnostics
            println(
                "Simulated $(idx) of $(repetitions) repetitions of ACES. The time elapsed since simulation started is $(round(time() - repetitions_start, digits = 3)) s.",
            )
        end
    end
    aces_data = ACESData(
        d,
        noise_est_coll,
        noise_error_coll,
        budget_set,
        shots_set,
        repetitions,
        seed,
        seeds,
        calculation_times,
        overall_time,
    )
    if save_data && repetitions > saved_repetitions
        save_aces(aces_data; clear_design = clear_design)
    end
    if diagnostics
        println(
            "Simulated all $(repetitions) repetitions of ACES. The time elapsed since simulation started is $(round(time() - repetitions_start, digits = 3)) s.",
        )
    end
    return aces_data::ACESData
end

"""
    get_noise_score(aces_data::ACESData, merit::Merit)

Returns the z-scores for the supplied normalised root-mean-square error (NRMSE) data in the noise errors stored in `aces_data`, given the merit `merit`.
"""
function get_noise_score(aces_data::ACESData, merit::Merit)
    noise_score_coll = get_noise_score(aces_data.noise_error_coll, merit)
    return noise_score_coll::Matrix{NoiseScore}
end

"""
    get_model_violation(aces_data::ACESData; projected::Bool = false)

Returns the noise model violation z-score for the generalised residual sum of squares corresponding to the noise estimates stored in `aces_data`, given the design also stored in `aces_data`, calculating with the projected gate eigenvalues if `projected` is `true`.
"""
function get_model_violation(aces_data::ACESData; projected::Bool = false)
    model_violation_coll =
        get_model_violation(aces_data.d, aces_data.noise_est_coll; projected = projected)
    return model_violation_coll::Matrix{Float64}
end

# Import Stim
const stim = PythonCall.pynew()
function __init__()
    return PythonCall.pycopy!(stim, pyimport("stim"))
end

"""
    ACESData

ACES noise characterisation experiment simulation results.

# Fields

  - `d::Design`: Experimental design.
  - `budget_set::Vector{Int}`: Measurement budgets.
  - `shots_set::Vector{Int}`: Measurement shots corresponding to the measurement budgets in `budget_set`.
  - `repetitions::Int`: Number of times to repeat the ACES estimation procedure.
  - `seeds::Vector{UInt64}`: Seeds for each of the repetitions.
  - `eigenvalues::Vector{Float64}`: Circuit eigenvalues.
  - `covariance::SparseMatrixCSC{Float64, Int32}`: Circuit eigenvalue estimator covariance matrix.
  - `est_eigenvalues_coll::Matrix{Vector{Float64}}`: Estimated circuit eigenvalues for each of the measurement budgets and repetitions.
  - `fgls_gate_eigenvalues_coll::Matrix{Vector{Float64}}`: FGLS estimated gate eigenvalues for each of the measurement budgets and repetitions.
  - `gls_gate_eigenvalues_coll::Matrix{Vector{Float64}}`: GLS estimated gate eigenvalues for each of the measurement budgets and repetitions, which uses the true circuit eigenvalue estimator covariance matrix.
  - `wls_gate_eigenvalues_coll::Matrix{Vector{Float64}}`: WLS estimated gate eigenvalues for each of the measurement budgets and repetitions.
  - `ols_gate_eigenvalues_coll::Matrix{Vector{Float64}}`: OLS estimated gate eigenvalues for each of the measurement budgets and repetitions.
  - `fgls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}`: FGLS estimated gate probability distributions for each of the measurement budgets and repetitions.
  - `gls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}`: GLS estimated gate probability distributions for each of the measurement budgets and repetitions.
  - `wls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}`: WLS estimated gate probability distributions for each of the measurement budgets and repetitions.
  - `ols_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}`: OLS estimated gate probability distributions for each of the measurement budgets and repetitions.
  - `fgls_gate_norm_coll::Matrix{Float64}`: Normalised RMS error between the FGLS estimated gate eigenvalues and the true gate eigenvalues for each of the measurement budgets and repetitions.
  - `gls_gate_norm_coll::Matrix{Float64}`: Normalised RMS error between the GLS estimated gate eigenvalues and the true gate eigenvalues for each of the measurement budgets and repetitions.
  - `wls_gate_norm_coll::Matrix{Float64}`: Normalised RMS error between the WLS estimated gate eigenvalues and the true gate eigenvalues for each of the measurement budgets and repetitions.
  - `ols_gate_norm_coll::Matrix{Float64}`: Normalised RMS error between the OLS estimated gate eigenvalues and the true gate eigenvalues for each of the measurement budgets and repetitions.
  - `calculation_times::Matrix{Float64}`: Time taken to simulate sampling and estimate the gate eigenvalues with FGLS, GLS, WLS, and OLS for each repetition.
  - `overall_time::Float64`: Overall time taken to simulate ACES across all repetitions.
"""
struct ACESData
    d::Design
    budget_set::Vector{Int}
    shots_set::Vector{Int}
    repetitions::Int
    seeds::Vector{UInt64}
    eigenvalues::Vector{Float64}
    covariance::SparseMatrixCSC{Float64, Int32}
    est_eigenvalues_coll::Matrix{Vector{Float64}}
    fgls_gate_eigenvalues_coll::Matrix{Vector{Float64}}
    gls_gate_eigenvalues_coll::Matrix{Vector{Float64}}
    wls_gate_eigenvalues_coll::Matrix{Vector{Float64}}
    ols_gate_eigenvalues_coll::Matrix{Vector{Float64}}
    fgls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}
    gls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}
    wls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}
    ols_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}
    fgls_gate_norm_coll::Matrix{Float64}
    gls_gate_norm_coll::Matrix{Float64}
    wls_gate_norm_coll::Matrix{Float64}
    ols_gate_norm_coll::Matrix{Float64}
    calculation_times::Matrix{Float64}
    overall_time::Float64
end

function Base.show(io::IO, a::ACESData)
    return print(
        io,
        "ACES data from a design for a $(a.d.c.circuit_param.circuit_name) circuit with $(length(a.d.tuple_set)) tuples and $(a.d.experiment_number) experiments.",
    )
end

@struct_hash_equal_isequal ACESData

"""
    get_stim_circuit_string(circuit::Vector{Layer}, gate_probabilities::Dict{Gate, Vector{Float64}}, add_prep::Bool, add_meas::Bool)

Returns a Stim string representation of the circuit `circuit` alongside error probabilities specified by `gate_probabilities`, as well as noisy preparations if `add_prep` is `true`, and noisy measurements if `add_meas` is `true`.
"""
function get_stim_circuit_string(
    circuit::Vector{Layer},
    gate_probabilities::Dict{Gate, Vector{Float64}},
    add_prep::Bool,
    add_meas::Bool,
)
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
    # Currently supported gate types
    preparation_gates = ["PZ+", "PZ-", "PX+", "PX-", "PY+", "PY-"]
    measurement_gates = ["MZ", "MX", "MY"]
    one_qubit_gates = ["H", "S", "X", "Y", "Z", "I"]
    two_qubit_gates = ["CX", "CZ"]
    # Generate the gate list for Stim
    string_vector = Vector{String}(undef, length(circuit))
    for (i, l) in enumerate(circuit)
        noise_string_vector = Vector{String}(undef, length(l.layer))
        gate_string_vector = Vector{String}(undef, length(l.layer))
        for (j, gate) in enumerate(l.layer)
            if gate.type in one_qubit_gates
                noise_string = "PAULI_CHANNEL_1($(join(gate_probabilities[gate][order_1], ", "))) $(gate.targets[1])\n"
                gate_string = "$(gate.type) $(gate.targets[1])\n"
                noise_string_vector[j] = noise_string
                gate_string_vector[j] = gate_string
            elseif gate.type in two_qubit_gates
                noise_string = "PAULI_CHANNEL_2($(join(gate_probabilities[gate][order_2], ", "))) $(join(gate.targets, " "))\n"
                gate_string = "$(gate.type) $(join(gate.targets, " "))\n"
                noise_string_vector[j] = noise_string
                gate_string_vector[j] = gate_string
            elseif gate.type in measurement_gates
                if add_meas
                    meas_string = "$(gate.type)($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                else
                    meas_string = "$(gate.type) $(gate.targets[1])\n"
                end
                noise_string_vector[j] = ""
                gate_string_vector[j] = meas_string
            elseif gate.type in preparation_gates
                gate_prep = gate.type[1:(end - 1)]
                gate_sign = gate.type[end]
                if add_prep
                    prep_error_string = "X_ERROR($(gate_probabilities[gate][2])) $(gate.targets[1])\n"
                else
                    prep_error_string = ""
                end
                if gate_sign == '-'
                    sign_string = "X $(gate.targets[1])\n"
                elseif gate_sign == '+'
                    sign_string = ""
                else
                    throw(error("Invalid preparation gate sign $(gate_sign)."))
                end
                if gate_prep == "PZ"
                    prep_string = ""
                elseif gate_prep == "PX"
                    prep_string = "H $(gate.targets[1])\n"
                elseif gate_prep == "PY"
                    prep_string = "H $(gate.targets[1])\nS $(gate.targets[1])\n"
                else
                    throw(error("Invalid preparation gate type $(gate_prep)."))
                end
                noise_string_vector[j] = prep_error_string
                gate_string_vector[j] = sign_string * prep_string
            else
                throw(error("Invalid gate type $(gate.type)."))
            end
        end
        string_vector[i] = join(noise_string_vector) * join(gate_string_vector)
    end
    stim_circuit_string = join(string_vector)
    return stim_circuit_string::String
end

"""
    stim_sample(stim_circuit_string::String, shots::Int; stim_seed::Union{UInt64, Nothing} = nothing, force_gc::Bool = false)

Returns bit-packaged measurement outcomes from simulating the circuit `stim_circuit_string` over `shots` measurement shots using the Python package `Stim`.
While the seed `stim_seed` can be fixed, `Stim` guarantees inconsistency when using the same seed on different versions.
"""
function stim_sample(
    stim_circuit_string::String,
    shots::Int;
    stim_seed::Union{UInt64, Nothing} = nothing,
    force_gc::Bool = false,
)
    # Generate the circuit from the string and sample shots
    # Disable PythonCall garbage collection and ensure Python is run on thread 1 to avoid crashes
    @assert Threads.threadid() == 1
    if force_gc
        GC.gc()
    end
    PythonCall.GC.disable()
    # Generate the Stim circuit and sampler
    if stim_seed !== nothing
        stim_circuit_sampler =
            stim.Circuit(stim_circuit_string).compile_sampler(; seed = stim_seed)
    else
        stim_circuit_sampler = stim.Circuit(stim_circuit_string).compile_sampler()
    end
    # Sample all of the shots at once
    stim_samples = pyconvert(
        Matrix{UInt8},
        stim_circuit_sampler.sample(; shots = shots, bit_packed = true),
    )
    # Re-enable PythonCall garbage collection
    PythonCall.GC.enable()
    if force_gc
        GC.gc()
    end
    return stim_samples::Matrix{UInt8}
end

"""
    batch_shots(shots::Int, measurements::Int, max_samples::Int)

Returns the shots divided into batches for sampling from Stim.
"""
function batch_shots(shots::Int, measurements::Int, max_samples::Int)
    # Divide the shots into batches
    # Stim samples in batches of 256, so we want batch sizes that are multiples of 256
    base = 256
    batch_num = cld(measurements * shots, max_samples)
    base_mult = fld(fld(shots, base) + 1, batch_num)
    remainder = shots - (batch_num * base_mult - 1) * base
    shot_batches = Vector{Int}(undef, batch_num)
    for b in 1:batch_num
        if b < batch_num
            shot_batches[b] = base_mult * base
        else
            shot_batches[b] = (base_mult - 1) * base + remainder
        end
    end
    return shot_batches::Vector{Int}
end

"""
    estimate_eigenvalues(d::Design, shots_set::Vector{Int}; kwargs...)

Returns simulated estimated circuit eigenvalues for the experimental design `d` across each of the measurement shots counts specified in `shots_set`.

# Keyword arguments

  - `seed::Union{UInt64, Nothing} = nothing`: Seed controlling the random seeds for Stim.
  - `epsilon::Float64 = 0.1`: Set estimated circuit eigenvalues below this threshold to this value, which sometimes occurs for experiments with small shot weights when the measurement budget is small.
  - `detailed_diagnostics::Bool = false`: Whether to print diagnostic information about the simulation.
  - `max_samples::Int = 10^10`: Maximum number of samples to take in a single Stim simulation.
  - `force_gc::Bool = false`: Whether to force garbage collection before and after each Stim simulation.
"""
function estimate_eigenvalues(
    d::Design,
    shots_set::Vector{Int};
    seed::Union{UInt64, Nothing} = nothing,
    epsilon::Float64 = 0.1,
    detailed_diagnostics::Bool = false,
    max_samples::Int = 10^10,
    force_gc::Bool = false,
)
    # Set up some variables describing organisation of the circuits
    n = d.c.qubit_num
    gate_probabilities = d.c.gate_probabilities
    add_prep = d.c.add_prep
    add_meas = d.c.add_meas
    budget_count = length(shots_set)
    tuple_number = length(d.tuple_set)
    start_time = time()
    # Divide the shots between the experiments
    @assert sum(d.shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
    @assert all(d.shot_weights .> 0.0) "The shot weights are not all positive."
    shots_divided_float = [
        shots_set[s] * d.shot_weights[t] / d.experiment_numbers[t] for
        t in 1:tuple_number, s in 1:budget_count
    ]
    @assert vec(sum(shots_divided_float .* d.experiment_numbers; dims = 1)) ≈ shots_set "The shots have not been divided correctly."
    shots_divided = ceil.(Int, shots_divided_float)
    shots_maximum = vec(maximum(shots_divided; dims = 2))
    # Generate the requisite random seeds for Stim
    if seed !== nothing
        Random.seed!(seed)
    end
    batched_shots = any(shots_maximum * n .> max_samples)
    if batched_shots
        stim_seeds = Vector{Vector{Vector{Vector{UInt64}}}}(undef, tuple_number)
        for i in 1:tuple_number
            batches = length(batch_shots(shots_maximum[i], n, max_samples))
            tuple_circuit_number = length(d.prep_ensemble[i])
            stim_seeds[i] = Vector{Vector{Vector{UInt64}}}(undef, tuple_circuit_number)
            for j in 1:tuple_circuit_number
                sign_circuit_number = length(d.prep_ensemble[i][j])
                stim_seeds[i][j] = Vector{Vector{UInt64}}(undef, sign_circuit_number)
                for k in 1:sign_circuit_number
                    stim_seeds[i][j][k] = rand(UInt64, batches)
                end
            end
        end
    else
        stim_seeds = Vector{Vector{Vector{UInt64}}}(undef, tuple_number)
        for i in 1:tuple_number
            tuple_circuit_number = length(d.prep_ensemble[i])
            stim_seeds[i] = Vector{Vector{UInt64}}(undef, tuple_circuit_number)
            for j in 1:tuple_circuit_number
                sign_circuit_number = length(d.prep_ensemble[i][j])
                stim_seeds[i][j] = rand(UInt64, sign_circuit_number)
            end
        end
    end
    if seed !== nothing
        Random.seed!()
    end
    # Determine and sample from the circuits and then process their contributions to the eigenvalues
    eigenvalues_coll = Vector{Vector{Vector{Float64}}}(undef, budget_count)
    for s in 1:budget_count
        eigenvalues_coll[s] = Vector{Vector{Float64}}(undef, tuple_number)
    end
    for i in 1:tuple_number
        # Initialise variables
        circuit_tuple = d.tuple_set[i]
        tuple_circuit = d.c.circuit[circuit_tuple]
        tuple_circuit_string =
            get_stim_circuit_string(tuple_circuit, gate_probabilities, add_prep, add_meas)
        # Initialise the Pauli mappings
        mappings = d.mapping_ensemble[i]
        L = length(mappings)
        initial_set = Vector{Pauli}(undef, L)
        initial_support_set = Vector{Vector{Int}}(undef, L)
        final_set = Vector{Pauli}(undef, L)
        final_support_set = Vector{Vector{Int}}(undef, L)
        for l in 1:L
            m = mappings[l]
            initial_set[l] = m.initial
            initial_support_set[l] = get_support(m.initial)
            final_set[l] = m.final
            final_support_set[l] = get_support(m.final)
        end
        # Initialise the eigenvalue estimator containers
        eigenvalues_pos = [zeros(Int, L) for _ in 1:budget_count]
        eigenvalues_neg = [zeros(Int, L) for _ in 1:budget_count]
        shots_pos = [zeros(Int, L) for _ in 1:budget_count]
        shots_neg = [zeros(Int, L) for _ in 1:budget_count]
        eigenvalues = [zeros(Float64, L) for _ in 1:budget_count]
        shots = [zeros(Int, L) for _ in 1:budget_count]
        # Generate the circuits and estimate the eigenvalues
        tuple_circuit_number = length(d.meas_ensemble[i])
        for j in 1:tuple_circuit_number
            # Initialise variables
            prep_layers = d.prep_ensemble[i][j]
            sign_circuit_number = length(prep_layers)
            meas_layer = d.meas_ensemble[i][j]
            meas_layer_string = get_stim_circuit_string(
                [meas_layer],
                gate_probabilities,
                add_prep,
                add_meas,
            )
            # Determine the qubits prepared and measured by the circuit
            prep_qubits = [
                [gate.targets[1] for gate in prep_layers[k].layer] for
                k in 1:sign_circuit_number
            ]
            meas_qubits = [gate.targets[1] for gate in d.meas_ensemble[i][j].layer]
            # Sample the data for all the different sign configurations
            sign_data = Vector{Matrix{UInt8}}(undef, sign_circuit_number)
            for k in 1:sign_circuit_number
                # Initialise variables
                prep_layer_string = get_stim_circuit_string(
                    [prep_layers[k]],
                    gate_probabilities,
                    add_prep,
                    add_meas,
                )
                # Determine the circuit
                stim_circuit_string =
                    prep_layer_string * tuple_circuit_string * meas_layer_string
                # Sample from the circuit
                if batched_shots
                    # If there are too many shots to sample at once, divide them into batches
                    shot_batches = batch_shots(shots_maximum[i], n, max_samples)
                    batches = length(shot_batches)
                    sign_data_coll = Vector{Matrix{UInt8}}(undef, batches)
                    for b in 1:batches
                        sign_data_coll[b] = stim_sample(
                            stim_circuit_string,
                            shot_batches[b];
                            stim_seed = stim_seeds[i][j][k][b],
                            force_gc = force_gc,
                        )
                    end
                    sign_data[k] = vcat(sign_data_coll...)
                else
                    # Sample all of the shots at once
                    sign_data[k] = stim_sample(
                        stim_circuit_string,
                        shots_maximum[i];
                        stim_seed = stim_seeds[i][j][k],
                        force_gc = force_gc,
                    )
                end
            end
            # Determine the circuit's contribution to the circuit eigenvalue estimators
            experiment = d.experiment_ensemble[i][j]
            for pauli_index in experiment
                # Determine the sign of the Pauli preparations
                prep_indices = [
                    [
                        findfirst(qubit .== prep_qubits[k]) for
                        qubit in initial_support_set[pauli_index]
                    ] for k in 1:sign_circuit_number
                ]
                prep_signs = [
                    iseven(
                        sum(
                            gate.type[end] == '-' for
                            gate in prep_layers[k].layer[prep_indices[k]]
                        ),
                    ) for k in 1:sign_circuit_number
                ]
                # Determine the measured qubits over which to marginalise
                meas_indices = [
                    findfirst(qubit .== meas_qubits) for
                    qubit in final_support_set[pauli_index]
                ]
                # Marginalise the sign data to estimate the eigenvalue
                for k in 1:sign_circuit_number
                    for idx in 1:budget_count
                        shots_value = shots_divided[i, idx]
                        pauli_shots = 0
                        for s in 1:shots_value
                            pauli_shot = 0
                            for m in meas_indices
                                # The measurement outcomes are packed into UInt8
                                pauli_shot +=
                                    (sign_data[k][s, 1 + fld(m - 1, 8)] >> ((m - 1) % 8)) &
                                    1
                            end
                            pauli_shots += 1 - 2 * (pauli_shot % 2)
                        end
                        if prep_signs[k]
                            eigenvalues_pos[idx][pauli_index] += pauli_shots
                            shots_pos[idx][pauli_index] += shots_value
                        else
                            eigenvalues_neg[idx][pauli_index] += pauli_shots
                            shots_neg[idx][pauli_index] += shots_value
                        end
                    end
                end
            end
            if detailed_diagnostics
                println(
                    "Simulated sampling experiment $(j) of $(tuple_circuit_number) for tuple $(i) of $(tuple_number). The time elapsed since simulation started is $(round(time() - start_time, digits = 3)) s.",
                )
            end
        end
        # Flip each eigenvalue if the sign of its final Pauli is negative
        for l in 1:L
            if final_set[l].pauli[2n + 1]
                for s in 1:budget_count
                    eigenvalues_pos[s][l] = -eigenvalues_pos[s][l]
                    eigenvalues_neg[s][l] = -eigenvalues_neg[s][l]
                end
            end
        end
        # Estimate and save the eigenvalues
        for s in 1:budget_count
            if hasproperty(d.c, :partition)
                @assert shots_pos[s] == shots_neg[s] "The number of positive and negative shots is not equal."
                shots[s] = shots_pos[s] + shots_neg[s]
                eigenvalues[s] =
                    (
                        (eigenvalues_pos[s] ./ shots_pos[s]) .-
                        (eigenvalues_neg[s] ./ shots_neg[s])
                    ) / 2
            else
                @assert all(shots_neg[s] .== 0) "There are negative shots when negative sign configurations should not have been prepared."
                shots[s] = shots_pos[s]
                eigenvalues[s] = eigenvalues_pos[s] ./ shots_pos[s]
            end
            eigenvalues_coll[s][i] = eigenvalues[s]
        end
    end
    # Force Julia to run garbage collection, which aims to prevent occasional segfaults
    # Garbage collection of PythonCall objects mixes badly with multithreaded Julia
    GC.gc()
    # Concatenate the results
    eigenvalues_set = Vector{Vector{Float64}}(undef, budget_count)
    for s in 1:budget_count
        # Collate the eigenvalues
        eigenvalues_set[s] = vcat(eigenvalues_coll[s]...)
        small_eigenvalues = (eigenvalues_set[s] .< epsilon)
        if sum(small_eigenvalues) > 0
            @debug "When sampling $(round(shots_set[s], sigdigits = 4)) shots, $(sum(small_eigenvalues)) of $(length(eigenvalues_set[s])) estimated eigenvalues are smaller than $(epsilon) and have been set to that value, the smallest being $(minimum(eigenvalues_set[s][small_eigenvalues])), and the largest being $(maximum(eigenvalues_set[s][small_eigenvalues]))."
            eigenvalues_set[s][small_eigenvalues] .= epsilon
        end
        large_eigenvalues = (eigenvalues_set[s] .> 1.0)
        if sum(large_eigenvalues) > 0
            @debug "When sampling $(round(shots_set[s], sigdigits = 4)) shots, $(sum(large_eigenvalues)) of $(length(eigenvalues_set[s])) estimated eigenvalues are larger than 1 and have been set to that value, the smallest being $(minimum(eigenvalues_set[s][large_eigenvalues])), and the largest being $(maximum(eigenvalues_set[s][large_eigenvalues]))."
            eigenvalues_set[s][large_eigenvalues] .= 1.0
        end
    end
    return eigenvalues_set::Vector{Vector{Float64}}
end

"""
    fgls_estimate_gate_eigenvalues(d::Design, est_eigenvalues::Vector{Float64}; kwargs...)

Returns the gate eigenvalues estimated from the estimated circuit eigenvalues `est_eigenvalues` corresponding to the design `d` with feasible generalised least squares.

# Keyword arguments

  - `epsilon::Float64 = 1e-10`: Threshold for convergence of the feasible generalised least squares algorithm.
  - `max_iterations::Int = 10`: Maximum number of iterations for the feasible generalised least squares algorithm.
  - `recalc_eigenvalues::Bool = true`: If `true`, the circuit eigenvalues are recalculated from the estimated gate eigenvalues at each iteration, which ensures that the estimated covariance matrix is positive-definite.
  - `constrain::Bool = true`: If `true`, the gate eigenvalues are constrained to be at most 1.
  - `diagnostics::Bool = false`: Whether to print diagnostic information about the feasible generalised least squares algorithm.
"""
function fgls_estimate_gate_eigenvalues(
    d::Design,
    est_eigenvalues::Vector{Float64};
    epsilon::Float64 = 1e-10,
    max_iterations::Int = 10,
    recalc_eigenvalues::Bool = true,
    constrain::Bool = true,
    diagnostics::Bool = false,
)
    # Initialise variables
    N = size(d.matrix, 2)
    # Initialise the gate eigenvalues using the diagonal covariance matrix estimator
    old_gate_eigenvalues = wls_estimate_gate_eigenvalues(d, est_eigenvalues)
    # Recursively perform generalised least squares until convergence
    iter = 1
    while true
        # Calculate the covariance matrix
        if recalc_eigenvalues
            old_eigenvalues = exp.(-(d.matrix * (-log.(old_gate_eigenvalues))))
            covariance = calc_covariance(d, old_eigenvalues, old_gate_eigenvalues)
        else
            covariance = calc_covariance(d, est_eigenvalues, old_gate_eigenvalues)
        end
        # Perform generalised least squares
        new_gate_eigenvalues = gls_estimate_gate_eigenvalues(
            d,
            est_eigenvalues,
            covariance;
            constrain = constrain,
        )
        # Check for convergence
        norm_difference = norm(new_gate_eigenvalues - old_gate_eigenvalues, 2) / sqrt(N)
        old_gate_eigenvalues = new_gate_eigenvalues
        # If the gate eigenvalues have converged, or the maximum iteration number has been reached, return the results
        if norm_difference < epsilon || iter >= max_iterations
            if diagnostics
                if norm_difference < epsilon
                    println(
                        "Feasible generalised least squares has converged after $(iter) iterations.",
                    )
                else
                    println(
                        "The maximum number of feasible generalised least squares iterations $(max_iterations) has been reached without convergence. The difference is $(round(norm_difference, sigdigits = 4)), whereas the threshold for convergence is $(epsilon).",
                    )
                end
            end
            fgls_gate_eigenvalues = new_gate_eigenvalues
            return fgls_gate_eigenvalues::Vector{Float64}
        else
            iter += 1
        end
    end
end

"""
    gls_estimate_gate_eigenvalues(d::Design, est_eigenvalues::Vector{Float64}, est_covariance::SparseMatrixCSC{Float64, Int32}; constrain::Bool = true)

Returns the gate eigenvalues estimated from the estimated circuit eigenvalues `est_eigenvalues` corresponding to the design `d` with generalised least squares, with estimated circuit eigenvalue covariance matrix `est_covariance`.
If `constrain` is `true`, the gate eigenvalues are constrained to be at most 1.
"""
function gls_estimate_gate_eigenvalues(
    d::Design,
    est_eigenvalues::Vector{Float64},
    est_covariance::SparseMatrixCSC{Float64, Int32};
    constrain::Bool = true,
)
    # Make sure that the supplied covariance matrix is positive-definite
    # Determine the scaling factor for generalised least squares from the covariance matrix
    # Use a first-order Taylor approximation to estimate the covariance matrix of the circuit log-eigenvalues
    est_eigenvalues_inv_diag = Diagonal(est_eigenvalues .^ (-1))
    est_covariance_log = sparse(
        Symmetric(est_eigenvalues_inv_diag * est_covariance * est_eigenvalues_inv_diag),
    )
    if any(isinf, est_eigenvalues_inv_diag)
        println("est_eigenvalues_inv_diag contains infs.")
        display(sort(est_eigenvalues))
    end
    # Use the block diagonal structure of the covariance matrix to speed up computation of the GLS factor
    mapping_lengths = length.(d.mapping_ensemble)
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    mapping_number = sum(mapping_lengths)
    covariance_factor = spzeros(Float64, Int, mapping_number, mapping_number)
    for i in eachindex(mapping_lengths)
        mapping_range = mapping_lower[i]:mapping_upper[i]
        mapping_length = mapping_lengths[i]
        # The sparse Cholesky decomposition of A actually computes a decomposition of a permuted matrix PAP'=LL'
        # Hence A=P'LL'P, and our GLS factor is the inverse of P'L, L^(-1)P
        # Note also that A^(-1) = P'L^(-1)'L^(-1)P
        block_chol = cholesky(est_covariance_log[mapping_range, mapping_range])
        block_L_inv = sparse(inv(LowerTriangular(sparse(block_chol.L))))
        block_perm = sparse(1:mapping_length, block_chol.p, ones(mapping_length))
        covariance_factor[mapping_range, mapping_range] = block_L_inv * block_perm
    end
    # Estimate the gate log-eigenvalues using generalised least squares
    design_matrix = convert(SparseMatrixCSC{Float64, Int}, d.matrix)
    gls_gate_log_eigenvalues =
        (covariance_factor * design_matrix) \ (covariance_factor * (-log.(est_eigenvalues)))
    # Ensure that the gate (negative) log-eigenvalues are not less than 0, which corresponds to the gate eigenvalues being greater than 1
    if constrain
        gls_gate_log_eigenvalues[(gls_gate_log_eigenvalues .< 0.0)] .= 0.0
    end
    gls_gate_eigenvalues = exp.(-gls_gate_log_eigenvalues)
    return gls_gate_eigenvalues::Vector{Float64}
end

"""
    wls_estimate_gate_eigenvalues(d::Design, est_eigenvalues::Vector{Float64}; constrain::Bool = true)

Returns the gate eigenvalues estimated from the estimated circuit eigenvalues `est_eigenvalues` corresponding to the design `d` with weighted least squares.
If `constrain` is `true`, the gate eigenvalues are constrained to be at most 1.
"""
function wls_estimate_gate_eigenvalues(
    d::Design,
    est_eigenvalues::Vector{Float64};
    constrain::Bool = true,
)
    # Generate the diagonal covariance matrix
    M = length(est_eigenvalues)
    tuple_number = length(d.tuple_set)
    mapping_lengths = length.(d.mapping_ensemble)
    index_lower = cumsum([0; mapping_lengths[1:(end - 1)]])
    tuple_times_factor = sum(d.shot_weights .* d.tuple_times)
    covariance_diag = Vector{Float64}(undef, M)
    for i in 1:tuple_number
        tuple_covariance_dict = d.covariance_dict_ensemble[i]
        tuple_experiment_number = d.experiment_numbers[i]
        shot_weight = d.shot_weights[i]
        for j in 1:mapping_lengths[i]
            # Calculate the variance scaling factor
            eigenvalue_experiments = tuple_covariance_dict[CartesianIndex(j, j)][2]
            scaling_var =
                (tuple_times_factor / shot_weight) *
                (tuple_experiment_number / eigenvalue_experiments)
            # Calculate the variance of the eigenvalue
            if est_eigenvalues[index_lower[i] + j] < 1.0
                eigenvalue_variance =
                    scaling_var * (1 - est_eigenvalues[index_lower[i] + j]^2)
            else
                eigenvalue_variance =
                    scaling_var * (1 - maximum(est_eigenvalues[est_eigenvalues .< 1.0])^2)
            end
            # Set the covariance matrix term
            covariance_diag[index_lower[i] + j] = eigenvalue_variance
        end
    end
    # Generate the covariance factor 
    covariance_factor = Diagonal(covariance_diag .^ (-1 / 2) .* est_eigenvalues)
    # Estimate the gate log-eigenvalues using generalised least squares
    wls_gate_log_eigenvalues =
        (covariance_factor * convert(SparseMatrixCSC{Float64, Int}, d.matrix)) \
        (covariance_factor * (-log.(est_eigenvalues)))
    # Ensure that the gate (negative) log-eigenvalues are not less than 0, which corresponds to the gate eigenvalues being greater than 1
    if constrain
        wls_gate_log_eigenvalues[(wls_gate_log_eigenvalues .< 0.0)] .= 0.0
    end
    wls_gate_eigenvalues = exp.(-wls_gate_log_eigenvalues)
    return wls_gate_eigenvalues::Vector{Float64}
end

"""
    ols_estimate_gate_eigenvalues(d::Design, est_eigenvalues::Vector{Float64}; constrain::Bool = true)

Returns the gate eigenvalues estimated from the estimated circuit eigenvalues `est_eigenvalues` corresponding to the design `d` with ordinary least squares.
If `constrain` is `true`, the gate eigenvalues are constrained to be at most 1.
"""
function ols_estimate_gate_eigenvalues(
    d::Design,
    est_eigenvalues::Vector{Float64};
    constrain::Bool = true,
)
    # Estimate the gate log-eigenvalues using ordinary least squares
    design_matrix = convert(SparseMatrixCSC{Float64, Int}, d.matrix)
    ols_gate_log_eigenvalues = design_matrix \ (-log.(est_eigenvalues))
    # Ensure that the gate (negative) log-eigenvalues are not less than 0, which corresponds to the gate eigenvalues being greater than 1
    if constrain
        ols_gate_log_eigenvalues[(ols_gate_log_eigenvalues .< 0.0)] .= 0.0
    end
    ols_gate_eigenvalues = exp.(-ols_gate_log_eigenvalues)
    return ols_gate_eigenvalues::Vector{Float64}
end

"""
    estimate_gate_probabilities(d::Design, est_gate_eigenvalues::Vector{Float64})

Returns the gate Pauli error probabilities estimated from the estimated gate eigenvalues `est_gate_eigenvalues` corresponding to the design `d`.
The estimated probability distributions are projected into the probability simplex.
"""
function estimate_gate_probabilities(d::Design, est_gate_eigenvalues::Vector{Float64})
    # Determine the transform matrices
    H₁ = [1 1; 1 -1]
    # The index orders 1-qubit Paulis as:
    # (I=0),  X=1,    Z=2,    Y=3
    # This is the natural bit string ordering.
    W_1 = wht_matrix(1)
    # The index orders 2-qubit Paulis as:
    # (II=0), XI=1,   IX=2,   XX=3
    # ZI=4,   YI=5,   ZX=6,   YX=7
    # IZ=8,   XZ=9,   IY=10,  XY=11
    # ZZ=12,  YZ=13,  ZY=14,  YY=15
    # This is the natural bit string ordering.
    W_2 = wht_matrix(2)
    # Estimate the gate error probabilities from the eigenvalues, projecting into the probability simplex
    gates = d.c.total_gates
    gate_index = d.c.gate_index
    est_gate_probabilities = Dict{Gate, Vector{Float64}}()
    for gate in gates
        # Determine the transformation matrix to turn eigenvalues into error probabilities
        if gate.type ∈ ["MZ", "MX", "MY"]
            transform_matrix = H₁
        elseif length(gate.targets) == 1
            transform_matrix = W_1
        elseif length(gate.targets) == 2
            transform_matrix = W_2
        else
            throw(error("Invalid gate $(gate)."))
        end
        # The transform matrices are all square
        transform_size = size(transform_matrix)[1]
        # Determine the eigenvalues
        est_gate = [
            1.0
            est_gate_eigenvalues[(gate_index[gate] + 1):(gate_index[gate] + transform_size - 1)]
        ]
        # Generate the probability distribution
        est_gate_probabilities[gate] =
            project_simplex(transform_matrix * est_gate / transform_size)
    end
    return est_gate_probabilities::Dict{Gate, Vector{Float64}}
end

"""
    simulate_aces(d::Design, budget_set::Vector{Int}; kwargs...)

Simulates ACES noise characterisation experiments for the experimental design `d` across each of the measurement budgets in `budget_set`.

WARNING: Seeding has the same features as in Stim.
The behaviour of the same random seed will differ across different versions of Stim.
Also, when measurement shots are sampled in batches, which occurs when `max_samples` is exceeded, the results will differ from when all shots are sampled at once.

# Arguments

  - `d::Design`: Experimental design.
  - `budget_set::Vector{Int}`: Measurement budgets for which to simulate ACES.

# Keyword arguments

  - `repetitions::Int = 1`: Number of simulation repetitions.
  - `seed::Union{UInt64, Nothing} = nothing`: the seed to use for the random number generator.
  - `N_warn::Int = 3 * 10^4`: Number of circuit eigenvalues above which to warn the user about certain keyword argument choices.
  - `max_samples::Int = 10^10`: Maximum number of Stim samples collected in a single simulation.
  - `force_gc::Bool = false`: Whether to force garbage collection before and after each Stim simulation; this was added to prevent occasional segfaults but massively slows down the simulation, and currently does not appear to be necessary.
  - `diagnostics::Bool = true`: Whether to print diagnostics.
  - `detailed_diagnostics::Bool = false`: Whether to print detailed diagnostics.
  - `save_data::Bool = false`: Whether to save the data.
  - `save_interval::Int = 50`: Repetition interval at which to save the data.
  - `clear_design::Bool = false`: Whether to clear the saved design data after saving the full simulation data.
"""
function simulate_aces(
    d::Design,
    budget_set::Vector{Int};
    repetitions::Int = 1,
    seed::Union{UInt64, Nothing} = nothing,
    N_warn::Int = 3 * 10^4,
    max_samples::Int = 10^10,
    force_gc::Bool = false,
    diagnostics::Bool = true,
    detailed_diagnostics::Bool = false,
    save_data::Bool = false,
    save_interval::Int = 50,
    clear_design::Bool = false,
)
    # Warn the user if they have unadvisable settings for a large circuit
    if d.c.N >= N_warn
        if ~diagnostics
            @warn "This ACES simulation is for a very large circuit: turning on diagnostics is advised."
        end
        if ~detailed_diagnostics
            @warn "This ACES simulation is for a very large circuit: turning on detailed diagnostics is advised."
        end
        if ~save_data
            @warn "This ACES simulation is for a very large circuit: saving the data is advised."
        end
    end
    # Generate synthetic ACES data
    N = size(d.matrix, 2)
    gate_eigenvalues = d.c.gate_eigenvalues
    (eigenvalues, covariance) = calc_eigenvalues_covariance(d)
    # Normalise the sampled shot count by the amount of time taken to perform the circuits
    budget_count = length(budget_set)
    tuple_times_factor = sum(d.shot_weights .* d.tuple_times)
    shots_set = round.(Int, budget_set / tuple_times_factor)
    @assert all(shots_set .<= budget_set) "The normalised shots are larger than the original shots."
    # Generate the requisite random seeds using the fixed seed
    if seed !== nothing
        Random.seed!(seed)
    end
    seeds = rand(UInt64, repetitions)
    if seed !== nothing
        Random.seed!()
    end
    # Initialise all data storage
    est_eigenvalues_coll = Matrix{Vector{Float64}}(undef, repetitions, budget_count)
    # Only estimate the covariance matrix and perform FGLS and GLS if the full covariance matrix is included in the design
    if d.full_covariance
        # FGLS estimates
        fgls_gate_eigenvalues_coll =
            Matrix{Vector{Float64}}(undef, repetitions, budget_count)
        fgls_gate_probabilities_coll =
            Matrix{Dict{Gate, Vector{Float64}}}(undef, repetitions, budget_count)
        fgls_gate_norm_coll = Matrix{Float64}(undef, repetitions, budget_count)
        # GLS estimates
        gls_gate_eigenvalues_coll =
            Matrix{Vector{Float64}}(undef, repetitions, budget_count)
        gls_gate_probabilities_coll =
            Matrix{Dict{Gate, Vector{Float64}}}(undef, repetitions, budget_count)
        gls_gate_norm_coll = Matrix{Float64}(undef, repetitions, budget_count)
    else
        # FGLS estimates
        fgls_gate_eigenvalues_coll = Matrix{Vector{Float64}}(undef, 0, 0)
        fgls_gate_probabilities_coll = Matrix{Dict{Gate, Vector{Float64}}}(undef, 0, 0)
        fgls_gate_norm_coll = Matrix{Float64}(undef, 0, 0)
        # GLS estimates
        gls_gate_eigenvalues_coll = Matrix{Vector{Float64}}(undef, 0, 0)
        gls_gate_probabilities_coll = Matrix{Dict{Gate, Vector{Float64}}}(undef, 0, 0)
        gls_gate_norm_coll = Matrix{Float64}(undef, 0, 0)
    end
    # WLS estimates
    wls_gate_eigenvalues_coll = Matrix{Vector{Float64}}(undef, repetitions, budget_count)
    wls_gate_probabilities_coll =
        Matrix{Dict{Gate, Vector{Float64}}}(undef, repetitions, budget_count)
    wls_gate_norm_coll = Matrix{Float64}(undef, repetitions, budget_count)
    # OLS estimates
    ols_gate_eigenvalues_coll = Matrix{Vector{Float64}}(undef, repetitions, budget_count)
    ols_gate_probabilities_coll =
        Matrix{Dict{Gate, Vector{Float64}}}(undef, repetitions, budget_count)
    ols_gate_norm_coll = Matrix{Float64}(undef, repetitions, budget_count)
    # Time data
    calculation_times = Matrix{Float64}(undef, repetitions, 5)
    overall_time = 0.0
    # Load saved data
    saved_repetitions = 0
    filename = aces_data_filename(d, budget_set)
    if isfile(pwd() * "/data/" * filename)
        saved_aces_data = load_aces(d, budget_set)
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
        same_data =
            (saved_aces_data.d.c.circuit_param == d.c.circuit_param) &&
            (saved_aces_data.d.c.noise_param == d.c.noise_param) &&
            (saved_aces_data.d.matrix == d.matrix) &&
            (saved_aces_data.d.tuple_set == d.tuple_set) &&
            (saved_aces_data.d.full_covariance == d.full_covariance) &&
            (saved_aces_data.d.shot_weights ≈ d.shot_weights) &&
            (saved_aces_data.budget_set ≈ budget_set) &&
            (saved_aces_data.shots_set ≈ shots_set)
        # Only used the saved data if it has the same parameters
        if same_data
            # Set parameters
            saved_repetitions = saved_aces_data.repetitions
            seeds[1:saved_repetitions] = saved_aces_data.seeds
            # Circuit eigenvalue estimates
            est_eigenvalues_coll[1:saved_repetitions, :] =
                saved_aces_data.est_eigenvalues_coll
            if d.full_covariance
                # FGLS estimates
                fgls_gate_eigenvalues_coll[1:saved_repetitions, :] =
                    saved_aces_data.fgls_gate_eigenvalues_coll
                fgls_gate_probabilities_coll[1:saved_repetitions, :] =
                    saved_aces_data.fgls_gate_probabilities_coll
                fgls_gate_norm_coll[1:saved_repetitions, :] =
                    saved_aces_data.fgls_gate_norm_coll
                # GLS estimates
                gls_gate_eigenvalues_coll[1:saved_repetitions, :] =
                    saved_aces_data.gls_gate_eigenvalues_coll
                gls_gate_probabilities_coll[1:saved_repetitions, :] =
                    saved_aces_data.gls_gate_probabilities_coll
                gls_gate_norm_coll[1:saved_repetitions, :] =
                    saved_aces_data.gls_gate_norm_coll
            end
            # WLS estimates
            wls_gate_eigenvalues_coll[1:saved_repetitions, :] =
                saved_aces_data.wls_gate_eigenvalues_coll
            wls_gate_probabilities_coll[1:saved_repetitions, :] =
                saved_aces_data.wls_gate_probabilities_coll
            wls_gate_norm_coll[1:saved_repetitions, :] = saved_aces_data.wls_gate_norm_coll
            # OLS estimates
            ols_gate_eigenvalues_coll[1:saved_repetitions, :] =
                saved_aces_data.ols_gate_eigenvalues_coll
            ols_gate_probabilities_coll[1:saved_repetitions, :] =
                saved_aces_data.ols_gate_probabilities_coll
            ols_gate_norm_coll[1:saved_repetitions, :] = saved_aces_data.ols_gate_norm_coll
            # Time data
            calculation_times[1:saved_repetitions, :] = saved_aces_data.calculation_times
            overall_time = saved_aces_data.overall_time
            if diagnostics
                println(
                    "Loading $(saved_repetitions) of $(repetitions) repetitions of simulated ACES data.",
                )
            end
        else
            @warn "Attempted and failed to load and reuse the saved data, as its parameters differ from the supplied data."
        end
    end
    repetitions_start = time()
    for idx in (saved_repetitions + 1):repetitions
        # Simulate ACES circuits and process the data to estimate the circuit eigenvalues
        simulate_start = time()
        est_eigenvalues_set = estimate_eigenvalues(
            d,
            shots_set;
            seed = seeds[idx],
            detailed_diagnostics = detailed_diagnostics,
            max_samples = max_samples,
            force_gc = force_gc,
        )
        simulate_time = time() - simulate_start
        if diagnostics
            println(
                "Sampling $(round(maximum(shots_set), sigdigits=4)) shots, normalised from $(round(maximum(budget_set), sigdigits=4)) shots, divided between $(d.experiment_number) experiments for the $(length(d.tuple_set)) tuples in the set, and then estimating the circuit eigenvalues for each of the $(budget_count) different shot counts, took $(round(simulate_time, digits = 3)) s.",
            )
        end
        # Estimate the gate eigenvalues and probabilities for each of the shots
        ls_times = zeros(Float64, 4, budget_count)
        for s in 1:budget_count
            # Estimate the gate eigenvalues and probabilities
            est_eigenvalues = est_eigenvalues_set[s]
            est_eigenvalues_coll[idx, s] = est_eigenvalues
            if d.full_covariance
                # FGLS estimates
                time_1 = time()
                fgls_gate_eigenvalues = fgls_estimate_gate_eigenvalues(
                    d,
                    est_eigenvalues;
                    diagnostics = diagnostics,
                )
                fgls_gate_probabilities =
                    estimate_gate_probabilities(d, fgls_gate_eigenvalues)
                fgls_gate_norm =
                    sqrt(budget_set[s] / N) *
                    norm(fgls_gate_eigenvalues - gate_eigenvalues, 2)
                fgls_gate_eigenvalues_coll[idx, s] = fgls_gate_eigenvalues
                fgls_gate_probabilities_coll[idx, s] = fgls_gate_probabilities
                fgls_gate_norm_coll[idx, s] = fgls_gate_norm
                ls_times[1, s] = time() - time_1
                # GLS estimates
                time_2 = time()
                gls_gate_eigenvalues =
                    gls_estimate_gate_eigenvalues(d, est_eigenvalues, covariance)
                gls_gate_probabilities =
                    estimate_gate_probabilities(d, gls_gate_eigenvalues)
                gls_gate_norm =
                    sqrt(budget_set[s] / N) *
                    norm(gls_gate_eigenvalues - gate_eigenvalues, 2)
                gls_gate_eigenvalues_coll[idx, s] = gls_gate_eigenvalues
                gls_gate_probabilities_coll[idx, s] = gls_gate_probabilities
                gls_gate_norm_coll[idx, s] = gls_gate_norm
                ls_times[2, s] = time() - time_2
            end
            # WLS estimates
            time_3 = time()
            wls_gate_eigenvalues = wls_estimate_gate_eigenvalues(d, est_eigenvalues)
            wls_gate_probabilities = estimate_gate_probabilities(d, wls_gate_eigenvalues)
            wls_gate_norm =
                sqrt(budget_set[s] / N) * norm(wls_gate_eigenvalues - gate_eigenvalues, 2)
            wls_gate_eigenvalues_coll[idx, s] = wls_gate_eigenvalues
            wls_gate_probabilities_coll[idx, s] = wls_gate_probabilities
            wls_gate_norm_coll[idx, s] = wls_gate_norm
            ls_times[3, s] = time() - time_3
            # OLS estimates
            time_4 = time()
            ols_gate_eigenvalues = ols_estimate_gate_eigenvalues(d, est_eigenvalues)
            ols_gate_probabilities = estimate_gate_probabilities(d, ols_gate_eigenvalues)
            ols_gate_norm =
                sqrt(budget_set[s] / N) * norm(ols_gate_eigenvalues - gate_eigenvalues, 2)
            ols_gate_eigenvalues_coll[idx, s] = ols_gate_eigenvalues
            ols_gate_probabilities_coll[idx, s] = ols_gate_probabilities
            ols_gate_norm_coll[idx, s] = ols_gate_norm
            ls_times[4, s] = time() - time_4
        end
        ls_times = vec(sum(ls_times; dims = 2))
        # Time data
        calculation_times[idx, :] = [simulate_time; ls_times]
        overall_time += simulate_time + sum(ls_times)
        if diagnostics
            println(
                "Estimating the gate eigenvalues for each of the $(length(budget_set)) different shot counts with $(d.full_covariance ? "FGLS, GLS, WLS, OLS" : "WLS, OLS") took $(join(round.((d.full_covariance ? ls_times : ls_times[3:4]), digits = 3), ", ")) s.",
            )
        end
        if (idx % save_interval == 0) && idx != repetitions
            if save_data
                if d.full_covariance
                    aces_data = ACESData(
                        d,
                        budget_set,
                        shots_set,
                        idx,
                        seeds[1:idx],
                        eigenvalues,
                        covariance,
                        est_eigenvalues_coll[1:idx, :],
                        fgls_gate_eigenvalues_coll[1:idx, :],
                        gls_gate_eigenvalues_coll[1:idx, :],
                        wls_gate_eigenvalues_coll[1:idx, :],
                        ols_gate_eigenvalues_coll[1:idx, :],
                        fgls_gate_probabilities_coll[1:idx, :],
                        gls_gate_probabilities_coll[1:idx, :],
                        wls_gate_probabilities_coll[1:idx, :],
                        ols_gate_probabilities_coll[1:idx, :],
                        fgls_gate_norm_coll[1:idx, :],
                        gls_gate_norm_coll[1:idx, :],
                        wls_gate_norm_coll[1:idx, :],
                        ols_gate_norm_coll[1:idx, :],
                        calculation_times[1:idx, :],
                        overall_time,
                    )
                else
                    aces_data = ACESData(
                        d,
                        budget_set,
                        shots_set,
                        idx,
                        seeds[1:idx],
                        eigenvalues,
                        covariance,
                        est_eigenvalues_coll[1:idx, :],
                        fgls_gate_eigenvalues_coll,
                        gls_gate_eigenvalues_coll,
                        wls_gate_eigenvalues_coll[1:idx, :],
                        ols_gate_eigenvalues_coll[1:idx, :],
                        fgls_gate_probabilities_coll,
                        gls_gate_probabilities_coll,
                        wls_gate_probabilities_coll[1:idx, :],
                        ols_gate_probabilities_coll[1:idx, :],
                        fgls_gate_norm_coll,
                        gls_gate_norm_coll,
                        wls_gate_norm_coll[1:idx, :],
                        ols_gate_norm_coll[1:idx, :],
                        calculation_times[1:idx, :],
                        overall_time,
                    )
                end
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
        budget_set,
        shots_set,
        repetitions,
        seeds,
        eigenvalues,
        covariance,
        est_eigenvalues_coll,
        fgls_gate_eigenvalues_coll,
        gls_gate_eigenvalues_coll,
        wls_gate_eigenvalues_coll,
        ols_gate_eigenvalues_coll,
        fgls_gate_probabilities_coll,
        gls_gate_probabilities_coll,
        wls_gate_probabilities_coll,
        ols_gate_probabilities_coll,
        fgls_gate_norm_coll,
        gls_gate_norm_coll,
        wls_gate_norm_coll,
        ols_gate_norm_coll,
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

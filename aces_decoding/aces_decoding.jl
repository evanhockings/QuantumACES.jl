using QuantumACES, JLD2
enter_folder("aces_decoding")
# Set up parameters
dist = 3
r_1 = 0.05 / 100
r_2 = 0.4 / 100
r_m = 0.8 / 100
r_im = r_m
r_r = 0.2 / 100
r_1_std_log = 0.5
r_2_std_log = r_1_std_log
r_m_std_log = 0.25
r_im_std_log = r_m_std_log
r_r_std_log = r_m_std_log
seed = UInt(0)
combined = true
# Generate the circuit
dep_param_load = get_dep_param(r_1, r_2, r_m)
dep_param = get_dep_param(;
    r_1 = r_1,
    r_2 = r_2,
    r_m = r_m,
    r_im = r_im,
    r_r = r_r,
    combined = combined,
)
log_param = get_log_param(;
    r_1 = r_1,
    r_2 = r_2,
    r_m = r_m,
    r_im = r_im,
    r_r = r_r,
    r_1_std_log = r_1_std_log,
    r_2_std_log = r_2_std_log,
    r_m_std_log = r_m_std_log,
    r_im_std_log = r_im_std_log,
    r_r_std_log = r_r_std_log,
    seed = seed,
    combined = combined,
)
rotated_param = get_rotated_param(dist)
# Load the design
tuple_number = 30
repeat_numbers = [49; 99; 25; 49; 49; 99]
d = load_design(rotated_param, dep_param_load, tuple_number, repeat_numbers, true, :none)
# Simulate decoding with ACES noise estimates across distances and noise model instances
dist_max = 13
dist_set = collect(dist:2:dist_max)
dist_num = length(dist_set)
rounds_set = [3; 5; 9; 17; 33]
rounds_num = length(rounds_set)
shots_set = [10^6; 10^7]
shots_num = length(shots_set)
times_factor = sum(d.shot_weights .* d.tuple_times)
budget_set = round.(Int, times_factor * shots_set)
shots = 10^5
seed_num_set = [1500; 300; 100; 80; 60; 50]
@assert length(seed_num_set) == dist_num
seed_coll = [collect(UInt64, 1:seed_num) for seed_num in seed_num_set]
print_repetitions = 10
save_between = true
function aces_decoding()
    start_time = time()
    memory_summary_coll = Vector{Vector{Vector{MemorySummary}}}(undef, dist_num)
    for (idx, dist_val) in pairs(dist_set)
        # Construct the design at the distance
        dist_start = time()
        rotated_param_dist = deepcopy(rotated_param)
        rotated_param_dist.params[:vertical_dist] = dist_val
        rotated_param_dist.params[:horizontal_dist] = dist_val
        circuit_dist = get_circuit(rotated_param_dist, dep_param)
        d_dist = generate_design(circuit_dist, d)
        # Simulate ACES and decoding across noise model instances
        memory_summary_coll[idx] = Vector{Vector{MemorySummary}}()
        for seed_val in seed_coll[idx]
            log_param_seed = deepcopy(log_param)
            log_param_seed.params[:seed] = seed_val
            d_dist_log = update_noise(d_dist, log_param_seed)
            circuit_dist_log = d_dist_log.c
            aces_data_dist_log = simulate_aces(
                d_dist_log,
                budget_set;
                split = true,
                seed = seed_val,
                diagnostics = false,
            )
            noise_est_set = aces_data_dist_log.noise_est_coll[1, :]
            decoder_gate_probabilities = [
                circuit_dist_log.gate_probabilities
                circuit_dist.gate_probabilities
                [noise_est.gls_gate_probabilities for noise_est in noise_est_set]
            ]
            decoder_labels = [
                "True"
                "Depolarising"
                ["ACES S'=$(shots_set[idx])" for idx in 1:shots_num]
            ]
            memory_summary_rounds = Vector{MemorySummary}()
            for rounds in rounds_set
                memory_data = simulate_memory(
                    circuit_dist_log,
                    rounds,
                    shots;
                    seed = seed_val,
                    decoder_gate_probabilities = decoder_gate_probabilities,
                    decoder_labels = decoder_labels,
                )
                memory_summary = get_memory_summary(memory_data)
                push!(memory_summary_rounds, memory_summary)
            end
            push!(memory_summary_coll[idx], memory_summary_rounds)
            if seed_val % print_repetitions == 0 || seed_val == seed_coll[idx][end]
                println(
                    "Decoding at distance $(dist_val) for seeds up to $(seed_val) has taken $(round(time() - dist_start, digits = 2)) s.",
                )
                if save_between
                    save(
                        pwd() * "/data/aces_decoding_data_dist_$(dist_val).jld2",
                        Dict("memory_summary_set" => memory_summary_coll[idx]),
                    )
                end
            end
        end
        println(
            "Decoding at distance $(dist_val) took $(round(time() - dist_start, digits = 2)) s. The total time elapsed is $(round(time() - start_time, digits = 2)) s.",
        )
    end
    save(
        pwd() * "/data/aces_decoding_data.jld2",
        Dict("memory_summary_coll" => memory_summary_coll),
    )
    return nothing
end
aces_decoding()
println("Complete: aces_decoding.jl")

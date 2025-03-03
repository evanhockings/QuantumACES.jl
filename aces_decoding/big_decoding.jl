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
merit = calc_merit(d)
# Generate the design at distance 25 and simulate ACES
dist_big = 25
rotated_param_big = get_rotated_param(dist_big)
rotated_planar_big = get_circuit(rotated_param_big, dep_param)
rotated_planar_big_log = get_circuit(rotated_param_big, log_param)
shots_set = [10^6; 10^7]
times_factor = sum(d.shot_weights .* d.tuple_times)
budget_set = round.(Int, times_factor * shots_set)
aces_data_big = load_aces(
    rotated_param_big,
    log_param,
    tuple_number,
    repeat_numbers,
    false,
    :none,
    budget_set,
    seed,
)
noise_est_set = aces_data_big.noise_est_coll[1, :]
shots_num = length(shots_set)
@assert aces_data_big.d.c == rotated_planar_big_log
# Simulate a memory experiment
big_rounds = dist_big
big_shots = 5 * 10^6
decoder_gate_probabilities = [
    rotated_planar_big_log.gate_probabilities
    rotated_planar_big.gate_probabilities
    [noise_est.wls_gate_probabilities for noise_est in noise_est_set]
]
decoder_labels = [
    "True"
    "Depolarising"
    ["ACES S'=$(shots_set[idx])" for idx in 1:shots_num]
]
big_memory_data = simulate_memory(
    rotated_planar_big_log,
    big_rounds,
    big_shots;
    seed = seed,
    decoder_gate_probabilities = decoder_gate_probabilities,
    decoder_labels = decoder_labels,
    diagnostics = true,
)
big_memory_summary = get_memory_summary(big_memory_data)
# Simulate memory experiments across rounds
rounds_set = [3; 5; 9; 17; 33]
shots = 10^6
memory_summary_rounds = Vector{MemorySummary}()
for rounds in rounds_set
    memory_data = simulate_memory(
        rotated_planar_big_log,
        rounds,
        shots;
        seed = seed,
        decoder_gate_probabilities = decoder_gate_probabilities,
        decoder_labels = decoder_labels,
        diagnostics = true,
    )
    memory_summary = get_memory_summary(memory_data)
    push!(memory_summary_rounds, memory_summary)
end
# Save the data
save(
    "data/big_decoding_data.jld2",
    Dict(
        "big_memory_data" => big_memory_data,
        "big_memory_summary" => big_memory_summary,
        "memory_summary_rounds" => memory_summary_rounds,
    ),
)
println("Complete: big_decoding.jl")

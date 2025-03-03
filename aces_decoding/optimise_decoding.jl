using QuantumACES
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
dep_param = get_dep_param(r_1, r_2, r_m)
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
rotated_planar = get_circuit(rotated_param, dep_param)
# Optimise the design
# The rotated planar syndrome extraction circuit has 7 unique orindary layers and 1 unique reset layer
# The basic tuple set consists of those 8 layers and the empty tuple
# The repeat tuple set consists of the 7 ordinary layers
# Target 30 tuples: 9 (basic tuple set), 2 x 7 (repeat tuple set), 7 extra
add_circuit = false
max_depth = 100
repeat_points = 2
extra_tuple_number = 7
excursion_length = 10
excursions_unchanged = 5
d_opt = optimise_design(
    rotated_planar;
    options = OptimOptions(;
        add_circuit = add_circuit,
        max_depth = max_depth,
        repeat_points = repeat_points,
        extra_tuple_number = extra_tuple_number,
        excursion_length = excursion_length,
        excursions_unchanged = excursions_unchanged,
        seed = seed,
        save_data = true,
    ),
)
# Set the shot weights as if it was a randomised design
min_randomisations = 512
target_shot_budget = 10^7
experiment_shots = 64
(randomisations, shot_budget) =
    get_randomisations(d_opt, min_randomisations, target_shot_budget, experiment_shots)
shot_weights =
    (randomisations .* d_opt.experiment_numbers) /
    sum((randomisations .* d_opt.experiment_numbers))
d = generate_design(
    rotated_planar,
    d_opt.tuple_set_data;
    shot_weights = shot_weights,
    save_data = true,
)
merit = calc_merit(d)
# Generate the design at distance 25 and simulate ACES
dist_big = 25
rotated_param_big = get_rotated_param(dist_big)
rotated_planar_big_log = get_circuit(rotated_param_big, log_param)
shots_set = [10^6; 10^7]
times_factor = sum(d.shot_weights .* d.tuple_times)
budget_set = round.(Int, times_factor * shots_set)
d_big = generate_design(
    rotated_planar_big_log,
    d;
    full_covariance = false,
    diagnostics = true,
    save_data = true,
)
aces_data_big = simulate_aces(
    d_big,
    budget_set;
    split = true,
    seed = seed,
    save_data = true,
    clear_design = true,
)
println("Complete: optimise_decoding.jl")

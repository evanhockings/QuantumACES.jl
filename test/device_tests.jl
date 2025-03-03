using QuantumACES, StatsBase, Test
# Make sure to delete the data folder even if previous tests failed
enter_folder("test")
rm("data"; force = true, recursive = true)
# Set up codes
dist = 3
r_1 = 0.1 / 100
r_2 = 0.5 / 100
r_m = 3.0 / 100
std_log = sqrt(log(10 / 9))
seed = UInt(0)
combined = true
hex_param = get_hex_param(dist)
dep_param = get_dep_param(r_1, r_2, r_m; combined = combined)
log_param = get_log_param(r_1, r_2, r_m, std_log; seed = seed, combined = combined)
heavy_hex = get_circuit(hex_param, dep_param)
# Set up optimisation parameters
ls_type = :wls
est_type = :sum
max_steps = 1
max_cycles = 0
repeat_points = 2
min_depth = 50
max_depth = 50
max_excursions = 0
# Optimise the design
d_hex_dep = optimise_design(
    heavy_hex;
    options = OptimOptions(;
        ls_type = ls_type,
        est_type = est_type,
        max_steps = max_steps,
        max_cycles = max_cycles,
        repeat_points = repeat_points,
        min_depth = min_depth,
        max_depth = max_depth,
        max_excursions = max_excursions,
    ),
)
d_hex = update_noise(d_hex_dep, log_param)
hex_tuple_number = length(d_hex.tuple_set)
gate_eigenvalues = d_hex.c.gate_eigenvalues
N = d_hex.c.gate_data.N
# Test that bit string conversion works properly
@testset "Conversion between UInt8 vectors and bitstrings" begin
    uint8_length = 13
    bit_num = 101
    @assert bit_num % 8 != 0
    rand_uint8_vector = rand(UInt8, uint8_length)
    rand_bitstring = join(rand(["0", "1"], bit_num))
    @test rand_uint8_vector ==
          parse_uint8_vector(parse_bitstring(rand_uint8_vector, 8 * uint8_length))
    @test rand_bitstring == parse_bitstring(parse_uint8_vector(rand_bitstring), bit_num)
end
# Generate a randomised design
target_shot_budget = 10^7
experiment_shots = 10^4
min_randomisations = 2
d_hex_rand = generate_rand_design(
    d_hex,
    min_randomisations,
    target_shot_budget,
    experiment_shots;
    seed = seed,
)
display(d_hex_rand)
total_randomisations = sum(d_hex_rand.randomisations)
meas_budget = round(Int, get_meas_budget(d_hex_rand))
# Test saving and loading a randomised design
@testset "Saving and loading randomised designs" begin
    save_rand_design(d_hex_rand)
    @test load_rand_design(d_hex, total_randomisations, seed) == d_hex_rand
    delete_rand_design(d_hex_rand)
end
# Generate an ensemble of Stim and Qiskit circuits
# Place the circuit on the layout of IBM Torino
# Circuit layout based on IBM Torino
torino_qubit_num = 49
torino_data_qubits = [5; 24; 43; 7; 26; 45; 9; 28; 47]
torino_ancilla_X_qubits = [4; 23; 6; 25; 44; 8; 27; 46; 29; 48]
torino_ancilla_Z_qubits = [16; 35; 17; 36]
torino_qubit_map = [torino_data_qubits; torino_ancilla_X_qubits; torino_ancilla_Z_qubits]
(stim_hex_ensemble, qiskit_hex_ensemble) =
    get_stim_qiskit_ensemble(d_hex_rand, torino_qubit_num, torino_qubit_map)
# Test saving and loading a randomised ensemble
# This test is disabled by default because it is very slow
long_test = false
if long_test
    using PythonCall
    @testset "Saving and loading randomised ensembles" begin
        # Saving and loading the Qiskit ensemble
        save_qiskit_ensemble(d_hex_rand, qiskit_hex_ensemble)
        @test pyconvert(Bool, load_qiskit_ensemble(d_hex_rand) == qiskit_hex_ensemble)
        delete_qiskit_ensemble(d_hex_rand)
    end
end
# Get the design corresponding to the randomised design, with the appropriate shot weights
d_hex_shot = get_design(d_hex_rand)
merit_hex_shot = calc_merit(d_hex_shot)
# Simulate the designs using the original, Stim randomised ensemble, and Qiskit randomised ensemble functions
z_score_cutoff_upper = 3.5
z_score_cutoff_lower = -3.5
@testset "Noise estimation" begin
    # Simulate the design using the original function
    aces_data = simulate_aces(d_hex_shot, [meas_budget]; seed = seed)
    noise_est = aces_data.noise_est_coll[1, 1]
    noise_error = aces_data.noise_error_coll[1, 1]
    noise_score = get_noise_score(noise_error, merit_hex_shot)
    model_violation = get_model_violation(d_hex_shot, noise_est)
    rss = get_rss(d_hex_shot, noise_est)
    aic = get_aic(d_hex_shot, noise_est)
    bic = get_bic(d_hex_shot, noise_est)
    @test is_score_expected(noise_score, z_score_cutoff_lower, z_score_cutoff_upper)
    @test z_score_cutoff_lower <= model_violation && model_violation <= z_score_cutoff_upper
end
@testset "Stim randomised ensemble" begin
    # Simulate and process the data in Stim
    simulate_stim_ensemble(d_hex_rand, stim_hex_ensemble, experiment_shots; seed = seed)
    stim_noise_est =
        estimate_stim_ensemble(d_hex_rand, experiment_shots; simulation_seed = seed)
    stim_eigenvalues = get_eigenvalues(stim_noise_est)
    stim_noise_error = get_noise_error(d_hex_rand, stim_noise_est)
    # Check agreement with the merits
    stim_noise_score = get_noise_score(stim_noise_error, merit_hex_shot)
    stim_model_violation = get_model_violation(d_hex_rand, stim_noise_est)
    stim_rss = get_rss(d_hex_rand, stim_noise_est)
    stim_aic = get_aic(d_hex_rand, stim_noise_est)
    stim_bic = get_bic(d_hex_rand, stim_noise_est)
    @test is_score_expected(stim_noise_score, z_score_cutoff_lower, z_score_cutoff_upper)
    @test z_score_cutoff_lower <= stim_model_violation &&
          stim_model_violation <= z_score_cutoff_upper
end
# Use a reduced number of shots because Qiskit simulations are extremely slow
reduced_shots = 64
@testset "Qiskit randomised ensemble" begin
    # Simulate and process the data in Qiskit
    simulate_qiskit_ensemble(d_hex_rand, qiskit_hex_ensemble, reduced_shots)
    process_qiskit_ensemble(d_hex_rand, torino_qubit_num, torino_qubit_map, reduced_shots)
    qiskit_noise_est = estimate_qiskit_ensemble(d_hex_rand, torino_qubit_map, reduced_shots)
    qiskit_eigenvalues = get_eigenvalues(qiskit_noise_est)
    qiskit_gate_eigenvalues = qiskit_noise_est.gls_unproj_gate_eigenvalues
    qiskit_gate_probabilities = qiskit_noise_est.gls_gate_probabilities
    # Test that the eigenvalues are all 1, as the Qiskit circuit simulations are noiseless
    @test qiskit_eigenvalues == ones(Float64, length(qiskit_eigenvalues))
    @test qiskit_gate_eigenvalues == ones(Float64, length(qiskit_gate_eigenvalues))
    no_gate_errors = true
    for value in values(qiskit_gate_probabilities)
        no_error_value = zeros(Float64, length(value))
        no_error_value[1] = 1.0
        if value != no_error_value
            no_gate_errors = false
            break
        end
    end
    @test no_gate_errors
end
# Remove the data folder generated for these tests
rm("data"; force = true, recursive = true)
exit_folder("test")

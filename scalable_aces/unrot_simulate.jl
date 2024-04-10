using AveragedCircuitEigenvalueSampling, JLD2
enter_folder("scalable_aces")
# Set up the parameters
dist = 3
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
shots_set = [10^6; 10^7; 10^8]
repetitions = 1000
unrotated_param = UnrotatedPlanarParameters(dist)
dep_param = DepolarisingParameters(r_1, r_2, r_m)
log_param = LogNormalParameters(r_1, r_2, r_m, total_std_log; seed = seed)
# Load the design
metadata_dict = load("data/design_metadata_$(code_filename(unrotated_param)).jld2")
@assert unrotated_param == metadata_dict["unrotated_param"]
@assert dep_param == metadata_dict["dep_param"]
dep_param_set = metadata_dict["dep_param_set"]
tuple_number_set = metadata_dict["tuple_number_set"]
repeat_numbers_set = metadata_dict["repeat_numbers_set"]
dep_idx = 14
@assert dep_param == dep_param_set[dep_idx]
d = load_design(
    unrotated_param,
    dep_param_set[dep_idx],
    tuple_number_set[dep_idx],
    repeat_numbers_set[dep_idx],
    true,
)
@assert d.code.noise_param == dep_param
d_log = Update(d, log_param)
# Generate the trivial design
code = Code(unrotated_param, dep_param)
trivial_tuple_set = TrivialTupleSet(code)
d_triv = GenerateDesign(code, trivial_tuple_set)
@assert d_triv.code.noise_param == dep_param
d_triv_log = Update(d_triv, log_param)
# Simualte ACES for the optimised design and depolarising noise
aces_data_dep =
    SimulateACES(d, shots_set; repetitions = repetitions, seed = seed, save_data = true)
# Simualte ACES for the optimised design and log-normal noise
aces_data_log =
    SimulateACES(d_log, shots_set; repetitions = repetitions, seed = seed, save_data = true)
# Simualte ACES for the trivial design and depolarising noise
aces_data_triv_dep = SimulateACES(
    d_triv,
    shots_set;
    repetitions = repetitions,
    seed = seed,
    save_data = true,
)
# Simualte ACES for the trivial design and log-normal noise
aces_data_triv_log = SimulateACES(
    d_triv_log,
    shots_set;
    repetitions = repetitions,
    seed = seed,
    save_data = true,
)

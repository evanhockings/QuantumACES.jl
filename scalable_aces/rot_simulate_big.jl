using AveragedCircuitEigenvalueSampling, JLD2
enter_folder("scalable_aces")
# Set up the parameters
dist = 3
dist_big = 25
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
shots_set = [10^6; 10^7; 10^8; 10^9]
rotated_param = RotatedPlanarParameters(dist)
rotated_param_big = RotatedPlanarParameters(dist_big)
dep_param = DepolarisingParameters(r_1, r_2, r_m)
log_param = LogNormalParameters(r_1, r_2, r_m, total_std_log; seed = seed)
# Load the design
metadata_dict = load("data/design_metadata_$(code_filename(rotated_param)).jld2")
@assert rotated_param == metadata_dict["rotated_param"]
@assert dep_param == metadata_dict["dep_param"]
dep_param_set = metadata_dict["dep_param_set"]
tuple_number_set = metadata_dict["tuple_number_set"]
repeat_numbers_set = metadata_dict["repeat_numbers_set"]
dep_idx = 14
@assert dep_param == dep_param_set[dep_idx]
tuple_number = tuple_number_set[dep_idx]
repeat_numbers = repeat_numbers_set[dep_idx]
d = load_design(rotated_param, dep_param, tuple_number, repeat_numbers, true)
@assert d.code.noise_param == dep_param
# Generate the design at a large code distance
if isfile(
    pwd() *
    "/data/" *
    design_filename(rotated_param_big, dep_param, tuple_number, repeat_numbers, false),
)
    println("Loading the design.")
    d_big = load_design(rotated_param_big, dep_param, tuple_number, repeat_numbers, false)
else
    println("Calculating the design.")
    code_big = Code(rotated_param_big, dep_param)
    d_big = GenerateDesign(
        code_big,
        d.tuple_set_data;
        shot_weights = d.shot_weights,
        full_covariance = false,
        diagnostics = true,
        save_data = true,
    )
end
# Simualte ACES for the optimised design and depolarising noise
aces_data_dep = SimulateACES(
    d_big,
    shots_set;
    seed = seed,
    detailed_diagnostics = true,
    save_data = true,
    force_gc = true,
)
aces_data_dep = nothing
# Simualte ACES for the optimised design and log-normal noise
d_big_log = Update(d_big, log_param)
d_big = nothing
aces_data_log = SimulateACES(
    d_big_log,
    shots_set;
    seed = seed,
    detailed_diagnostics = true,
    save_data = true,
    force_gc = true,
)

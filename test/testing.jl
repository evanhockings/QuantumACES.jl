using AveragedCircuitEigenvalueSampling
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
dep_param = DepolarisingParameters(r_1, r_2, r_m)
dep_name = noise_filename(dep_param)
log_param_seeded = LognormalParameters(r_1, r_2, r_m, total_std_log; seed = seed)
log_name_seeded = noise_filename(log_param_seeded)
log_param_unseeded = LognormalParameters(r_1, r_2, r_m, total_std_log)
log_name_unseeded = noise_filename(log_param_unseeded)

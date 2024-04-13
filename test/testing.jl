using AveragedCircuitEigenvalueSampling
#
dist = 3
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
shots_set = [10^7; 2 * 10^7]
#
dep_param = DepolarisingParameters(r_1, r_2, r_m)
rotated_param = RotatedPlanarParameters(dist)
rotated_planar = Code(rotated_param, dep_param)
#
tuple_set_data = get_tuple_set_data(rotated_planar)
tuple_set_data_opt = optimise_repetitions(rotated_planar, tuple_set_data)
#
# d_rot = generate_design(rotated_planar)
# aces_data_rot = simulate_aces(d_rot, shots_set)

using AveragedCircuitEigenvalueSampling, JLD2
enter_folder("scalable_aces")
# Set up the parameters
dist = 3
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
seed = UInt(0)
rotated_param = RotatedPlanarParameters(dist)
dep_param = DepolarisingParameters(r_1, r_2, r_m)
ratio_set = [1 / 2, 1, 2]
dep_param_set = [
    DepolarisingParameters(r_1 * f_1, r_2 * f_2, r_m * f_m) for f_1 in ratio_set for
    f_2 in ratio_set for f_m in ratio_set
]
dep_param_num = length(dep_param_set)
dep_idx = 14
@assert dep_param == dep_param_set[dep_idx]
# Optimise designs for a range of noise strengths
tuple_number_set = Vector{Int}(undef, dep_param_num)
repeat_numbers_set = Vector{Vector{Int}}(undef, dep_param_num)
for (idx, param) in enumerate(dep_param_set)
    code = Code(rotated_param, param)
    d = OptimiseDesign(code; ls_type = :wls, seed = seed, save_data = true)
    tuple_number_set[idx] = length(d.tuple_set)
    repeat_numbers_set[idx] = d.tuple_set_data.repeat_numbers
end
# Optimise designs for the other least squares estimators
code = Code(rotated_param, dep_param)
d_gls = OptimiseDesign(code; ls_type = :gls, seed = seed, save_data = true)
gls_tuple_number = length(d_gls.tuple_set)
gls_repeat_numbers = d_gls.tuple_set_data.repeat_numbers
d_ols = OptimiseDesign(code; ls_type = :ols, seed = seed, save_data = true)
ols_tuple_number = length(d_ols.tuple_set)
ols_repeat_numbers = d_ols.tuple_set_data.repeat_numbers
# Save the file data
jldsave(
    pwd() * "/data/design_metadata_$(code_filename(rotated_param)).jld2";
    rotated_param,
    dep_param,
    dep_param_set,
    tuple_number_set,
    repeat_numbers_set,
    gls_tuple_number,
    gls_repeat_numbers,
    ols_tuple_number,
    ols_repeat_numbers,
)

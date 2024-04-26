using QuantumACES, JLD2
enter_folder("scalable_aces")
# Set up the parameters
dist = 3
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
seed = UInt(0)
unrotated_param = get_unrotated_param(dist)
dep_param = get_dep_param(r_1, r_2, r_m)
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
    c = get_circuit(unrotated_param, param)
    d = optimise_design(
        c;
        options = OptimOptions(; ls_type = :wls, save_data = true, seed = seed),
    )
    tuple_number_set[idx] = length(d.tuple_set)
    repeat_numbers_set[idx] = d.tuple_set_data.repeat_numbers
end
# Optimise designs for the other least squares estimators
c = get_circuit(unrotated_param, dep_param)
d_gls = optimise_design(
    c;
    options = OptimOptions(; ls_type = :gls, save_data = true, seed = seed),
)
gls_tuple_number = length(d_gls.tuple_set)
gls_repeat_numbers = d_gls.tuple_set_data.repeat_numbers
d_ols = optimise_design(
    c;
    options = OptimOptions(; ls_type = :ols, save_data = true, seed = seed),
)
ols_tuple_number = length(d_ols.tuple_set)
ols_repeat_numbers = d_ols.tuple_set_data.repeat_numbers
# Save the file data
jldsave(
    pwd() * "/data/design_metadata_$(unrotated_param.circuit_name).jld2";
    unrotated_param,
    dep_param,
    dep_param_set,
    tuple_number_set,
    repeat_numbers_set,
    gls_tuple_number,
    gls_repeat_numbers,
    ols_tuple_number,
    ols_repeat_numbers,
)

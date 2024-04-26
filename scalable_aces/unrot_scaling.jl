using QuantumACES, JLD2, Random
enter_folder("scalable_aces")
# Set up the parameters
dist = 3
dist_max = 9
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
ls_type = :wls
unrotated_param = get_unrotated_param(dist)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
unrotated_planar = get_circuit(unrotated_param, dep_param)
# Load the designs
metadata_dict = load("data/design_metadata_$(unrotated_param.circuit_name).jld2")
@assert unrotated_param == metadata_dict["unrotated_param"]
@assert dep_param == metadata_dict["dep_param"]
dep_param_set = metadata_dict["dep_param_set"]
tuple_number_set = metadata_dict["tuple_number_set"]
repeat_numbers_set = metadata_dict["repeat_numbers_set"]
gls_tuple_number = metadata_dict["gls_tuple_number"]
gls_repeat_numbers = metadata_dict["gls_repeat_numbers"]
ols_tuple_number = metadata_dict["ols_tuple_number"]
ols_repeat_numbers = metadata_dict["ols_repeat_numbers"]
dep_param_num = length(dep_param_set)
dep_idx = 14
@assert dep_param == dep_param_set[dep_idx]
dep_worst_idx = 3
repetitions = 400
d_gls = load_design(
    unrotated_param,
    dep_param,
    gls_tuple_number,
    gls_repeat_numbers,
    true,
    :gls,
)
d_wls = load_design(
    unrotated_param,
    dep_param_set[dep_idx],
    tuple_number_set[dep_idx],
    repeat_numbers_set[dep_idx],
    true,
    ls_type,
)
d_wls_worst = load_design(
    unrotated_param,
    dep_param_set[dep_worst_idx],
    tuple_number_set[dep_worst_idx],
    repeat_numbers_set[dep_worst_idx],
    true,
    ls_type,
)
d_ols = load_design(
    unrotated_param,
    dep_param,
    ols_tuple_number,
    ols_repeat_numbers,
    true,
    :ols,
)
d_basic = generate_design(unrotated_planar, get_basic_tuple_set(unrotated_planar))
# Calculate the merits over a range of random instances of log-normal noise
expectation_array = Matrix{Float64}(undef, dep_param_num, repetitions)
variance_array = Matrix{Float64}(undef, dep_param_num, repetitions)
gls_expectation_set = Vector{Float64}(undef, repetitions)
gls_variance_set = Vector{Float64}(undef, repetitions)
ols_expectation_set = Vector{Float64}(undef, repetitions)
ols_variance_set = Vector{Float64}(undef, repetitions)
# Set the random seeds
Random.seed!(seed)
seeds = rand(UInt64, repetitions)
Random.seed!()
start_time = time()
for idx in 1:dep_param_num
    d = load_design(
        unrotated_param,
        dep_param_set[idx],
        tuple_number_set[idx],
        repeat_numbers_set[idx],
        true,
        ls_type,
    )
    for rep in 1:repetitions
        log_param_rep = get_log_param(
            log_param.r_1,
            log_param.r_2,
            log_param.r_m,
            log_param.total_std_log;
            seed = seeds[rep],
        )
        d_log = update_noise(d, log_param_rep)
        covariance_log = calc_covariance_log(d_log)
        (expectation, variance) = calc_ls_moments(d_log, covariance_log, :wls)
        expectation_array[idx, rep] = expectation
        variance_array[idx, rep] = variance
    end
    println(
        "Calculated the merits for design $(idx) of $(dep_param_num). The time elapsed since starting is $(round(time() - start_time, digits = 3)) s.",
    )
end
for rep in 1:repetitions
    log_param_rep = get_log_param(
        log_param.r_1,
        log_param.r_2,
        log_param.r_m,
        log_param.total_std_log;
        seed = seeds[rep],
    )
    d_gls_log = update_noise(d_gls, log_param_rep)
    gls_covariance_log = calc_covariance_log(d_gls_log)
    (gls_expectation, gls_variance) = calc_ls_moments(d_gls_log, gls_covariance_log, :gls)
    gls_expectation_set[rep] = gls_expectation
    gls_variance_set[rep] = gls_variance
    d_ols_log = update_noise(d_ols, log_param_rep)
    ols_covariance_log = calc_covariance_log(d_ols_log)
    (ols_expectation, ols_variance) = calc_ls_moments(d_ols_log, ols_covariance_log, :ols)
    ols_expectation_set[rep] = ols_expectation
    ols_variance_set[rep] = ols_variance
end
# Save the merit data
jldsave(
    pwd() *
    "/data/design_merit_data_$(unrotated_param.circuit_name)_$(log_param.noise_name).jld2";
    expectation_array,
    variance_array,
    gls_expectation_set,
    gls_variance_set,
    ols_expectation_set,
    ols_variance_set,
)
# Load the design and calculate the depolarising and log-normal noise scaling data
# Optimised WLS design
@assert d_wls.c.noise_param == dep_param
dep_planar_scaling_wls =
    calc_depolarising_planar_scaling(d_wls, dist_max; ls_type = :wls, save_data = true)
d_wls_log = update_noise(d_wls, log_param)
log_planar_scaling_wls = calc_lognormal_planar_scaling(
    d_wls_log,
    dist_max;
    ls_type = :wls,
    seed = seed,
    save_data = true,
)
# Load the design and calculate the depolarising noise scaling data
# Optimised GLS design
@assert d_gls.c.noise_param == dep_param
dep_planar_scaling_gls =
    calc_depolarising_planar_scaling(d_gls, dist_max; ls_type = :gls, save_data = true)
# Optimised OLS design
@assert d_ols.c.noise_param == dep_param
dep_planar_scaling_ols =
    calc_depolarising_planar_scaling(d_ols, dist_max; ls_type = :ols, save_data = true)
# Badly optimised WLS design
d_wls_worst = update_noise(d_wls_worst, dep_param)
@assert d_wls_worst.c.noise_param == dep_param
dep_planar_scaling_wls_worst = calc_depolarising_planar_scaling(
    d_wls_worst,
    dist_max;
    ls_type = :wls,
    save_data = true,
)
# Trivial design
@assert d_basic.c.noise_param == dep_param
dep_planar_scaling_basic =
    calc_depolarising_planar_scaling(d_basic, dist_max; ls_type = :wls, save_data = true)

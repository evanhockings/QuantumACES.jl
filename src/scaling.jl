"""
    DepolarisingScaling(d::Design, dist_range::Vector{Int}; ls_type::Symbol = :none, save_data::Bool = false, diagnostics::Bool = true)

Calculate the merit of the design for the supplied code distances to examine its scaling.
"""
function DepolarisingScaling(
    d::Design,
    dist_range::Vector{Int};
    ls_type::Symbol = :none,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    # Check the parameters
    @assert typeof(d.code.noise_param) == DepolarisingParameters "This function requires depolarising noise."
    @assert minimum(dist_range) >= 3 "The supplied distances must all be at least 3."
    # Set some variables
    code_param = d.code.code_param
    noise_param = d.code.noise_param
    tuple_set_data = d.tuple_set_data
    tuple_set = TupleSet(tuple_set_data)
    shot_weights = d.shot_weights
    if ls_type == :none
        if d.ls_type == :none
            ls_type = :wls
        else
            ls_type = d.ls_type
        end
    end
    # Initialise the merit scaling
    merit_scaling = Vector{Merit}(undef, length(dist_range))
    # Calculate the figure of merit scaling with code distance
    start_time = time()
    calculation_times = Matrix{Float64}(undef, length(dist_range), 2)
    for (idx, dist) in enumerate(dist_range)
        # Initialise up the code
        code_param_dist = deepcopy(code_param)
        if typeof(code_param_dist) == RotatedPlanarParameters
            @reset code_param_dist.vertical_dist = dist
            @reset code_param_dist.horizontal_dist = dist
        elseif typeof(code_param_dist) == UnrotatedPlanarParameters
            @reset code_param_dist.vertical_dist = dist
            @reset code_param_dist.horizontal_dist = dist
        else
            throw(error("Unsupported code type $(code_type)."))
        end
        code = Code(code_param_dist, noise_param)
        # Generate the design
        time_1 = time()
        d_dist = GenerateDesign(code, tuple_set_data; shot_weights = shot_weights)
        time_2 = time()
        design_time = time_2 - time_1
        if diagnostics
            println(
                "Generating the design at distance $(dist) took $(round(design_time, digits = 3)) s.",
            )
        end
        # Calculate the merit
        merit_dist = LSMerit(d_dist, ls_type)
        merit_scaling[idx] = merit_dist
        time_3 = time()
        merit_time = time_3 - time_2
        calculation_times[idx, :] = [design_time; merit_time]
        if diagnostics
            println(
                "Calculating the merit at distance $(dist) took $(round(merit_time, digits = 3)) s.",
            )
        end
    end
    # Save and return the results
    overall_time = time() - start_time
    dep_scaling_data = DepolarisingScalingData(
        merit_scaling,
        dist_range,
        code_param,
        noise_param,
        tuple_set,
        tuple_set_data,
        shot_weights,
        ls_type,
        calculation_times,
        overall_time,
    )
    if save_data
        save_scaling(dep_scaling_data)
    end
    if diagnostics
        println(
            "Finished calculating merit scaling with distance for depolarising noise. The time elapsed since calculations started is $(round(overall_time, digits = 3)) s.",
        )
    end
    return dep_scaling_data::DepolarisingScalingData
end

"""
    DepolarisingScaling(d::Design, dist_max::Int; ls_type::Symbol = :none, save_data::Bool = false, diagnostics::Bool = true)

Calculate the merit of the design for code distances from 3 to the supplied maximum to examine its scaling.
"""
function DepolarisingScaling(
    d::Design,
    dist_max::Int;
    ls_type::Symbol = :none,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    # Generate the distance range
    dist_range = collect(3:dist_max)
    # Calculate the scaling data
    dep_scaling_data = DepolarisingScaling(
        d,
        dist_range;
        ls_type = ls_type,
        save_data = save_data,
        diagnostics = diagnostics,
    )
    return dep_scaling_data::DepolarisingScalingData
end

#
function DepolarisingFits(dep_scaling_data::DepolarisingScalingData)
    # Set up variables
    dist_range = dep_scaling_data.dist_range
    G_scaling = [merit.G for merit in dep_scaling_data.merit_scaling]
    N_scaling = [merit.N for merit in dep_scaling_data.merit_scaling]
    trace_scaling = [sum(merit.eigenvalues) for merit in dep_scaling_data.merit_scaling]
    trace_sq_scaling =
        [sum(merit.eigenvalues .^ 2) for merit in dep_scaling_data.merit_scaling]
    expectation_scaling = [merit.expectation for merit in dep_scaling_data.merit_scaling]
    variance_scaling = [merit.variance for merit in dep_scaling_data.merit_scaling]
    # Initialise the quadratic model
    @. quadratic(d, c) = c[1] + c[2] * d + c[3] * d^2
    # Fit the gate number
    gate_model = lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = G_scaling, x = dist_range))
    gate_params = round.(Int, coef(gate_model))
    gate_number(d) = quadratic(d, gate_params)
    # Test the fit
    @assert gate_params ≈ coef(gate_model) "The coefficients of the quadratic fit of the gate numbers are not integers."
    @assert gate_number(dist_range) ≈ G_scaling "The gate numbers are not well-fit by a quadratic."
    # Fit the gate eigenvalue number
    gate_eigenvalue_model =
        lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = N_scaling, x = dist_range))
    gate_eigenvalue_params = round.(Int, coef(gate_eigenvalue_model))
    gate_eigenvalue_number(d) = quadratic(d, gate_eigenvalue_params)
    # Test the fit
    @assert gate_eigenvalue_params ≈ coef(gate_eigenvalue_model) "The coefficients of the quadratic fit of the gate eigenvalue numbers are not integers."
    @assert gate_eigenvalue_number(dist_range) ≈ N_scaling "The gate eigenvalue numbers are not well-fit by a quadratic."
    # Fit the trace of the gate eigenvalue estimator covariance matrix
    trace_model =
        lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = trace_scaling, x = dist_range))
    trace_params = coef(trace_model)
    trace_fit(d) = quadratic(d, trace_params)
    # Test the fit
    if ~(isapprox(trace_fit(dist_range), trace_scaling; rtol = 1e-3))
        @warn "The traces are not well-fit by a quadratic."
    end
    # Fit the trace of the square of the gate eigenvalue estimator covariance matrix
    trace_sq_model =
        lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = trace_sq_scaling, x = dist_range))
    trace_sq_params = coef(trace_sq_model)
    trace_sq_fit(d) = quadratic(d, trace_sq_params)
    # Test the fit
    if ~(isapprox(trace_sq_fit(dist_range), trace_sq_scaling; rtol = 1e-3))
        @warn "The traces of the square are not well-fit by a quadratic."
    end
    # Fit the NRMSE expectation and variance
    @. expectation_fit(d) =
        sqrt(trace_fit(d) / gate_eigenvalue_number(d)) *
        (1 - (trace_sq_fit(d) / (4 * trace_fit(d)^2)))
    @. variance_fit(d) =
        (trace_sq_fit(d) / (2 * gate_eigenvalue_number(d) * trace_fit(d))) *
        (1 - (trace_sq_fit(d) / (8 * trace_fit(d)^2)))
    # Test the fits
    if ~(isapprox(expectation_fit(dist_range), expectation_scaling; rtol = 1e-3))
        @warn "The NRMSE expectations are not well-fit."
    end
    if ~(isapprox(variance_fit(dist_range), variance_scaling; rtol = 1e-3))
        @warn "The NRMSE variances are not well-fit."
    end
    # Return the functions
    return (
        gate_number::Function,
        gate_params::Vector{Int},
        gate_eigenvalue_number::Function,
        gate_eigenvalue_params::Vector{Int},
        trace_fit::Function,
        trace_params::Vector{Float64},
        trace_sq_fit::Function,
        trace_sq_params::Vector{Float64},
        expectation_fit::Function,
        variance_fit::Function,
    )
end

"""
    LogNormalScaling(d::Design, dist_range::Vector{Int}; ls_type::Symbol = :none, precision::Float64 = 1e-3, max_repetitions::Int = 10000, min_repetitions::Int = 100, print_repetitions::Int = 100, seed::Union{UInt64, Nothing} = nothing, save_data::Bool = false, diagnostics::Bool = true)

Calculate the merit of the design for the supplied code distances to examine its scaling.
"""
function LogNormalScaling(
    d::Design,
    dist_range::Vector{Int};
    ls_type::Symbol = :none,
    precision::Float64 = 1e-3,
    max_repetitions::Int = 10000,
    min_repetitions::Int = 100,
    print_repetitions::Int = 100,
    seed::Union{UInt64, Nothing} = nothing,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    # Check the parameters
    @assert typeof(d.code.noise_param) == LogNormalParameters "This function requires log-normal Pauli noise."
    @assert minimum(dist_range) >= 3 "The supplied distances must all be at least 3."
    @assert precision > 0 "The precision must be positive."
    @assert max_repetitions > 0 "The maximum number of repetitions must be positive."
    @assert min_repetitions > 0 "The minimum number of repetitions must be positive."
    @assert max_repetitions >= min_repetitions "The maximum number of repetitions must be greater than or equal to the minimum number of repetitions."
    @assert print_repetitions > 0 "The number of repetitions between printing must be positive."
    # Set some variables
    code_param = d.code.code_param
    noise_param = d.code.noise_param
    r_1 = noise_param.r_1
    r_2 = noise_param.r_2
    r_m = noise_param.r_m
    total_std_log = noise_param.total_std_log
    tuple_set_data = d.tuple_set_data
    tuple_set = TupleSet(tuple_set_data)
    shot_weights = d.shot_weights
    if ls_type == :none
        if d.ls_type == :none
            ls_type = :wls
        else
            ls_type = d.ls_type
        end
    end
    # Set the random seeds
    if seed !== nothing
        Random.seed!(seed)
    end
    seeds = rand(UInt64, max_repetitions)
    if seed !== nothing
        Random.seed!()
    end
    # Initialise the variable scalings
    expectation_scaling = Vector{Vector{Float64}}(undef, length(dist_range))
    variance_scaling = Vector{Vector{Float64}}(undef, length(dist_range))
    eigenvalues_scaling = Vector{Vector{Vector{Float64}}}(undef, length(dist_range))
    # Calculate the figure of merit scaling with code distance
    start_time = time()
    calculation_times = Matrix{Float64}(undef, length(dist_range), 2)
    for (idx, dist) in enumerate(dist_range)
        # Initialise up the code
        code_param_dist = deepcopy(code_param)
        if typeof(code_param_dist) == RotatedPlanarParameters
            @reset code_param_dist.vertical_dist = dist
            @reset code_param_dist.horizontal_dist = dist
        elseif typeof(code_param_dist) == UnrotatedPlanarParameters
            @reset code_param_dist.vertical_dist = dist
            @reset code_param_dist.horizontal_dist = dist
        else
            throw(error("Unsupported code type $(code_type)."))
        end
        code = Code(code_param_dist, noise_param)
        # Generate the design
        time_1 = time()
        d_dist = GenerateDesign(code, tuple_set_data; shot_weights = shot_weights)
        time_2 = time()
        design_time = time_2 - time_1
        if diagnostics
            println(
                "Generating the design at distance $(dist) took $(round(design_time, digits = 3)) s.",
            )
        end
        # Calculate the merit
        expectation_scaling[idx] = Vector{Float64}(undef, 0)
        variance_scaling[idx] = Vector{Float64}(undef, 0)
        eigenvalues_scaling[idx] = Vector{Vector{Float64}}(undef, 0)
        rep = 1
        generating = true
        while generating
            # Generate the noise and design
            noise_param_rep =
                LogNormalParameters(r_1, r_2, r_m, total_std_log; seed = seeds[rep])
            d_rep = Update(d_dist, noise_param_rep)
            # Calculate the variables
            covariance_log_rep = MeritData(d_rep)
            gate_eigenvalues_cov = LSCovariance(d_rep, covariance_log_rep, ls_type)
            eigenvalues_rep = eigvals(gate_eigenvalues_cov)
            (expectation_rep, variance_rep) = NRMSEMoments(eigenvalues_rep)
            push!(expectation_scaling[idx], expectation_rep)
            push!(variance_scaling[idx], variance_rep)
            push!(eigenvalues_scaling[idx], eigenvalues_rep)
            if diagnostics && rep % print_repetitions == 0
                println(
                    "Calculated $(rep) merit variables at distance $(dist). The time elapsed is $(round(time() - time_2, digits = 3)) s.",
                )
            end
            # Check if the scaling has converged
            if (
                rep >= min_repetitions &&
                std(expectation_scaling[idx]) / sqrt(rep) < precision
            ) || rep >= max_repetitions
                generating = false
            else
                rep += 1
            end
        end
        time_3 = time()
        merit_time = time_3 - time_2
        calculation_times[idx, :] = [design_time; merit_time]
        if diagnostics
            println(
                "Calculating the $(rep) merit variables at distance $(dist) took $(round(merit_time, digits = 3)) s.",
            )
        end
    end
    # Save and return the results
    overall_time = time() - start_time
    log_scaling_data = LogNormalScalingData(
        expectation_scaling,
        variance_scaling,
        eigenvalues_scaling,
        dist_range,
        code_param,
        noise_param,
        seeds,
        tuple_set,
        tuple_set_data,
        shot_weights,
        ls_type,
        calculation_times,
        overall_time,
    )
    if save_data
        save_scaling(log_scaling_data)
    end
    if diagnostics
        println(
            "Finished calculating merit variable scaling with distance for log-normal Pauli noise. The time elapsed since calculations started is $(round(overall_time, digits = 3)) s.",
        )
    end
    return log_scaling_data::LogNormalScalingData
end

function LogNormalScaling(
    d::Design,
    dist_max::Int;
    ls_type::Symbol = :none,
    precision::Float64 = 1e-3,
    max_repetitions::Int = 10000,
    min_repetitions::Int = 100,
    seed::Union{UInt64, Nothing} = nothing,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    # Generate the distance range
    dist_range = collect(3:dist_max)
    # Calculate the scaling data
    log_scaling_data = LogNormalScaling(
        d,
        dist_range;
        ls_type = ls_type,
        precision = precision,
        max_repetitions = max_repetitions,
        min_repetitions = min_repetitions,
        seed = seed,
        save_data = save_data,
        diagnostics = diagnostics,
    )
    return log_scaling_data::LogNormalScalingData
end

#
function LogNormalFits(
    log_scaling_data::LogNormalScalingData,
    gate_eigenvalue_number::Function,
    trace_params_init::Vector{Float64},
    trace_sq_params_init::Vector{Float64},
)
    # Set up variables
    dist_range = log_scaling_data.dist_range
    expectation_scaling = log_scaling_data.expectation_scaling
    variance_scaling = log_scaling_data.variance_scaling
    repetitions = length.(expectation_scaling)
    @assert repetitions == length.(variance_scaling) "The number of repetitions for the expectation and variance scaling data do not match."
    # Create the quadratic models for the NRMSE expectation and variance
    @. expectation_model(d, c) =
        sqrt((c[1] + c[2] * d + c[3] * d^2) / gate_eigenvalue_number(d)) *
        (1 - ((c[4] + c[5] * d + c[6] * d^2) / (4 * (c[1] + c[2] * d + c[3] * d^2)^2)))
    @. variance_model(d, c) =
        (
            (c[4] + c[5] * d + c[6] * d^2) /
            (2 * gate_eigenvalue_number(d) * (c[1] + c[2] * d + c[3] * d^2))
        ) * (1 - ((c[4] + c[5] * d + c[6] * d^2) / (8 * (c[1] + c[2] * d + c[3] * d^2)^2)))
    # Simultaneously fit the mean of the NRMSE expectation and variance across the instances of log-normal noise
    mean_expectation_scaling = mean.(expectation_scaling)
    mean_variance_scaling = mean.(variance_scaling)
    # Rescale the variance so it is considered appropriately when fitting
    pair_rescale = mean_expectation_scaling[1] / sqrt(mean_variance_scaling[1])
    pair_scaling = [mean_expectation_scaling; pair_rescale * sqrt.(mean_variance_scaling)]
    pair_model(d, c) = [expectation_model(d, c); pair_rescale * sqrt.(variance_model(d, c))]
    param_init = [trace_params_init; trace_sq_params_init]
    pair_fit = curve_fit(pair_model, dist_range, pair_scaling, param_init)
    # Fit the NRMSE expectation
    expectation_fit_log(d) = expectation_model(d, pair_fit.param)
    if ~(isapprox(expectation_fit_log(dist_range), mean_expectation_scaling; rtol = 1e-3))
        @warn "The mean NRMSE expectations are not well-fit."
    end
    # Fit the NRMSE variance
    variance_fit_log(d) = variance_model(d, pair_fit.param)
    if ~(isapprox(variance_fit_log(dist_range), mean_variance_scaling; rtol = 1e-3))
        @warn "The mean NRMSE variances are not well-fit."
    end
    # Return the functions
    return (expectation_fit_log::Function, variance_fit_log::Function)
end

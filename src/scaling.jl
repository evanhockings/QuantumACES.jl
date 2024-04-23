struct DepolarisingPlanarScaling <: AbstractScalingData
    # Code distances
    dist_range::Vector{Int}
    # Merit of the design for a range of code distances
    merit_scaling::Vector{Merit}
    # Gate number fit
    G_fit::Function
    # Gate number fit parameters
    # a + bd + cd^2
    G_params::Vector{Int}
    # Gate eigenvalue number fit
    N_fit::Function
    # Gate eigenvalue number fit parameters
    # a + bd + cd^2
    N_params::Vector{Int}
    # Trace of the gate eigenvalue estimator covariance matrix fit
    trace_fit::Function
    # Trace of the gate eigenvalue estimator covariance matrix fit parameters 
    # a + bd + cd^2
    trace_params::Vector{Float64}
    # Trace of the gate eigenvalue estimator covariance matrix squared fit
    trace_sq_fit::Function
    # Trace of the gate eigenvalue estimator covariance matrix squared fit parameters
    # a + bd + cd^2
    trace_sq_params::Vector{Float64}
    # NRMSE expectation fit
    expectation_fit::Function
    # NRMSE variance fit
    variance_fit::Function
    # Code parameters
    circuit_param::AbstractCircuitParameters
    # Depolarising noise parameters
    noise_param::DepolarisingParameters
    # Circuit rearrangements used to generate the design matrix
    tuple_set::Vector{Vector{Int}}
    # Data used to generate the tuple set
    tuple_set_data::TupleSetData
    # Weighting of the shots allocated to each tuple
    shot_weights::Vector{Float64}
    # Type of least squares estimator for which the merits were calculated
    ls_type::Symbol
    # The time taken to generate the design and calculate the merit for each distance
    # (design_time, merit_time)
    calculation_times::Matrix{Float64}
    # The overall time taken to calculate the merit scaling for depolarising noise
    overall_time::Float64
end

function Base.show(io::IO, s::DepolarisingPlanarScaling)
    return print(
        io,
        "Merit scaling data with depolarising noise of a design for a $(s.circuit_param.circuit_name) code with $(length(s.tuple_set)) tuples.",
    )
end

struct LognormalPlanarScaling <: AbstractScalingData
    # Code distances
    dist_range::Vector{Int}
    # Gate eigenvalue number fit
    N_fit::Function
    # Gate eigenvalue number fit parameters
    # a + bd + cd^2
    N_params::Vector{Int}
    # Expected NRMSE for a range of code distances
    expectation_scaling::Vector{Vector{Float64}}
    # Average expected NRMSE fit
    expectation_fit::Function
    # NRMSE variance for a range of code distances
    variance_scaling::Vector{Vector{Float64}}
    # Average NRMSE variance fit
    variance_fit::Function
    # Eigenvalues of the gate log-eigenvalue estimator covariance matrix for a range of code distances
    eigenvalues_scaling::Vector{Vector{Vector{Float64}}}
    # Code parameters
    circuit_param::AbstractCircuitParameters
    # Log-normal random noise parameters
    noise_param::LognormalParameters
    # Random seeds for the noise parameters
    seeds::Vector{UInt64}
    # Circuit rearrangements used to generate the design matrix
    tuple_set::Vector{Vector{Int}}
    # Data used to generate the tuple set
    tuple_set_data::TupleSetData
    # Weighting of the shots allocated to each tuple
    shot_weights::Vector{Float64}
    # Type of least squares estimator for which the merits were calculated
    ls_type::Symbol
    # The time taken to generate the design and calculate the merits for each distance
    # (design_time, log_merit_time, dep_merit_time)
    calculation_times::Matrix{Float64}
    # The overall time taken to calculate the merit scaling for log-normal random noise
    overall_time::Float64
end

function Base.show(io::IO, s::LognormalPlanarScaling)
    return print(
        io,
        "Merit scaling data with log-normal random noise of a design for a $(s.circuit_param.circuit_name) code with $(length(s.tuple_set)) tuples.",
    )
end

"""
    calc_depolarising_planar_scaling(d::Design, dist_range::Vector{Int}; ls_type::Symbol = :none, save_data::Bool = false, diagnostics::Bool = true)

Calculate the merit of the design for the supplied code distances to examine its scaling.
"""
function calc_depolarising_planar_scaling(
    d::Design,
    dist_range::Vector{Int};
    ls_type::Symbol = :none,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    # Check the parameters
    @assert typeof(d.c.noise_param) == DepolarisingParameters "This function requires depolarising noise."
    @assert (
        typeof(d.c.circuit_param) == RotatedPlanarParameters &&
        typeof(d.c) == RotatedPlanarCircuit
    ) || (
        typeof(d.c.circuit_param) == UnrotatedPlanarParameters &&
        typeof(d.c) == UnrotatedPlanarCircuit
    ) "This function requires planar codes."
    @assert minimum(dist_range) >= 3 "The supplied distances must all be at least 3."
    # Set some variables
    circuit_param = d.c.circuit_param
    noise_param = d.c.noise_param
    tuple_set_data = d.tuple_set_data
    tuple_set = get_tuple_set(tuple_set_data)
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
        # Initialise the circuit
        if typeof(circuit_param) == RotatedPlanarParameters
            circuit_param_dist = get_rotated_param(
                dist;
                check_type = circuit_param.check_type,
                gate_type = circuit_param.gate_type,
                dynamically_decouple = circuit_param.dynamically_decouple,
                pad_identity = circuit_param.pad_identity,
                single_qubit_time = circuit_param.layer_time_dict[:single_qubit],
                two_qubit_time = circuit_param.layer_time_dict[:two_qubit],
                dynamical_decoupling_time = circuit_param.layer_time_dict[:dynamical],
                meas_reset_time = circuit_param.layer_time_dict[:meas_reset],
            )
        elseif typeof(circuit_param) == UnrotatedPlanarParameters
            circuit_param_dist = get_unrotated_param(
                dist;
                gate_type = circuit_param.gate_type,
                pad_identity = circuit_param.pad_identity,
                single_qubit_time = circuit_param.layer_time_dict[:single_qubit],
                two_qubit_time = circuit_param.layer_time_dict[:two_qubit],
                meas_reset_time = circuit_param.layer_time_dict[:meas_reset],
            )
        else
            throw(error("Unsupported circuit type $(typeof(circuit_param))."))
        end
        c = get_circuit(circuit_param_dist, noise_param)
        # Generate the design
        time_1 = time()
        d_dist = generate_design(c, tuple_set_data; shot_weights = shot_weights)
        time_2 = time()
        design_time = time_2 - time_1
        if diagnostics
            println(
                "Generating the design at distance $(dist) took $(round(design_time, digits = 3)) s.",
            )
        end
        # Calculate the merit
        merit_dist = calc_ls_merit(d_dist, ls_type)
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
    # Fit the trends
    # Set up variables
    G_scaling = [merit.G for merit in merit_scaling]
    N_scaling = [merit.N for merit in merit_scaling]
    trace_scaling = [sum(merit.eigenvalues) for merit in merit_scaling]
    trace_sq_scaling = [sum(merit.eigenvalues .^ 2) for merit in merit_scaling]
    expectation_scaling = [merit.expectation for merit in merit_scaling]
    variance_scaling = [merit.variance for merit in merit_scaling]
    # Initialise the quadratic model
    @. quadratic(d, c) = c[1] + c[2] * d + c[3] * d^2
    # Fit the gate number
    G_model = lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = G_scaling, x = dist_range))
    G_params = round.(Int, coef(G_model))
    G_fit(d) = quadratic(d, G_params)
    @assert G_params ≈ coef(G_model) "The coefficients of the quadratic fit of the gate numbers are not integers."
    @assert G_fit(dist_range) ≈ G_scaling "The gate numbers are not well-fit by a quadratic."
    # Fit the gate eigenvalue number
    N_model = lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = N_scaling, x = dist_range))
    N_params = round.(Int, coef(N_model))
    N_fit(d) = quadratic(d, N_params)
    @assert N_params ≈ coef(N_model) "The coefficients of the quadratic fit of the gate eigenvalue numbers are not integers."
    @assert N_fit(dist_range) ≈ N_scaling "The gate eigenvalue numbers are not well-fit by a quadratic."
    # Fit the trace of the gate eigenvalue estimator covariance matrix
    trace_model =
        lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = trace_scaling, x = dist_range))
    trace_params = coef(trace_model)
    trace_fit(d) = quadratic(d, trace_params)
    if ~(isapprox(trace_fit(dist_range), trace_scaling; rtol = 1e-3))
        @warn "The traces are not well-fit by a quadratic."
    end
    # Fit the trace of the square of the gate eigenvalue estimator covariance matrix
    trace_sq_model =
        lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = trace_sq_scaling, x = dist_range))
    trace_sq_params = coef(trace_sq_model)
    trace_sq_fit(d) = quadratic(d, trace_sq_params)
    if ~(isapprox(trace_sq_fit(dist_range), trace_sq_scaling; rtol = 1e-3))
        @warn "The traces of the square are not well-fit by a quadratic."
    end
    # Fit the NRMSE expectation and variance
    @. expectation_fit(d) =
        sqrt(trace_fit(d) / N_fit(d)) * (1 - (trace_sq_fit(d) / (4 * trace_fit(d)^2)))
    @. variance_fit(d) =
        (trace_sq_fit(d) / (2 * N_fit(d) * trace_fit(d))) *
        (1 - (trace_sq_fit(d) / (8 * trace_fit(d)^2)))
    if ~(isapprox(expectation_fit(dist_range), expectation_scaling; rtol = 1e-3))
        @warn "The NRMSE expectations are not well-fit."
    end
    if ~(isapprox(variance_fit(dist_range), variance_scaling; rtol = 1e-3))
        @warn "The NRMSE variances are not well-fit."
    end
    # Save and return the results
    overall_time = time() - start_time
    dep_planar_scaling = DepolarisingPlanarScaling(
        dist_range,
        merit_scaling,
        G_fit,
        G_params,
        N_fit,
        N_params,
        trace_fit,
        trace_params,
        trace_sq_fit,
        trace_sq_params,
        expectation_fit,
        variance_fit,
        circuit_param,
        noise_param,
        tuple_set,
        tuple_set_data,
        shot_weights,
        ls_type,
        calculation_times,
        overall_time,
    )
    if save_data
        save_scaling(dep_planar_scaling)
    end
    if diagnostics
        println(
            "Finished calculating merit scaling with distance for depolarising noise. The time elapsed since calculations started is $(round(overall_time, digits = 3)) s.",
        )
    end
    return dep_planar_scaling::DepolarisingPlanarScaling
end

"""
    calc_depolarising_planar_scaling(d::Design, dist_max::Int; ls_type::Symbol = :none, save_data::Bool = false, diagnostics::Bool = true)

Calculate the merit of the design for code distances from 3 to the supplied maximum to examine its scaling.
"""
function calc_depolarising_planar_scaling(
    d::Design,
    dist_max::Int;
    ls_type::Symbol = :none,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    # Generate the distance range
    if typeof(d.c.circuit_param) == RotatedPlanarParameters
        dist_range = collect(3:dist_max)
    elseif typeof(d.c.circuit_param) == UnrotatedPlanarParameters
        dist_range = collect(3:dist_max)
    else
        throw(error("Unsupported circuit type $(typeof(d.c.circuit_param))."))
    end
    # Calculate the scaling data
    dep_planar_scaling = calc_depolarising_planar_scaling(
        d,
        dist_range;
        ls_type = ls_type,
        save_data = save_data,
        diagnostics = diagnostics,
    )
    return dep_planar_scaling::DepolarisingPlanarScaling
end

"""
    calc_lognormal_planar_scaling(d::Design, dist_range::Vector{Int}; ls_type::Symbol = :none, precision::Float64 = 1e-3, max_repetitions::Int = 10000, min_repetitions::Int = 100, print_repetitions::Int = 100, seed::Union{UInt64, Nothing} = nothing, save_data::Bool = false, diagnostics::Bool = true)

Calculate the merit of the design for the supplied code distances to examine its scaling.
"""
function calc_lognormal_planar_scaling(
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
    @assert typeof(d.c.noise_param) == LognormalParameters "This function requires log-normal Pauli noise."
    @assert (
        typeof(d.c.circuit_param) == RotatedPlanarParameters &&
        typeof(d.c) == RotatedPlanarCircuit
    ) || (
        typeof(d.c.circuit_param) == UnrotatedPlanarParameters &&
        typeof(d.c) == UnrotatedPlanarCircuit
    ) "This function requires planar codes."
    @assert minimum(dist_range) >= 3 "The supplied distances must all be at least 3."
    @assert precision > 0 "The precision must be positive."
    @assert max_repetitions > 0 "The maximum number of repetitions must be positive."
    @assert min_repetitions > 0 "The minimum number of repetitions must be positive."
    @assert max_repetitions >= min_repetitions "The maximum number of repetitions must be greater than or equal to the minimum number of repetitions."
    @assert print_repetitions > 0 "The number of repetitions between printing must be positive."
    # Set some variables
    circuit_param = d.c.circuit_param
    noise_param = d.c.noise_param
    r_1 = noise_param.r_1
    r_2 = noise_param.r_2
    r_m = noise_param.r_m
    total_std_log = noise_param.total_std_log
    dep_noise_param = get_dep_param(r_1, r_2, r_m)
    tuple_set_data = d.tuple_set_data
    tuple_set = get_tuple_set(tuple_set_data)
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
    N_scaling = Vector{Int}(undef, length(dist_range))
    expectation_scaling = Vector{Vector{Float64}}(undef, length(dist_range))
    variance_scaling = Vector{Vector{Float64}}(undef, length(dist_range))
    eigenvalues_scaling = Vector{Vector{Vector{Float64}}}(undef, length(dist_range))
    dep_trace_scaling = Vector{Float64}(undef, length(dist_range))
    dep_trace_sq_scaling = Vector{Float64}(undef, length(dist_range))
    # Calculate the figure of merit scaling with code distance
    start_time = time()
    calculation_times = Matrix{Float64}(undef, length(dist_range), 3)
    for (idx, dist) in enumerate(dist_range)
        # Initialise the circuit
        if typeof(circuit_param) == RotatedPlanarParameters
            circuit_param_dist = get_rotated_param(
                dist;
                check_type = circuit_param.check_type,
                gate_type = circuit_param.gate_type,
                dynamically_decouple = circuit_param.dynamically_decouple,
                pad_identity = circuit_param.pad_identity,
                single_qubit_time = circuit_param.layer_time_dict[:single_qubit],
                two_qubit_time = circuit_param.layer_time_dict[:two_qubit],
                dynamical_decoupling_time = circuit_param.layer_time_dict[:dynamical],
                meas_reset_time = circuit_param.layer_time_dict[:meas_reset],
            )
        elseif typeof(circuit_param) == UnrotatedPlanarParameters
            circuit_param_dist = get_unrotated_param(
                dist;
                gate_type = circuit_param.gate_type,
                pad_identity = circuit_param.pad_identity,
                single_qubit_time = circuit_param.layer_time_dict[:single_qubit],
                two_qubit_time = circuit_param.layer_time_dict[:two_qubit],
                meas_reset_time = circuit_param.layer_time_dict[:meas_reset],
            )
        else
            throw(error("Unsupported circuit type $(typeof(circuit_param))."))
        end
        c = get_circuit(circuit_param_dist, noise_param)
        N_scaling[idx] = c.N
        # Generate the design
        time_1 = time()
        d_dist = generate_design(c, tuple_set_data; shot_weights = shot_weights)
        time_2 = time()
        design_time = time_2 - time_1
        if diagnostics
            println(
                "Generating the design at distance $(dist) took $(round(design_time, digits = 3)) s.",
            )
        end
        # Calculate the lognormal noise merit
        expectation_scaling[idx] = Vector{Float64}(undef, 0)
        variance_scaling[idx] = Vector{Float64}(undef, 0)
        eigenvalues_scaling[idx] = Vector{Vector{Float64}}(undef, 0)
        rep = 1
        generating = true
        while generating
            # Generate the noise and design
            noise_param_rep = get_log_param(r_1, r_2, r_m, total_std_log; seed = seeds[rep])
            d_rep = update_noise(d_dist, noise_param_rep)
            # Calculate the variables
            covariance_log_rep = calc_covariance_log(d_rep)
            gate_eigenvalues_cov = calc_ls_covariance(d_rep, covariance_log_rep, ls_type)
            eigenvalues_rep = eigvals(gate_eigenvalues_cov)
            (expectation_rep, variance_rep) = nrmse_moments(eigenvalues_rep)
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
        log_merit_time = time_3 - time_2
        if diagnostics
            println(
                "Calculating the $(rep) merit variables at distance $(dist) took $(round(log_merit_time, digits = 3)) s.",
            )
        end
        # Calculate the depolarising noise trace scaling
        d_dep = update_noise(d_dist, dep_noise_param)
        covariance_log_dep = calc_covariance_log(d_dep)
        gate_eigenvalues_cov_dep = calc_ls_covariance(d_dep, covariance_log_dep, ls_type)
        eigenvalues_dep = eigvals(gate_eigenvalues_cov_dep)
        dep_trace_scaling[idx] = sum(eigenvalues_dep)
        dep_trace_sq_scaling[idx] = sum(eigenvalues_dep .^ 2)
        time_4 = time()
        dep_merit_time = time_4 - time_3
        calculation_times[idx, :] = [design_time; log_merit_time; dep_merit_time]
    end
    # Fit the trends
    # Set up variables
    repetitions = length.(expectation_scaling)
    @assert repetitions == length.(variance_scaling) "The number of repetitions for the expectation and variance scaling data do not match."
    # Initialise the quadratic model
    @. quadratic(d, c) = c[1] + c[2] * d + c[3] * d^2
    # Fit the gate eigenvalue number
    N_model = lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = N_scaling, x = dist_range))
    N_params = round.(Int, coef(N_model))
    N_fit(d) = quadratic(d, N_params)
    @assert N_params ≈ coef(N_model) "The coefficients of the quadratic fit of the gate eigenvalue numbers are not integers."
    @assert N_fit(dist_range) ≈ N_scaling "The gate eigenvalue numbers are not well-fit by a quadratic."
    # Fit the trace of the gate eigenvalue estimator covariance matrix
    dep_trace_model =
        lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = dep_trace_scaling, x = dist_range))
    dep_trace_params = coef(dep_trace_model)
    dep_trace_fit(d) = quadratic(d, dep_trace_params)
    if ~(isapprox(dep_trace_fit(dist_range), dep_trace_scaling; rtol = 1e-3))
        @warn "The traces are not well-fit by a quadratic."
    end
    # Fit the trace of the square of the gate eigenvalue estimator covariance matrix
    dep_trace_sq_model =
        lm(@formula(y ~ 1 + x + x^2), DataFrame(; y = dep_trace_sq_scaling, x = dist_range))
    dep_trace_sq_params = coef(dep_trace_sq_model)
    dep_trace_sq_fit(d) = quadratic(d, dep_trace_sq_params)
    if ~(isapprox(dep_trace_sq_fit(dist_range), dep_trace_sq_scaling; rtol = 1e-3))
        @warn "The traces of the square are not well-fit by a quadratic."
    end
    # Create the quadratic-based models for the NRMSE expectation and variance
    @. expectation_model(d, c) =
        sqrt((c[1] + c[2] * d + c[3] * d^2) / N_fit(d)) *
        (1 - ((c[4] + c[5] * d + c[6] * d^2) / (4 * (c[1] + c[2] * d + c[3] * d^2)^2)))
    @. variance_model(d, c) =
        ((c[4] + c[5] * d + c[6] * d^2) / (2 * N_fit(d) * (c[1] + c[2] * d + c[3] * d^2))) *
        (1 - ((c[4] + c[5] * d + c[6] * d^2) / (8 * (c[1] + c[2] * d + c[3] * d^2)^2)))
    # Simultaneously fit the mean of the NRMSE expectation and variance across the instances of log-normal noise
    mean_expectation_scaling = mean.(expectation_scaling)
    mean_variance_scaling = mean.(variance_scaling)
    # Rescale the variance so it is considered appropriately when fitting
    pair_rescale = mean_expectation_scaling[1] / sqrt(mean_variance_scaling[1])
    pair_scaling = [mean_expectation_scaling; pair_rescale * sqrt.(mean_variance_scaling)]
    pair_model(d, c) = [expectation_model(d, c); pair_rescale * sqrt.(variance_model(d, c))]
    param_init = [dep_trace_params; dep_trace_sq_params]
    pair_fit = curve_fit(pair_model, dist_range, pair_scaling, param_init)
    # Fit the NRMSE expectation
    expectation_fit(d) = expectation_model(d, pair_fit.param)
    if ~(isapprox(expectation_fit(dist_range), mean_expectation_scaling; rtol = 1e-3))
        @warn "The mean NRMSE expectations are not well-fit."
    end
    # Fit the NRMSE variance
    variance_fit(d) = variance_model(d, pair_fit.param)
    if ~(isapprox(variance_fit(dist_range), mean_variance_scaling; rtol = 1e-3))
        @warn "The mean NRMSE variances are not well-fit."
    end
    # Save and return the results
    overall_time = time() - start_time
    log_planar_scaling = LognormalPlanarScaling(
        dist_range,
        N_fit,
        N_params,
        expectation_scaling,
        expectation_fit,
        variance_scaling,
        variance_fit,
        eigenvalues_scaling,
        circuit_param,
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
        save_scaling(log_planar_scaling)
    end
    if diagnostics
        println(
            "Finished calculating merit variable scaling with distance for log-normal Pauli noise. The time elapsed since calculations started is $(round(overall_time, digits = 3)) s.",
        )
    end
    return log_planar_scaling::LognormalPlanarScaling
end

function calc_lognormal_planar_scaling(
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
    if typeof(d.c.circuit_param) == RotatedPlanarParameters
        dist_range = collect(3:dist_max)
    elseif typeof(d.c.circuit_param) == UnrotatedPlanarParameters
        dist_range = collect(3:dist_max)
    else
        throw(error("Unsupported circuit type $(typeof(d.c.circuit_param))."))
    end
    # Calculate the scaling data
    log_planar_scaling = calc_lognormal_planar_scaling(
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
    return log_planar_scaling::LognormalPlanarScaling
end

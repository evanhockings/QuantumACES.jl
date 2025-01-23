"""
    MeritScaling

Scaling data for an experimental design for a circuit with vertical and horizontal distance parameters, conventionally under depolarising noise.

# Fields

  - `d::Design`: Design for which the scaling data is calculated.
  - `dist_range::Vector{Int}`: Circuit distances.
  - `merits::Vector{Merit}`: Merits of the design for a range of distances.
  - `calculation_times::Matrix{Float64}`: Time taken to generate the design, and calculate the merit, for each distance.
  - `overall_time::Float64`: The overall time taken to calculate the merit scaling data for depolarising noise.
"""
struct MeritScaling <: AbstractScalingData
    d::Design
    dist_range::Vector{Int}
    merits::Vector{Merit}
    calculation_times::Matrix{Float64}
    overall_time::Float64
end

function Base.show(io::IO, s::MeritScaling)
    print(
        io,
        "Merit scaling data for a design for a $(s.d.c.circuit_param.circuit_name) code with $(length(s.d.tuple_set)) tuples.",
    )
    return nothing
end

@struct_hash_equal_isequal MeritScaling

"""
    EnsembleScaling

Ensemble scaling data for an experimental design for a circuit with vertical and horizontal distance parameters under a random noise model.

# Fields

  - `d::Design`: Design for which the scaling data is calculated.
  - `dist_range::Vector{Int}`: Circuit distances.
  - `merits::Vector{Merit}`: Merits of the design for each distance.
  - `merit_ensemble::Vector{Vector{Merit}}`: Merits of the design across random noise models for each distance.
  - `seeds::Vector{UInt64}`: Seeds for the random noise parameters.
  - `calculation_times::Matrix{Float64}`: The time to generate the design, calculate the depolarising merit, and calculate the merit ensemble, for each distance.
  - `overall_time::Float64`: The overall time taken to calculate the merit ensemble scaling data.
"""
struct EnsembleScaling <: AbstractScalingData
    d::Design
    dist_range::Vector{Int}
    merits::Vector{Merit}
    merit_ensemble::Vector{Vector{Merit}}
    seeds::Vector{UInt64}
    calculation_times::Matrix{Float64}
    overall_time::Float64
end

function Base.show(io::IO, s::EnsembleScaling)
    print(
        io,
        "Merit ensemble scaling data over random noise models for a design for a $(s.d.c.circuit_param.circuit_name) code with $(length(s.d.tuple_set)) tuples.",
    )
    return nothing
end

@struct_hash_equal_isequal EnsembleScaling

"""
    ScalingFit

Scaling fit data for an experimental design for a circuit with vertical and horizontal distance parameters, conventionally under depolarising noise.

All fit parameters describe a quadratic as ``p[1] + p[2] * dist + p[3] * dist^2``.

# Fields

  - `N::Vector{Int}`: Gate eigenvalues quadratic fit parameters.
  - `N_model::Function`: Quadratic model for the gate eigenvalues.
  - `N_marginal::Vector{Int}`: Marginal gate eigenvalues quadratic fit parameters.
  - `N_marginal_model::Function`: Quadratic model for the marginal gate eigenvalues.
  - `N_relative::Vector{Int}`: Relative gate eigenvalues quadratic fit parameters.
  - `N_relative_model::Function`: Quadratic model for the relative gate eigenvalues.
  - `G::Vector{Int}`: Gate number quadratic fit parameters.
  - `G_model::Function`: Quadratic model for the gate number.
  - `gls_tr::Vector{Float64}`: Ordinary GLS trace quadratic fit parameters.
  - `gls_tr_model::Function`: Quadratic model for the ordinary GLS trace.
  - `gls_tr_sq::Vector{Float64}`: Ordinary GLS trace squared quadratic fit parameters.
  - `gls_tr_sq_model::Function`: Quadratic model for the ordinary GLS trace squared.
  - `gls_expectation_model::Function`: Model for the ordinary GLS expectation.
  - `gls_variance_model::Function`: Model for the ordinary GLS variance.
  - `gls_marginal_tr::Vector{Float64}`: Marginal GLS trace quadratic fit parameters.
  - `gls_marginal_tr_model::Function`: Quadratic model for the marginal GLS trace.
  - `gls_marginal_tr_sq::Vector{Float64}`: Marginal GLS trace squared quadratic fit parameters.
  - `gls_marginal_tr_sq_model::Function`: Quadratic model for the marginal GLS trace squared.
  - `gls_marginal_expectation_model::Function`: Model for the marginal GLS expectation.
  - `gls_marginal_variance_model::Function`: Model for the marginal GLS variance.
  - `gls_relative_tr::Vector{Float64}`: Relative GLS trace quadratic fit parameters.
  - `gls_relative_tr_model::Function`: Quadratic model for the relative GLS trace.
  - `gls_relative_tr_sq::Vector{Float64}`: Relative GLS trace squared quadratic fit parameters.
  - `gls_relative_tr_sq_model::Function`: Quadratic model for the relative GLS trace squared.
  - `gls_relative_expectation_model::Function`: Model for the relative GLS expectation.
  - `gls_relative_variance_model::Function`: Model for the relative GLS variance.
  - `wls_tr::Vector{Float64}`: Ordinary WLS trace quadratic fit parameters.
  - `wls_tr_model::Function`: Quadratic model for the ordinary WLS trace.
  - `wls_tr_sq::Vector{Float64}`: Ordinary WLS trace squared quadratic fit parameters.
  - `wls_tr_sq_model::Function`: Quadratic model for the ordinary WLS trace squared.
  - `wls_expectation_model::Function`: Model for the ordinary WLS expectation.
  - `wls_variance_model::Function`: Model for the ordinary WLS variance.
  - `wls_marginal_tr::Vector{Float64}`: Marginal WLS trace quadratic fit parameters.
  - `wls_marginal_tr_model::Function`: Quadratic model for the marginal WLS trace.
  - `wls_marginal_tr_sq::Vector{Float64}`: Marginal WLS trace squared quadratic fit parameters.
  - `wls_marginal_tr_sq_model::Function`: Quadratic model for the marginal WLS trace squared.
  - `wls_marginal_expectation_model::Function`: Model for the marginal WLS expectation.
  - `wls_marginal_variance_model::Function`: Model for the marginal WLS variance.
  - `wls_relative_tr::Vector{Float64}`: Relative WLS trace quadratic fit parameters.
  - `wls_relative_tr_model::Function`: Quadratic model for the relative WLS trace.
  - `wls_relative_tr_sq::Vector{Float64}`: Relative WLS trace squared quadratic fit parameters.
  - `wls_relative_tr_sq_model::Function`: Quadratic model for the relative WLS trace squared.
  - `wls_relative_expectation_model::Function`: Model for the relative WLS expectation.
  - `wls_relative_variance_model::Function`: Model for the relative WLS variance.
  - `ols_tr::Vector{Float64}`: Ordinary OLS trace quadratic fit parameters.
  - `ols_tr_model::Function`: Quadratic model for the ordinary OLS trace.
  - `ols_tr_sq::Vector{Float64}`: Ordinary OLS trace squared quadratic fit parameters.
  - `ols_tr_sq_model::Function`: Quadratic model for the ordinary OLS trace squared.
  - `ols_expectation_model::Function`: Model for the ordinary OLS expectation.
  - `ols_variance_model::Function`: Model for the ordinary OLS variance.
  - `ols_marginal_tr::Vector{Float64}`: Marginal OLS trace quadratic fit parameters.
  - `ols_marginal_tr_model::Function`: Quadratic model for the marginal OLS trace.
  - `ols_marginal_tr_sq::Vector{Float64}`: Marginal OLS trace squared quadratic fit parameters.
  - `ols_marginal_tr_sq_model::Function`: Quadratic model for the marginal OLS trace squared.
  - `ols_marginal_expectation_model::Function`: Model for the marginal OLS expectation.
  - `ols_marginal_variance_model::Function`: Model for the marginal OLS variance.
  - `ols_relative_tr::Vector{Float64}`: Relative OLS trace quadratic fit parameters.
  - `ols_relative_tr_model::Function`: Quadratic model for the relative OLS trace.
  - `ols_relative_tr_sq::Vector{Float64}`: Relative OLS trace squared quadratic fit parameters.
  - `ols_relative_tr_sq_model::Function`: Quadratic model for the relative OLS trace squared.
  - `ols_relative_expectation_model::Function`: Model for the relative OLS expectation.
  - `ols_relative_variance_model::Function`: Model for the relative OLS variance.
"""
struct ScalingFit
    N::Vector{Int}
    N_model::Function
    N_marginal::Vector{Int}
    N_marginal_model::Function
    N_relative::Vector{Int}
    N_relative_model::Function
    G::Vector{Int}
    G_model::Function
    gls_tr::Vector{Float64}
    gls_tr_model::Function
    gls_tr_sq::Vector{Float64}
    gls_tr_sq_model::Function
    gls_expectation_model::Function
    gls_variance_model::Function
    gls_marginal_tr::Vector{Float64}
    gls_marginal_tr_model::Function
    gls_marginal_tr_sq::Vector{Float64}
    gls_marginal_tr_sq_model::Function
    gls_marginal_expectation_model::Function
    gls_marginal_variance_model::Function
    gls_relative_tr::Vector{Float64}
    gls_relative_tr_model::Function
    gls_relative_tr_sq::Vector{Float64}
    gls_relative_tr_sq_model::Function
    gls_relative_expectation_model::Function
    gls_relative_variance_model::Function
    wls_tr::Vector{Float64}
    wls_tr_model::Function
    wls_tr_sq::Vector{Float64}
    wls_tr_sq_model::Function
    wls_expectation_model::Function
    wls_variance_model::Function
    wls_marginal_tr::Vector{Float64}
    wls_marginal_tr_model::Function
    wls_marginal_tr_sq::Vector{Float64}
    wls_marginal_tr_sq_model::Function
    wls_marginal_expectation_model::Function
    wls_marginal_variance_model::Function
    wls_relative_tr::Vector{Float64}
    wls_relative_tr_model::Function
    wls_relative_tr_sq::Vector{Float64}
    wls_relative_tr_sq_model::Function
    wls_relative_expectation_model::Function
    wls_relative_variance_model::Function
    ols_tr::Vector{Float64}
    ols_tr_model::Function
    ols_tr_sq::Vector{Float64}
    ols_tr_sq_model::Function
    ols_expectation_model::Function
    ols_variance_model::Function
    ols_marginal_tr::Vector{Float64}
    ols_marginal_tr_model::Function
    ols_marginal_tr_sq::Vector{Float64}
    ols_marginal_tr_sq_model::Function
    ols_marginal_expectation_model::Function
    ols_marginal_variance_model::Function
    ols_relative_tr::Vector{Float64}
    ols_relative_tr_model::Function
    ols_relative_tr_sq::Vector{Float64}
    ols_relative_tr_sq_model::Function
    ols_relative_expectation_model::Function
    ols_relative_variance_model::Function
end

function Base.show(io::IO, scaling_fit::ScalingFit)
    print(io, "Scaling fit data.")
    return nothing
end

@struct_hash_equal_isequal ScalingFit

"""
    EnsembleFit

Ensemble scaling fit data for an experimental design for a circuit with vertical and horizontal distance parameters under a random noise model.

All pair fit parameters correspond to the trace and trace squared quadratic terms appearing in the expressions for the expectation and variance.

# Fields

  - `gls_pair::Vector{Float64}`: Ordinary GLS pair quadratic fit parameters.
  - `gls_expectation_model::Function`: Model for the ordinary GLS expectation.
  - `gls_variance_model::Function`: Model for the ordinary GLS variance.
  - `gls_marginal_pair::Vector{Float64}`: Marginal GLS pair quadratic fit parameters.
  - `gls_marginal_expectation_model::Function`: Model for the marginal GLS expectation.
  - `gls_marginal_variance_model::Function`: Model for the marginal GLS variance.
  - `gls_relative_pair::Vector{Float64}`: Relative GLS pair quadratic fit parameters.
  - `gls_relative_expectation_model::Function`: Model for the relative GLS expectation.
  - `gls_relative_variance_model::Function`: Model for the relative GLS variance.
  - `wls_pair::Vector{Float64}`: Ordinary WLS pair quadratic fit parameters.
  - `wls_expectation_model::Function`: Model for the ordinary WLS expectation.
  - `wls_variance_model::Function`: Model for the ordinary WLS variance.
  - `wls_marginal_pair::Vector{Float64}`: Marginal WLS pair quadratic fit parameters.
  - `wls_marginal_expectation_model::Function`: Model for the marginal WLS expectation.
  - `wls_marginal_variance_model::Function`: Model for the marginal WLS variance.
  - `wls_relative_pair::Vector{Float64}`: Relative WLS pair quadratic fit parameters.
  - `wls_relative_expectation_model::Function`: Model for the relative WLS expectation.
  - `wls_relative_variance_model::Function`: Model for the relative WLS variance.
  - `ols_pair::Vector{Float64}`: Ordinary OLS pair quadratic fit parameters.
  - `ols_expectation_model::Function`: Model for the ordinary OLS expectation.
  - `ols_variance_model::Function`: Model for the ordinary OLS variance.
  - `ols_marginal_pair::Vector{Float64}`: Marginal OLS pair quadratic fit parameters.
  - `ols_marginal_expectation_model::Function`: Model for the marginal OLS expectation.
  - `ols_marginal_variance_model::Function`: Model for the marginal OLS variance.
  - `ols_relative_pair::Vector{Float64}`: Relative OLS pair quadratic fit parameters.
  - `ols_relative_expectation_model::Function`: Model for the relative OLS expectation.
  - `ols_relative_variance_model::Function`: Model for the relative OLS variance.
"""
struct EnsembleFit
    gls_pair::Vector{Float64}
    gls_expectation_model::Function
    gls_variance_model::Function
    gls_marginal_pair::Vector{Float64}
    gls_marginal_expectation_model::Function
    gls_marginal_variance_model::Function
    gls_relative_pair::Vector{Float64}
    gls_relative_expectation_model::Function
    gls_relative_variance_model::Function
    wls_pair::Vector{Float64}
    wls_expectation_model::Function
    wls_variance_model::Function
    wls_marginal_pair::Vector{Float64}
    wls_marginal_expectation_model::Function
    wls_marginal_variance_model::Function
    wls_relative_pair::Vector{Float64}
    wls_relative_expectation_model::Function
    wls_relative_variance_model::Function
    ols_pair::Vector{Float64}
    ols_expectation_model::Function
    ols_variance_model::Function
    ols_marginal_pair::Vector{Float64}
    ols_marginal_expectation_model::Function
    ols_marginal_variance_model::Function
    ols_relative_pair::Vector{Float64}
    ols_relative_expectation_model::Function
    ols_relative_variance_model::Function
end

function Base.show(io::IO, ensemble_fit::EnsembleFit)
    print(io, "Ensemble scaling fit data.")
    return nothing
end

@struct_hash_equal_isequal EnsembleFit

"""
    calc_merit_scaling(d::Design, dist_max::Integer; kwargs...)
    calc_merit_scaling(d::Design, dist_range::Vector{Int}; kwargs...)

Returns the scaling data as a [`MeritScaling`](@ref) object for the merit of the design `d`, for a circuit with `vertical_dist` and `horizontal_dist` parameters, as a function of these distances.

# Arguments

  - `d::Design`: Design for which the merit scaling is calculated.
  - `dist_max::Int`: Maximum distance for which the merit scaling is calculated.
  - `dist_range::Vector{Int}`: Distances for which the merit scaling is calculated.

# Keyword arguments

  - `warning::Bool = true`: Whether to print a warning if the noise model is not depolarising.
  - `diagnostics::Bool = true`: Whether to print diagnostic information.
  - `save_data::Bool = false`: Whether to save the merit scaling data.
"""
function calc_merit_scaling(
    d::Design,
    dist_range::Vector{Int};
    warning::Bool = true,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    # Initialise parameters
    circuit_param = d.c.circuit_param
    circuit_params = circuit_param.params
    @assert haskey(circuit_params, :vertical_dist) "The circuit parameters lack a `vertical_dist` field."
    @assert haskey(circuit_params, :horizontal_dist) "The circuit parameters lack  a `horizontal_dist` field."
    noise_param = d.c.noise_param
    if warning && typeof(noise_param) != DepolarisingParameters
        println(
            "WARNING: This function works best with deterministic noise models such as depolarising noise. Disable this warning by setting `warning` to be `false`.",
        )
    end
    @assert minimum(dist_range) >= 3 "The supplied distances must all be at least 3."
    if length(dist_range) < 3
        @warn "Fewer than 3 distances were supplied, which is insufficient for quadratically fitting the data."
    end
    tuple_set_data = d.tuple_set_data
    shot_weights = d.shot_weights
    # Calculate the merit scaling with distance
    start_time = time()
    merits = Vector{Merit}(undef, length(dist_range))
    calculation_times = Matrix{Float64}(undef, length(dist_range), 2)
    for (idx, dist) in pairs(dist_range)
        # Initialise the circuit parameters
        circuit_param_dist = deepcopy(circuit_param)
        circuit_param_dist.params[:vertical_dist] = dist
        circuit_param_dist.params[:horizontal_dist] = dist
        # Generate the design
        time_1 = time()
        c = get_circuit(circuit_param_dist, noise_param)
        d_dist = generate_design(
            c,
            tuple_set_data;
            shot_weights = shot_weights,
            full_covariance = true,
        )
        time_2 = time()
        design_time = time_2 - time_1
        if diagnostics
            println(
                "Generating the design at distance $(dist) took $(round(design_time, digits = 3)) s.",
            )
        end
        # Calculate the merit
        merit_dist = calc_merit(d_dist)
        merits[idx] = merit_dist
        time_3 = time()
        merit_time = time_3 - time_2
        calculation_times[idx, :] = [design_time; merit_time]
        if diagnostics
            println(
                "Calculating the merit at distance $(dist) took $(round(merit_time, digits = 3)) s.",
            )
        end
        GC.gc()
    end
    # Save and return the results
    overall_time = time() - start_time
    merit_scaling = MeritScaling(d, dist_range, merits, calculation_times, overall_time)
    if save_data
        save_scaling(merit_scaling)
    end
    if diagnostics
        println(
            "Finished calculating merit scaling with distance. The time elapsed since calculations started is $(round(overall_time, digits = 3)) s.",
        )
    end
    return merit_scaling::MeritScaling
end
function calc_merit_scaling(
    d::Design,
    dist_max::Integer;
    warning::Bool = true,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    dist_range = collect(3:dist_max)
    merit_scaling = calc_merit_scaling(
        d,
        dist_range;
        warning = warning,
        save_data = save_data,
        diagnostics = diagnostics,
    )
    return merit_scaling::MeritScaling
end

"""
    calc_ensemble_scaling(d::Design, dist_max::Integer; kwargs...)
    calc_ensemble_scaling(d::Design, dist_range::Vector{Int}; kwargs...)

Returns ensemble scaling data as a [`EnsembleScaling`](@ref) object for the merit of the design `d` over random instances of noise models, for a circuit with `vertical_dist` and `horizontal_dist` parameters, as a function of these distances.

# Arguments

  - `d::Design`: Design for which the merit scaling is calculated.
  - `dist_max::Int`: Maximum code distance for which the merit scaling is calculated.
  - `dist_range::Vector{Int}`: Code distances for which the merit scaling is calculated.

# Keyword arguments

  - `precision::Real = 2e-3`: Precision to which the figure of merit is estimated, corresponding to the target standard error of the mean.
  - `max_repetitions::Integer = 10000`: Maximum number of random instances of log-normal Pauli noise over which the figure of merit is calculated.
  - `min_repetitions::Integer = 50`: Minimum number of random instances of log-normal Pauli noise over which the figure of merit is calculated.
  - `print_repetitions::Integer = 50`: Number of random instances of log-normal Pauli noise between printing diagnostics.
  - `seed::Union{UInt64, Nothing} = nothing`: Seeds used to generate instances of log-normal Pauli noise.
  - `diagnostics::Bool = true`: Whether to print diagnostic information.
  - `save_data::Bool = false`: Whether to save the merit scaling data.
"""
function calc_ensemble_scaling(
    d::Design,
    dist_range::Vector{Int};
    precision::Real = 2e-3,
    max_repetitions::Integer = 10000,
    min_repetitions::Integer = 50,
    print_repetitions::Integer = 50,
    seed::Union{UInt64, Nothing} = nothing,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    # Initialise parameters
    @assert precision > 0 "The precision must be positive."
    @assert max_repetitions >= 2 "The maximum number of repetitions must be at least 2."
    @assert min_repetitions >= 2 "The minimum number of repetitions must be at least 2."
    @assert max_repetitions >= min_repetitions "The maximum number of repetitions must be greater than or equal to the minimum number of repetitions."
    @assert print_repetitions > 0 "The number of repetitions between printing must be positive."
    circuit_param = d.c.circuit_param
    circuit_params = circuit_param.params
    @assert haskey(circuit_params, :vertical_dist) "The circuit parameters lack a `vertical_dist` field."
    @assert haskey(circuit_params, :horizontal_dist) "The circuit parameters lack  a `horizontal_dist` field."
    noise_param = d.c.noise_param
    noise_params = noise_param.params
    @assert haskey(noise_params, :seed) "The noise parameters lack a `seed` field."
    @assert haskey(noise_params, :r_1) "The noise parameters lack a `r_1` field."
    @assert haskey(noise_params, :r_2) "The noise parameters lack a `r_2` field."
    @assert haskey(noise_params, :r_m) "The noise parameters lack a `r_m` field."
    @assert haskey(noise_params, :r_im) "The noise parameters lack a `r_im` field."
    @assert haskey(noise_params, :r_r) "The noise parameters lack a `r_r` field."
    @assert haskey(noise_params, :combined) "The noise parameters lack a `combined` field."
    @assert minimum(dist_range) >= 3 "The supplied distances must all be at least 3."
    if length(dist_range) < 3
        @warn "Fewer than 3 distances were supplied, which is insufficient for quadratically fitting the data."
    end
    tuple_set_data = d.tuple_set_data
    shot_weights = d.shot_weights
    dep_noise_param = get_dep_param(;
        r_1 = noise_params[:r_1],
        r_2 = noise_params[:r_2],
        r_m = noise_params[:r_m],
        r_im = noise_params[:r_im],
        r_r = noise_params[:r_r],
        combined = noise_params[:combined],
    )
    # Set the random seeds
    if seed !== nothing
        Random.seed!(seed)
    end
    seeds = rand(UInt64, max_repetitions)
    if seed !== nothing
        Random.seed!()
    end
    # Calculate the merit ensemble scaling with distance
    start_time = time()
    merits = Vector{Merit}(undef, length(dist_range))
    merit_ensemble = Vector{Vector{Merit}}(undef, length(dist_range))
    calculation_times = Matrix{Float64}(undef, length(dist_range), 3)
    for (idx, dist) in pairs(dist_range)
        # Initialise the circuit parameters
        circuit_param_dist = deepcopy(circuit_param)
        circuit_param_dist.params[:vertical_dist] = dist
        circuit_param_dist.params[:horizontal_dist] = dist
        # Generate the design
        time_1 = time()
        c = get_circuit(circuit_param_dist, dep_noise_param)
        d_dist = generate_design(
            c,
            tuple_set_data;
            shot_weights = shot_weights,
            full_covariance = true,
        )
        time_2 = time()
        design_time = time_2 - time_1
        if diagnostics
            println(
                "Generating the design at distance $(dist) took $(round(design_time, digits = 3)) s.",
            )
        end
        # Calculate the merit
        merit_dist = calc_merit(d_dist)
        merits[idx] = merit_dist
        time_3 = time()
        merit_time = time_3 - time_2
        if diagnostics
            println(
                "Calculating the depolarising merit at distance $(dist) took $(round(merit_time, digits = 3)) s.",
            )
        end
        # Calculate the merit ensemble
        rep = 1
        generating = true
        merit_ensemble[idx] = Vector{Merit}()
        while generating
            # Generate the noise and design
            noise_param_rep = deepcopy(noise_param)
            noise_param_rep.params[:seed] = seeds[rep]
            d_rep = update_noise(d_dist, noise_param_rep)
            # Calculate the merit
            merit_rep = calc_merit(d_rep)
            push!(merit_ensemble[idx], merit_rep)
            # Check if the scaling has converged, only examining GLS and WLS
            gls_expectations = [merit.gls_expectation for merit in merit_ensemble[idx]]
            wls_expectations = [merit.wls_expectation for merit in merit_ensemble[idx]]
            gls_rel_std = std(gls_expectations) / (mean(gls_expectations) * sqrt(rep))
            wls_rel_std = std(wls_expectations) / (mean(wls_expectations) * sqrt(rep))
            if (
                rep >= min_repetitions &&
                gls_rel_std < precision &&
                wls_rel_std < precision
            ) || rep >= max_repetitions
                generating = false
            else
                # Display a diagnostic if appropriate
                if diagnostics && rep % print_repetitions == 0
                    println(
                        "Calculated $(rep) merits at distance $(dist). The GLS and WLS relative standard errors of the mean are $(round(gls_rel_std, sigdigits = 3)) and $(round(wls_rel_std, sigdigits = 3)), respectively, compared to the desired precision $(round(precision, sigdigits = 3)). The time elapsed for this ensemble is $(round(time() - time_3, digits = 3)) s.",
                    )
                end
                rep += 1
            end
        end
        time_4 = time()
        rep_time = time_4 - time_3
        calculation_times[idx, :] = [design_time; merit_time; rep_time]
        if diagnostics
            println(
                "Calculating the ensemble of $(rep) merits at distance $(dist) took $(round(rep_time, digits = 3)) s.",
            )
        end
        GC.gc()
    end
    # Save and return the results
    overall_time = time() - start_time
    ensemble_scaling = EnsembleScaling(
        d,
        dist_range,
        merits,
        merit_ensemble,
        seeds,
        calculation_times,
        overall_time,
    )
    if save_data
        save_scaling(ensemble_scaling)
    end
    if diagnostics
        println(
            "Finished calculating the merit ensemble scaling with distance. The time elapsed since calculations started is $(round(overall_time, digits = 3)) s.",
        )
    end
    return ensemble_scaling::EnsembleScaling
end
function calc_ensemble_scaling(
    d::Design,
    dist_max::Integer;
    precision::Real = 2e-3,
    max_repetitions::Integer = 10000,
    min_repetitions::Integer = 50,
    print_repetitions::Integer = 50,
    seed::Union{UInt64, Nothing} = nothing,
    diagnostics::Bool = true,
    save_data::Bool = false,
)
    dist_range = collect(3:dist_max)
    ensemble_scaling = calc_ensemble_scaling(
        d,
        dist_range;
        precision = precision,
        max_repetitions = max_repetitions,
        min_repetitions = min_repetitions,
        print_repetitions = print_repetitions,
        seed = seed,
        save_data = save_data,
        diagnostics = diagnostics,
    )
    return ensemble_scaling::EnsembleScaling
end

"""
    quadratic_model(dist, params)

Quadratic model as a function of the distance `dist` with parameters `params`.
"""
function quadratic_model(dist, params)
    return params[1] .+ params[2] * dist .+ params[3] * dist .^ 2
end

"""
    expectation_model(dist, tr_pair_params, N_params)

Model of the expectation as a function of the distance `dist` with trace pair parameters `tr_pair_params`, consisting of trace and trace squared quadratic parameters, and gate eigenvalue quadratic parameters `N_params`.
"""
function expectation_model(dist, tr_pair_params, N_params)
    return sqrt.(
        quadratic_model(dist, tr_pair_params[1:3]) ./ quadratic_model(dist, N_params)
    ) .* (
        1 .- (
            quadratic_model(dist, tr_pair_params[4:6]) ./
            (4 * quadratic_model(dist, tr_pair_params[1:3]) .^ 2)
        )
    )
end

"""
    variance_model(dist, tr_pair_params, N_params)

Model of the variance as a function of the distance `dist` with trace pair parameters `tr_pair_params`, consisting of trace and trace squared quadratic parameters, and gate eigenvalue quadratic parameters `N_params`.
"""
function variance_model(dist, tr_pair_params, N_params)
    return (
        quadratic_model(dist, tr_pair_params[4:6]) ./
        (2 * quadratic_model(dist, N_params) .* quadratic_model(dist, tr_pair_params[1:3]))
    ) .* (
        1 .- (
            quadratic_model(dist, tr_pair_params[4:6]) ./
            (8 * quadratic_model(dist, tr_pair_params[1:3]) .^ 2)
        )
    )
end

"""
    check_model(dist_range, model_range, model; precision::Real = 1e-2)

Displays a warning if the model predictions do not match the data to relative precision `precision`.
"""
function check_model(dist_range, model_range, model; precision::Real = 1e-2)
    if ~all(isapprox.(model_range, model(dist_range); atol = 0.0, rtol = precision))
        println(
            "WARNING: model predictions do not match the data to relative precision $(precision); respectively:\n$(model(dist_range))\n$(model_range).",
        )
    end
    return nothing
end

"""
    fit_quadratic_model(dist_range::Vector{Int}, quadratic_range::Vector{Float64}; precision::Real = 1e-2)
    fit_quadratic_model(dist_range::Vector{Int}, quadratic_range::Vector{Int})

Returns the quadratic fit parameters for the data `quadratic_range` as a function of the distance `dist_range`, requiring exact fit for integer data, and displaying a warning if non-integer data is not fit to relative precision `precision`.
"""
function fit_quadratic_model(
    dist_range::Vector{Int},
    quadratic_range::Vector{Float64};
    precision::Real = 1e-2,
)
    quadratic_fit =
        LsqFit.curve_fit(quadratic_model, dist_range, quadratic_range, zeros(Float64, 3))
    quadratic_params = quadratic_fit.param
    check_quadratic_model(dist) = quadratic_model(dist, quadratic_params)
    check_model(dist_range, quadratic_range, check_quadratic_model; precision = precision)
    return quadratic_params::Vector{Float64}
end
function fit_quadratic_model(dist_range::Vector{Int}, quadratic_range::Vector{Int})
    quadratic_fit =
        LsqFit.curve_fit(quadratic_model, dist_range, quadratic_range, zeros(Float64, 3))
    quadratic_params = round.(Int, quadratic_fit.param)
    @assert quadratic_params ≈ quadratic_fit.param
    @assert quadratic_model(dist_range, quadratic_params) == quadratic_range
    return quadratic_params::Vector{Int}
end

"""
    get_scaling_fit(merit_scaling::MeritScaling; precision::Real = 1e-2)
    get_scaling_fit(merits::Vector{Merit}, dist_range::Vector{Int}; precision::Real = 1e-2)

Returns merit scaling fit data as a [`ScalingFit`](@ref) object for the merit scaling data `merits` fit against the distances `dist_range`, contained in the scaling data `merit_scaling`, displaying a warning if data is not fit to relative precision `precision`.
"""
function get_scaling_fit(
    merits::Vector{Merit},
    dist_range::Vector{Int};
    precision::Real = 1e-2,
)
    # Fit the integer quantities
    N_range = [merit.N for merit in merits]
    N_params = fit_quadratic_model(dist_range, N_range)
    N_model(dist) = quadratic_model(dist, N_params)
    N_marginal_range = [merit.N_marginal for merit in merits]
    N_marginal_params = fit_quadratic_model(dist_range, N_marginal_range)
    N_marginal_model(dist) = quadratic_model(dist, N_marginal_params)
    N_relative_range = [merit.N_relative for merit in merits]
    N_relative_params = fit_quadratic_model(dist_range, N_relative_range)
    N_relative_model(dist) = quadratic_model(dist, N_relative_params)
    G_range = [merit.G for merit in merits]
    G_params = fit_quadratic_model(dist_range, G_range)
    G_model(dist) = quadratic_model(dist, G_params)
    # Fit the ordinary GLS merit
    gls_tr_range = [sum(merit.gls_cov_eigenvalues) for merit in merits]
    gls_tr_params = fit_quadratic_model(dist_range, gls_tr_range; precision = precision)
    gls_tr_model(dist) = quadratic_model(dist, gls_tr_params)
    gls_tr_sq_range = [sum(merit.gls_cov_eigenvalues .^ 2) for merit in merits]
    gls_tr_sq_params =
        fit_quadratic_model(dist_range, gls_tr_sq_range; precision = precision)
    gls_tr_sq_model(dist) = quadratic_model(dist, gls_tr_sq_params)
    gls_expectation_range = [merit.gls_expectation for merit in merits]
    gls_expectation_model(dist) =
        expectation_model(dist, [gls_tr_params; gls_tr_sq_params], N_params)
    check_model(
        dist_range,
        gls_expectation_range,
        gls_expectation_model;
        precision = precision,
    )
    gls_variance_range = [merit.gls_variance for merit in merits]
    gls_variance_model(dist) =
        variance_model(dist, [gls_tr_params; gls_tr_sq_params], N_params)
    check_model(dist_range, gls_variance_range, gls_variance_model; precision = precision)
    # Fit the marginal GLS merit
    gls_marginal_tr_range = [sum(merit.gls_marginal_cov_eigenvalues) for merit in merits]
    gls_marginal_tr_params =
        fit_quadratic_model(dist_range, gls_marginal_tr_range; precision = precision)
    gls_marginal_tr_model(dist) = quadratic_model(dist, gls_marginal_tr_params)
    gls_marginal_tr_sq_range =
        [sum(merit.gls_marginal_cov_eigenvalues .^ 2) for merit in merits]
    gls_marginal_tr_sq_params =
        fit_quadratic_model(dist_range, gls_marginal_tr_sq_range; precision = precision)
    gls_marginal_tr_sq_model(dist) = quadratic_model(dist, gls_marginal_tr_sq_params)
    gls_marginal_expectation_range = [merit.gls_marginal_expectation for merit in merits]
    gls_marginal_expectation_model(dist) = expectation_model(
        dist,
        [gls_marginal_tr_params; gls_marginal_tr_sq_params],
        N_marginal_params,
    )
    check_model(
        dist_range,
        gls_marginal_expectation_range,
        gls_marginal_expectation_model;
        precision = precision,
    )
    gls_marginal_variance_range = [merit.gls_marginal_variance for merit in merits]
    gls_marginal_variance_model(dist) = variance_model(
        dist,
        [gls_marginal_tr_params; gls_marginal_tr_sq_params],
        N_marginal_params,
    )
    check_model(
        dist_range,
        gls_marginal_variance_range,
        gls_marginal_variance_model;
        precision = precision,
    )
    # Fit the relative GLS merit
    gls_relative_tr_range = [sum(merit.gls_relative_cov_eigenvalues) for merit in merits]
    gls_relative_tr_params =
        fit_quadratic_model(dist_range, gls_relative_tr_range; precision = precision)
    gls_relative_tr_model(dist) = quadratic_model(dist, gls_relative_tr_params)
    gls_relative_tr_sq_range =
        [sum(merit.gls_relative_cov_eigenvalues .^ 2) for merit in merits]
    gls_relative_tr_sq_params =
        fit_quadratic_model(dist_range, gls_relative_tr_sq_range; precision = precision)
    gls_relative_tr_sq_model(dist) = quadratic_model(dist, gls_relative_tr_sq_params)
    gls_relative_expectation_range = [merit.gls_relative_expectation for merit in merits]
    gls_relative_expectation_model(dist) = expectation_model(
        dist,
        [gls_relative_tr_params; gls_relative_tr_sq_params],
        N_relative_params,
    )
    check_model(
        dist_range,
        gls_relative_expectation_range,
        gls_relative_expectation_model;
        precision = precision,
    )
    gls_relative_variance_range = [merit.gls_relative_variance for merit in merits]
    gls_relative_variance_model(dist) = variance_model(
        dist,
        [gls_relative_tr_params; gls_relative_tr_sq_params],
        N_relative_params,
    )
    check_model(
        dist_range,
        gls_relative_variance_range,
        gls_relative_variance_model;
        precision = precision,
    )
    # Fit the ordinary WLS merit
    wls_tr_range = [sum(merit.wls_cov_eigenvalues) for merit in merits]
    wls_tr_params = fit_quadratic_model(dist_range, wls_tr_range; precision = precision)
    wls_tr_model(dist) = quadratic_model(dist, wls_tr_params)
    wls_tr_sq_range = [sum(merit.wls_cov_eigenvalues .^ 2) for merit in merits]
    wls_tr_sq_params =
        fit_quadratic_model(dist_range, wls_tr_sq_range; precision = precision)
    wls_tr_sq_model(dist) = quadratic_model(dist, wls_tr_sq_params)
    wls_expectation_range = [merit.wls_expectation for merit in merits]
    wls_expectation_model(dist) =
        expectation_model(dist, [wls_tr_params; wls_tr_sq_params], N_params)
    check_model(
        dist_range,
        wls_expectation_range,
        wls_expectation_model;
        precision = precision,
    )
    wls_variance_range = [merit.wls_variance for merit in merits]
    wls_variance_model(dist) =
        variance_model(dist, [wls_tr_params; wls_tr_sq_params], N_params)
    check_model(dist_range, wls_variance_range, wls_variance_model; precision = precision)
    # Fit the marginal WLS merit
    wls_marginal_tr_range = [sum(merit.wls_marginal_cov_eigenvalues) for merit in merits]
    wls_marginal_tr_params =
        fit_quadratic_model(dist_range, wls_marginal_tr_range; precision = precision)
    wls_marginal_tr_model(dist) = quadratic_model(dist, wls_marginal_tr_params)
    wls_marginal_tr_sq_range =
        [sum(merit.wls_marginal_cov_eigenvalues .^ 2) for merit in merits]
    wls_marginal_tr_sq_params =
        fit_quadratic_model(dist_range, wls_marginal_tr_sq_range; precision = precision)
    wls_marginal_tr_sq_model(dist) = quadratic_model(dist, wls_marginal_tr_sq_params)
    wls_marginal_expectation_range = [merit.wls_marginal_expectation for merit in merits]
    wls_marginal_expectation_model(dist) = expectation_model(
        dist,
        [wls_marginal_tr_params; wls_marginal_tr_sq_params],
        N_marginal_params,
    )
    check_model(
        dist_range,
        wls_marginal_expectation_range,
        wls_marginal_expectation_model;
        precision = precision,
    )
    wls_marginal_variance_range = [merit.wls_marginal_variance for merit in merits]
    wls_marginal_variance_model(dist) = variance_model(
        dist,
        [wls_marginal_tr_params; wls_marginal_tr_sq_params],
        N_marginal_params,
    )
    check_model(
        dist_range,
        wls_marginal_variance_range,
        wls_marginal_variance_model;
        precision = precision,
    )
    # Fit the relative WLS merit
    wls_relative_tr_range = [sum(merit.wls_relative_cov_eigenvalues) for merit in merits]
    wls_relative_tr_params =
        fit_quadratic_model(dist_range, wls_relative_tr_range; precision = precision)
    wls_relative_tr_model(dist) = quadratic_model(dist, wls_relative_tr_params)
    wls_relative_tr_sq_range =
        [sum(merit.wls_relative_cov_eigenvalues .^ 2) for merit in merits]
    wls_relative_tr_sq_params =
        fit_quadratic_model(dist_range, wls_relative_tr_sq_range; precision = precision)
    wls_relative_tr_sq_model(dist) = quadratic_model(dist, wls_relative_tr_sq_params)
    wls_relative_expectation_range = [merit.wls_relative_expectation for merit in merits]
    wls_relative_expectation_model(dist) = expectation_model(
        dist,
        [wls_relative_tr_params; wls_relative_tr_sq_params],
        N_relative_params,
    )
    check_model(
        dist_range,
        wls_relative_expectation_range,
        wls_relative_expectation_model;
        precision = precision,
    )
    wls_relative_variance_range = [merit.wls_relative_variance for merit in merits]
    wls_relative_variance_model(dist) = variance_model(
        dist,
        [wls_relative_tr_params; wls_relative_tr_sq_params],
        N_relative_params,
    )
    check_model(
        dist_range,
        wls_relative_variance_range,
        wls_relative_variance_model;
        precision = precision,
    )
    # Fit the ordinary OLS merit
    ols_tr_range = [sum(merit.ols_cov_eigenvalues) for merit in merits]
    ols_tr_params = fit_quadratic_model(dist_range, ols_tr_range; precision = precision)
    ols_tr_model(dist) = quadratic_model(dist, ols_tr_params)
    ols_tr_sq_range = [sum(merit.ols_cov_eigenvalues .^ 2) for merit in merits]
    ols_tr_sq_params =
        fit_quadratic_model(dist_range, ols_tr_sq_range; precision = precision)
    ols_tr_sq_model(dist) = quadratic_model(dist, ols_tr_sq_params)
    ols_expectation_range = [merit.ols_expectation for merit in merits]
    ols_expectation_model(dist) =
        expectation_model(dist, [ols_tr_params; ols_tr_sq_params], N_params)
    check_model(
        dist_range,
        ols_expectation_range,
        ols_expectation_model;
        precision = precision,
    )
    ols_variance_range = [merit.ols_variance for merit in merits]
    ols_variance_model(dist) =
        variance_model(dist, [ols_tr_params; ols_tr_sq_params], N_params)
    check_model(dist_range, ols_variance_range, ols_variance_model; precision = precision)
    # Fit the marginal OLS merit
    ols_marginal_tr_range = [sum(merit.ols_marginal_cov_eigenvalues) for merit in merits]
    ols_marginal_tr_params =
        fit_quadratic_model(dist_range, ols_marginal_tr_range; precision = precision)
    ols_marginal_tr_model(dist) = quadratic_model(dist, ols_marginal_tr_params)
    ols_marginal_tr_sq_range =
        [sum(merit.ols_marginal_cov_eigenvalues .^ 2) for merit in merits]
    ols_marginal_tr_sq_params =
        fit_quadratic_model(dist_range, ols_marginal_tr_sq_range; precision = precision)
    ols_marginal_tr_sq_model(dist) = quadratic_model(dist, ols_marginal_tr_sq_params)
    ols_marginal_expectation_range = [merit.ols_marginal_expectation for merit in merits]
    ols_marginal_expectation_model(dist) = expectation_model(
        dist,
        [ols_marginal_tr_params; ols_marginal_tr_sq_params],
        N_marginal_params,
    )
    check_model(
        dist_range,
        ols_marginal_expectation_range,
        ols_marginal_expectation_model;
        precision = precision,
    )
    ols_marginal_variance_range = [merit.ols_marginal_variance for merit in merits]
    ols_marginal_variance_model(dist) = variance_model(
        dist,
        [ols_marginal_tr_params; ols_marginal_tr_sq_params],
        N_marginal_params,
    )
    check_model(
        dist_range,
        ols_marginal_variance_range,
        ols_marginal_variance_model;
        precision = precision,
    )
    # Fit the relative OLS merit
    ols_relative_tr_range = [sum(merit.ols_relative_cov_eigenvalues) for merit in merits]
    ols_relative_tr_params =
        fit_quadratic_model(dist_range, ols_relative_tr_range; precision = precision)
    ols_relative_tr_model(dist) = quadratic_model(dist, ols_relative_tr_params)
    ols_relative_tr_sq_range =
        [sum(merit.ols_relative_cov_eigenvalues .^ 2) for merit in merits]
    ols_relative_tr_sq_params =
        fit_quadratic_model(dist_range, ols_relative_tr_sq_range; precision = precision)
    ols_relative_tr_sq_model(dist) = quadratic_model(dist, ols_relative_tr_sq_params)
    ols_relative_expectation_range = [merit.ols_relative_expectation for merit in merits]
    ols_relative_expectation_model(dist) = expectation_model(
        dist,
        [ols_relative_tr_params; ols_relative_tr_sq_params],
        N_relative_params,
    )
    check_model(
        dist_range,
        ols_relative_expectation_range,
        ols_relative_expectation_model;
        precision = precision,
    )
    ols_relative_variance_range = [merit.ols_relative_variance for merit in merits]
    ols_relative_variance_model(dist) = variance_model(
        dist,
        [ols_relative_tr_params; ols_relative_tr_sq_params],
        N_relative_params,
    )
    check_model(
        dist_range,
        ols_relative_variance_range,
        ols_relative_variance_model;
        precision = precision,
    )
    # Return the results
    scaling_fit = ScalingFit(
        N_params,
        N_model,
        N_marginal_params,
        N_marginal_model,
        N_relative_params,
        N_relative_model,
        G_params,
        G_model,
        gls_tr_params,
        gls_tr_model,
        gls_tr_sq_params,
        gls_tr_sq_model,
        gls_expectation_model,
        gls_variance_model,
        gls_marginal_tr_params,
        gls_marginal_tr_model,
        gls_marginal_tr_sq_params,
        gls_marginal_tr_sq_model,
        gls_marginal_expectation_model,
        gls_marginal_variance_model,
        gls_relative_tr_params,
        gls_relative_tr_model,
        gls_relative_tr_sq_params,
        gls_relative_tr_sq_model,
        gls_relative_expectation_model,
        gls_relative_variance_model,
        wls_tr_params,
        wls_tr_model,
        wls_tr_sq_params,
        wls_tr_sq_model,
        wls_expectation_model,
        wls_variance_model,
        wls_marginal_tr_params,
        wls_marginal_tr_model,
        wls_marginal_tr_sq_params,
        wls_marginal_tr_sq_model,
        wls_marginal_expectation_model,
        wls_marginal_variance_model,
        wls_relative_tr_params,
        wls_relative_tr_model,
        wls_relative_tr_sq_params,
        wls_relative_tr_sq_model,
        wls_relative_expectation_model,
        wls_relative_variance_model,
        ols_tr_params,
        ols_tr_model,
        ols_tr_sq_params,
        ols_tr_sq_model,
        ols_expectation_model,
        ols_variance_model,
        ols_marginal_tr_params,
        ols_marginal_tr_model,
        ols_marginal_tr_sq_params,
        ols_marginal_tr_sq_model,
        ols_marginal_expectation_model,
        ols_marginal_variance_model,
        ols_relative_tr_params,
        ols_relative_tr_model,
        ols_relative_tr_sq_params,
        ols_relative_tr_sq_model,
        ols_relative_expectation_model,
        ols_relative_variance_model,
    )
    return scaling_fit::ScalingFit
end
function get_scaling_fit(merit_scaling::MeritScaling; precision::Real = 1e-2)
    scaling_fit = get_scaling_fit(
        merit_scaling.merits,
        merit_scaling.dist_range;
        precision = precision,
    )
    return scaling_fit::ScalingFit
end

"""
    fit_pair_model(expectation_model::Function, variance_model::Function, dist_range::Vector{Int}, expectation_range::Vector{Float64}, variance_range::Vector{Float64}, param_init::Vector{Float64}; precision::Real = 1e-2)

Returns the pair fit parameters for the expectation and variance models `expectation_model` and `variance_model` as functions of the distance `dist_range`, for the expectation and variance data `expectation_range` and `variance_range`, respectively, and initial parameters `param_init`, displaying a warning if the model predictions do not match the data to relative precision `precision`.
"""
function fit_pair_model(
    expectation_model::Function,
    variance_model::Function,
    dist_range::Vector{Int},
    expectation_range::Vector{Float64},
    variance_range::Vector{Float64},
    param_init::Vector{Float64};
    precision::Real = 1e-2,
)
    try
        pair_scale = expectation_range[1] / sqrt(variance_range[1])
        pair_range = [expectation_range; pair_scale * sqrt.(variance_range)]
        pair_model(dist, params) = [
            expectation_model(dist, params)
            pair_scale * sqrt.(variance_model(dist, params))
        ]
        pair_fit = LsqFit.curve_fit(pair_model, dist_range, pair_range, param_init)
        pair_params = pair_fit.param
        check_expectation_model(dist) = expectation_model(dist, pair_params)
        check_model(
            dist_range,
            expectation_range,
            check_expectation_model;
            precision = precision,
        )
        check_variance_model(dist) = variance_model(dist, pair_params)
        check_model(dist_range, variance_range, check_variance_model; precision = precision)
        return pair_params::Vector{Float64}
    catch
        println(
            "WARNING: Ensemble fit failed. Falling back to depolarising initialisation parameters.",
        )
        return param_init::Vector{Float64}
    end
end

"""
    get_ensemble_fit(ensemble_scaling::EnsembleScaling; precision::Real = 1e-1)

Returns ensemble scaling fit data as a [`EnsembleFit`](@ref) object for the ensemble scaling data `ensemble_scaling`, displaying a warning if data is not fit to relative precision `precision`.
"""
function get_ensemble_fit(ensemble_scaling::EnsembleScaling; precision::Real = 1e-1)
    # Initialise quantities
    scaling_fit = get_scaling_fit(
        ensemble_scaling.merits,
        ensemble_scaling.dist_range;
        precision = precision,
    )
    merit_ensemble = ensemble_scaling.merit_ensemble
    dist_range = ensemble_scaling.dist_range
    # Generate the expectation and variance models
    N_params = scaling_fit.N
    ordinary_expectation_model(dist, params) = expectation_model(dist, params, N_params)
    ordinary_variance_model(dist, params) = variance_model(dist, params, N_params)
    N_marginal_params = scaling_fit.N_marginal
    marginal_expectation_model(dist, params) =
        expectation_model(dist, params, N_marginal_params)
    marginal_variance_model(dist, params) = variance_model(dist, params, N_marginal_params)
    N_relative_params = scaling_fit.N_relative
    relative_expectation_model(dist, params) =
        expectation_model(dist, params, N_relative_params)
    relative_variance_model(dist, params) = variance_model(dist, params, N_relative_params)
    # Fit the ordinary GLS merit
    gls_mean_expectation_scaling =
        [mean(merit.gls_expectation for merit in merit_set) for merit_set in merit_ensemble]
    gls_mean_variance_scaling =
        [mean(merit.gls_variance for merit in merit_set) for merit_set in merit_ensemble]
    gls_params_init = [scaling_fit.gls_tr; scaling_fit.gls_tr_sq]
    gls_pair_params = fit_pair_model(
        ordinary_expectation_model,
        ordinary_variance_model,
        dist_range,
        gls_mean_expectation_scaling,
        gls_mean_variance_scaling,
        gls_params_init;
        precision = precision,
    )
    gls_expectation_model(dist) = ordinary_expectation_model(dist, gls_pair_params)
    gls_variance_model(dist) = ordinary_variance_model(dist, gls_pair_params)
    # Fit the marginal GLS merit
    gls_marginal_mean_expectation_scaling = [
        mean(merit.gls_marginal_expectation for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    gls_marginal_mean_variance_scaling = [
        mean(merit.gls_marginal_variance for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    gls_marginal_params_init = [scaling_fit.gls_marginal_tr; scaling_fit.gls_marginal_tr_sq]
    gls_marginal_pair_params = fit_pair_model(
        marginal_expectation_model,
        marginal_variance_model,
        dist_range,
        gls_marginal_mean_expectation_scaling,
        gls_marginal_mean_variance_scaling,
        gls_marginal_params_init;
        precision = precision,
    )
    gls_marginal_expectation_model(dist) =
        marginal_expectation_model(dist, gls_marginal_pair_params)
    gls_marginal_variance_model(dist) =
        marginal_variance_model(dist, gls_marginal_pair_params)
    # Fit the relative GLS merit
    gls_relative_mean_expectation_scaling = [
        mean(merit.gls_relative_expectation for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    gls_relative_mean_variance_scaling = [
        mean(merit.gls_relative_variance for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    gls_relative_params_init = [scaling_fit.gls_relative_tr; scaling_fit.gls_relative_tr_sq]
    gls_relative_pair_params = fit_pair_model(
        relative_expectation_model,
        relative_variance_model,
        dist_range,
        gls_relative_mean_expectation_scaling,
        gls_relative_mean_variance_scaling,
        gls_relative_params_init;
        precision = precision,
    )
    gls_relative_expectation_model(dist) =
        relative_expectation_model(dist, gls_relative_pair_params)
    gls_relative_variance_model(dist) =
        relative_variance_model(dist, gls_relative_pair_params)
    # Fit the ordinary WLS merit
    wls_mean_expectation_scaling =
        [mean(merit.wls_expectation for merit in merit_set) for merit_set in merit_ensemble]
    wls_mean_variance_scaling =
        [mean(merit.wls_variance for merit in merit_set) for merit_set in merit_ensemble]
    wls_params_init = [scaling_fit.wls_tr; scaling_fit.wls_tr_sq]
    wls_pair_params = fit_pair_model(
        ordinary_expectation_model,
        ordinary_variance_model,
        dist_range,
        wls_mean_expectation_scaling,
        wls_mean_variance_scaling,
        wls_params_init;
        precision = precision,
    )
    wls_expectation_model(dist) = ordinary_expectation_model(dist, wls_pair_params)
    wls_variance_model(dist) = ordinary_variance_model(dist, wls_pair_params)
    # Fit the marginal WLS merit
    wls_marginal_mean_expectation_scaling = [
        mean(merit.wls_marginal_expectation for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    wls_marginal_mean_variance_scaling = [
        mean(merit.wls_marginal_variance for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    wls_marginal_params_init = [scaling_fit.wls_marginal_tr; scaling_fit.wls_marginal_tr_sq]
    wls_marginal_pair_params = fit_pair_model(
        marginal_expectation_model,
        marginal_variance_model,
        dist_range,
        wls_marginal_mean_expectation_scaling,
        wls_marginal_mean_variance_scaling,
        wls_marginal_params_init;
        precision = precision,
    )
    wls_marginal_expectation_model(dist) =
        marginal_expectation_model(dist, wls_marginal_pair_params)
    wls_marginal_variance_model(dist) =
        marginal_variance_model(dist, wls_marginal_pair_params)
    # Fit the relative WLS merit
    wls_relative_mean_expectation_scaling = [
        mean(merit.wls_relative_expectation for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    wls_relative_mean_variance_scaling = [
        mean(merit.wls_relative_variance for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    wls_relative_params_init = [scaling_fit.wls_relative_tr; scaling_fit.wls_relative_tr_sq]
    wls_relative_pair_params = fit_pair_model(
        relative_expectation_model,
        relative_variance_model,
        dist_range,
        wls_relative_mean_expectation_scaling,
        wls_relative_mean_variance_scaling,
        wls_relative_params_init;
        precision = precision,
    )
    wls_relative_expectation_model(dist) =
        relative_expectation_model(dist, wls_relative_pair_params)
    wls_relative_variance_model(dist) =
        relative_variance_model(dist, wls_relative_pair_params)
    # Fit the ordinary OLS merit
    ols_mean_expectation_scaling =
        [mean(merit.ols_expectation for merit in merit_set) for merit_set in merit_ensemble]
    ols_mean_variance_scaling =
        [mean(merit.ols_variance for merit in merit_set) for merit_set in merit_ensemble]
    ols_params_init = [scaling_fit.ols_tr; scaling_fit.ols_tr_sq]
    ols_pair_params = fit_pair_model(
        ordinary_expectation_model,
        ordinary_variance_model,
        dist_range,
        ols_mean_expectation_scaling,
        ols_mean_variance_scaling,
        ols_params_init;
        precision = precision,
    )
    ols_expectation_model(dist) = ordinary_expectation_model(dist, ols_pair_params)
    ols_variance_model(dist) = ordinary_variance_model(dist, ols_pair_params)
    # Fit the marginal OLS merit
    ols_marginal_mean_expectation_scaling = [
        mean(merit.ols_marginal_expectation for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    ols_marginal_mean_variance_scaling = [
        mean(merit.ols_marginal_variance for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    ols_marginal_params_init = [scaling_fit.ols_marginal_tr; scaling_fit.ols_marginal_tr_sq]
    ols_marginal_pair_params = fit_pair_model(
        marginal_expectation_model,
        marginal_variance_model,
        dist_range,
        ols_marginal_mean_expectation_scaling,
        ols_marginal_mean_variance_scaling,
        ols_marginal_params_init;
        precision = precision,
    )
    ols_marginal_expectation_model(dist) =
        marginal_expectation_model(dist, ols_marginal_pair_params)
    ols_marginal_variance_model(dist) =
        marginal_variance_model(dist, ols_marginal_pair_params)
    # Fit the relative OLS merit
    ols_relative_mean_expectation_scaling = [
        mean(merit.ols_relative_expectation for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    ols_relative_mean_variance_scaling = [
        mean(merit.ols_relative_variance for merit in merit_set) for
        merit_set in merit_ensemble
    ]
    ols_relative_params_init = [scaling_fit.ols_relative_tr; scaling_fit.ols_relative_tr_sq]
    ols_relative_pair_params = fit_pair_model(
        relative_expectation_model,
        relative_variance_model,
        dist_range,
        ols_relative_mean_expectation_scaling,
        ols_relative_mean_variance_scaling,
        ols_relative_params_init;
        precision = precision,
    )
    ols_relative_expectation_model(dist) =
        relative_expectation_model(dist, ols_relative_pair_params)
    ols_relative_variance_model(dist) =
        relative_variance_model(dist, ols_relative_pair_params)
    # Return the results
    ensemble_fit = EnsembleFit(
        gls_pair_params,
        gls_expectation_model,
        gls_variance_model,
        gls_marginal_pair_params,
        gls_marginal_expectation_model,
        gls_marginal_variance_model,
        gls_relative_pair_params,
        gls_relative_expectation_model,
        gls_relative_variance_model,
        wls_pair_params,
        wls_expectation_model,
        wls_variance_model,
        wls_marginal_pair_params,
        wls_marginal_expectation_model,
        wls_marginal_variance_model,
        wls_relative_pair_params,
        wls_relative_expectation_model,
        wls_relative_variance_model,
        ols_pair_params,
        ols_expectation_model,
        ols_variance_model,
        ols_marginal_pair_params,
        ols_marginal_expectation_model,
        ols_marginal_variance_model,
        ols_relative_pair_params,
        ols_relative_expectation_model,
        ols_relative_variance_model,
    )
    return ensemble_fit::EnsembleFit
end

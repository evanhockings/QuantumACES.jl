#
function enter_folder(folder::String)
    # Check if the folder starts with `/`
    if folder[1] != '/'
        folder = "/" * folder
    end
    # Enter the folder if not currently in the folder
    if pwd()[(end - length(folder) + 1):end] != folder
        cd(pwd() * folder)
    end
    return nothing
end

#
function exit_folder(folder::String)
    # Check if the folder starts with `/`
    if folder[1] != '/'
        folder = "/" * folder
    end
    # Exit the folder if currently in the folder
    if pwd()[(end - length(folder) + 1):end] == folder
        cd(dirname(pwd()))
    end
    return nothing
end

#
function code_filename(circuit_param::AbstractCircuitParameters)
    filename = "$(circuit_param.code_name)"
    if typeof(circuit_param) == RotatedPlanarParameters ||
       typeof(circuit_param) == UnrotatedPlanarParameters
        filename *= "_$(circuit_param.vertical_dist)_$(circuit_param.horizontal_dist)"
    else
        throw(error("Unsupported code type $(typeof(circuit_param))."))
    end
    return filename::String
end

#
function noise_filename(noise_param::AbstractNoiseParameters)
    if typeof(noise_param) == DepolarisingParameters
        filename = "depolarising"
    elseif typeof(noise_param) == LognormalParameters
        filename = "lognormal"
    else
        throw(error("Unsupported noise type $(typeof(noise_param))."))
    end
    filename *= "_$(join(round.([getfield(noise_param, field_name) for field_name in fieldnames(typeof(noise_param))], sigdigits=4),"_"))"
    return filename::String
end

#
function tuples_filename(tuple_number::Int, repeat_numbers::Vector{Int})
    if length(repeat_numbers) > 0
        filename = "$(tuple_number)_$(join(repeat_numbers, "_"))"
    else
        filename = "$(tuple_number)"
    end
    return filename::String
end

#
function design_filename(
    circuit_param::AbstractCircuitParameters,
    noise_param::AbstractNoiseParameters,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
)
    filename = "design_$(code_filename(circuit_param))_$(noise_filename(noise_param))_$(tuples_filename(tuple_number, repeat_numbers))_$(full_covariance).jld2"
    return filename::String
end

#
function design_filename(d::Design)
    filename = design_filename(
        d.code.circuit_param,
        d.code.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        d.full_covariance,
    )
    return filename::String
end

#
function dep_scaling_filename(
    circuit_param::AbstractCircuitParameters,
    noise_param::AbstractNoiseParameters,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    ls_type::Symbol,
)
    filename = "depolarising_scaling_$(code_filename(circuit_param))_$(noise_filename(noise_param))_$(tuples_filename(tuple_number, repeat_numbers))_$(ls_type).jld2"
    return filename::String
end

#
function dep_scaling_filename(d::Design, ls_type::Symbol)
    filename = dep_scaling_filename(
        d.code.circuit_param,
        d.code.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        ls_type,
    )
    return filename::String
end

#
function dep_scaling_filename(dep_scaling_data::DepolarisingScalingData)
    filename = dep_scaling_filename(
        dep_scaling_data.circuit_param,
        dep_scaling_data.noise_param,
        length(dep_scaling_data.tuple_set),
        dep_scaling_data.tuple_set_data.repeat_numbers,
        dep_scaling_data.ls_type,
    )
    return filename::String
end

#
function log_scaling_filename(
    circuit_param::AbstractCircuitParameters,
    noise_param::AbstractNoiseParameters,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    ls_type::Symbol,
)
    filename = "lognormal_scaling_$(code_filename(circuit_param))_$(noise_filename(noise_param))_$(tuples_filename(tuple_number, repeat_numbers))_$(ls_type).jld2"
    return filename::String
end

#
function log_scaling_filename(d::Design, ls_type::Symbol)
    filename = log_scaling_filename(
        d.code.circuit_param,
        d.code.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        ls_type,
    )
    return filename::String
end

#
function log_scaling_filename(log_scaling_data::LognormalScalingData)
    filename = log_scaling_filename(
        log_scaling_data.circuit_param,
        log_scaling_data.noise_param,
        length(log_scaling_data.tuple_set),
        log_scaling_data.tuple_set_data.repeat_numbers,
        log_scaling_data.ls_type,
    )
    return filename::String
end

#
function aces_data_filename(
    circuit_param::AbstractCircuitParameters,
    noise_param::AbstractNoiseParameters,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    shots_set::Vector{Int},
)
    filename = "aces_data_$(code_filename(circuit_param))_$(noise_filename(noise_param))_$(tuples_filename(tuple_number, repeat_numbers))_$(full_covariance)_$(join(round.(shots_set, sigdigits=4), "_")).jld2"
    return filename::String
end

#
function aces_data_filename(d::Design, shots_set::Vector{Int})
    filename = aces_data_filename(
        d.code.circuit_param,
        d.code.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        d.full_covariance,
        shots_set,
    )
    return filename::String
end

#
function aces_data_filename(aces_data::ACESData)
    filename = aces_data_filename(aces_data.d, aces_data.shots_set)
    return filename::String
end

"""
    Save(d::Design)

Save the design to a file specified by the design itself.
"""
function save_design(d::Design)
    # Check for the data directory
    if !isdir("data")
        mkdir("data")
    end
    # Save the design
    filename = design_filename(d)
    save_object(pwd() * "/data/" * filename, d)
    return nothing
end

"""
    Load(vertical_dist::Int, horizontal_dist::Int, circuit_number::Int, full_covariance::Bool)

Load the design from a file specified by the supplied variables.
"""
function load_design(
    circuit_param::AbstractCircuitParameters,
    noise_param::AbstractNoiseParameters,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
)
    # Load the design
    filename = design_filename(
        circuit_param,
        noise_param,
        tuple_number,
        repeat_numbers,
        full_covariance,
    )
    d = load_object(pwd() * "/data/" * filename)
    return d::Design
end

#
function delete_design(d::Design)
    # Delete the saved design data
    filename = design_filename(d)
    if isfile(pwd() * "/data/" * filename)
        rm(pwd() * "/data/" * filename)
    else
        @warn "Unable to find and delete the file."
    end
    return nothing
end

#
function save_scaling(dep_scaling_data::DepolarisingScalingData)
    # Check for the data directory
    if !isdir("data")
        mkdir("data")
    end
    # Save the design
    filename = dep_scaling_filename(dep_scaling_data)
    save_object(pwd() * "/data/" * filename, dep_scaling_data)
    return nothing
end

#
function load_scaling(
    circuit_param::AbstractCircuitParameters,
    noise_param::AbstractNoiseParameters,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    ls_type::Symbol,
    scaling_type::Symbol,
)
    # Load the design
    if scaling_type == :dep
        filename = dep_scaling_filename(
            circuit_param,
            noise_param,
            tuple_number,
            repeat_numbers,
            ls_type,
        )
        dep_scaling_data = load_object(pwd() * "/data/" * filename)
        return dep_scaling_data::DepolarisingScalingData
    elseif scaling_type == :log
        filename = log_scaling_filename(
            circuit_param,
            noise_param,
            tuple_number,
            repeat_numbers,
            ls_type,
        )
        log_scaling_data = load_object(pwd() * "/data/" * filename)
        return log_scaling_data::LognormalScalingData
    else
        throw(
            error(
                "Unsupported scaling type $(scaling_type); supported types include :dep and :log.",
            ),
        )
    end
end

#
function load_scaling(d::Design, ls_type::Symbol, scaling_type::Symbol)
    # Load the design
    if scaling_type == :dep
        filename = dep_scaling_filename(d, ls_type)
        dep_scaling_data = load_object(pwd() * "/data/" * filename)
        return dep_scaling_data::DepolarisingScalingData
    elseif scaling_type == :log
        filename = log_scaling_filename(d, ls_type)
        log_scaling_data = load_object(pwd() * "/data/" * filename)
        return log_scaling_data::LognormalScalingData
    else
        throw(
            error(
                "Unsupported scaling type $(scaling_type); supported types are :dep and :log.",
            ),
        )
    end
end

#
function delete_scaling(dep_scaling_data::DepolarisingScalingData)
    # Delete the saved design data
    filename = dep_scaling_filename(dep_scaling_data)
    if isfile(pwd() * "/data/" * filename)
        rm(pwd() * "/data/" * filename)
    else
        @warn "Unable to find and delete the file."
    end
    return nothing
end

#
function save_scaling(log_scaling_data::LognormalScalingData)
    # Check for the data directory
    if !isdir("data")
        mkdir("data")
    end
    # Save the design
    filename = log_scaling_filename(log_scaling_data)
    save_object(pwd() * "/data/" * filename, log_scaling_data)
    return nothing
end

#
function delete_scaling(log_scaling_data::LognormalScalingData)
    # Delete the saved design data
    filename = log_scaling_filename(log_scaling_data)
    if isfile(pwd() * "/data/" * filename)
        rm(pwd() * "/data/" * filename)
    else
        @warn "Unable to find and delete the file."
    end
    return nothing
end

"""
    save(aces_data::ACESData; clear_design::Bool = false)

Save the processed ACES data to a file specified by the design and the noise index.
"""
function save_aces(aces_data::ACESData; clear_design::Bool = false)
    # Check for the data directory
    if !isdir("data")
        mkdir("data")
    end
    # Save the ACES data
    filename = aces_data_filename(aces_data)
    save_object(pwd() * "/data/" * filename, aces_data)
    # Delete the design file so it isn't saved twice
    if clear_design && isfile(pwd() * "/data/" * design_filename(aces_data.d))
        delete_design(aces_data.d)
    end
    return nothing
end

"""
    Load(code_type::String, vertical_dist::Int, horizontal_dist::Int, circuit_number::Int, full_covariance::Bool, p_idx::Int)

Load the processed ACES data from a file specified by the design and the shots.
"""
function load_aces(
    circuit_param::AbstractCircuitParameters,
    noise_param::AbstractNoiseParameters,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    shots_set::Vector{Int},
)
    # Load the ACES data
    filename = aces_data_filename(
        circuit_param,
        noise_param,
        tuple_number,
        repeat_numbers,
        full_covariance,
        shots_set,
    )
    aces_data = load_object(pwd() * "/data/" * filename)
    return aces_data::ACESData
end

#
function load_aces(d::Design, shots_set::Vector{Int})
    # Load the ACES data
    filename = aces_data_filename(d, shots_set)
    aces_data = load_object(pwd() * "/data/" * filename)
    return aces_data::ACESData
end

#
function delete_aces(aces_data::ACESData)
    # Delete the saved ACES data
    filename = aces_data_filename(aces_data)
    if isfile(pwd() * "/data/" * filename)
        rm(pwd() * "/data/" * filename)
    else
        @warn "Unable to find and delete the file."
    end
    return nothing
end

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
    circuit_param::T,
    noise_param::U,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    filename = "design_$(circuit_param.circuit_name)_$(noise_param.noise_name)_$(tuples_filename(tuple_number, repeat_numbers))_$(full_covariance)"
    if ls_type == :none
        filename *= ".jld2"
    else
        filename *= "_$(ls_type).jld2"
    end
    return filename::String
end

#
function design_filename(d::Design)
    filename = design_filename(
        d.c.circuit_param,
        d.c.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        d.full_covariance,
        d.ls_type,
    )
    return filename::String
end

#
function scaling_filename(
    circuit_param::T,
    noise_param::U,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    ls_type::Symbol,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    filename = "scaling_$(circuit_param.circuit_name)_$(noise_param.noise_name)_$(tuples_filename(tuple_number, repeat_numbers))"
    if ls_type == :none
        filename *= ".jld2"
    else
        filename *= "_$(ls_type).jld2"
    end
    return filename::String
end

#
function scaling_filename(d::Design, ls_type::Symbol)
    filename = scaling_filename(
        d.c.circuit_param,
        d.c.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        ls_type,
    )
    return filename::String
end

#
function scaling_filename(scaling_data::T) where {T <: AbstractScalingData}
    filename = scaling_filename(
        scaling_data.circuit_param,
        scaling_data.noise_param,
        length(scaling_data.tuple_set),
        scaling_data.tuple_set_data.repeat_numbers,
        scaling_data.ls_type,
    )
    return filename::String
end

#
function aces_data_filename(
    circuit_param::T,
    noise_param::U,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
    shots_set::Vector{Int},
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    filename = "aces_data_$(circuit_param.circuit_name)_$(noise_param.noise_name)_$(tuples_filename(tuple_number, repeat_numbers))_$(full_covariance)_$(join(round.(shots_set, sigdigits=4), "_"))"
    if ls_type == :none
        filename *= ".jld2"
    else
        filename *= "_$(ls_type).jld2"
    end
    return filename::String
end

#
function aces_data_filename(d::Design, shots_set::Vector{Int})
    filename = aces_data_filename(
        d.c.circuit_param,
        d.c.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        d.full_covariance,
        d.ls_type,
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
    circuit_param::T,
    noise_param::U,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    # Load the design
    filename = design_filename(
        circuit_param,
        noise_param,
        tuple_number,
        repeat_numbers,
        full_covariance,
        ls_type,
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
function save_scaling(scaling_data::T) where {T <: AbstractScalingData}
    # Check for the data directory
    if !isdir("data")
        mkdir("data")
    end
    # Save the scaling data
    filename = scaling_filename(scaling_data)
    save_object(pwd() * "/data/" * filename, scaling_data)
    return nothing
end

#
function load_scaling(
    circuit_param::T,
    noise_param::U,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    ls_type::Symbol,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    # Load the scaling data
    filename =
        scaling_filename(circuit_param, noise_param, tuple_number, repeat_numbers, ls_type)
    scaling_data = load_object(pwd() * "/data/" * filename)
    return scaling_data::T where {T <: AbstractScalingData}
end

#
function load_scaling(d::Design, ls_type::Symbol)
    # Load the scaling data
    filename = scaling_filename(d, ls_type)
    scaling_data = load_object(pwd() * "/data/" * filename)
    return scaling_data::T where {T <: AbstractScalingData}
end

#
function delete_scaling(scaling_data::T) where {T <: AbstractScalingData}
    # Delete the saved design data
    filename = scaling_filename(scaling_data)
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
    Load(circuit_type::String, vertical_dist::Int, horizontal_dist::Int, circuit_number::Int, full_covariance::Bool, p_idx::Int)

Load the processed ACES data from a file specified by the design and the shots.
"""
function load_aces(
    circuit_param::T,
    noise_param::U,
    tuple_number::Int,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
    shots_set::Vector{Int},
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    # Load the ACES data
    filename = aces_data_filename(
        circuit_param,
        noise_param,
        tuple_number,
        repeat_numbers,
        full_covariance,
        ls_type,
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

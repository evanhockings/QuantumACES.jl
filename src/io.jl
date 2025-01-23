"""
    enter_folder(folder::String)

If not currently in the folder `folder`, enter it.
"""
function enter_folder(folder::String)
    # Enter the folder if not currently in the folder
    if pwd()[(end - length(folder) + 1):end] != folder
        cd(folder)
    end
    return nothing
end

"""
    exit_folder(folder::String)

If currently in the folder `folder`, exit it.
"""
function exit_folder(folder::String)
    # Exit the folder if currently in the folder
    if pwd()[(end - length(folder) + 1):end] == folder
        cd(dirname(pwd()))
    end
    return nothing
end

"""
    tuples_filename(tuple_number::Integer, repeat_numbers::Vector{Int})

Returns a string describing the filename for the supplied tuple set data.
"""
function tuples_filename(tuple_number::Integer, repeat_numbers::Vector{Int})
    if length(repeat_numbers) > 0
        filename = "$(tuple_number)_$(join(repeat_numbers, "_"))"
    else
        filename = "$(tuple_number)"
    end
    return filename::String
end

"""
    design_filename(d::Design)
    design_filename(circuit_param::AbstractCircuitParameters, noise_param::AbstractNoiseParameters, tuple_number::Integer, repeat_numbers::Vector{Int}, full_covariance::Bool, ls_type::Symbol)

Returns a string describing the filename corresponding to the supplied design data.
"""
function design_filename(
    circuit_param::T,
    noise_param::U,
    tuple_number::Integer,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    filename = "design_$(circuit_param.circuit_name)_$(noise_param.noise_name)_$(tuples_filename(tuple_number, repeat_numbers))_$(full_covariance)"
    if ls_type != :none
        filename *= "_$(ls_type)"
    end
    filename *= ".jld2"
    return filename::String
end
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

"""
    scaling_filename(scaling_data::AbstractScalingData)
    scaling_filename(d::Design, ls_type::Symbol)
    scaling_filename(circuit_param::AbstractCircuitParameters, noise_param::AbstractNoiseParameters, tuple_number::Integer, repeat_numbers::Vector{Int}, ls_type::Symbol)

Returns a string describing the filename for the scaling data corresponding to the supplied design data.
"""
function scaling_filename(
    circuit_param::T,
    noise_param::U,
    tuple_number::Integer,
    repeat_numbers::Vector{Int},
    ls_type::Symbol,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    filename = "scaling_$(circuit_param.circuit_name)_$(noise_param.noise_name)_$(tuples_filename(tuple_number, repeat_numbers))"
    if ls_type != :none
        filename *= "_$(ls_type)"
    end
    filename *= ".jld2"
    return filename::String
end
function scaling_filename(d::Design)
    filename = scaling_filename(
        d.c.circuit_param,
        d.c.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        d.ls_type,
    )
    return filename::String
end
function scaling_filename(scaling_data::T) where {T <: AbstractScalingData}
    filename = scaling_filename(scaling_data.d)
    return filename::String
end

"""
    aces_data_filename(aces_data::ACESData)
    aces_data_filename(d::Design, budget_set::Vector{Int})
    aces_data_filename(circuit_param::AbstractCircuitParameters, noise_param::AbstractNoiseParameters, tuple_number::Integer, repeat_numbers::Vector{Int}, full_covariance::Bool, ls_type::Symbol, budget_set::Vector{Int})

Returns a string describing the filename corresponding to the ACES data.
"""
function aces_data_filename(
    circuit_param::T,
    noise_param::U,
    tuple_number::Integer,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
    budget_set::Vector{Int},
    seed::UInt64,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    filename = "aces_data_$(circuit_param.circuit_name)_$(noise_param.noise_name)_$(tuples_filename(tuple_number, repeat_numbers))_$(full_covariance)_$(join(round.(budget_set, sigdigits = 4), "_"))"
    if ls_type != :none
        filename *= "_$(ls_type)"
    end
    filename *= "_$(seed).jld2"
    return filename::String
end
function aces_data_filename(d::Design, budget_set::Vector{Int}, seed::UInt64)
    filename = aces_data_filename(
        d.c.circuit_param,
        d.c.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        d.full_covariance,
        d.ls_type,
        budget_set,
        seed,
    )
    return filename::String
end
function aces_data_filename(aces_data::ACESData)
    filename = aces_data_filename(aces_data.d, aces_data.budget_set, aces_data.seed)
    return filename::String
end

"""
    rand_design_filename(d_rand::RandDesign)
    rand_design_filename(d::Design, total_randomisations::Integer, seed::UInt64)
    rand_design_filename(circuit_param::AbstractCircuitParameters, noise_param::AbstractNoiseParameters, tuple_number::Integer, repeat_numbers::Vector{Int}, full_covariance::Bool, ls_type::Symbol, total_randomisations::Integer, seed::UInt64)

Returns a string describing the filename corresponding to the randomised design data.
"""
function rand_design_filename(
    circuit_param::T,
    noise_param::U,
    tuple_number::Integer,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
    total_randomisations::Integer,
    seed::UInt64,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    # Get the 
    filename = design_filename(
        circuit_param,
        noise_param,
        tuple_number,
        repeat_numbers,
        full_covariance,
        ls_type,
    )
    @assert filename[(end - 4):end] == ".jld2"
    filename = filename[1:(end - 5)] * "_$(total_randomisations)_$(seed).jld2"
    return filename::String
end
function rand_design_filename(d::Design, total_randomisations::Integer, seed::UInt64)
    filename = rand_design_filename(
        d.c.circuit_param,
        d.c.noise_param,
        length(d.tuple_set),
        d.tuple_set_data.repeat_numbers,
        d.full_covariance,
        d.ls_type,
        total_randomisations,
        seed,
    )
    return filename::String
end
function rand_design_filename(d_rand::RandDesign)
    filename = rand_design_filename(d_rand.d, sum(d_rand.randomisations), d_rand.seed)
    return filename::String
end

"""
    save_design(d::Design)

Saves the design `d` with the appropriate filename.
"""
function save_design(d::Design)
    # Check for the data directory
    if ~isdir("data")
        mkdir("data")
    end
    # Save the design
    filename = design_filename(d)
    save_object(pwd() * "/data/" * filename, d)
    return nothing
end

"""
    load_design(circuit_param::AbstractCircuitParameters, noise_param::AbstractNoiseParameters, tuple_number::Integer, repeat_numbers::Vector{Int}, full_covariance::Bool, ls_type::Symbol)

Loads the design whose filename is specified by the supplied variables.
"""
function load_design(
    circuit_param::T,
    noise_param::U,
    tuple_number::Integer,
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

"""
    delete_design(d::Design)

Deletes the file corresponding to the design `d`.
"""
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

"""
    save_scaling(scaling_data::AbstractScalingData)

Saves the scaling data `scaling_data` with the appropriate filename.
"""
function save_scaling(scaling_data::T) where {T <: AbstractScalingData}
    # Check for the data directory
    if ~isdir("data")
        mkdir("data")
    end
    # Save the scaling data
    filename = scaling_filename(scaling_data)
    save_object(pwd() * "/data/" * filename, scaling_data)
    return nothing
end

"""
    load_scaling(d::Design, ls_type::Symbol)
    load_scaling(circuit_param::AbstractCircuitParameters, noise_param::AbstractNoiseParameters, tuple_number::Integer, repeat_numbers::Vector{Int}, ls_type::Symbol)

Loads the scaling data whose filename is specified by the supplied variables.
"""
function load_scaling(
    circuit_param::T,
    noise_param::U,
    tuple_number::Integer,
    repeat_numbers::Vector{Int},
    ls_type::Symbol,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    # Load the scaling data
    filename =
        scaling_filename(circuit_param, noise_param, tuple_number, repeat_numbers, ls_type)
    scaling_data = load_object(pwd() * "/data/" * filename)
    return scaling_data::T where {T <: AbstractScalingData}
end
function load_scaling(d::Design)
    # Load the scaling data
    filename = scaling_filename(d)
    scaling_data = load_object(pwd() * "/data/" * filename)
    return scaling_data::T where {T <: AbstractScalingData}
end

"""
    delete_scaling(scaling_data::AbstractScalingData)

Deletes the file corresponding to the scaling data `scaling_data`.
"""
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
    save_aces(aces_data::ACESData; clear_design::Bool = false)

Saves the ACES data `aces_data` with the appropriate filename, and deletes the design file if it exists and `clear_design` is `true`.
"""
function save_aces(aces_data::ACESData; clear_design::Bool = false)
    # Check for the data directory
    if ~isdir("data")
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
    load_aces(d::Design, budget_set::Vector{Int})
    load_aces(circuit_param::AbstractCircuitParameters, noise_param::AbstractNoiseParameters, tuple_number::Integer, repeat_numbers::Vector{Int}, full_covariance::Bool, ls_type::Symbol, budget_set::Vector{Int})

Loads the ACES data whose filename is specified by the supplied variables.
"""
function load_aces(
    circuit_param::T,
    noise_param::U,
    tuple_number::Integer,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
    budget_set::Vector{Int},
    seed::UInt64,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    # Load the ACES data
    filename = aces_data_filename(
        circuit_param,
        noise_param,
        tuple_number,
        repeat_numbers,
        full_covariance,
        ls_type,
        budget_set,
        seed,
    )
    aces_data = load_object(pwd() * "/data/" * filename)
    return aces_data::ACESData
end
function load_aces(d::Design, budget_set::Vector{Int}, seed::UInt64)
    # Load the ACES data
    filename = aces_data_filename(d, budget_set, seed)
    aces_data = load_object(pwd() * "/data/" * filename)
    return aces_data::ACESData
end

"""
    delete_aces(aces_data::ACESData)

Deletes the file corresponding to the ACES data `aces_data`.
"""
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

"""
    save_rand_design(d_rand::RandDesign)

Saves the randomised experimental design `d_rand` with the appropriate filename.
"""
function save_rand_design(d_rand::RandDesign)
    # Check for the data directory
    if ~isdir("data")
        mkdir("data")
    end
    # Save the randomised ensemble
    filename = rand_design_filename(d_rand)
    save_object(pwd() * "/data/" * filename, d_rand)
    return nothing
end

"""
    load_rand_design(d::Design, total_randomisations::Integer, seed::UInt64)
    load_rand_design(circuit_param::AbstractCircuitParameters, noise_param::AbstractNoiseParameters, tuple_number::Integer, repeat_numbers::Vector{Int}, full_covariance::Bool, ls_type::Symbol, total_randomisations::Integer, seed::UInt64)

Loads the randomised experimental design whose filename is specified by the supplied variables.
"""
function load_rand_design(
    circuit_param::T,
    noise_param::U,
    tuple_number::Integer,
    repeat_numbers::Vector{Int},
    full_covariance::Bool,
    ls_type::Symbol,
    total_randomisations::Integer,
    seed::UInt64,
) where {T <: AbstractCircuitParameters, U <: AbstractNoiseParameters}
    # Load the randomised ensemble
    filename = rand_design_filename(
        circuit_param,
        noise_param,
        tuple_number,
        repeat_numbers,
        full_covariance,
        ls_type,
        total_randomisations,
        seed,
    )
    d_rand = load_object(pwd() * "/data/" * filename)
    return d_rand::RandDesign
end
function load_rand_design(d::Design, total_randomisations::Integer, seed::UInt64)
    # Load the randomised ensemble
    filename = rand_design_filename(d, total_randomisations, seed)
    d_rand = load_object(pwd() * "/data/" * filename)
    return d_rand::RandDesign
end

"""
    delete_rand_design(d_rand::RandDesign)

Deletes the file corresponding to the randomised experimental design `d_rand`.
"""
function delete_rand_design(d_rand::RandDesign)
    # Delete the saved randomised ensemble
    filename = rand_design_filename(d_rand)
    if isfile(pwd() * "/data/" * filename)
        rm(pwd() * "/data/" * filename)
    else
        @warn "Unable to find and delete the file."
    end
    return nothing
end

"""
    save_qiskit_ensemble(d_rand::RandDesign, qiskit_ensemble::Py)

Saves the Qiskit ensemble `qiskit_ensemble` with the appropriate filename.
"""
function save_qiskit_ensemble(d_rand::RandDesign, qiskit_ensemble::Py)
    # Check for the data directory
    if ~isdir("data")
        mkdir("data")
    end
    # Save the Qiskit ensemble
    filename = rand_design_filename(d_rand)
    @assert filename[1:6] == "design"
    filename = "qiskit_ensemble" * filename[7:end]
    @assert filename[(end - 4):end] == ".jld2"
    filename = filename[1:(end - 5)] * ".pickle"
    open(pwd() * "/data/" * filename, "w") do f
        pickle.dump(qiskit_ensemble, f)
    end
    return nothing
end

"""
    load_qiskit_ensemble(d_rand::RandDesign)

Loads the Qiskit ensemble whose filename is specified by the supplied randomised experimental design.
"""
function load_qiskit_ensemble(d_rand::RandDesign)
    # Load the Qiskit ensemble
    filename = rand_design_filename(d_rand)
    @assert filename[1:6] == "design"
    filename = "qiskit_ensemble" * filename[7:end]
    @assert filename[(end - 4):end] == ".jld2"
    filename = filename[1:(end - 5)] * ".pickle"
    qiskit_ensemble = open(pwd() * "/data/" * filename, "r") do f
        pickle.load(f)
    end
    return qiskit_ensemble::Py
end

"""
    delete_qiskit_ensemble(d_rand::RandDesign)

Deletes the Qiskit ensemble whose filename is specified by the supplied randomised experimental design.
"""
function delete_qiskit_ensemble(d_rand::RandDesign)
    # Delete the saved Qiskit ensemble
    filename = rand_design_filename(d_rand)
    @assert filename[1:6] == "design"
    filename = "qiskit_ensemble" * filename[7:end]
    @assert filename[(end - 4):end] == ".jld2"
    filename = filename[1:(end - 5)] * ".pickle"
    if isfile(pwd() * "/data/" * filename)
        rm(pwd() * "/data/" * filename)
    else
        @warn "Unable to find and delete the file."
    end
    return nothing
end

"""
    save_rand_design_job(d_rand::RandDesign, backend::String, job_counts::Union{Vector{Matrix{UInt8}}, Vector{Dict{String, Int}}}, job_idx::Integer)

Saves the job counts data `job_counts` with index `job_idx` for the randomised experimental design `d_rand` with the filename determined by the supplied `backend`.
"""
function save_rand_design_job(
    d_rand::RandDesign,
    backend::String,
    job_counts::Union{Vector{Matrix{UInt8}}, Vector{Dict{String, Int}}},
    job_idx::Integer,
)
    # Check for the data directory
    if ~isdir("data")
        mkdir("data")
    end
    # Check for the job directory
    filename = backend * "_" * rand_design_filename(d_rand)
    @assert filename[(end - 4):end] == ".jld2"
    filename = filename[1:(end - 5)]
    if ~isdir("data/" * filename)
        mkdir("data/" * filename)
    end
    # Save the job data
    filename *= "/job_$(job_idx).jld2"
    save_object(pwd() * "/data/" * filename, job_counts)
    return nothing
end

"""
    load_rand_design_job(d_rand::RandDesign, backend::String, job_idx::Integer)

Loads the job counts data for the randomised experimental design `d_rand` whose filename is specified by the supplied variables.
"""
function load_rand_design_job(d_rand::RandDesign, backend::String, job_idx::Integer)
    # Load the job data
    filename = backend * "_" * rand_design_filename(d_rand)
    @assert filename[(end - 4):end] == ".jld2"
    filename = filename[1:(end - 5)] * "/job_$(job_idx).jld2"
    job_counts = load_object(pwd() * "/data/" * filename)
    return job_counts::Union{Vector{Matrix{UInt8}}, Vector{Dict{String, Int}}}
end

"""
    delete_rand_design_job(d_rand::RandDesign, backend::String, job_idx::Integer)

Deletes the job counts data for the randomised experimental design `d_rand` whose filename is specified by the supplied variables.
"""
function delete_rand_design_job(d_rand::RandDesign, backend::String, job_idx::Integer)
    # Delete the saved Qiskit ensemble
    filename = backend * "_" * rand_design_filename(d_rand)
    @assert filename[(end - 4):end] == ".jld2"
    filename = filename[1:(end - 5)] * "/job_$(job_idx).jld2"
    if isfile(pwd() * "/data/" * filename)
        rm(pwd() * "/data/" * filename)
    else
        @warn "Unable to find and delete the file."
    end
    return nothing
end

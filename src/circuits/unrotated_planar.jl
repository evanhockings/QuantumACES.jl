"""
    UnrotatedPlanarParameters

Parameters for the syndrome extraction circuit of an unrotated planar code.

# Fields

  - `params::Dict{Symbol, Any}`: Dictionary of the circuit parameters described below.
  - `circuit_name::String`: Name of the circuit used for saving data.

# Parameters

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `gate_type::Symbol`: Type of two-qubit gate used in the circuit, which must be `:cx`.
  - `ancilla_measurement::Bool`: Whether to include mid-circuit ancilla measurements.
  - `pad_identity::Bool`: Whether to pad layers with single-qubit identity gates.
  - `layer_time_dict::Dict{Symbol, Float64}`: Dictionary of layer times.
"""
struct UnrotatedPlanarParameters <: AbstractCircuitParameters
    params::Dict{Symbol, Any}
    circuit_name::String
    # Default constructor
    function UnrotatedPlanarParameters(params::Dict{Symbol, Any}, circuit_name::String)
        # Check circuit parameters are present
        @assert haskey(params, :vertical_dist) "The vertical distance is missing."
        @assert haskey(params, :horizontal_dist) "The horizontal distance is missing."
        @assert haskey(params, :gate_type) "The gate type is missing."
        @assert haskey(params, :ancilla_measurement) "The ancilla measurement flag is missing."
        @assert haskey(params, :pad_identity) "The pad identity flag is missing."
        @assert haskey(params, :layer_time_dict) "The layer time dictionary is missing."
        vertical_dist = params[:vertical_dist]
        horizontal_dist = params[:horizontal_dist]
        gate_type = params[:gate_type]
        ancilla_measurement = params[:ancilla_measurement]
        pad_identity = params[:pad_identity]
        layer_time_dict = params[:layer_time_dict]
        # Check some conditions
        @assert (vertical_dist >= 3 && horizontal_dist >= 3) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 3 x 3."
        @assert gate_type == :cx "Invalid gate type $(gate_type). Must be :cx."
        @assert haskey(layer_time_dict, :single_qubit) "The layer time dictionary must contain the key :single_qubit."
        @assert haskey(layer_time_dict, :two_qubit) "The layer time dictionary must contain the key :two_qubit."
        @assert haskey(layer_time_dict, :meas_reset) "The layer time dictionary must contain the key :meas_reset."
        @assert haskey(layer_time_dict, :mid_reset) "The layer time dictionary must contain the key :mid_reset."
        @assert layer_time_dict[:single_qubit] > 0.0 "The single-qubit layer time must be positive."
        @assert layer_time_dict[:two_qubit] > 0.0 "The two-qubit layer time must be positive."
        @assert layer_time_dict[:meas_reset] > 0.0 "The measurement and reset time must be positive."
        @assert layer_time_dict[:mid_reset] > 0.0 "The mid-circuit reset layer time must be positive."
        # Set the circuit name
        new_circuit_name = "unrotated_planar_$(vertical_dist)_$(horizontal_dist)"
        if gate_type != :cx
            new_circuit_name *= "_gate_type_$(gate_type)"
        end
        if ancilla_measurement != true
            new_circuit_name *= "_no_ancilla_measurement"
        end
        if pad_identity != true
            new_circuit_name *= "_no_pad_identity"
        end
        # Return parameters
        return new(params, new_circuit_name)::UnrotatedPlanarParameters
    end
end

@struct_hash_equal_isequal UnrotatedPlanarParameters

"""
    get_unrotated_param(vertical_dist::Integer, horizontal_dist::Integer; kwargs...)
    get_unrotated_param(dist::Integer; kwargs...)

Returns an [`UnrotatedPlanarParameters`](@ref) object that parameterises the syndrome extraction circuit of an unrotated planar code.

Default gate layer times are estimated from `Suppressing quantum errors by scaling a surface code logical qubit` by Google Quantum AI (2023).

# Arguments

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `dist::Int`: Distance of the code; this is equivalent to setting `vertical_dist = dist` and `horizontal_dist = dist`.

# Keyword arguments

  - `gate_type::Symbol = :cx`: Type of two-qubit gate used in the circuit, which must be `:cx`.
  - `ancilla_measurement::Bool = true`: Whether to include mid-circuit reset.
  - `pad_identity::Bool = true`: Whether to pad layers with single-qubit identity gates.
  - `single_qubit_time::Real = 29`: Time taken to implement a single-qubit gate in nanoseconds.
  - `two_qubit_time::Real = 29`: Time taken to implement a two-qubit gate in nanoseconds.
  - `meas_reset_time::Real = 660`: Time taken to perform measurement and reset at the end of the circuit in nanoseconds.
  - `mid_reset_time::Real = 660`: Time taken to perform mid-circuit reset in nanoseconds.
"""
function get_unrotated_param(
    vertical_dist::Integer,
    horizontal_dist::Integer;
    gate_type::Symbol = :cx,
    ancilla_measurement::Bool = true,
    pad_identity::Bool = true,
    single_qubit_time::Real = 29,
    two_qubit_time::Real = 29,
    meas_reset_time::Real = 660,
    mid_reset_time::Real = 660,
)
    # Construct the layer time dictionary
    layer_time_dict = Dict{Symbol, Float64}(
        :single_qubit => single_qubit_time,
        :two_qubit => two_qubit_time,
        :meas_reset => meas_reset_time,
        :mid_reset => mid_reset_time,
    )
    # Construct the circuit parameters
    params = Dict{Symbol, Any}(
        :vertical_dist => vertical_dist,
        :horizontal_dist => horizontal_dist,
        :gate_type => gate_type,
        :ancilla_measurement => ancilla_measurement,
        :pad_identity => pad_identity,
        :layer_time_dict => layer_time_dict,
    )
    # Return parameters
    unrotated_param = UnrotatedPlanarParameters(params, "unrotated_planar")
    return unrotated_param::UnrotatedPlanarParameters
end
function get_unrotated_param(
    dist::Integer;
    gate_type::Symbol = :cx,
    ancilla_measurement::Bool = true,
    pad_identity::Bool = true,
    single_qubit_time::Real = 29,
    two_qubit_time::Real = 29,
    meas_reset_time::Real = 660,
    mid_reset_time::Real = 660,
)
    # Return parameters
    unrotated_param = get_unrotated_param(
        dist,
        dist;
        gate_type = gate_type,
        ancilla_measurement = ancilla_measurement,
        pad_identity = pad_identity,
        single_qubit_time = single_qubit_time,
        two_qubit_time = two_qubit_time,
        meas_reset_time = meas_reset_time,
        mid_reset_time = mid_reset_time,
    )
    return unrotated_param::UnrotatedPlanarParameters
end

"""
    unrotated_planar_circuit(unrotated_param::UnrotatedPlanarParameters)

Returns fields used to construct the syndrome extraction circuit of an unrotated planar code based on the supplied parameters `unrotated_param`.
"""
function unrotated_planar_circuit(unrotated_param::UnrotatedPlanarParameters)
    # Set up variables
    v = unrotated_param.params[:vertical_dist]
    h = unrotated_param.params[:horizontal_dist]
    gate_type = unrotated_param.params[:gate_type]
    ancilla_measurement = unrotated_param.params[:ancilla_measurement]
    pad_identity = unrotated_param.params[:pad_identity]
    layer_time_dict = unrotated_param.params[:layer_time_dict]
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    mid_reset_type = :mid_reset
    # Generate the qubits
    data_qubits = vcat([(i, j) for i in 1:(2v - 1), j in 1:(2h - 1) if (i + j) % 2 == 0]...)
    data_num = length(data_qubits)
    ancilla_qubits =
        vcat([(i, j) for i in 1:(2v - 1), j in 1:(2h - 1) if (i + j) % 2 == 1]...)
    ancilla_num = length(ancilla_qubits)
    qubits = vcat(data_qubits, ancilla_qubits)
    qubit_num = length(qubits)
    # Generate the qubit indices
    data_indices = collect(1:data_num)
    ancilla_indices = collect((data_num + 1):qubit_num)
    ancilla_x_indices =
        ancilla_indices[([qubits[ancilla_indices][i][1] for i in 1:ancilla_num] .% 2) .== 0]
    ancilla_x_num = length(ancilla_x_indices)
    ancilla_z_indices =
        ancilla_indices[([qubits[ancilla_indices][i][2] for i in 1:ancilla_num] .% 2) .== 0]
    ancilla_z_num = length(ancilla_z_indices)
    # Generate the qubit layout
    inverse_indices = Dict(qubits[i] => i for i in 1:qubit_num)
    qubit_layout = [" " for i in 1:(2v - 1), j in 1:(2h - 1)]
    for (i, j) in qubits[data_indices]
        qubit_layout[i, j] = "o($(inverse_indices[(i, j)]))"
    end
    for (i, j) in qubits[ancilla_x_indices]
        qubit_layout[i, j] = "x($(inverse_indices[(i, j)]))"
    end
    for (i, j) in qubits[ancilla_z_indices]
        qubit_layout[i, j] = "z($(inverse_indices[(i, j)]))"
    end
    # Check qubit numbers
    @assert data_num == v * h + (v - 1) * (h - 1)
    @assert ancilla_num == v * (h - 1) + (v - 1) * h
    @assert qubit_num == (2v - 1) * (2h - 1)
    @assert qubit_num == data_num + ancilla_num
    # Generate the X ancilla gate layers
    layers_x = Vector{Vector{Int}}[[], [], [], []]
    check_x_indices = Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, ancilla_x_num)
    for (idx, ancilla_x_qubit) in pairs(qubits[ancilla_x_indices])
        # Get adjacent qubits in a reverse N order, viewed after rotating the code 45 degrees clockwise
        adj_n_qubits = [
            (ancilla_x_qubit[1] - 1, ancilla_x_qubit[2]),
            (ancilla_x_qubit[1], ancilla_x_qubit[2] + 1),
            (ancilla_x_qubit[1], ancilla_x_qubit[2] - 1),
            (ancilla_x_qubit[1] + 1, ancilla_x_qubit[2]),
        ]
        check_x_indices[idx] = ([inverse_indices[ancilla_x_qubit]], Int[])
        for (layer_idx, adj_qubit) in pairs(adj_n_qubits)
            if haskey(inverse_indices, adj_qubit)
                push!(
                    layers_x[layer_idx],
                    [inverse_indices[ancilla_x_qubit], inverse_indices[adj_qubit]],
                )
                push!(check_x_indices[idx][2], inverse_indices[adj_qubit])
            end
        end
    end
    # Generate the Z ancilla gate layers
    layers_z = Vector{Vector{Int}}[[], [], [], []]
    check_z_indices = Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, ancilla_z_num)
    for (idx, ancilla_z_qubit) in pairs(qubits[ancilla_z_indices])
        # Get adjacent qubits in a horizontally-mirrored Z order, viewed after rotating the code 45 degrees clockwise
        adj_z_qubits = [
            (ancilla_z_qubit[1] - 1, ancilla_z_qubit[2]),
            (ancilla_z_qubit[1], ancilla_z_qubit[2] - 1),
            (ancilla_z_qubit[1], ancilla_z_qubit[2] + 1),
            (ancilla_z_qubit[1] + 1, ancilla_z_qubit[2]),
        ]
        check_z_indices[idx] = ([inverse_indices[ancilla_z_qubit]], Int[])
        for (layer_idx, adj_qubit) in pairs(adj_z_qubits)
            if haskey(inverse_indices, adj_qubit)
                push!(
                    layers_z[layer_idx],
                    [inverse_indices[adj_qubit], inverse_indices[ancilla_z_qubit]],
                )
                push!(check_z_indices[idx][2], inverse_indices[adj_qubit])
            end
        end
    end
    # Construct the circuit
    if gate_type == :cx
        circuit = [
            make_layer("H", ancilla_x_indices, qubit_num),
            make_layer("CX", vcat(layers_z[1], layers_x[1]), qubit_num),
            make_layer("CX", vcat(layers_z[2], layers_x[2]), qubit_num),
            make_layer("CX", vcat(layers_z[3], layers_x[3]), qubit_num),
            make_layer("CX", vcat(layers_z[4], layers_x[4]), qubit_num),
            make_layer("H", ancilla_x_indices, qubit_num),
        ]
        layer_types = [
            single_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            single_qubit_type,
        ]
    else
        throw(error("Unsupported gate type $(gate_type)."))
    end
    # Add the ancilla measurements if appropriate
    if ancilla_measurement
        ancilla_meas_layer = make_layer("R", ancilla_indices, qubit_num)
        push!(circuit, ancilla_meas_layer)
        push!(layer_types, mid_reset_type)
    end
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [pad_layer(l) for l in circuit]
    end
    # Get the layer times
    layer_times = get_layer_times(layer_types, layer_time_dict)
    # Generate the code parameters
    init_x_indices = data_indices
    init_z_indices = Int[]
    logical_x_indices = [inverse_indices[(1, 2j - 1)] for j in 1:h]
    logical_z_indices = [inverse_indices[(2i - 1, 1)] for i in 1:v]
    code_param = CodeParameters(
        qubits,
        qubit_layout,
        inverse_indices,
        data_indices,
        ancilla_indices,
        ancilla_x_indices,
        ancilla_z_indices,
        check_x_indices,
        check_z_indices,
        init_x_indices,
        init_z_indices,
        logical_x_indices,
        logical_z_indices,
    )
    extra_fields = Dict{Symbol, Any}(:code_param => code_param)
    # Return circuit data
    return (
        circuit::Vector{Layer},
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
        extra_fields::Dict{Symbol, Any},
    )
end

"""
    get_circuit(unrotated_param::UnrotatedPlanarParameters, noise_param::AbstractNoiseParameters; kwargs...)

Returns an unrotated planar code syndrome extraction circuit in the form of a [`Circuit`](@ref) object parameterised by the supplied circuit and noise parameters.

# Arguments

  - `unrotated_param::UnrotatedPlanarParameters`: Parameters for an unrotated planar code.
  - `noise_param::AbstractNoiseParameters`: Noise parameters for the circuit.

# Keyword arguments

  - `noisy_prep::Bool = false`: Whether to treat preparations as noisy and characterise the associated noise, defaulting to `false`; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `noisy_meas::Bool = true`: Whether to treat measurements as noisy and characterise the associated noise, defaulting to `true`; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] : false,`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
  - `strict::Bool = false`: Whether to be strict about which gates count as estimable to relative precision.
"""
function get_circuit(
    unrotated_param::UnrotatedPlanarParameters,
    noise_param::T;
    noisy_prep::Bool = false,
    noisy_meas::Bool = true,
    combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] :
                     false,
    strict::Bool = false,
) where {T <: AbstractNoiseParameters}
    # Construct the circuit
    (circuit, layer_types, layer_times, extra_fields) =
        unrotated_planar_circuit(unrotated_param)
    c = get_circuit(
        circuit,
        layer_types,
        layer_times,
        noise_param;
        circuit_param = unrotated_param,
        extra_fields = extra_fields,
        noisy_prep = noisy_prep,
        noisy_meas = noisy_meas,
        combined = combined,
        strict = strict,
    )
    return c::Circuit
end

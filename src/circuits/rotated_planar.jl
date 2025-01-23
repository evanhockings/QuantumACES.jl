"""
    RotatedPlanarParameters

Parameters for the syndrome extraction circuit of a rotated surface code.

# Fields

  - `params::Dict{Symbol, Any}`: Dictionary of the circuit parameters described below.
  - `circuit_name::String`: Name of the circuit used for saving data.

# Parameters

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `check_type::Symbol`: Type of stabiliser used in the circuit, either `:xzzx` or `:standard`.
  - `gate_type::Symbol`: Type of two-qubit gate used in the circuit, either `:cx` or `:cz`.
  - `dynamically_decouple::Bool`: Whether to dynamically decouple the circuit; `true` is currently only supported for `:xzzx` and `:cz`.
  - `ancilla_measurement::Bool`: Whether to include mid-circuit ancilla measurements.
  - `pad_identity::Bool`: Whether to pad layers with single-qubit identity gates.
  - `layer_time_dict::Dict{Symbol, Float64}`: Dictionary of layer times.
"""
struct RotatedPlanarParameters <: AbstractCircuitParameters
    params::Dict{Symbol, Any}
    circuit_name::String
    # Default constructor
    function RotatedPlanarParameters(params::Dict{Symbol, Any}, circuit_name::String)
        # Check circuit parameters are present
        @assert haskey(params, :vertical_dist) "The vertical distance is missing."
        @assert haskey(params, :horizontal_dist) "The horizontal distance is missing."
        @assert haskey(params, :check_type) "The check type is missing."
        @assert haskey(params, :gate_type) "The gate type is missing."
        @assert haskey(params, :dynamically_decouple) "The dynamical decoupling flag is missing."
        @assert haskey(params, :ancilla_measurement) "The ancilla measurement flag is missing."
        @assert haskey(params, :pad_identity) "The pad identity flag is missing."
        @assert haskey(params, :layer_time_dict) "The layer time dictionary is missing."
        vertical_dist = params[:vertical_dist]
        horizontal_dist = params[:horizontal_dist]
        check_type = params[:check_type]
        gate_type = params[:gate_type]
        dynamically_decouple = params[:dynamically_decouple]
        ancilla_measurement = params[:ancilla_measurement]
        pad_identity = params[:pad_identity]
        layer_time_dict = params[:layer_time_dict]
        # Check some conditions
        @assert (vertical_dist >= 3 && horizontal_dist >= 3) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 3 x 3."
        @assert (check_type == :xzzx || check_type == :standard) "Invalid check type $(check_type). Must be either :xzzx or :standard."
        @assert (gate_type == :cx || gate_type == :cz) "Invalid gate type $(gate_type). Must be either :cx or :cz."
        @assert (check_type == :xzzx && gate_type == :cz) ||
                (check_type == :standard && gate_type == :cx) "Unsupported pairing of check type $(check_type) and gate type $(gate_type)."
        if dynamically_decouple && ~(check_type == :xzzx && gate_type == :cz)
            @warn "Dynamical decoupling is currently only supported for check type :xzzx and gate type :cz."
        end
        @assert haskey(layer_time_dict, :single_qubit) "The layer time dictionary must contain the key :single_qubit."
        @assert haskey(layer_time_dict, :two_qubit) "The layer time dictionary must contain the key :two_qubit."
        @assert haskey(layer_time_dict, :meas_reset) "The layer time dictionary must contain the key :meas_reset."
        @assert haskey(layer_time_dict, :dynamical) "The layer time dictionary must contain the key :dynamical."
        @assert haskey(layer_time_dict, :mid_reset) "The layer time dictionary must contain the key :mid_reset."
        @assert layer_time_dict[:single_qubit] > 0.0 "The single-qubit layer time must be positive."
        @assert layer_time_dict[:two_qubit] > 0.0 "The two-qubit layer time must be positive."
        @assert layer_time_dict[:meas_reset] > 0.0 "The measurement and reset time must be positive."
        @assert layer_time_dict[:dynamical] > 0.0 "The dynamical decoupling layer time must be positive."
        @assert layer_time_dict[:mid_reset] > 0.0 "The mid-circuit reset layer time must be positive."
        # Set the circuit name
        new_circuit_name = "rotated_planar_$(vertical_dist)_$(horizontal_dist)"
        if check_type != :xzzx
            new_circuit_name *= "_check_type_$(check_type)"
        end
        if gate_type != :cz
            new_circuit_name *= "_gate_type_$(gate_type)"
        end
        if dynamically_decouple != true
            new_circuit_name *= "_no_dynamically_decouple"
        end
        if ancilla_measurement != true
            new_circuit_name *= "_no_ancilla_measurement"
        end
        if pad_identity != true
            new_circuit_name *= "_no_pad_identity"
        end
        # Return parameters
        return new(params, new_circuit_name)::RotatedPlanarParameters
    end
end

@struct_hash_equal_isequal RotatedPlanarParameters

"""
    get_rotated_param(vertical_dist::Integer, horizontal_dist::Integer; kwargs...)
    get_rotated_param(dist::Integer; kwargs...)

Returns a [`RotatedPlanarParameters`](@ref) object that parameterises the syndrome extraction circuit of a rotated surface code.

Default gate layer times are estimated from `Suppressing quantum errors by scaling a surface code logical qubit` by Google Quantum AI (2023).

# Arguments

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `dist::Int`: Distance of the code; this is equivalent to setting `vertical_dist = dist` and `horizontal_dist = dist`.

# Keyword arguments

  - `check_type::Symbol = :xzzx`: Type of stabiliser used in the circuit, either `:xzzx` or `:standard`.
  - `gate_type::Symbol = :cz`: Type of two-qubit gate used in the circuit, either `:cx` or `:cz`.
  - `dynamically_decouple::Bool = true`: Whether to dynamically decouple the circuit; `true` is currently only supported for `:xzzx` and `:cz`.
  - `ancilla_measurement::Bool = true`: Whether to include mid-circuit reset.
  - `pad_identity::Bool = true`: Whether to pad layers with single-qubit identity gates.
  - `single_qubit_time::Real = 29`: Time taken to implement a single-qubit gate in nanoseconds.
  - `two_qubit_time::Real = 29`: Time taken to implement a two-qubit gate in nanoseconds.
  - `meas_reset_time::Real = 660`: Time taken to perform measurement and reset at the end of the circuit in nanoseconds.
  - `dynamical_decoupling_time::Real = 29`: Time taken to implement a dynamical decoupling layer in nanoseconds.
  - `mid_reset_time::Real = 660`: Time taken to perform mid-circuit reset in nanoseconds.
"""
function get_rotated_param(
    vertical_dist::Integer,
    horizontal_dist::Integer;
    check_type::Symbol = :xzzx,
    gate_type::Symbol = :cz,
    dynamically_decouple::Bool = true,
    ancilla_measurement::Bool = true,
    pad_identity::Bool = true,
    single_qubit_time::Real = 29,
    two_qubit_time::Real = 29,
    meas_reset_time::Real = 660,
    dynamical_decoupling_time::Real = 29,
    mid_reset_time::Real = 660,
)
    # Construct the layer time dictionary
    layer_time_dict = Dict{Symbol, Float64}(
        :single_qubit => single_qubit_time,
        :two_qubit => two_qubit_time,
        :meas_reset => meas_reset_time,
        :dynamical => dynamical_decoupling_time,
        :mid_reset => mid_reset_time,
    )
    # Construct the circuit parameters
    params = Dict{Symbol, Any}(
        :vertical_dist => vertical_dist,
        :horizontal_dist => horizontal_dist,
        :check_type => check_type,
        :gate_type => gate_type,
        :dynamically_decouple => dynamically_decouple,
        :ancilla_measurement => ancilla_measurement,
        :pad_identity => pad_identity,
        :layer_time_dict => layer_time_dict,
    )
    # Return parameters
    rotated_param = RotatedPlanarParameters(params, "rotated_planar")
    return rotated_param::RotatedPlanarParameters
end
function get_rotated_param(
    dist::Integer;
    check_type::Symbol = :xzzx,
    gate_type::Symbol = :cz,
    dynamically_decouple::Bool = true,
    ancilla_measurement::Bool = true,
    pad_identity::Bool = true,
    single_qubit_time::Real = 29,
    two_qubit_time::Real = 29,
    meas_reset_time::Real = 660,
    dynamical_decoupling_time::Real = 29,
    mid_reset_time::Real = 660,
)
    # Return parameters
    rotated_param = get_rotated_param(
        dist,
        dist;
        check_type = check_type,
        gate_type = gate_type,
        dynamically_decouple = dynamically_decouple,
        ancilla_measurement = ancilla_measurement,
        pad_identity = pad_identity,
        single_qubit_time = single_qubit_time,
        two_qubit_time = two_qubit_time,
        meas_reset_time = meas_reset_time,
        dynamical_decoupling_time = dynamical_decoupling_time,
        mid_reset_time = mid_reset_time,
    )
    return rotated_param::RotatedPlanarParameters
end

"""
    rotated_planar_circuit(rotated_param::RotatedPlanarParameters)

Returns fields used to construct the syndrome extraction circuit of a rotated surface code in the form of a [`RotatedPlanarCircuit`](@ref) object, based on the supplied parameters `rotated_param`.
"""
function rotated_planar_circuit(rotated_param::RotatedPlanarParameters)
    # Set up variables
    v = rotated_param.params[:vertical_dist]
    h = rotated_param.params[:horizontal_dist]
    check_type = rotated_param.params[:check_type]
    gate_type = rotated_param.params[:gate_type]
    dynamically_decouple = rotated_param.params[:dynamically_decouple]
    ancilla_measurement = rotated_param.params[:ancilla_measurement]
    pad_identity = rotated_param.params[:pad_identity]
    layer_time_dict = rotated_param.params[:layer_time_dict]
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    dynamical_decoupling_type = :dynamical
    mid_reset_type = :mid_reset
    # Generate the qubits
    data_qubits = vcat([(2i, 2j) for i in 1:v, j in 1:h]...)
    data_num = length(data_qubits)
    inner_ancilla_qubits = vcat([(2i + 1, 2j + 1) for i in 1:(v - 1), j in 1:(h - 1)]...)
    inner_num = length(inner_ancilla_qubits)
    boundary_locations = vcat(
        [(2i + 1, 1) for i in 1:(v - 1)],
        [(2v + 1, 2j + 1) for j in 1:(h - 1)],
        [(2v + 1 - 2i, 2h + 1) for i in 1:(v - 1)],
        [(1, 2h + 1 - 2j) for j in 1:(h - 1)],
    )
    boundary_ancilla_qubits = boundary_locations[1:2:end]
    boundary_num = length(boundary_ancilla_qubits)
    ancilla_qubits = vcat(inner_ancilla_qubits, boundary_ancilla_qubits)
    ancilla_num = length(ancilla_qubits)
    qubits = vcat(data_qubits, ancilla_qubits)
    qubit_num = length(qubits)
    # Generate the qubit indices
    data_indices = collect(1:data_num)
    ancilla_indices = collect((data_num + 1):qubit_num)
    ancilla_x_indices = ancilla_indices[(sum.(ancilla_qubits) .% 4) .== 0]
    ancilla_x_num = length(ancilla_x_indices)
    ancilla_z_indices = ancilla_indices[(sum.(ancilla_qubits) .% 4) .== 2]
    ancilla_z_num = length(ancilla_z_indices)
    # Generate the qubit layout
    inverse_indices = Dict(qubits[i] => i for i in 1:qubit_num)
    qubit_layout = [" " for i in 1:(2v + 1), j in 1:(2h + 1)]
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
    @assert data_num == v * h
    @assert inner_num == (v - 1) * (h - 1)
    @assert boundary_num == v + h - 2
    @assert ancilla_num == v * h - 1
    @assert qubit_num == 2v * h - 1
    @assert ancilla_num == inner_num + boundary_num
    @assert qubit_num == data_num + ancilla_num
    # Generate the X ancilla gate layers
    layers_x = Vector{Vector{Int}}[[], [], [], []]
    check_x_indices = Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, ancilla_x_num)
    for (idx, ancilla_x_qubit) in pairs(qubits[ancilla_x_indices])
        # Get adjacent qubits in a vertically-mirrored N order
        adj_n_qubits = [
            (ancilla_x_qubit[1] - 1, ancilla_x_qubit[2] - 1),
            (ancilla_x_qubit[1] + 1, ancilla_x_qubit[2] - 1),
            (ancilla_x_qubit[1] - 1, ancilla_x_qubit[2] + 1),
            (ancilla_x_qubit[1] + 1, ancilla_x_qubit[2] + 1),
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
        # Get adjacent qubits in a Z order
        adj_z_qubits = [
            (ancilla_z_qubit[1] - 1, ancilla_z_qubit[2] - 1),
            (ancilla_z_qubit[1] - 1, ancilla_z_qubit[2] + 1),
            (ancilla_z_qubit[1] + 1, ancilla_z_qubit[2] - 1),
            (ancilla_z_qubit[1] + 1, ancilla_z_qubit[2] + 1),
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
    if check_type == :xzzx && gate_type == :cz
        if dynamically_decouple
            circuit = [
                make_layer(["X", "H"], [data_indices, ancilla_indices], qubit_num),
                make_layer("CZ", vcat(layers_z[1], layers_x[1]), qubit_num),
                make_layer(["H", "X"], [data_indices, ancilla_indices], qubit_num),
                make_layer("CZ", vcat(layers_z[2], layers_x[2]), qubit_num),
                make_layer(["X", "X"], [data_indices, ancilla_indices], qubit_num),
                make_layer("CZ", vcat(layers_z[3], layers_x[3]), qubit_num),
                make_layer(["H", "X"], [data_indices, ancilla_indices], qubit_num),
                make_layer("CZ", vcat(layers_z[4], layers_x[4]), qubit_num),
                make_layer(["X", "H"], [data_indices, ancilla_indices], qubit_num),
            ]
            layer_types = [
                single_qubit_type,
                two_qubit_type,
                single_qubit_type,
                two_qubit_type,
                dynamical_decoupling_type,
                two_qubit_type,
                single_qubit_type,
                two_qubit_type,
                single_qubit_type,
            ]
        else
            circuit = [
                make_layer("H", ancilla_indices, qubit_num),
                make_layer("CZ", vcat(layers_z[1], layers_x[1]), qubit_num),
                make_layer("H", data_indices, qubit_num),
                make_layer("CZ", vcat(layers_z[2], layers_x[2]), qubit_num),
                make_layer("CZ", vcat(layers_z[3], layers_x[3]), qubit_num),
                make_layer("H", data_indices, qubit_num),
                make_layer("CZ", vcat(layers_z[4], layers_x[4]), qubit_num),
                make_layer("H", ancilla_indices, qubit_num),
            ]
            layer_types = [
                single_qubit_type,
                two_qubit_type,
                single_qubit_type,
                two_qubit_type,
                two_qubit_type,
                single_qubit_type,
                two_qubit_type,
                single_qubit_type,
            ]
        end
    elseif check_type == :standard && gate_type == :cx
        if dynamically_decouple
            @warn "Dynamical decoupling is currently only supported for check type :xzzx and gate type :cz."
        end
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
        throw(
            error(
                "Unsupported pairing of check type $(check_type) and gate type $(gate_type).",
            ),
        )
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
    if check_type == :xzzx
        init_x_indices = data_indices[(sum.(data_qubits) .% 4) .== 0]
        init_z_indices = data_indices[(sum.(data_qubits) .% 4) .== 2]
    elseif check_type == :standard
        init_x_indices = data_indices
        init_z_indices = Int[]
    end
    logical_x_indices = [inverse_indices[(2, 2j)] for j in 1:h]
    logical_z_indices = [inverse_indices[(2i, 2)] for i in 1:v]
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
    get_circuit(rotated_param::RotatedPlanarParameters, noise_param::AbstractNoiseParameters; kwargs...)

Returns a [`RotatedPlanarCircuit`](@ref) circuit object parameterised by the supplied circuit and noise parameters.

# Arguments

  - `rotated_param::RotatedPlanarParameters`: Parameters for a rotated surface code.
  - `noise_param::AbstractNoiseParameters`: Noise parameters for the circuit.

# Keyword arguments

  - `noisy_prep::Bool = false`: Whether to treat preparations as noisy and characterise the associated noise, defaulting to `false`; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `noisy_meas::Bool = true`: Whether to treat measurements as noisy and characterise the associated noise, defaulting to `true`; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] : false,`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
  - `strict::Bool = false`: Whether to be strict about which gates count as estimable to relative precision.
"""
function get_circuit(
    rotated_param::RotatedPlanarParameters,
    noise_param::T;
    noisy_prep::Bool = false,
    noisy_meas::Bool = true,
    combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] :
                     false,
    strict::Bool = false,
) where {T <: AbstractNoiseParameters}
    # Construct the circuit
    (circuit, layer_types, layer_times, extra_fields) =
        rotated_planar_circuit(rotated_param)
    c = get_circuit(
        circuit,
        layer_types,
        layer_times,
        noise_param;
        circuit_param = rotated_param,
        extra_fields = extra_fields,
        noisy_prep = noisy_prep,
        noisy_meas = noisy_meas,
        combined = combined,
        strict = strict,
    )
    return c::Circuit
end

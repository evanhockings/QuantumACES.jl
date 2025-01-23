"""
    HeavyHexParameters

Parameters for the syndrome extraction circuit of a heavy hex code.

# Fields

  - `params::Dict{Symbol, Any}`: Dictionary of the circuit parameters described below.
  - `circuit_name::String`: Name of the circuit used for saving data.

# Parameters

  - `vertical_dist::Int`: Vertical (Z) distance of the code.
  - `horizontal_dist::Int`: Horizontal (X) distance of the code.
  - `flipped::Bool`: Whether ancilla qubit layout is flipped left to right.
  - `gate_type::Symbol`: Type of two-qubit gate used in the circuit, which must be `:cx`.
  - `ancilla_measurement::Bool`: Whether to include mid-circuit ancilla measurements.
  - `pad_identity::Bool`: Whether to pad layers with single-qubit identity gates.
  - `layer_time_dict::Dict{Symbol, Float64}`: Dictionary of layer times.
"""
struct HeavyHexParameters <: AbstractCircuitParameters
    params::Dict{Symbol, Any}
    circuit_name::String
    # Default constructor
    function HeavyHexParameters(params::Dict{Symbol, Any}, circuit_name::String)
        # Check circuit parameters are present
        @assert haskey(params, :vertical_dist) "The vertical distance is missing."
        @assert haskey(params, :horizontal_dist) "The horizontal distance is missing."
        @assert haskey(params, :flipped) "The flipped flag is missing."
        @assert haskey(params, :gate_type) "The gate type is missing."
        @assert haskey(params, :ancilla_measurement) "The ancilla measurement flag is missing."
        @assert haskey(params, :pad_identity) "The pad identity flag is missing."
        @assert haskey(params, :layer_time_dict) "The layer time dictionary is missing."
        vertical_dist = params[:vertical_dist]
        horizontal_dist = params[:horizontal_dist]
        flipped = params[:flipped]
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
        new_circuit_name = "heavy_hex_$(vertical_dist)_$(horizontal_dist)"
        if flipped != false
            new_circuit_name *= "_flipped"
        end
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
        return new(params, new_circuit_name)::HeavyHexParameters
    end
end

@struct_hash_equal_isequal HeavyHexParameters

"""
    get_hex_param(vertical_dist::Integer, horizontal_dist::Integer; kwargs...)
    get_hex_param(dist::Integer; kwargs...)

Returns a [`HeavyHexParameters`](@ref) object that parameterises the syndrome extraction circuit of a heavy hex code.

Default gate layer times are estimated from IBM Torino in 2024.

# Arguments

  - `vertical_dist::Int`: Vertical distance of the code.
  - `horizontal_dist::Int`: Horizontal distance of the code.
  - `dist::Int`: Distance of the code; this is equivalent to setting `vertical_dist = dist` and `horizontal_dist = dist`.

# Keyword arguments

  - `flipped::Bool = false`: Whether to flip the ancilla qubit layout left to right.
  - `gate_type::Symbol = :cx`: Type of two-qubit gate used in the circuit.
  - `ancilla_measurement::Bool = true`: Whether to include mid-circuit reset.
  - `pad_identity::Bool = true`: Whether to pad layers with single-qubit identity gates.
  - `single_qubit_time::Real = 32`: Time taken to implement a single-qubit gate in nanoseconds.
  - `two_qubit_time::Real = 200`: Time taken to implement a two-qubit gate in nanoseconds.
  - `meas_reset_time::Real = 2500`: Time taken to perform measurement and reset at the end of the circuit in nanoseconds.
  - `meas_reset_time::Real = 2500`: Time taken to perform mid-circuit reset in nanoseconds.
"""
function get_hex_param(
    vertical_dist::Integer,
    horizontal_dist::Integer;
    flipped::Bool = false,
    gate_type::Symbol = :cx,
    ancilla_measurement::Bool = true,
    pad_identity = true,
    single_qubit_time::Real = 32,
    two_qubit_time::Real = 200,
    meas_reset_time::Real = 2500,
    mid_reset_time::Real = 2500,
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
        :flipped => flipped,
        :gate_type => gate_type,
        :ancilla_measurement => ancilla_measurement,
        :pad_identity => pad_identity,
        :layer_time_dict => layer_time_dict,
    )
    # Return parameters
    hex_param = HeavyHexParameters(params, "heavy_hex")
    return hex_param::HeavyHexParameters
end
function get_hex_param(
    dist::Integer;
    flipped::Bool = false,
    gate_type::Symbol = :cx,
    ancilla_measurement::Bool = true,
    pad_identity = true,
    single_qubit_time::Real = 32,
    two_qubit_time::Real = 200,
    meas_reset_time::Real = 2500,
    mid_reset_time::Real = 2500,
)
    # Return parameters
    hex_param = get_hex_param(
        dist,
        dist;
        flipped = flipped,
        gate_type = gate_type,
        ancilla_measurement = ancilla_measurement,
        pad_identity = pad_identity,
        single_qubit_time = single_qubit_time,
        two_qubit_time = two_qubit_time,
        meas_reset_time = meas_reset_time,
        mid_reset_time = mid_reset_time,
    )
    return hex_param::HeavyHexParameters
end

"""
    heavy_hex_circuit(hex_param::HeavyHexParameters)

Returns fields used to construct the syndrome extraction circuit of a heavy hex code in the form of a [`HeavyHexCircuit`](@ref) object, based on the supplied parameters `hex_param`.
"""
function heavy_hex_circuit(hex_param::HeavyHexParameters)
    # Set up variables
    v = hex_param.params[:vertical_dist]
    h = hex_param.params[:horizontal_dist]
    flipped = hex_param.params[:flipped]
    gate_type = hex_param.params[:gate_type]
    ancilla_measurement = hex_param.params[:ancilla_measurement]
    pad_identity = hex_param.params[:pad_identity]
    layer_time_dict = hex_param.params[:layer_time_dict]
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    mid_reset_type = :mid_reset
    # Generate the qubits
    data_qubits = vcat([(2i, 2j) for i in 1:v, j in 1:h]...)
    data_num = length(data_qubits)
    if flipped && isodd(h)
        isparity = isodd
    else
        isparity = iseven
    end
    ancilla_x_qubits =
        vcat([(2i + 1, 2j - 1) for i in 1:(v - 1), j in 1:(h + 1) if isparity(i + j)]...)
    ancilla_x_num = length(ancilla_x_qubits)
    ancilla_z_qubits = vcat([(2i, 2j - 1) for i in 1:v, j in 1:(h + 1)]...)
    if flipped && isodd(h)
        deleteat!(ancilla_z_qubits, findfirst(x -> x == (2, 1), ancilla_z_qubits))
        if isodd(v)
            deleteat!(ancilla_z_qubits, findfirst(x -> x == (2v, 2h + 1), ancilla_z_qubits))
        else
            deleteat!(ancilla_z_qubits, findfirst(x -> x == (2v, 1), ancilla_z_qubits))
        end
    else
        if isodd(h)
            deleteat!(ancilla_z_qubits, findfirst(x -> x == (2, 2h + 1), ancilla_z_qubits))
        end
        if isodd(v)
            deleteat!(ancilla_z_qubits, findfirst(x -> x == (2v, 1), ancilla_z_qubits))
        end
        if isodd(v + h)
            deleteat!(ancilla_z_qubits, findfirst(x -> x == (2v, 2h + 1), ancilla_z_qubits))
        end
    end
    ancilla_z_num = length(ancilla_z_qubits)
    ancilla_qubits = vcat(ancilla_x_qubits, ancilla_z_qubits)
    ancilla_num = length(ancilla_qubits)
    qubits = vcat(data_qubits, ancilla_qubits)
    qubit_num = length(qubits)
    # Generate the qubit indices
    data_indices = collect(1:data_num)
    ancilla_indices = collect((data_num + 1):qubit_num)
    ancilla_x_indices = collect((data_num + 1):(data_num + ancilla_x_num))
    ancilla_z_indices = collect((data_num + ancilla_x_num + 1):qubit_num)
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
    if iseven(v) && iseven(h)
        @assert ancilla_x_num == (v * h + v - h) // 2
        @assert ancilla_z_num == v * h + v
        @assert ancilla_num == (3v * h + 3v - h) // 2
        @assert qubit_num == (5v * h + 3v - h) // 2
    else
        @assert ancilla_x_num == (v * h + v - h - 1) // 2
        @assert ancilla_z_num == v * h + v - 2
        @assert ancilla_num == (3v * h + 3v - h - 5) // 2
        @assert qubit_num == (5v * h + 3v - h - 5) // 2
    end
    @assert ancilla_num == ancilla_x_num + ancilla_z_num
    @assert qubit_num == data_num + ancilla_num
    # Generate the X check indices
    # First generate the boundary Z ancillas
    left_ancilla_z_qubits =
        ancilla_z_qubits[[ancilla_z_qubit[2] == 1 for ancilla_z_qubit in ancilla_z_qubits]]
    right_ancilla_z_qubits = ancilla_z_qubits[[
        ancilla_z_qubit[2] == 2h + 1 for ancilla_z_qubit in ancilla_z_qubits
    ]]
    side_ancilla_z_qubits = vcat(left_ancilla_z_qubits, right_ancilla_z_qubits)
    side_ancilla_z_num = length(side_ancilla_z_qubits)
    if iseven(v) && iseven(h)
        @assert side_ancilla_z_num == 2 * v
    else
        @assert side_ancilla_z_num == 2 * v - 2
    end
    check_x_indices =
        Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, v - 1 + side_ancilla_z_num)
    for i in 1:(v - 1)
        # In the heavy hex code, the X stabilisers form horizontal strips across the code
        # This is why the code does not have a threshold
        check_x_ancilla_qubits =
            vcat([(2i + 1, 2j - 1) for j in 1:(h + 1) if isparity(i + j)]...)
        check_x_ancilla_indices =
            [inverse_indices[ancilla_x_qubit] for ancilla_x_qubit in check_x_ancilla_qubits]
        # The included data qubits are those adjacent to the strip
        check_x_data_qubits = vcat(vcat([[(2i, 2j); (2i + 2, 2j)] for j in 1:h]...)...)
        check_x_data_indices =
            [inverse_indices[data_qubit] for data_qubit in check_x_data_qubits]
        check_x_indices[i] = (check_x_ancilla_indices, check_x_data_indices)
    end
    # Append the boundary X checks (which do not include any data qubits)
    for (i, ancilla_z_qubit) in pairs(side_ancilla_z_qubits)
        check_x_indices[i + v - 1] = ([inverse_indices[ancilla_z_qubit]], Int[])
    end
    # Generate the Z check indices
    check_z_virtual_qubits =
        vcat([(2i - 1, 2j + 1) for i in 1:(v + 1), j in 1:(h - 1) if ~isparity(i + j)]...)
    check_z_num = length(check_z_virtual_qubits)
    if iseven(v) && iseven(h)
        @assert check_z_num == (v * h - v + h - 2) // 2
    else
        @assert check_z_num == (v * h - v + h - 1) // 2
    end
    check_z_indices = Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, check_z_num)
    for idx in 1:check_z_num
        check_z_virtual_qubit = check_z_virtual_qubits[idx]
        check_z_ancilla_indices = Int[]
        check_z_data_indices = Int[]
        # In the heavy hex code, the Z stabilisers are only pairs of qubits
        # They anticommute with the X stabilisers, which must then become strips
        check_z_above_qubit = (check_z_virtual_qubit[1] - 1, check_z_virtual_qubit[2])
        check_z_above_data_qubits = [
            (check_z_virtual_qubit[1] - 1, check_z_virtual_qubit[2] - 1)
            (check_z_virtual_qubit[1] - 1, check_z_virtual_qubit[2] + 1)
        ]
        if haskey(inverse_indices, check_z_above_qubit)
            push!(check_z_ancilla_indices, inverse_indices[check_z_above_qubit])
            append!(
                check_z_data_indices,
                [
                    inverse_indices[check_z_data_qubit] for
                    check_z_data_qubit in check_z_above_data_qubits
                ],
            )
        end
        check_z_below_qubit = (check_z_virtual_qubit[1] + 1, check_z_virtual_qubit[2])
        check_z_below_data_qubits = [
            (check_z_virtual_qubit[1] + 1, check_z_virtual_qubit[2] - 1)
            (check_z_virtual_qubit[1] + 1, check_z_virtual_qubit[2] + 1)
        ]
        if haskey(inverse_indices, check_z_below_qubit)
            push!(check_z_ancilla_indices, inverse_indices[check_z_below_qubit])
            append!(
                check_z_data_indices,
                [
                    inverse_indices[check_z_data_qubit] for
                    check_z_data_qubit in check_z_below_data_qubits
                ],
            )
        end
        check_z_indices[idx] = (check_z_ancilla_indices, check_z_data_indices)
    end
    # Generate the circuit gate layers
    layer_num = 14
    layers_cx = [Vector{Int}[] for idx in 1:layer_num]
    # Construct the X stabiliser gates
    for ancilla_x_qubit in qubits[ancilla_x_indices]
        above_ancilla_z_qubit = (ancilla_x_qubit[1] - 1, ancilla_x_qubit[2])
        below_ancilla_z_qubit = (ancilla_x_qubit[1] + 1, ancilla_x_qubit[2])
        # Construct layers 1-3
        above_right_data_qubit = (ancilla_x_qubit[1] - 1, ancilla_x_qubit[2] + 1)
        if haskey(inverse_indices, above_right_data_qubit)
            push!(
                layers_cx[1],
                [
                    inverse_indices[above_ancilla_z_qubit],
                    inverse_indices[above_right_data_qubit],
                ],
            )
            push!(
                layers_cx[2],
                [inverse_indices[ancilla_x_qubit], inverse_indices[above_ancilla_z_qubit]],
            )
            push!(
                layers_cx[3],
                [
                    inverse_indices[above_ancilla_z_qubit],
                    inverse_indices[above_right_data_qubit],
                ],
            )
        end
        # Construct layers 4-6
        below_right_data_qubit = (ancilla_x_qubit[1] + 1, ancilla_x_qubit[2] + 1)
        if haskey(inverse_indices, below_right_data_qubit)
            push!(
                layers_cx[4],
                [
                    inverse_indices[below_ancilla_z_qubit],
                    inverse_indices[below_right_data_qubit],
                ],
            )
            push!(
                layers_cx[5],
                [inverse_indices[ancilla_x_qubit], inverse_indices[below_ancilla_z_qubit]],
            )
            push!(
                layers_cx[6],
                [
                    inverse_indices[below_ancilla_z_qubit],
                    inverse_indices[below_right_data_qubit],
                ],
            )
        end
        # Construct layers 7-9, which resemble layers 6-4
        above_left_data_qubit = (ancilla_x_qubit[1] - 1, ancilla_x_qubit[2] - 1)
        if haskey(inverse_indices, above_left_data_qubit)
            push!(
                layers_cx[7],
                [
                    inverse_indices[above_ancilla_z_qubit],
                    inverse_indices[above_left_data_qubit],
                ],
            )
            push!(
                layers_cx[8],
                [inverse_indices[ancilla_x_qubit], inverse_indices[above_ancilla_z_qubit]],
            )
            push!(
                layers_cx[9],
                [
                    inverse_indices[above_ancilla_z_qubit],
                    inverse_indices[above_left_data_qubit],
                ],
            )
        end
        # Construct layers 10-12, which resemble layers 3-1
        below_left_data_qubit = (ancilla_x_qubit[1] + 1, ancilla_x_qubit[2] - 1)
        if haskey(inverse_indices, below_left_data_qubit)
            push!(
                layers_cx[10],
                [
                    inverse_indices[below_ancilla_z_qubit],
                    inverse_indices[below_left_data_qubit],
                ],
            )
            push!(
                layers_cx[11],
                [inverse_indices[ancilla_x_qubit], inverse_indices[below_ancilla_z_qubit]],
            )
            push!(
                layers_cx[12],
                [
                    inverse_indices[below_ancilla_z_qubit],
                    inverse_indices[below_left_data_qubit],
                ],
            )
        end
    end
    # Get the boundary X ancilla qubits
    left_ancilla_x_qubits = [(2i + 1, 1) for i in 1:(v - 1) if ~isparity(i)]
    right_ancilla_x_qubits = [(2i + 1, 2h + 1) for i in 1:(v - 1) if ~isparity(i + h)]
    @assert left_ancilla_x_qubits ⊆ ancilla_x_qubits
    @assert right_ancilla_x_qubits ⊆ ancilla_x_qubits
    # Add the left boundary gates to layers 4 and 7
    for ancilla_x_qubit in left_ancilla_x_qubits
        above_ancilla_z_qubit = (ancilla_x_qubit[1] - 1, ancilla_x_qubit[2])
        below_ancilla_z_qubit = (ancilla_x_qubit[1] + 1, ancilla_x_qubit[2])
        push!(
            layers_cx[4],
            [inverse_indices[ancilla_x_qubit], inverse_indices[above_ancilla_z_qubit]],
        )
        push!(
            layers_cx[7],
            [inverse_indices[ancilla_x_qubit], inverse_indices[below_ancilla_z_qubit]],
        )
    end
    # Add the right boundary gates to layers 6 and 9
    for ancilla_x_qubit in right_ancilla_x_qubits
        above_ancilla_z_qubit = (ancilla_x_qubit[1] - 1, ancilla_x_qubit[2])
        below_ancilla_z_qubit = (ancilla_x_qubit[1] + 1, ancilla_x_qubit[2])
        push!(
            layers_cx[6],
            [inverse_indices[ancilla_x_qubit], inverse_indices[above_ancilla_z_qubit]],
        )
        push!(
            layers_cx[9],
            [inverse_indices[ancilla_x_qubit], inverse_indices[below_ancilla_z_qubit]],
        )
    end
    # Construct the Z stabiliser gates
    for ancilla_z_qubit in qubits[ancilla_z_indices]
        # Construct layer 13, which implements Z stabiliser gates
        right_data_qubit = (ancilla_z_qubit[1], ancilla_z_qubit[2] + 1)
        if haskey(inverse_indices, right_data_qubit) && (ancilla_z_qubit[2] != 1)
            push!(
                layers_cx[13],
                [inverse_indices[right_data_qubit], inverse_indices[ancilla_z_qubit]],
            )
        end
        # Construct layer 14, which implements the other Z stabiliser gates
        left_data_qubit = (ancilla_z_qubit[1], ancilla_z_qubit[2] - 1)
        if haskey(inverse_indices, left_data_qubit) && (ancilla_z_qubit[2] != 2 * h + 1)
            push!(
                layers_cx[14],
                [inverse_indices[left_data_qubit], inverse_indices[ancilla_z_qubit]],
            )
        end
    end
    # Construct the circuit, adding the Hadamard gates to layers 1 and 12
    circuit = [make_layer("CX", layers_cx[idx], qubit_num) for idx in 1:layer_num]
    ancilla_x_hadamard_gates = make_layer("H", ancilla_x_indices, qubit_num).layer
    circuit[1] = Layer(vcat(circuit[1].layer, ancilla_x_hadamard_gates), qubit_num)
    circuit[12] = Layer(vcat(circuit[12].layer, ancilla_x_hadamard_gates), qubit_num)
    layer_types = [two_qubit_type for idx in 1:layer_num]
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
    return (
        circuit::Vector{Layer},
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
        extra_fields::Dict{Symbol, Any},
    )
end

"""
    get_circuit(hex_param::HeavyHexParameters, noise_param::AbstractNoiseParameters; kwargs...)

Returns a `HeavyHexCircuit` circuit object parameterised by the supplied circuit and noise parameters.

# Arguments

  - `hex_param::HeavyHexParameters`: Parameters for a heavy hex code.
  - `noise_param::AbstractNoiseParameters`: Noise parameters for the circuit.

# Keyword arguments

  - `noisy_prep::Bool = false`: Whether to treat preparations as noisy and characterise the associated noise, defaulting to `false`; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `noisy_meas::Bool = true`: Whether to treat measurements as noisy and characterise the associated noise, defaulting to `true`; a full-rank design cannot be produced if both `noisy_prep` and `noisy_meas` are `true`.
  - `combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] : false,`: Whether to treat Pauli X, Y, and Z basis SPAM noise as the same.
  - `strict::Bool = false`: Whether to be strict about which gates count as estimable to relative precision.
"""
function get_circuit(
    hex_param::HeavyHexParameters,
    noise_param::T;
    noisy_prep::Bool = false,
    noisy_meas::Bool = true,
    combined::Bool = haskey(noise_param.params, :combined) ? noise_param.params[:combined] :
                     false,
    strict::Bool = false,
) where {T <: AbstractNoiseParameters}
    # Construct the circuit
    (circuit, layer_types, layer_times, extra_fields) = heavy_hex_circuit(hex_param)
    c = get_circuit(
        circuit,
        layer_types,
        layer_times,
        noise_param;
        circuit_param = hex_param,
        extra_fields = extra_fields,
        noisy_prep = noisy_prep,
        noisy_meas = noisy_meas,
        combined = combined,
        strict = strict,
    )
    return c::Circuit
end

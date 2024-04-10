#
function LayerTimes(
    layer_types::Vector{Symbol};
    single_qubit_time::Float64 = 0.0,
    two_qubit_time::Float64 = 0.0,
    dynamical_decoupling_time::Float64 = 0.0,
    meas_reset_time::Float64 = 0.0,
)
    # Append the layer times
    layer_times = Vector{Float64}(undef, 0)
    for type in layer_types
        if type == :single_qubit
            push!(layer_times, single_qubit_time)
        elseif type == :two_qubit
            push!(layer_times, two_qubit_time)
        elseif type == :dynamical
            push!(layer_times, dynamical_decoupling_time)
        else
            throw(error("Unsupported layer type $(type)."))
        end
    end
    push!(layer_times, meas_reset_time)
    @assert all(layer_times .> 0.0) "The layer times must be positive; check that layer times have been supplied for all of the relevant layer types."
    return layer_times::Vector{Float64}
end

"""
    RotatedPlanar(vertical_dist::Int, horizontal_dist::Int, check_type::Symbol, gate_type::Symbol, dynamically_decouple::Bool)

Generate the error dection circuit for a rotated planar surface code with the specified vertical (Z) and horizontal (X) distances. The checks can be either "xzzx" or "standard" and the gates can be either "cz" or "cx". Dynamical decoupling currently only works for "xzzx" and "cz". Note that the most natural pairings are "xzzx" and "cz", and "standard" and "cx".
"""
function RotatedPlanar(
    vertical_dist::Int,
    horizontal_dist::Int,
    check_type::Symbol,
    gate_type::Symbol,
    dynamically_decouple::Bool,
    pad_identity::Bool,
    single_qubit_time::Float64,
    two_qubit_time::Float64,
    dynamical_decoupling_time::Float64,
    meas_reset_time::Float64,
)
    # Set up variables for notational convenience
    v = vertical_dist
    h = horizontal_dist
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    dynamical_decoupling_type = :dynamical
    # Generate the qubit indices
    data_qubits = vcat([(2i, 2j) for i in 1:v, j in 1:h]...)
    data_num = length(data_qubits)
    @assert data_num == v * h
    inner_ancilla_qubits = vcat([(2i + 1, 2j + 1) for i in 1:(v - 1), j in 1:(h - 1)]...)
    inner_num = length(inner_ancilla_qubits)
    @assert inner_num == (v - 1) * (h - 1)
    boundary_locations = vcat(
        [(2i + 1, 1) for i in 1:(v - 1)],
        [(2v + 1, 2j + 1) for j in 1:(h - 1)],
        [(2v + 1 - 2i, 2h + 1) for i in 1:(v - 1)],
        [(1, 2h + 1 - 2j) for j in 1:(h - 1)],
    )
    boundary_ancilla_qubits = boundary_locations[1:2:end]
    boundary_num = length(boundary_ancilla_qubits)
    @assert boundary_num == v + h - 2
    ancilla_qubits = vcat(inner_ancilla_qubits, boundary_ancilla_qubits)
    ancilla_num = length(ancilla_qubits)
    @assert ancilla_num == inner_num + boundary_num
    qubits = vcat(data_qubits, ancilla_qubits)
    qubit_num = length(qubits)
    @assert qubit_num == data_num + ancilla_num
    @assert qubit_num == 2v * h - 1
    data_indices = collect(1:data_num)
    ancilla_indices = collect((data_num + 1):qubit_num)
    ancilla_X_indices = ancilla_indices[(sum.(ancilla_qubits) .% 4) .== 0]
    ancilla_Z_indices = ancilla_indices[(sum.(ancilla_qubits) .% 4) .== 2]
    inverse_indices = Dict(qubits[i] => i for i in 1:qubit_num)
    # Generate the qubit layout
    qubit_layout = [" " for i in 1:(2v + 1), j in 1:(2h + 1)]
    for (i, j) in qubits[data_indices]
        qubit_layout[i, j] = "o"
    end
    for (i, j) in qubits[ancilla_X_indices]
        qubit_layout[i, j] = "x"
    end
    for (i, j) in qubits[ancilla_Z_indices]
        qubit_layout[i, j] = "z"
    end
    # Generate the X-type ancilla gate layers
    layers_X = Vector{Vector{Int}}[[], [], [], []]
    for ancilla_X_qubit in qubits[ancilla_X_indices]
        # Get indices of adjacent qubits in a vertically-mirrored N order
        adj_N_indices = [
            (ancilla_X_qubit[1] - 1, ancilla_X_qubit[2] - 1),
            (ancilla_X_qubit[1] + 1, ancilla_X_qubit[2] - 1),
            (ancilla_X_qubit[1] - 1, ancilla_X_qubit[2] + 1),
            (ancilla_X_qubit[1] + 1, ancilla_X_qubit[2] + 1),
        ]
        for (idx, adj_idx) in enumerate(adj_N_indices)
            if haskey(inverse_indices, adj_idx)
                push!(
                    layers_X[idx],
                    [inverse_indices[ancilla_X_qubit], inverse_indices[adj_idx]],
                )
            end
        end
    end
    # Generate the Z-type ancilla gate layers
    layers_Z = Vector{Vector{Int}}[[], [], [], []]
    for ancilla_Z_qubit in qubits[ancilla_Z_indices]
        # Get indices of adjacent qubits in a Z order
        adj_Z_indices = [
            (ancilla_Z_qubit[1] - 1, ancilla_Z_qubit[2] - 1),
            (ancilla_Z_qubit[1] - 1, ancilla_Z_qubit[2] + 1),
            (ancilla_Z_qubit[1] + 1, ancilla_Z_qubit[2] - 1),
            (ancilla_Z_qubit[1] + 1, ancilla_Z_qubit[2] + 1),
        ]
        for (idx, adj_idx) in enumerate(adj_Z_indices)
            if haskey(inverse_indices, adj_idx)
                push!(
                    layers_Z[idx],
                    [inverse_indices[adj_idx], inverse_indices[ancilla_Z_qubit]],
                )
            end
        end
    end
    # Construct the circuit
    if check_type == :xzzx && gate_type == :cz
        if dynamically_decouple
            circuit = [
                Layer(["X", "H"], [data_indices, ancilla_indices], qubit_num),
                Layer("CZ", vcat(layers_Z[1], layers_X[1]), qubit_num),
                Layer(["H", "X"], [data_indices, ancilla_indices], qubit_num),
                Layer("CZ", vcat(layers_Z[2], layers_X[2]), qubit_num),
                Layer(["X", "X"], [data_indices, ancilla_indices], qubit_num),
                Layer("CZ", vcat(layers_Z[3], layers_X[3]), qubit_num),
                Layer(["H", "X"], [data_indices, ancilla_indices], qubit_num),
                Layer("CZ", vcat(layers_Z[4], layers_X[4]), qubit_num),
                Layer(["X", "H"], [data_indices, ancilla_indices], qubit_num),
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
                Layer("H", ancilla_indices, qubit_num),
                Layer("CZ", vcat(layers_Z[1], layers_X[1]), qubit_num),
                Layer("H", data_indices, qubit_num),
                Layer("CZ", vcat(layers_Z[2], layers_X[2]), qubit_num),
                Layer("CZ", vcat(layers_Z[3], layers_X[3]), qubit_num),
                Layer("H", data_indices, qubit_num),
                Layer("CZ", vcat(layers_Z[4], layers_X[4]), qubit_num),
                Layer("H", ancilla_indices, qubit_num),
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
        circuit = [
            Layer("H", ancilla_X_indices, qubit_num),
            Layer("CX", vcat(layers_Z[1], layers_X[1]), qubit_num),
            Layer("CX", vcat(layers_Z[2], layers_X[2]), qubit_num),
            Layer("CX", vcat(layers_Z[3], layers_X[3]), qubit_num),
            Layer("CX", vcat(layers_Z[4], layers_X[4]), qubit_num),
            Layer("H", ancilla_X_indices, qubit_num),
        ]
        layer_types = [
            single_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            single_qubit_type,
        ]
    end
    layer_times = LayerTimes(
        layer_types;
        single_qubit_time = single_qubit_time,
        two_qubit_time = two_qubit_time,
        dynamical_decoupling_time = dynamical_decoupling_time,
        meas_reset_time = meas_reset_time,
    )
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [PadIdentity(l) for l in circuit]
    end
    # Return the code data
    return (
        circuit::Vector{Layer},
        qubit_num::Int,
        qubits::Vector{Tuple{Int, Int}},
        inverse_indices::Dict{Tuple{Int, Int}, Int},
        data_indices::Vector{Int},
        ancilla_indices::Vector{Int},
        ancilla_X_indices::Vector{Int},
        ancilla_Z_indices::Vector{Int},
        qubit_layout::Matrix{String},
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
    )
end

"""
    UnrotatedPlanar(vertical_dist::Int, horizontal_dist::Int, gate_type::Symbol, pad_identity::Bool, single_qubit_time::Float64, two_qubit_time::Float64, meas_reset_time::Float64)

Generate the error dection circuit for an unrotated planar surface code with the specified vertical and horizontal distances. The gates can be either "cx" or "cz".
"""
function UnrotatedPlanar(
    vertical_dist::Int,
    horizontal_dist::Int,
    gate_type::Symbol,
    pad_identity::Bool,
    single_qubit_time::Float64,
    two_qubit_time::Float64,
    meas_reset_time::Float64,
)
    # Set up variables for notational convenience
    v = vertical_dist
    h = horizontal_dist
    single_qubit_type = :single_qubit
    two_qubit_type = :two_qubit
    # Generate the qubit indices
    qubits = vcat([(i, j) for i in 1:(2v - 1), j in 1:(2h - 1)]...)
    qubit_num = length(qubits)
    qubit_indices = collect(1:qubit_num)
    @assert qubit_num == (2v - 1) * (2h - 1)
    data_indices = qubit_indices[(sum.(qubits) .% 2) .== 0]
    data_num = length(data_indices)
    @assert data_num == v * h + (v - 1) * (h - 1)
    ancilla_indices = qubit_indices[(sum.(qubits) .% 2) .== 1]
    ancilla_num = length(ancilla_indices)
    @assert ancilla_num == v * (h - 1) + (v - 1) * h
    ancilla_X_indices =
        ancilla_indices[([qubits[ancilla_indices][i][1] for i in 1:ancilla_num] .% 2) .== 0]
    ancilla_Z_indices =
        ancilla_indices[([qubits[ancilla_indices][i][2] for i in 1:ancilla_num] .% 2) .== 0]
    inverse_indices = Dict(qubits[i] => i for i in 1:qubit_num)
    # Generate the qubit layout
    qubit_layout = [" " for i in 1:(2v - 1), j in 1:(2h - 1)]
    for (i, j) in qubits[data_indices]
        qubit_layout[i, j] = "o"
    end
    for (i, j) in qubits[ancilla_X_indices]
        qubit_layout[i, j] = "x"
    end
    for (i, j) in qubits[ancilla_Z_indices]
        qubit_layout[i, j] = "z"
    end
    # Generate the X-type ancilla gate layers
    layers_X = Vector{Vector{Int}}[[], [], [], []]
    for ancilla_X_qubit in qubits[ancilla_X_indices]
        # Get indices of adjacent qubits in a reverse N order, viewed after rotating the code 45 degrees clockwise
        adj_N_indices = [
            (ancilla_X_qubit[1] - 1, ancilla_X_qubit[2]),
            (ancilla_X_qubit[1], ancilla_X_qubit[2] + 1),
            (ancilla_X_qubit[1], ancilla_X_qubit[2] - 1),
            (ancilla_X_qubit[1] + 1, ancilla_X_qubit[2]),
        ]
        for (idx, adj_idx) in enumerate(adj_N_indices)
            if haskey(inverse_indices, adj_idx)
                push!(
                    layers_X[idx],
                    [inverse_indices[ancilla_X_qubit], inverse_indices[adj_idx]],
                )
            end
        end
    end
    # Generate the Z-type ancilla gate layers
    layers_Z = Vector{Vector{Int}}[[], [], [], []]
    for ancilla_Z_qubit in qubits[ancilla_Z_indices]
        # Get indices of adjacent qubits in a horizontally-mirrored Z order, viewed after rotating the code 45 degrees clockwise
        adj_Z_indices = [
            (ancilla_Z_qubit[1] - 1, ancilla_Z_qubit[2]),
            (ancilla_Z_qubit[1], ancilla_Z_qubit[2] - 1),
            (ancilla_Z_qubit[1], ancilla_Z_qubit[2] + 1),
            (ancilla_Z_qubit[1] + 1, ancilla_Z_qubit[2]),
        ]
        for (idx, adj_idx) in enumerate(adj_Z_indices)
            if haskey(inverse_indices, adj_idx)
                push!(
                    layers_Z[idx],
                    [inverse_indices[adj_idx], inverse_indices[ancilla_Z_qubit]],
                )
            end
        end
    end
    # Construct the circuit
    if gate_type == :cx
        circuit = [
            Layer("H", ancilla_X_indices, qubit_num),
            Layer("CX", vcat(layers_Z[1], layers_X[1]), qubit_num),
            Layer("CX", vcat(layers_Z[2], layers_X[2]), qubit_num),
            Layer("CX", vcat(layers_Z[3], layers_X[3]), qubit_num),
            Layer("CX", vcat(layers_Z[4], layers_X[4]), qubit_num),
            Layer("H", ancilla_X_indices, qubit_num),
        ]
        layer_types = [
            single_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            two_qubit_type,
            single_qubit_type,
        ]
    end
    layer_times = LayerTimes(
        layer_types;
        single_qubit_time = single_qubit_time,
        two_qubit_time = two_qubit_time,
        meas_reset_time = meas_reset_time,
    )
    # Pad each layer with identity gates if appropriate
    if pad_identity
        circuit = [PadIdentity(l) for l in circuit]
    end
    # Return the code data
    return (
        circuit::Vector{Layer},
        qubit_num::Int,
        qubits::Vector{Tuple{Int, Int}},
        inverse_indices::Dict{Tuple{Int, Int}, Int},
        data_indices::Vector{Int},
        ancilla_indices::Vector{Int},
        ancilla_X_indices::Vector{Int},
        ancilla_Z_indices::Vector{Int},
        qubit_layout::Matrix{String},
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
    )
end

# Generate a code
function GenerateCode(code_param::AbstractCodeParameters)
    if typeof(code_param) == RotatedPlanarParameters
        return RotatedPlanar(
            code_param.vertical_dist,
            code_param.horizontal_dist,
            code_param.check_type,
            code_param.gate_type,
            code_param.dynamically_decouple,
            code_param.pad_identity,
            code_param.single_qubit_time,
            code_param.two_qubit_time,
            code_param.dynamical_decoupling_time,
            code_param.meas_reset_time,
        )
    elseif typeof(code_param) == UnrotatedPlanarParameters
        return UnrotatedPlanar(
            code_param.vertical_dist,
            code_param.horizontal_dist,
            code_param.gate_type,
            code_param.pad_identity,
            code_param.single_qubit_time,
            code_param.two_qubit_time,
            code_param.meas_reset_time,
        )
    else
        throw(error("Unsupported code parameter type $(typeof(code_param))."))
    end
end

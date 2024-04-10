module AveragedCircuitEigenvalueSampling

@warn "This package is currently in a prerelease state. The API is NOT STABLE or appropriately documented. The next release will change the entire API. Please wait for it before building on this package."

using PythonCall
using LinearAlgebra, SparseArrays
using Random, Distributions, Combinatorics, StatsBase
using GLM, LsqFit, QuadGK, FiniteDifferences, Optim, DataFrames
using Base.Threads
using Accessors, StructEquality, PrettyTables, FileIO, JLD2

mutable struct Tableau
    # Tableau
    tableau::Matrix{Bool}
    # Qubit number
    qubit_num::Int16
    # Defaunt constructor
    function Tableau(n::Integer)
        tableau = zeros(Bool, 2n + 1, 2n + 1)
        for i in 1:(2n)
            tableau[i, i] = 1
        end
        return new(tableau, n)::Tableau
    end
end

function Base.show(io::IO, t::Tableau)
    for i in 1:(t.qubit_num)
        println(io, PauliString(Pauli(t.tableau[t.qubit_num + i, :], t.qubit_num)))
    end
end

struct Gate
    # Gate type
    type::String
    # Index labelling the unique occurrences of the gate in a circuit
    # Two gates are considered unique if they occur in non-identical layers
    index::Int32
    # Gate targets
    targets::Vector{Int16}
end

Base.show(io::IO, g::Gate) = print(io, "($(g.type)-$(Int(g.index)):$(Int.(g.targets)))")

function Base.isless(g₁::Gate, g₂::Gate)
    return isless([g₁.type; g₁.index; g₁.targets], [g₂.type; g₂.index; g₂.targets])
end

struct Layer
    # Layer of gates
    layer::Vector{Gate}
    # Qubit number
    qubit_num::Int16
    # Default constructor
    function Layer(layer::Vector{Gate}, qubit_num::Integer)
        # Parameter checks
        target_set = sort(vcat([gate.targets for gate in layer]...))
        @assert all(target_set .>= 1) && all(target_set .<= qubit_num) "The layer $(layer) acts on a qubit that is out of bounds on $(qubit_num) qubits."
        @assert target_set == unique(target_set) "The layer $(layer) acts on the same target more than once."
        return new(layer, qubit_num)::Layer
    end
end

Base.show(io::IO, l::Layer) = print(io, [gate for gate in l.layer])

abstract type AbstractNoiseParameters end

struct DepolarisingParameters <: AbstractNoiseParameters
    # Single-qubit gate entanglement infidelity
    # This is the sum of all 3 non-identity Pauli error probabilities
    r_1::Float64
    # Two-qubit gate entanglement infidelity
    # This is the sum of all 15 non-identity Pauli error probabilities
    r_2::Float64
    # Measurement entanglement infidelity
    # This is the measurement error probability
    r_m::Float64
    # Default constructor
    function DepolarisingParameters(r_1::Float64, r_2::Float64, r_m::Float64)
        # Check the nosie parameters
        @assert (r_1 >= 0) && (r_1 <= 3 / 4) "The single-qubit gate entanglement infidelity $(r_1) is out of bounds."
        @assert (r_2 >= 0) && (r_2 <= 15 / 16) "The two-qubit gate entanglement infidelity $(r_2) is out of bounds."
        @assert (r_m >= 0) && (r_m <= 1 / 2) "The measurement entanglement infidelity $(r_m) is out of bounds."
        # Return parameters
        return new(r_1, r_2, r_m)::DepolarisingParameters
    end
end

struct LogNormalParameters <: AbstractNoiseParameters
    # Mean of the single-qubit gate entanglement infidelity
    r_1::Float64
    # Mean of the two-qubit gate entanglement infidelity
    r_2::Float64
    # Mean of the measurement entanglement infidelity
    r_m::Float64
    # Approximate standard deviation of the logarithm of the entanglement infidelity
    total_std_log::Float64
    # Random seed
    seed::UInt64
    # Constructor
    function LogNormalParameters(
        r_1::Float64,
        r_2::Float64,
        r_m::Float64,
        total_std_log::Float64;
        seed::Union{UInt64, Nothing} = nothing,
    )
        # Check the nosie parameters
        @assert (r_1 >= 0) && (r_1 <= 3 / 4) "The single-qubit gate entanglement infidelity $(r_1) is out of bounds."
        @assert (r_2 >= 0) && (r_2 <= 15 / 16) "The two-qubit gate entanglement infidelity $(r_2) is out of bounds."
        @assert (r_m >= 0) && (r_m <= 1 / 2) "The measurement entanglement infidelity $(r_m) is out of bounds."
        # Randomly set the seed if one isn't supplied
        if seed === nothing
            seed = rand(UInt64)
        end
        # Return parameters
        return new(
            r_1::Float64,
            r_2::Float64,
            r_m::Float64,
            total_std_log::Float64,
            seed::UInt64,
        )::LogNormalParameters
    end
end

abstract type AbstractCodeParameters end

struct RotatedPlanarParameters <: AbstractCodeParameters
    # Vertical and Z distance of the code
    vertical_dist::Int
    # Horizontal and X distance of the code
    horizontal_dist::Int
    # Type of stabiliser used in the circuit
    check_type::Symbol
    # Type of two-qubit gate used in the circuit
    gate_type::Symbol
    # Whether to dynamically decouple the circuit
    dynamically_decouple::Bool
    # Whether to pad layers with single-qubit identity gates
    pad_identity::Bool
    # Single-qubit gate layer time
    single_qubit_time::Float64
    # Two-qubit gate layer time
    two_qubit_time::Float64
    # Dynamical decoupling gate layer time
    dynamical_decoupling_time::Float64
    # Measurement and reset time
    meas_reset_time::Float64
    # Name of the code for saving data
    code_name::String
    # Default constructor
    function RotatedPlanarParameters(
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
        code_name::String,
    )
        # Check some conditions
        @assert (vertical_dist >= 2 && horizontal_dist >= 2) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 2 x 2."
        @assert (check_type == :xzzx || check_type == :standard) "Invalid check type $(check_type). Must be either :xzzx or :standard."
        @assert (gate_type == :cx || gate_type == :cz) "Invalid gate type $(gate_type). Must be either :cx or :cz."
        @assert (check_type == :xzzx && gate_type == :cz) ||
                (check_type == :standard && gate_type == :cx) "Unsupported pairing of check type $(check_type) and gate type $(gate_type)."
        if dynamically_decouple && ~(check_type == :xzzx && gate_type == :cz)
            @warn "Dynamical decoupling is only supported for check type :xzzx and gate type :cz."
        end
        test_code_name = "rotated_planar"
        if check_type != :xzzx
            test_code_name *= "_check_type_$(check_type)"
        end
        if gate_type != :cz
            test_code_name *= "_gate_type_$(gate_type)"
        end
        if dynamically_decouple != true
            test_code_name *= "_dynamically_decouple_$(dynamically_decouple)"
        end
        if pad_identity != true
            test_code_name *= "_pad_identity_$(pad_identity)"
        end
        @assert code_name == test_code_name "The code name $(code_name) does not match the code name generated by the supplied parameters $(test_code_name)."
        # Return parameters
        return new(
            vertical_dist,
            horizontal_dist,
            check_type,
            gate_type,
            dynamically_decouple,
            pad_identity,
            single_qubit_time,
            two_qubit_time,
            dynamical_decoupling_time,
            meas_reset_time,
            code_name,
        )::RotatedPlanarParameters
    end
    # Constructor
    # The default gate times are specified in nanoseconds (though units ultimately don't matter) and estimated from Google device data
    # In `Suppressing quantum errors by scaling a surface code logical qubit`, they specify measurement takes 500 ns and reset takes 160 ns
    # They also specify that the overall circuit, including measurement and reset, takes 921 ns
    # They say they achieve similar or improved results as `Exponential suppression of bit or phase errors with cyclic error correction`
    # This specifies 26 ns CZ gates, and 80 ns for two layers of H and two layers of CZ, implying 14 nz H gates
    # This implies Hadamard gates take 14 ns, but the single-qubit gate layers in the original paper take on average 31.4 ns
    # If we assume the dynamical decoupling X gates are decomposed into 2 H gates and a Z rotation, then we can imagine the single-qubit gate layers taking 31.4 ns
    # However, there's sufficient ambiguity that we'll simply treat all layers as taking the same amount of time, namely 29 ns
    function RotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int;
        check_type::Symbol = :xzzx,
        gate_type::Symbol = :cz,
        dynamically_decouple::Bool = true,
        pad_identity::Bool = true,
        single_qubit_time::Float64 = 29.0,
        two_qubit_time::Float64 = 29.0,
        dynamical_decoupling_time::Float64 = 29.0,
        meas_reset_time::Float64 = 660.0,
    )
        # Check some conditions
        @assert (vertical_dist >= 2 && horizontal_dist >= 2) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 2 x 2."
        @assert (check_type == :xzzx || check_type == :standard) "Invalid check type $(check_type). Must be either :xzzx or :standard."
        @assert (gate_type == :cx || gate_type == :cz) "Invalid gate type $(gate_type). Must be either :cx or :cz."
        @assert (check_type == :xzzx && gate_type == :cz) ||
                (check_type == :standard && gate_type == :cx) "Unsupported pairing of check type $(check_type) and gate type $(gate_type)."
        if dynamically_decouple && ~(check_type == :xzzx && gate_type == :cz)
            @warn "Dynamical decoupling is only supported for check type :xzzx and gate type :cz."
        end
        # Generate the code name
        code_name = "rotated_planar"
        if check_type != :xzzx
            code_name *= "_check_type_$(check_type)"
        end
        if gate_type != :cz
            code_name *= "_gate_type_$(gate_type)"
        end
        if dynamically_decouple != true
            code_name *= "_dynamically_decouple_$(dynamically_decouple)"
        end
        if pad_identity != true
            code_name *= "_pad_identity_$(pad_identity)"
        end
        # Return parameters
        return new(
            vertical_dist,
            horizontal_dist,
            check_type,
            gate_type,
            dynamically_decouple,
            pad_identity,
            single_qubit_time,
            two_qubit_time,
            dynamical_decoupling_time,
            meas_reset_time,
            code_name,
        )::RotatedPlanarParameters
    end
    # Square constructor
    function RotatedPlanarParameters(
        dist::Int;
        check_type::Symbol = :xzzx,
        gate_type::Symbol = :cz,
        dynamically_decouple::Bool = true,
        pad_identity::Bool = true,
        single_qubit_time::Float64 = 29.0,
        two_qubit_time::Float64 = 29.0,
        dynamical_decoupling_time::Float64 = 29.0,
        meas_reset_time::Float64 = 660.0,
    )
        # Return parameters
        return RotatedPlanarParameters(
            dist,
            dist;
            check_type = check_type,
            gate_type = gate_type,
            dynamically_decouple = dynamically_decouple,
            pad_identity = pad_identity,
            single_qubit_time = single_qubit_time,
            two_qubit_time = two_qubit_time,
            dynamical_decoupling_time = dynamical_decoupling_time,
            meas_reset_time = meas_reset_time,
        )::RotatedPlanarParameters
    end
end

struct UnrotatedPlanarParameters <: AbstractCodeParameters
    # Vertical and Z distance of the code
    vertical_dist::Int
    # Horizontal and X distance of the code
    horizontal_dist::Int
    # Type of two-qubit gate used in the circuit
    gate_type::Symbol
    # Whether to pad layers with single-qubit identity gates
    pad_identity::Bool
    # Single-qubit gate layer time
    single_qubit_time::Float64
    # Two-qubit gate layer time
    two_qubit_time::Float64
    # Measurement and reset time
    meas_reset_time::Float64
    # Name of the code for saving data
    code_name::String
    # Default constructor
    function UnrotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int,
        gate_type::Symbol,
        pad_identity::Bool,
        single_qubit_time::Float64,
        two_qubit_time::Float64,
        meas_reset_time::Float64,
        code_name::String,
    )
        # Check some conditions
        @assert (vertical_dist >= 2 && horizontal_dist >= 2) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 2 x 2."
        @assert gate_type == :cx "Invalid gate type $(gate_type). Must be :cx."
        test_code_name = "unrotated_planar"
        if gate_type != :cx
            test_code_name *= "_gate_type_$(gate_type)"
        end
        if pad_identity != true
            test_code_name *= "_pad_identity_$(pad_identity)"
        end
        @assert code_name == test_code_name "The code name $(code_name) does not match the code name generated by the supplied parameters $(test_code_name)."
        # Return parameters
        return new(
            vertical_dist,
            horizontal_dist,
            gate_type,
            pad_identity,
            single_qubit_time,
            two_qubit_time,
            meas_reset_time,
            code_name,
        )::UnrotatedPlanarParameters
    end
    # Constructor
    # The default gate times are specified in nanoseconds (though units ultimately don't matter) and estimated from Google device data
    # In `Suppressing quantum errors by scaling a surface code logical qubit`, they specify measurement takes 500 ns and reset takes 160 ns
    # They also specify that the overall circuit, including measurement and reset, takes 921 ns
    # They say they achieve similar or improved results as `Exponential suppression of bit or phase errors with cyclic error correction`
    # This specifies 26 ns CZ gates, and 80 ns for two layers of H and two layers of CZ, implying 14 nz H gates
    # This implies Hadamard gates take 14 ns, but the single-qubit gate layers in the original paper take on average 31.4 ns
    # If we assume the dynamical decoupling X gates are decomposed into 2 H gates and a Z rotation, then we can imagine the single-qubit gate layers taking 31.4 ns
    # However, there's sufficient ambiguity that we'll simply treat all layers as taking the same amount of time, namely 29 ns
    function UnrotatedPlanarParameters(
        vertical_dist::Int,
        horizontal_dist::Int;
        gate_type::Symbol = :cx,
        pad_identity::Bool = true,
        single_qubit_time::Float64 = 29.0,
        two_qubit_time::Float64 = 29.0,
        meas_reset_time::Float64 = 660.0,
    )
        # Check some conditions
        @assert (vertical_dist >= 2 && horizontal_dist >= 2) "Invalid distance $(vertical_dist) x $(horizontal_dist). Must be at least 2 x 2."
        @assert gate_type == :cx "Invalid gate type $(gate_type). Must be :cx."
        # Generate the code name
        code_name = "unrotated_planar"
        if gate_type != :cx
            code_name *= "_gate_type_$(gate_type)"
        end
        if pad_identity != true
            code_name *= "_pad_identity_$(pad_identity)"
        end
        # Return parameters
        return new(
            vertical_dist,
            horizontal_dist,
            gate_type,
            pad_identity,
            single_qubit_time,
            two_qubit_time,
            meas_reset_time,
            code_name,
        )::UnrotatedPlanarParameters
    end
    # Square constructor
    function UnrotatedPlanarParameters(
        dist::Int;
        gate_type::Symbol = :cx,
        pad_identity::Bool = true,
        single_qubit_time::Float64 = 29.0,
        two_qubit_time::Float64 = 29.0,
        meas_reset_time::Float64 = 660.0,
    )
        # Return parameters
        return UnrotatedPlanarParameters(
            dist,
            dist;
            gate_type = gate_type,
            pad_identity = pad_identity,
            single_qubit_time = single_qubit_time,
            two_qubit_time = two_qubit_time,
            meas_reset_time = meas_reset_time,
        )::UnrotatedPlanarParameters
    end
end

abstract type AbstractCircuit end

struct Circuit <: AbstractCircuit
    # Circuit
    circuit::Vector{Layer}
    # Tuple indexing the order of the circuit layers
    circuit_tuple::Vector{Int}
    # Qubit number
    qubit_num::Int
    # Whether to treat preparations as noisy and aim to characterise them
    add_prep::Bool
    # Whether to treat preparations as noisy and aim to characterise them
    add_meas::Bool
    # Indices of the unique layers in the original circuit, which become meaningless and are removed when a tuple is applied
    unique_layer_indices::Vector{Int}
    # Gates in the circuit
    gates::Vector{Gate}
    # Total gates in the circuit
    # Includes preparations if add_prep and measurements if add_meas
    total_gates::Vector{Gate}
    # Gate index labelling the ordering of the gate eigenvalues
    gate_index::Dict{Gate, Int}
    # Total number of gate eigenvalues
    N::Int
    # Default constructor
    function Circuit(
        circuit::Vector{Layer},
        circuit_tuple::Vector{Int},
        qubit_num::Int,
        add_prep::Bool,
        add_meas::Bool,
        unique_layer_indices::Vector{Int},
        gates::Vector{Gate},
        total_gates::Vector{Gate},
        gate_index::Dict{Gate, Int},
        N::Int,
    )
        # Return the circuit
        return new(
            circuit,
            circuit_tuple,
            qubit_num,
            add_prep,
            add_meas,
            unique_layer_indices,
            gates,
            total_gates,
            gate_index,
            N,
        )::Circuit
    end
    # Constructor
    function Circuit(
        circuit::Vector{Layer},
        qubit_num::Int;
        add_prep::Bool = false,
        add_meas::Bool = true,
    )
        # Generate the code
        circuit_tuple = collect(1:length(circuit))
        # Label the circuit
        (circuit, unique_layer_indices) = Label(circuit, qubit_num)
        # Generate the gates, total gates, and noise
        gates = Gates(circuit)
        (total_gates, gate_index, N) = Index(gates, qubit_num, add_prep, add_meas)
        # Return the circuit
        return new(
            circuit,
            circuit_tuple,
            qubit_num,
            add_prep,
            add_meas,
            unique_layer_indices,
            gates,
            total_gates,
            gate_index,
            N,
        )::Circuit
    end
end

struct Code <: AbstractCircuit
    # Code parameters
    code_param::AbstractCodeParameters
    # Circuit
    circuit::Vector{Layer}
    # Tuple indexing the order of the circuit layers
    circuit_tuple::Vector{Int}
    # Qubit number
    qubit_num::Int
    # Code qubits
    qubits::Vector{Tuple{Int, Int}}
    # Inverse of the code qubit indices
    inverse_indices::Dict{Tuple{Int, Int}, Int}
    # Data qubit indices
    data_indices::Vector{Int}
    # Ancilla qubit indices
    ancilla_indices::Vector{Int}
    # Ancilla X-check qubit indices
    ancilla_X_indices::Vector{Int}
    # Ancilla Z-check qubit indices
    ancilla_Z_indices::Vector{Int}
    # Code qubit layout
    qubit_layout::Matrix{String}
    # Whether to treat preparations as noisy and aim to characterise them
    add_prep::Bool
    # Whether to treat preparations as noisy and aim to characterise them
    add_meas::Bool
    # Indices of the unique layers in the original circuit, which become meaningless and are removed when a tuple is applied
    unique_layer_indices::Vector{Int}
    # Type of each layer
    layer_types::Vector{Symbol}
    # Time taken to perform each layer, including measurement and reset at the end
    layer_times::Vector{Float64}
    # Gates in the circuit
    gates::Vector{Gate}
    # Total gates in the circuit
    # Includes preparations if add_prep and measurements if add_meas
    total_gates::Vector{Gate}
    # Gate index labelling the ordering of the gate eigenvalues
    gate_index::Dict{Gate, Int}
    # Total number of gate eigenvalues
    N::Int
    # Noise parameters
    noise_param::AbstractNoiseParameters
    # Gate probabilities
    gate_probabilities::Dict{Gate, Vector{Float64}}
    # Gate eigenvalues
    gate_eigenvalues::Vector{Float64}
    # Default constructor
    function Code(
        code_param::AbstractCodeParameters,
        circuit::Vector{Layer},
        circuit_tuple::Vector{Int},
        qubit_num::Int,
        qubits::Vector{Tuple{Int, Int}},
        inverse_indices::Dict{Tuple{Int, Int}, Int},
        data_indices::Vector{Int},
        ancilla_indices::Vector{Int},
        ancilla_X_indices::Vector{Int},
        ancilla_Z_indices::Vector{Int},
        qubit_layout::Matrix{String},
        add_prep::Bool,
        add_meas::Bool,
        unique_layer_indices::Vector{Int},
        layer_types::Vector{Symbol},
        layer_times::Vector{Float64},
        gates::Vector{Gate},
        total_gates::Vector{Gate},
        gate_index::Dict{Gate, Int},
        N::Int,
        noise_param::AbstractNoiseParameters,
        gate_probabilities::Dict{Gate, Vector{Float64}},
        gate_eigenvalues::Vector{Float64},
    )
        # Return the code
        return new(
            code_param,
            circuit,
            circuit_tuple,
            qubit_num,
            qubits,
            inverse_indices,
            data_indices,
            ancilla_indices,
            ancilla_X_indices,
            ancilla_Z_indices,
            qubit_layout,
            add_prep,
            add_meas,
            unique_layer_indices,
            layer_types,
            layer_times,
            gates,
            total_gates,
            gate_index,
            N,
            noise_param,
            gate_probabilities,
            gate_eigenvalues,
        )::Code
    end
    # Constructor
    function Code(
        code_param::AbstractCodeParameters,
        noise_param::AbstractNoiseParameters;
        add_prep::Bool = false,
        add_meas::Bool = true,
    )
        # Generate the code
        (
            circuit,
            qubit_num,
            qubits,
            inverse_indices,
            data_indices,
            ancilla_indices,
            ancilla_X_indices,
            ancilla_Z_indices,
            qubit_layout,
            layer_types,
            layer_times,
        ) = GenerateCode(code_param)
        circuit_tuple = collect(1:length(circuit))
        # Check the parameters
        @assert length(layer_times) == length(circuit) + 1 "The layer times correspond to the times taken for the circuit layers, alongside measurement and reset at the end."
        for (idx, type) in enumerate(layer_types)
            if type == :single_qubit
                @assert layer_times[idx] == code_param.single_qubit_time "The layer time $(layer_times[idx]) does not match the single-qubit layer gate time $(code_param.single_qubit_time)."
                @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 1 "The single-qubit layer $(circuit[idx]) does not contain only single-qubit gates."
            elseif type == :two_qubit
                @assert layer_times[idx] == code_param.two_qubit_time "The layer time $(layer_times[idx]) does not match the two-qubit layer gate time $(code_param.two_qubit_time)."
                @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 2 "The two-qubit layer $(circuit[idx]) does not contain two-qubit gates."
            elseif type == :dynamical
                @assert layer_times[idx] == code_param.dynamical_decoupling_time "The layer time $(layer_times[idx]) does not match the dynamical decoupling layer gate time $(code_param.dynamical_decoupling_time)."
                @assert maximum(length(gate.targets) for gate in circuit[idx].layer) == 1 "The dynamical decoupling layer $(circuit[idx]) does not contain only single-qubit gates."
            else
                throw(error("Unsupported layer type $(type)."))
            end
        end
        @assert layer_times[end] == code_param.meas_reset_time "The layer time $(layer_times[end]) does not match the measurement and reset time $(code_param.meas_reset_time)."
        # Label the circuit
        (circuit, unique_layer_indices) = Label(circuit, qubit_num)
        # Generate the gates, total gates, and noise
        gates = Gates(circuit)
        (total_gates, gate_index, N) = Index(gates, qubit_num, add_prep, add_meas)
        gate_probabilities = GenerateNoise(total_gates, noise_param)
        gate_eigenvalues = GateEigenvalues(gate_probabilities, total_gates, gate_index, N)
        # Return the code
        return new(
            code_param,
            circuit,
            circuit_tuple,
            qubit_num,
            qubits,
            inverse_indices,
            data_indices,
            ancilla_indices,
            ancilla_X_indices,
            ancilla_Z_indices,
            qubit_layout,
            add_prep,
            add_meas,
            unique_layer_indices,
            layer_types,
            layer_times,
            gates,
            total_gates,
            gate_index,
            N,
            noise_param,
            gate_probabilities,
            gate_eigenvalues,
        )::Code
    end
end

Base.show(io::IO, code::Code) = show(io, MIME("text/plain"), code.qubit_layout)

struct Pauli
    # Boolean storage of a Pauli
    pauli::Vector{Bool}
    # Qubit number
    qubit_num::Int16
    # Constructor
    function Pauli(pauli::Vector{Bool}, n::Integer)
        @assert length(pauli) == 2n + 1 "The Pauli $(pauli) does not have a length commensurate with the qubit number $(n)."
        return new(pauli, n)::Pauli
    end
end

Base.show(io::IO, p::Pauli) = print(io, PauliString(p))

function Base.:(+)(p₁::Pauli, p₂::Pauli)
    @assert p₁.qubit_num == p₂.qubit_num "The Paulis $(p₁) and $(p₂) do not have the same number of qubits."
    return Pauli(convert(Vector{Bool}, (p₁.pauli .+ p₂.pauli) .% 2), p₁.qubit_num)
end

struct Mapping
    # Initial Pauli after preparation
    initial::Pauli
    # Final Pauli after the action of the circuit but before measurement
    final::Pauli
    # Row of the design matrix corresponding to the initial Pauli and the circuit which acts on the initial Pauli to produce the final Pauli
    design_row::SparseVector{Int32, Int32}
    # Track the support of the Pauli as it is acted upon by the layers of the circuit
    spread_track::Vector{Vector{Int16}}
end

function Base.show(io::IO, m::Mapping)
    return print(io, PauliString(m.initial) * " => " * PauliString(m.final))
end

mutable struct TupleSetData
    # The main tuple set
    tuple_set::Vector{Vector{Int}}
    # The tuple set whose tuples are repeated
    repeat_tuple_set::Vector{Vector{Int}}
    # The number of repetitions for each tuple
    repeat_numbers::Vector{Int}
    # Maps repeat numbers to tuples in the repeated tuple set
    repeat_indices::Vector{Int}
    # Default constructor
    function TupleSetData(
        tuple_set::Vector{Vector{Int}},
        repeat_tuple_set::Vector{Vector{Int}},
        repeat_numbers::Vector{Int},
        repeat_indices::Vector{Int},
    )
        @assert all(repeat_numbers .>= 0) "The repeat numbers must be non-negative."
        return new(
            tuple_set,
            repeat_tuple_set,
            repeat_numbers,
            repeat_indices,
        )::TupleSetData
    end
end

struct Design
    # The code circuit the design aims to characterise
    code::Code
    # Whether to generate the full covariance or just the terms on the diagonal
    full_covariance::Bool
    # The M x N design matrix (M circuit eigenvalues, N gate eigenvalues)
    matrix::SparseMatrixCSC{Int32, Int32}
    # Circuit rearrangements used to generate the design matrix
    tuple_set::Vector{Vector{Int}}
    # Data used to generate the tuple set
    tuple_set_data::TupleSetData
    # For each tuple in the set, a mapping set or vector of Mapping objects for each Pauli supported on some gate in the rearranged circuit, corresponding to the circuit eigenvalues
    mapping_ensemble::Vector{Vector{Mapping}}
    # For each tuple in the set, an experiment set or vector of experiments, which themselves contain indices of simultaneously preparable and measurable Mapping objects 
    experiment_ensemble::Vector{Vector{Vector{Int}}}
    # For each tuple in the set, a dictionary describing the non-zero entries of the covariance matrix of the circuit eigenvalue estimators produced by the corresponding experiment set
    covariance_dict_ensemble::Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}}
    # For each tuple in the set, a vector Layers implementing the appropriate preparations in all necessary sign configurations for each experiment in the corresponding experiment set
    prep_ensemble::Vector{Vector{Vector{Layer}}}
    # For each tuple in the set, a Layer implementing the appropriate measurement for each experiment in the corresponding experiment set
    meas_ensemble::Vector{Vector{Layer}}
    # Time taken to implement each tuple's circuit, normalised according to the trivial tuple set
    tuple_times::Vector{Float64}
    # Weighting of the shots allocated to each tuple
    shot_weights::Vector{Float64}
    # The number of experiments for each tuple in the set
    experiment_numbers::Vector{Int}
    # The total number of experiments
    experiment_number::Int
    # The time taken to generate components of the design for each tuple
    # (mapping_time, consistency_time, pauli_time, covariance_time, circuit_time)
    calculation_times::Matrix{Float64}
    # The overall time taken to generate the design
    overall_time::Float64
    # The time taken to optimise the design
    optimisation_time::Float64
    # Type of least squares for which the shot weights were optimised
    ls_type::Symbol
    # Default constructor
    function Design(
        code::Code,
        full_covariance::Bool,
        matrix::SparseMatrixCSC{Int32, Int32},
        tuple_set::Vector{Vector{Int}},
        tuple_set_data::TupleSetData,
        mapping_ensemble::Vector{Vector{Mapping}},
        experiment_ensemble::Vector{Vector{Vector{Int}}},
        covariance_dict_ensemble::Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}},
        prep_ensemble::Vector{Vector{Vector{Layer}}},
        meas_ensemble::Vector{Vector{Layer}},
        tuple_times::Vector{Float64},
        shot_weights::Vector{Float64},
        experiment_numbers::Vector{Int},
        experiment_number::Int,
        calculation_times::Matrix{Float64},
        overall_time::Float64,
        optimisation_time::Float64,
        ls_type::Symbol,
    )
        # Check parameters
        T = length(tuple_set)
        @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
        @assert tuple_set == TupleSet(tuple_set_data) "The tuple set doesn't align with the tuple set data."
        @assert length(mapping_ensemble) == T "The size of the mapping ensemble does not match the tuple set."
        @assert length(experiment_ensemble) == T "The size of the experiment ensemble does not match the tuple set."
        @assert length(covariance_dict_ensemble) == T "The size of the covariance dictionary ensemble does not match the tuple set."
        @assert length(prep_ensemble) == T "The size of the preparation ensemble does not match the tuple set."
        @assert length(meas_ensemble) == T "The size of the measurement ensemble does not match the tuple set."
        @assert length(tuple_times) == T "The number of tuple times does not match the tuple set."
        @assert length(shot_weights) == T "The number of shot weights does not match the tuple set."
        @assert sum(shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
        @assert all(shot_weights .> 0.0) "The shot weights are not all positive."
        @assert length(experiment_numbers) == T "The number of experiment numbers does not match the tuple set."
        @assert experiment_number == sum(experiment_numbers) "The experiment number $(experiment_number) does not match the sum of the experiment numbers $(sum(experiment_numbers))."
        @assert size(calculation_times) == (T, 5) "The calculation times do not match the tuple set."
        @assert ls_type ∈ [:none, :gls, :wls, :ols] "The least squares type $(ls_type) is not supported."
        # Return the design
        return new(
            code,
            full_covariance,
            matrix,
            tuple_set,
            tuple_set_data,
            mapping_ensemble,
            experiment_ensemble,
            covariance_dict_ensemble,
            prep_ensemble,
            meas_ensemble,
            tuple_times,
            shot_weights,
            experiment_numbers,
            experiment_number,
            calculation_times,
            overall_time,
            optimisation_time,
            ls_type,
        )::Design
    end
    # Constructor
    function Design(
        code::Code,
        full_covariance::Bool,
        matrix::SparseMatrixCSC{Int32, Int32},
        tuple_set::Vector{Vector{Int}},
        mapping_ensemble::Vector{Vector{Mapping}},
        experiment_ensemble::Vector{Vector{Vector{Int}}},
        covariance_dict_ensemble::Vector{Dict{CartesianIndex{2}, Tuple{Mapping, Int}}},
        prep_ensemble::Vector{Vector{Vector{Layer}}},
        meas_ensemble::Vector{Vector{Layer}},
        tuple_times::Vector{Float64},
        shot_weights::Vector{Float64},
        calculation_times::Matrix{Float64},
        overall_time::Float64,
    )
        # Initialise parameters
        T = length(tuple_set)
        tuple_set_data = TupleSetData(tuple_set, Vector{Int}[], Int[], Int[])
        experiment_numbers = length.([vcat(prep_set...) for prep_set in prep_ensemble])
        experiment_number = length(vcat(vcat(prep_ensemble...)...))
        optimisation_time = 0.0
        ls_type = :none
        # Check parameters
        @assert tuple_set == unique(tuple_set) "The tuple set contains repeated tuples."
        @assert length(mapping_ensemble) == T "The size of the mapping ensemble does not match the tuple set."
        @assert length(experiment_ensemble) == T "The size of the experiment ensemble does not match the tuple set."
        @assert length(covariance_dict_ensemble) == T "The size of the covariance dictionary ensemble does not match the tuple set."
        @assert length(prep_ensemble) == T "The size of the preparation ensemble does not match the tuple set."
        @assert length(meas_ensemble) == T "The size of the measurement ensemble does not match the tuple set."
        @assert length(tuple_times) == T "The number of tuple times does not match the tuple set."
        @assert length(shot_weights) == T "The number of shot weights does not match the tuple set."
        @assert sum(shot_weights) ≈ 1.0 "The shot weights are not appropriately normalised."
        @assert all(shot_weights .> 0.0) "The shot weights are not all positive."
        @assert length(experiment_numbers) == T "The number of experiment numbers does not match the tuple set."
        @assert experiment_number == sum(experiment_numbers) "The experiment number $(experiment_number) does not match the sum of the experiment numbers $(sum(experiment_numbers))."
        @assert size(calculation_times) == (T, 5) "The calculation times do not match the tuple set."
        @assert ls_type ∈ [:none, :gls, :wls, :ols] "The least squares type $(ls_type) is not supported."
        # Return the design
        return new(
            code,
            full_covariance,
            matrix,
            tuple_set,
            tuple_set_data,
            mapping_ensemble,
            experiment_ensemble,
            covariance_dict_ensemble,
            prep_ensemble,
            meas_ensemble,
            tuple_times,
            shot_weights,
            experiment_numbers,
            experiment_number,
            calculation_times,
            overall_time,
            optimisation_time,
            ls_type,
        )::Design
    end
end

function Base.show(io::IO, d::Design)
    return print(
        io,
        "Design for a $(d.code.code_param.code_name) code with $(length(d.tuple_set)) tuples and $(d.experiment_number) experiments.",
    )
end

struct Merit
    # Tuple set
    tuple_set::Vector{Vector{Int}}
    # Tuple set data
    tuple_set_data::TupleSetData
    # Code parameters
    code_param::AbstractCodeParameters
    # Noise parameters
    noise_param::AbstractNoiseParameters
    # The total number of gates
    G::Int
    # The number of gate eigenvalues
    N::Int
    # Type of least squares estimator used
    ls_type::Symbol
    # Expectation of the NRMSE
    expectation::Float64
    # Variance of the NRMSE
    variance::Float64
    # Eigenvalues of the gate log-eigenvalue estimator covariance matrix
    # These allows us to calculate the distribution of the NRMSE
    eigenvalues::Vector{Float64}
    # Condition number of the design matrix, the ratio of the largest and smallest singular values
    cond_num::Float64
    # Pseudoinverse norm of the design matrix, the inverse of the smallest singular value
    pinv_norm::Float64
    # Time taken to implement each tuple's circuit, normalised according to the trivial tuple set
    tuple_times::Vector{Float64}
    # Weighting of the shots allocated to each tuple
    shot_weights::Vector{Float64}
    # The number of experiments for each tuple in the set
    experiment_numbers::Vector{Int}
    # The total number of experiments
    experiment_number::Int
end

function Base.show(io::IO, merit::Merit)
    return show(io, round.((merit.expectation, sqrt(merit.variance)), digits = 5))
end

struct DepolarisingScalingData
    # Merit of the design for a range of code distances
    merit_scaling::Vector{Merit}
    # Code distances
    dist_range::Vector{Int}
    # Code parameters
    code_param::AbstractCodeParameters
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

function Base.show(io::IO, s::DepolarisingScalingData)
    return print(
        io,
        "Merit scaling data with depolarising noise of a design for a $(s.code_param.code_name) code with $(length(s.tuple_set)) tuples.",
    )
end

struct LogNormalScalingData
    # Expected NRMSE for a range of code distances
    expectation_scaling::Vector{Vector{Float64}}
    # NRMSE variance for a range of code distances
    variance_scaling::Vector{Vector{Float64}}
    # Eigenvalues of the gate log-eigenvalue estimator covariance matrix for a range of code distances
    eigenvalues_scaling::Vector{Vector{Vector{Float64}}}
    # Code distances
    dist_range::Vector{Int}
    # Code parameters
    code_param::AbstractCodeParameters
    # Log-normal random noise parameters
    noise_param::LogNormalParameters
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
    # (design_time, merit_time)
    calculation_times::Matrix{Float64}
    # The overall time taken to calculate the merit scaling for log-normal random noise
    overall_time::Float64
end

function Base.show(io::IO, s::LogNormalScalingData)
    return print(
        io,
        "Merit scaling data with log-normal random noise of a design for a $(s.code_param.code_name) code with $(length(s.tuple_set)) tuples.",
    )
end

struct ACESData
    # The design
    d::Design
    # The set of shots to repeatedly sample from the probability distribution
    shots_set::Vector{Int}
    # The set of shots to sample 
    shots_set_norm::Vector{Int}
    # The number of times to repeat the ACES estimation procedure
    repetitions::Int
    # Seeds for each of the repetitions
    seeds::Vector{UInt64}
    # Circuit eigenvalues
    eigenvalues::Vector{Float64}
    # Circuit eigenvalue estimator covariance matrix
    covariance::SparseMatrixCSC{Float64, Int32}
    # The estimated circuit eigenvalues for each of the shots in the set
    est_eigenvalues_coll::Matrix{Vector{Float64}}
    # The FGLS estimated gate eigenvalues for each of the shots in the set
    fgls_gate_eigenvalues_coll::Matrix{Vector{Float64}}
    # The GLS estimated gate eigenvalues for each of the shots in the set
    # This uses the true covariance matrix
    gls_gate_eigenvalues_coll::Matrix{Vector{Float64}}
    # The WLS estimated gate eigenvalues for each of the shots in the set
    wls_gate_eigenvalues_coll::Matrix{Vector{Float64}}
    # The OLS estimated gate eigenvalues for each of the shots in the set
    ols_gate_eigenvalues_coll::Matrix{Vector{Float64}}
    # The FGLS estimated gate probability distributions for each of the shots in the set
    fgls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}
    # The GLS estimated gate probability distributions for each of the shots in the set
    gls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}
    # The WLS estimated gate probability distributions for each of the shots in the set
    wls_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}
    # The OLS estimated gate probability distributions for each of the shots in the set
    ols_gate_probabilities_coll::Matrix{Dict{Gate, Vector{Float64}}}
    # The 2-norm between the FGLS estimated gate eigenvalues and the synthetic gate eigenvalues for each of the shots in the set
    fgls_gate_norm_coll::Matrix{Float64}
    # The 2-norm between the GLS estimated gate eigenvalues and the synthetic gate eigenvalues for each of the shots in the set
    gls_gate_norm_coll::Matrix{Float64}
    # The 2-norm between the WLS estimated gate eigenvalues and the synthetic gate eigenvalues for each of the shots in the set
    wls_gate_norm_coll::Matrix{Float64}
    # The 2-norm between the OLS estimated gate eigenvalues and the synthetic gate eigenvalues for each of the shots in the set
    ols_gate_norm_coll::Matrix{Float64}
    # The time taken to simulate sampling and estimate the gate eigenvalues for each repetition
    # (simulate_time, fgls_time, gls_time, wls_time, ols_time)
    calculation_times::Matrix{Float64}
    # The overall time taken to simulate ACES across all repetitions
    overall_time::Float64
end

function Base.show(io::IO, a::ACESData)
    return print(
        io,
        "ACES data from a design for a $(a.d.code.code_param.code_name) code with $(length(a.d.tuple_set)) tuples and $(a.d.experiment_number) experiments.",
    )
end

# Use StructEquality.jl to define hash, ==, and isequal for structs
@struct_hash_equal_isequal Tableau
@struct_hash_equal_isequal Gate
@struct_hash_equal_isequal Layer
@struct_hash_equal_isequal DepolarisingParameters
@struct_hash_equal_isequal LogNormalParameters
@struct_hash_equal_isequal RotatedPlanarParameters
@struct_hash_equal_isequal UnrotatedPlanarParameters
@struct_hash_equal_isequal Circuit
@struct_hash_equal_isequal Code
@struct_hash_equal_isequal Pauli
@struct_hash_equal_isequal Mapping
@struct_hash_equal_isequal TupleSetData
@struct_hash_equal_isequal Design
@struct_hash_equal_isequal Merit
@struct_hash_equal_isequal DepolarisingScalingData
@struct_hash_equal_isequal LogNormalScalingData
@struct_hash_equal_isequal ACESData

# ACES.jl
export AbstractNoiseParameters,
    AbstractCodeParameters,
    AbstractCircuit,
    Tableau,
    Gate,
    Layer,
    DepolarisingParameters,
    LogNormalParameters,
    RotatedPlanarParameters,
    UnrotatedPlanarParameters,
    Circuit,
    Code,
    Pauli,
    Mapping,
    TupleSetData,
    Design,
    Merit,
    DepolarisingScalingData,
    LogNormalScalingData,
    ACESData

# tableau.jl
export CX!,
    Hadamard!,
    Phase!,
    X!,
    Z!,
    Y!,
    CZ!,
    SQRT_ZZ!,
    SQRT_ZZ_DAG!,
    RowSum!,
    Measure!,
    Reset!,
    Apply!
include("tableau.jl")

# circuit.jl
export Layer, PadIdentity, Unwrap, Gates, Label, Index, ApplyTuple
include("circuit.jl")

# code.jl
export LayerTimes, RotatedPlanar, UnrotatedPlanar, GenerateCode
include("code.jl")

# noise.jl
export DepolarisingNoise, LogNormalNoise, GenerateNoise, GateEigenvalues, Update
include("noise.jl")

# tuples.jl
export TrivialTupleSet,
    TrivialExperimentNumbers, TrivialHarmMean, TupleSetParameters, TupleSetData, TupleSet
include("tuples.jl")

# design.jl
export PauliPreparationSet,
    PrepLayer,
    MeasLayer,
    Mapping,
    MappingSet,
    ConsistencySet,
    PackPaulis,
    CovarianceDict,
    PackCircuits,
    GenerateDesign,
    CompleteDesign
include("design.jl")

# merit.jl
export CalculateCovariance,
    SyntheticEigenvalues,
    MeritData,
    SparseCovarianceInv,
    GLSCovariance,
    WLSCovariance,
    OLSCovariance,
    LSCovariance,
    NRMSEMoments,
    GLSMoments,
    WLSMoments,
    OLSMoments,
    LSMoments,
    GLSMerit,
    WLSMerit,
    OLSMerit,
    LSMerit,
    MeritSet,
    ImhofIntegrand,
    nrmse_pdf
include("merit.jl")

# gradient.jl
export ShotFactor,
    ShotFactorInv,
    ShotGradient,
    NRMSEGradient,
    GLSGradient,
    GLSDescent,
    WLSGradient,
    WLSDescent,
    OLSGradient,
    OLSDescent,
    CompareDescents,
    OptimiseShotWeights
include("gradient.jl")

# optimise.jl
export OptimalExpectation,
    StepRepetitions,
    OptimiseRepetitions,
    SampleZipf,
    RandomTuple,
    Grow,
    Prune,
    GrowDesign,
    PruneDesign,
    OptimiseTupleSet,
    OptimiseDesign
include("optimise.jl")

# scaling.jl
export DepolarisingScaling, DepolarisingFits, LogNormalScaling, LogNormalFits
include("scaling.jl")

# simulate.jl
export StimCircuitString,
    StimSample,
    batch_shots,
    SampleEigenvalues,
    FGLSGateEigenvalues,
    GLSGateEigenvalues,
    WLSGateEigenvalues,
    OLSGateEigenvalues,
    EstimateGateProbabilities,
    SimulateACES
include("simulate.jl")

# utils.jl
export SimplexProject, Support, WHTMatrix, PauliString, MappingString, PrettyPrint
include("utils.jl")

# io.jl
export enter_folder,
    exit_folder,
    code_filename,
    noise_filename,
    tuples_filename,
    design_filename,
    dep_scaling_filename,
    log_scaling_filename,
    aces_data_filename,
    save_design,
    load_design,
    delete_design,
    save_scaling,
    load_scaling,
    delete_scaling,
    save_aces,
    load_aces,
    delete_aces
include("io.jl")

# kwargs.jl
include("kwargs.jl")

# IntelliSence for Julia VSCode is terrible, but this hacky trick makes it work properly
# It convinces the LSP that the following files are part of src
# Source: https://discourse.julialang.org/t/lsp-missing-reference-woes/98231/16
@static if false
    include("../test/runtests.jl")
    include("../test/merit_tests.jl")
    include("../test/design_tests.jl")
    include("../scalable_aces/rot_optimise.jl")
    include("../scalable_aces/rot_scaling.jl")
    include("../scalable_aces/rot_simulate.jl")
    include("../scalable_aces/rot_simulate_big.jl")
    include("../scalable_aces/unrot_optimise.jl")
    include("../scalable_aces/unrot_scaling.jl")
    include("../scalable_aces/unrot_simulate.jl")
    include("../scalable_aces/unrot_simulate_big.jl")
end

end

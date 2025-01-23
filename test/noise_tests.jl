using QuantumACES, Test
# Define the gate types
one_qubit_gates = get_one_qubit_gates()
two_qubit_gates = get_two_qubit_gates()
# Get the hardcoded gate orbits
orbit_indices_dict = get_orbit_indices_dict()
# Test single-qubit gates
@testset "Single qubit gate orbits" begin
    for type in one_qubit_gates
        if type ∈ ["M"; "R"]
            @test ~haskey(orbit_indices_dict, type)
        else
            gate = Gate(type, 0, [1])
            gate_orbits = QuantumACES.calc_gate_orbits(gate)
            orbit_indices = QuantumACES.get_orbit_indices(gate_orbits)
            @test orbit_indices_dict[type] == orbit_indices
        end
    end
end
# Test two-qubit gates
@testset "Two qubit gate orbits" begin
    for type in two_qubit_gates
        gate = Gate(type, 0, [1; 2])
        gate_orbits = QuantumACES.calc_gate_orbits(gate)
        orbit_indices = QuantumACES.get_orbit_indices(gate_orbits)
        @test orbit_indices_dict[type] == orbit_indices
    end
end
# Set up code
dist = 3
r_1 = 0.075 / 100
r_2 = 0.5 / 100
r_m = 2.0 / 100
total_std_log = sqrt(log(10 / 9))
seed = UInt(0)
rotated_param = get_rotated_param(dist)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
rotated_planar = get_circuit(rotated_param, log_param)
total_gates = rotated_planar.total_gates
gate_data = rotated_planar.gate_data
gate_probabilities = rotated_planar.gate_probabilities
gate_eigenvalues = rotated_planar.gate_eigenvalues
# Test conversion and marginalisation between gate eigenvalues and probabilities
@testset "Conversion and marginalisation" begin
    # Test conversion between gate probabilities and eigenvalues
    test_gate_eigenvalues = get_gate_eigenvalues(gate_probabilities, gate_data)
    @test test_gate_eigenvalues ≈ gate_eigenvalues
    test_gate_probabilities = get_gate_probabilities(gate_eigenvalues, gate_data)
    @test all(
        test_gate_probabilities[gate] ≈ gate_probabilities[gate] for gate in total_gates
    )
    # Generate marginal gate probabilities
    average_gate_probabilities =
        get_average_gate_probabilities(gate_probabilities, gate_data)
    full_average_gate_probabilities =
        get_full_average_gate_probabilities(gate_probabilities, gate_data)
    marginal_gate_probabilities =
        get_marginal_gate_probabilities(gate_probabilities, gate_data)
    relative_gate_probabilities =
        get_relative_gate_probabilities(gate_probabilities, gate_data)
    # Generate gate probabilities vectors
    gate_probabilities_vec =
        QuantumACES.get_gate_probabilities_vec(gate_probabilities, gate_data)
    marginal_gate_probabilities_vec = QuantumACES.get_marginal_gate_probabilities_vec(
        marginal_gate_probabilities,
        gate_data,
    )
    relative_gate_probabilities_vec = QuantumACES.get_relative_gate_probabilities_vec(
        relative_gate_probabilities,
        gate_data,
    )
    # Test conversion of gate probabilities from dictionary to vector
    test_gate_probabilities =
        QuantumACES.get_gate_probabilities_dict(gate_probabilities_vec, gate_data)
    test_marginal_gate_probabilities = QuantumACES.get_marginal_gate_probabilities_dict(
        marginal_gate_probabilities_vec,
        gate_data,
    )
    test_relative_gate_probabilities = QuantumACES.get_relative_gate_probabilities_dict(
        relative_gate_probabilities_vec,
        gate_data,
    )
    @test all(
        test_gate_probabilities[gate] ≈ gate_probabilities[gate] for gate in total_gates
    )
    @test all(
        test_marginal_gate_probabilities[gate] ≈ marginal_gate_probabilities[gate] for
        gate in total_gates
    )
    @test all(
        test_relative_gate_probabilities[gate] ≈ relative_gate_probabilities[gate] for
        gate in keys(relative_gate_probabilities)
    )
    # Generate marginal gate eigenvalues
    marginal_gate_eigenvalues = get_marginal_gate_eigenvalues(gate_eigenvalues, gate_data)
    relative_gate_eigenvalues = get_relative_gate_eigenvalues(gate_eigenvalues, gate_data)
    # Generate transform matrices
    marginal_transform = QuantumACES.get_marginal_transform(gate_data)
    relative_transform = QuantumACES.get_relative_transform(gate_data)
    pad_transform = QuantumACES.get_pad_transform(gate_data)
    pad_transform_inv = QuantumACES.get_pad_transform(gate_data; inverse = true)
    pad_mask = QuantumACES.get_pad_mask(gate_data)
    marginal_pad_transform = QuantumACES.get_marginal_pad_transform(gate_data)
    marginal_pad_transform_inv =
        QuantumACES.get_marginal_pad_transform(gate_data; inverse = true)
    marginal_pad_mask = QuantumACES.get_marginal_pad_mask(gate_data)
    relative_pad_transform = QuantumACES.get_relative_pad_transform(gate_data)
    relative_pad_transform_inv =
        QuantumACES.get_relative_pad_transform(gate_data; inverse = true)
    relative_pad_mask = QuantumACES.get_relative_pad_mask(gate_data)
    wht_transform = QuantumACES.get_wht_transform(gate_data)
    wht_transform_inv = QuantumACES.get_wht_transform(gate_data; inverse = true)
    marginal_wht_transform = QuantumACES.get_marginal_wht_transform(gate_data)
    marginal_wht_transform_inv =
        QuantumACES.get_marginal_wht_transform(gate_data; inverse = true)
    relative_wht_transform = QuantumACES.get_relative_wht_transform(gate_data)
    relative_wht_transform_inv =
        QuantumACES.get_relative_wht_transform(gate_data; inverse = true)
    # Test marginalisation of gate eigenvalues
    test_marginal_gate_eigenvalues = marginal_transform * gate_eigenvalues
    test_relative_gate_eigenvalues = relative_transform * gate_eigenvalues
    @test test_marginal_gate_eigenvalues ≈ marginal_gate_eigenvalues
    @test test_relative_gate_eigenvalues ≈ relative_gate_eigenvalues
    # Test transforming from gate eigenvalues to probabilities
    test_gate_probabilities_vec =
        wht_transform_inv * (pad_transform * gate_eigenvalues + pad_mask)
    test_marginal_gate_probabilities_vec =
        marginal_wht_transform_inv *
        (marginal_pad_transform * marginal_gate_eigenvalues + marginal_pad_mask)
    test_relative_gate_probabilities_vec =
        relative_wht_transform_inv *
        (relative_pad_transform * relative_gate_eigenvalues + relative_pad_mask)
    @test test_gate_probabilities_vec ≈ gate_probabilities_vec
    @test test_marginal_gate_probabilities_vec ≈ marginal_gate_probabilities_vec
    @test test_relative_gate_probabilities_vec ≈ relative_gate_probabilities_vec
    # Test transforming from gate probabilities to eigenvalues
    test_gate_eigenvalues =
        pad_transform_inv * (wht_transform * test_gate_probabilities_vec)
    test_marginal_gate_eigenvalues =
        marginal_pad_transform_inv *
        (marginal_wht_transform * test_marginal_gate_probabilities_vec)
    test_relative_gate_eigenvalues =
        relative_pad_transform_inv *
        (relative_wht_transform * test_relative_gate_probabilities_vec)
    @test test_gate_eigenvalues ≈ gate_eigenvalues
    @test test_marginal_gate_eigenvalues ≈ marginal_gate_eigenvalues
    @test test_relative_gate_eigenvalues ≈ relative_gate_eigenvalues
end

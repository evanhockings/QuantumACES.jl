using QuantumACES, Test
# Set up noise models
r_1 = 0.05 / 100
r_2 = 0.4 / 100
r_m = 0.8 / 100
total_std_log = 0.5
seed = UInt(0)
dep_param = get_dep_param(r_1, r_2, r_m)
log_param = get_log_param(r_1, r_2, r_m, total_std_log; seed = seed)
# Set up circuits
vertical_distances = [3; 4; 5]
horizontal_distances = [3; 4; 5]
time_multiplier = 100
primary_param_list = AbstractCircuitParameters[]
circuit_param_list = AbstractCircuitParameters[]
for vertical_dist in vertical_distances
    for horizontal_dist in horizontal_distances
        # Add the rotated planar circuits
        push!(primary_param_list, get_rotated_param(vertical_dist, horizontal_dist))
        push!(circuit_param_list, get_rotated_param(vertical_dist, horizontal_dist))
        push!(
            circuit_param_list,
            get_rotated_param(
                vertical_dist,
                horizontal_dist;
                dynamically_decouple = false,
                pad_identity = false,
                single_qubit_time = time_multiplier * rand(),
                two_qubit_time = time_multiplier * rand(),
                meas_reset_time = time_multiplier * rand(),
                mid_reset_time = time_multiplier * rand(),
            ),
        )
        push!(
            circuit_param_list,
            get_rotated_param(
                vertical_dist,
                horizontal_dist;
                gate_type = :cx,
                check_type = :standard,
                dynamically_decouple = false,
            ),
        )
        # Add the unrotated planar circuits
        push!(primary_param_list, get_unrotated_param(vertical_dist, horizontal_dist))
        push!(circuit_param_list, get_unrotated_param(vertical_dist, horizontal_dist))
        push!(
            circuit_param_list,
            get_unrotated_param(
                vertical_dist,
                horizontal_dist;
                pad_identity = false,
                single_qubit_time = time_multiplier * rand(),
                two_qubit_time = time_multiplier * rand(),
                meas_reset_time = time_multiplier * rand(),
                mid_reset_time = time_multiplier * rand(),
            ),
        )
        # Add the heavy hex circuits
        push!(primary_param_list, get_hex_param(vertical_dist, horizontal_dist))
        push!(circuit_param_list, get_hex_param(vertical_dist, horizontal_dist))
        push!(
            circuit_param_list,
            get_hex_param(
                vertical_dist,
                horizontal_dist;
                flipped = true,
                pad_identity = false,
                single_qubit_time = time_multiplier * rand(),
                two_qubit_time = time_multiplier * rand(),
                meas_reset_time = time_multiplier * rand(),
                mid_reset_time = time_multiplier * rand(),
            ),
        )
    end
end
# Test the memory circuits
shots = 256
reset_type_list = [:meas_reset; :meas]
@testset "Memory circuits" begin
    for circuit_param in circuit_param_list
        # Get parameters
        vertical_dist = circuit_param.params[:vertical_dist]
        horizontal_dist = circuit_param.params[:horizontal_dist]
        # Generate the circuits
        circuit_log = get_circuit(circuit_param, log_param)
        circuit_dep = get_circuit(circuit_param, dep_param)
        # Simulate the memory experiments
        decoder_gate_probabilities =
            [circuit_log.gate_probabilities, circuit_dep.gate_probabilities]
        for reset_type in reset_type_list
            if reset_type == :meas_reset
                rounds_list = [0; max(vertical_dist, horizontal_dist)]
            elseif reset_type == :meas
                rounds_list = [0; 1; max(vertical_dist, horizontal_dist)]
            else
                throw(error("Unsupported reset type $(reset_type)."))
            end
            # Simulate memory experiments and decoding
            for rounds in rounds_list
                if circuit_param == circuit_param_list[end] &&
                   reset_type == reset_type_list[end] &&
                   rounds == rounds_list[end]
                    memory_data = simulate_memory(
                        circuit_log,
                        rounds,
                        2 * shots;
                        seed = seed,
                        reset_type = reset_type,
                        decoder_gate_probabilities = decoder_gate_probabilities,
                        max_samples = shots,
                        diagnostics = true,
                    )
                else
                    memory_data = simulate_memory(
                        circuit_log,
                        rounds,
                        shots;
                        seed = seed,
                        reset_type = reset_type,
                        decoder_gate_probabilities = decoder_gate_probabilities,
                    )
                end
            end
            # Test the code distances
            (z_dist, x_dist) = calc_memory_distances(circuit_dep; reset_type = reset_type)
            @test vertical_dist == z_dist
            @test horizontal_dist == x_dist
        end
    end
end
# Test the primary circuits with SPAM noise in the preparations
@testset "Preparation noise" begin
    for circuit_param in primary_param_list
        # Get parameters
        vertical_dist = circuit_param.params[:vertical_dist]
        horizontal_dist = circuit_param.params[:horizontal_dist]
        # Generate the circuits
        circuit_log = get_circuit(
            circuit_param,
            log_param;
            noisy_prep = true,
            noisy_meas = false,
            combined = (vertical_dist == 5 || horizontal_dist == 5 ? true : false),
            strict = (vertical_dist == 3 || horizontal_dist == 3 ? true : false),
        )
        circuit_dep = get_circuit(
            circuit_param,
            dep_param;
            noisy_prep = true,
            noisy_meas = false,
            combined = (vertical_dist == 5 || horizontal_dist == 5 ? true : false),
            strict = (vertical_dist == 3 || horizontal_dist == 3 ? true : false),
        )
        # Simulate the memory experiments
        decoder_gate_probabilities =
            [circuit_log.gate_probabilities, circuit_dep.gate_probabilities]
        for reset_type in reset_type_list
            if reset_type == :meas_reset
                rounds_list = [0; max(vertical_dist, horizontal_dist)]
            elseif reset_type == :meas
                rounds_list = [0; 1; max(vertical_dist, horizontal_dist)]
            else
                throw(error("Unsupported reset type $(reset_type)."))
            end
            # Simulate memory experiments and decoding
            for rounds in rounds_list
                memory_data = simulate_memory(
                    circuit_log,
                    rounds,
                    shots;
                    seed = seed,
                    reset_type = reset_type,
                    decoder_type = (rounds == 0 ? :beliefmatching : :pymatching),
                    decoder_gate_probabilities = decoder_gate_probabilities,
                )
            end
            # Test the code distances
            (z_dist, x_dist) = calc_memory_distances(circuit_dep; reset_type = reset_type)
            @assert vertical_dist == z_dist
            @assert horizontal_dist == x_dist
        end
    end
end

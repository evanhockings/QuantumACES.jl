using QuantumACES, Test
# Run the tests
start_time = time()
@testset "QuantumACES.jl" begin
    @testset "Aqua" begin
        include("aqua_tests.jl")
    end
    println(
        "\nCompleted Aqua tests. Time elapsed: $(round(time() - start_time, digits = 1)) seconds.\n",
    )
    @testset "Noise tests" begin
        include("noise_tests.jl")
    end
    println(
        "\nCompleted noise tests. Time elapsed: $(round(time() - start_time, digits = 1)) seconds.\n",
    )
    @testset "Stim circuits" begin
        include("stim_tests.jl")
    end
    println(
        "\nCompleted Stim circuit tests. Time elapsed: $(round(time() - start_time, digits = 1)) seconds.\n",
    )
    @testset "Design merits and gradients" begin
        include("merit_tests.jl")
    end
    println(
        "\nCompleted design merit tests. Time elapsed: $(round(time() - start_time, digits = 1)) seconds.\n",
    )
    @testset "Optimising and simulating designs" begin
        include("design_tests.jl")
    end
    println(
        "\nCompleted design optimisation and simulation tests. Time elapsed: $(round(time() - start_time, digits = 1)) seconds.\n",
    )
    @testset "Creating circuits and noise models" begin
        include("creation_tests.jl")
    end
    println(
        "\nCompleted circuit and noise model creation tests. Time elapsed: $(round(time() - start_time, digits = 1)) seconds.\n",
    )
    @testset "Randomised compiling and Qiskit" begin
        include("device_tests.jl")
    end
    println(
        "\nCompleted randomised compiling and Qiskit tests. Time elapsed: $(round(time() - start_time, digits = 1)) seconds.\n",
    )
end

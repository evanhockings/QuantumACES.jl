using QuantumACES, Test

@testset "QuantumACES.jl" begin
    @testset "Aqua" begin
        include("aqua_tests.jl")
    end
    @testset "Design merits and gradients" begin
        include("merit_tests.jl")
    end
    @testset "Optimising and simulating designs" begin
        include("design_tests.jl")
    end
end

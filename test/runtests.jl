using AveragedCircuitEigenvalueSampling, Test

@testset "AveragedCircuitEigenvalueSampling.jl" begin
    @testset "Design merits and gradients" begin
        include("merit_tests.jl")
    end
    @testset "Optimising and simulating designs" begin
        include("design_tests.jl")
    end
end

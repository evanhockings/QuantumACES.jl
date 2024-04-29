using QuantumACES, Aqua
# Testing for ambiguities currently throws lots of warnings for dependencies, see:
# https://discourse.julialang.org/t/checking-ambiguities-in-aqua-jl/96908
Aqua.test_all(QuantumACES; ambiguities = false)

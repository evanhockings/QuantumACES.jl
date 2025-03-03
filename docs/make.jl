using QuantumACES, Documenter

DocMeta.setdocmeta!(QuantumACES, :DocTestSetup, :(using QuantumACES); recursive = true)

makedocs(;
    modules = [QuantumACES],
    authors = "Evan Hockings",
    sitename = "QuantumACES.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://evanhockings.github.io/QuantumACES.jl",
        edit_link = "main",
        assets = String[],
        size_threshold = 2^20,
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => Any["Guide" => "guide.md", "Performance" => "performance.md"],
        "Reference" => Any[
            "Public API" => "public.md",
            "Internal API" => Any[
                "QuantumACES.jl" => "internal/QuantumACES.md",
                "tableau.jl" => "internal/tableau.md",
                "noise.jl" => "internal/noise.md",
                "noises/depolarising.jl" => "internal/noises/depolarising.md",
                "noises/lognormal.jl" => "internal/noises/lognormal.md",
                "circuit.jl" => "internal/circuit.md",
                "circuits/rotated_planar.jl" => "internal/circuits/rotated_planar.md",
                "circuits/unrotated_planar.jl" => "internal/circuits/unrotated_planar.md",
                "circuits/heavy_hex.jl" => "internal/circuits/heavy_hex.md",
                "stim.jl" => "internal/stim.md",
                "tuples.jl" => "internal/tuples.md",
                "design.jl" => "internal/design.md",
                "rand_design.jl" => "internal/rand_design.md",
                "merit.jl" => "internal/merit.md",
                "optimise_weights.jl" => "internal/optimise_weights.md",
                "optimise_tuples.jl" => "internal/optimise_tuples.md",
                "estimate.jl" => "internal/estimate.md",
                "simulate.jl" => "internal/simulate.md",
                "device.jl" => "internal/device.md",
                "scaling.jl" => "internal/scaling.md",
                "kwargs.jl" => "internal/kwargs.md",
                "utils.jl" => "internal/utils.md",
                "io.jl" => "internal/io.md",
            ],
        ],
    ],
)

deploydocs(; repo = "github.com/evanhockings/QuantumACES.jl", devbranch = "main")

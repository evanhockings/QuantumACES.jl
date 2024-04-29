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
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Reference" => Any[
            "Public API" => "public.md",
            "Internal API" => Any[
                "QuantumACES.jl" => "internal/QuantumACES.md",
                "tableau.jl" => "internal/tableau.md",
                "noise.jl" => "internal/noise.md",
                "circuit.jl" => "internal/circuit.md",
                "tuples.jl" => "internal/tuples.md",
                "design.jl" => "internal/design.md",
                "merit.jl" => "internal/merit.md",
                "weights.jl" => "internal/weights.md",
                "optimise.jl" => "internal/optimise.md",
                "scaling.jl" => "internal/scaling.md",
                "simulate.jl" => "internal/simulate.md",
                "kwargs.jl" => "internal/kwargs.md",
                "utils.jl" => "internal/utils.md",
                "io.jl" => "internal/io.md",
            ],
        ],
    ],
)

deploydocs(; repo = "github.com/evanhockings/QuantumACES.jl", devbranch = "main")

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
            "Internal API" => map(
                s -> "internal/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/internal"))),
            ),
        ],
    ],
)

deploydocs(; repo = "github.com/evanhockings/QuantumACES.jl", devbranch = "main")

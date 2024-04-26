using Pkg
# pkg"activate .."
# push!(LOAD_PATH, "../src/")
script_dir = @__DIR__
Pkg.activate(script_dir)
parent_dir = dirname(script_dir)
Pkg.develop(PackageSpec(; path = parent_dir))
using Documenter, QuantumACES
#
DocMeta.setdocmeta!(QuantumACES, :DocTestSetup, :(using QuantumACES); recursive = true)

makedocs(;
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        # canonical = "https://evanhockings.github.io/QuantumACES.jl",
    ),
    modules = [QuantumACES],
    sitename = "QuantumACES.jl",
    authors = "Evan Hockings",
    repo = "https://github.com/evanhockings/QuantumACES.jl/blob/{commit}{path}#{line}",
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

#=

makedocs(;
    modules = [QuantumACES],
    authors = "Evan Hockings",
    repo = "https://github.com/evanhockings/QuantumACES.jl/blob/{commit}{path}#{line}",
    sitename = "QuantumACES.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://evanhockings.github.io/QuantumACES.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(;
    repo = "github.com/evanhockings/QuantumACES.jl",
    devbranch = "main",
)
=#

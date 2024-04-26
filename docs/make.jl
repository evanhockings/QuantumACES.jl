using Pkg
# pkg"activate .."
# push!(LOAD_PATH, "../src/")
script_dir = @__DIR__
Pkg.activate(script_dir)
parent_dir = dirname(script_dir)
Pkg.develop(PackageSpec(; path = parent_dir))
using Documenter, ACES
#
DocMeta.setdocmeta!(ACES, :DocTestSetup, :(using ACES); recursive = true)

makedocs(;
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        # canonical = "https://evanhockings.github.io/ACES.jl",
    ),
    modules = [ACES],
    sitename = "ACES.jl",
    authors = "Evan Hockings",
    repo = "https://github.com/evanhockings/ACES.jl/blob/{commit}{path}#{line}",
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
    modules = [ACES],
    authors = "Evan Hockings",
    repo = "https://github.com/evanhockings/ACES.jl/blob/{commit}{path}#{line}",
    sitename = "ACES.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://evanhockings.github.io/ACES.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(;
    repo = "github.com/evanhockings/ACES.jl",
    devbranch = "main",
)
=#

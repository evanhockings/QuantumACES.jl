using ACES
using Documenter

DocMeta.setdocmeta!(ACES, :DocTestSetup, :(using ACES); recursive = true)

makedocs(;
         modules = [ACES],
         authors = "Evan Hockings",
         repo = "https://github.com/EvanHockings/ACES.jl/blob/{commit}{path}#{line}",
         sitename = "ACES.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://EvanHockings.github.io/ACES.jl",
                                  assets = String[]),
         pages = [
             "Home" => "index.md",
         ])

deploydocs(;
           repo = "github.com/EvanHockings/ACES.jl",
           devbranch = "main")

using Documenter
using SpheroidalWaveFunctions

makedocs(
    sitename = "SpheroidalWaveFunctions.jl",
    modules = [SpheroidalWaveFunctions],
    authors = "SpheroidalWaveFunctions contributors",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Math and Usage" => "math-and-usage.md",
        "Backend Overrides" => "backend-overrides.md",
    ],
)

deploydocs(
    repo = "github.com/Brandyn/SpheroidalWaveFunctions.jl.git",
    devbranch = "master",
)

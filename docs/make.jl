using Documenter
using SpheroidalWaves

makedocs(
    sitename = "SpheroidalWaves.jl",
    modules = [SpheroidalWaves],
    authors = "SpheroidalWaves contributors",
    repo = "https://github.com/brandynlucca/SpheroidalWaves.jl/blob/{commit}{path}#{line}",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        repolink = "https://github.com/brandynlucca/SpheroidalWaves.jl",
        edit_link = "main",
        assets = ["assets/branding.css"],
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Math and Usage" => "math-and-usage.md",
        "Backend Overrides" => "backend-overrides.md",
    ],
)

deploydocs(
    repo = "github.com/brandynlucca/SpheroidalWaves.jl.git",
    devbranch = "main",
)


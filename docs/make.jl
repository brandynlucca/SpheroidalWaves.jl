using Documenter
using SpheroidalWaves

makedocs(
    sitename = "SpheroidalWaves.jl",
    modules = [SpheroidalWaves],
    checkdocs = :exports,
    authors = "SpheroidalWaves contributors",
    remotes = nothing,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        repolink = nothing,
        edit_link = nothing,
        assets = ["assets/branding.css"],
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Math and Usage" => "math-and-usage.md",
        "Mathematical Tools" => "mathematical-tools.md",
    ],
)

deploydocs(
    repo = "github.com/brandynlucca/SpheroidalWaves.jl.git",
    devbranch = "main",
)


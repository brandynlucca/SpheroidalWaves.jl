using Documenter
using SpheroidalWaves

makedocs(
    sitename = "SpheroidalWaves.jl",
    modules = [SpheroidalWaves],
    checkdocs = :exports,
    pagesonly = true,
    doctest = true,
    warnonly = false,
    authors = "SpheroidalWaves contributors",
    remotes = nothing,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = nothing,
        assets = ["assets/branding.css"],
        sidebar_sitename = false,
        footer = "Powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). " *
                 "[Logo attribution and license]("*
                 "https://github.com/brandynlucca/SpheroidalWaves.jl/blob/main/" *
                 "docs/src/assets/LICENSE).",
        repolink = "https://github.com/brandynlucca/SpheroidalWaves.jl"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Math and Usage" => "math-and-usage.md",
        "Mathematical Tools" => "mathematical-tools.md"
    ]
)

deploydocs(
    repo = "github.com/brandynlucca/SpheroidalWaves.jl.git",
    devbranch = "main"
)

using Documenter
using PoissonGrids

DocMeta.setdocmeta!(PoissonGrids, :DocTestSetup, :(using PoissonGrids); recursive = true)

makedocs(;
    sitename = "PoissonGrids.jl",
    modules = [PoissonGrids],
    checkdocs = :exports,
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/albert-de-montserrat/PoissonGrids.jl.git",
    branch = "gh-pages",
    devbranch = "main",
)

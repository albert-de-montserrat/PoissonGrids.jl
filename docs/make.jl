using Documenter
using PoissonGrids

makedocs(;
    sitename = "PoissonGrids.jl",
    modules = [PoissonGrids],
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

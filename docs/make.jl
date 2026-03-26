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

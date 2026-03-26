# PoissonGrids.jl

[![CI](https://github.com/albert-de-montserrat/PoissonGrids.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/albert-de-montserrat/PoissonGrids.jl/actions/workflows/ci.yml)
[![Docs](https://github.com/albert-de-montserrat/PoissonGrids.jl/actions/workflows/docs.yml/badge.svg)](https://albert-de-montserrat.github.io/PoissonGrids.jl/dev/)

`PoissonGrids.jl` generates one-dimensional adaptive grids from scalar monitor
functions.

## Quick Start

```julia
using PoissonGrids

M = gaussian_monitor(5.0, 0.0, 0.2)
xmin, xmax = -1.0, 1.0 # domain limits
nc = 32 # number of cells
solve_grid(xmin, xmax, M, nc; verbose = false)

```

The returned vector `u` contains the grid vertices.

![Gaussian refinement example](docs/src/assets/gaussian_refinement.png)

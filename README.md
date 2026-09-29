# PoissonGrids.jl

[![CI](https://github.com/albert-de-montserrat/PoissonGrids.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/albert-de-montserrat/PoissonGrids.jl/actions/workflows/ci.yml)
[![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://albert-de-montserrat.github.io/PoissonGrids.jl/stable/)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://albert-de-montserrat.github.io/PoissonGrids.jl/dev/)
[![codecov](https://codecov.io/gh/albert-de-montserrat/PoissonGrids.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/albert-de-montserrat/PoissonGrids.jl)
[![version](https://juliahub.com/docs/General/PoissonGrids/stable/version.svg)](https://juliahub.com/ui/Packages/General/PoissonGrids)

`PoissonGrids.jl` generates one-dimensional adaptive grids from scalar monitor
functions.

A *monitor function* `M(x) > 0` states where resolution is wanted: the solver
places grid vertices so that cell width is inversely proportional to `M`, so
regions where `M` is large end up finely resolved. Only the relative variation
of `M` matters — scaling it by a constant leaves the grid unchanged.

## Installation

`PoissonGrids.jl` is not yet registered in the General registry; install it
from the repository URL:

```julia
using Pkg
Pkg.add(url = "https://github.com/albert-de-montserrat/PoissonGrids.jl")
```

## Quick Start

```julia
using PoissonGrids

M = gaussian_monitor(5.0, 0.0, 0.2)
xmin, xmax = -1.0, 1.0 # domain limits
nc = 32 # number of cells
u = solve_grid(xmin, xmax, M, nc; verbose = false)
```

The returned vector `u` contains the `nc + 1` grid vertices, with the endpoints
fixed at `xmin` and `xmax`.

![Gaussian refinement example](docs/src/assets/gaussian_refinement.png)

## Monitor Functions

Three constructors cover the common refinement shapes. Any positive, finite
callable works too — the constructors are conveniences, not a required type.

```julia
# localized bump of amplitude α at xc, width σ
M = gaussian_monitor(5.0, 0.0, 0.2)
u = solve_grid(-1.0, 1.0, M, 32)

# one-sided interface: 1 → 1 + α across c, with sharpness κ
M = tanh_monitor(5.0, 20.0, 0.0)
u = solve_grid(-1.0, 1.0, M, 32)

# refined window of half-width b centered on c
M = window_monitor(2.0, 2.5, 3.0, 2.5)
u = solve_grid(-10.0, 10.0, M, 32)
```

`tanh_monitor` accepts `direction = :left` to refine toward smaller `x` instead.

## Convergence

`solve_grid` relaxes the grid until the normalized residual falls below `tol`,
and throws `ConvergenceError` if that has not happened within `maxiter`
iterations:

```julia
julia> solve_grid(-1.0, 1.0, gaussian_monitor(5.0, 0.0, 0.2), 32; maxiter = 1)
ERROR: solve_grid did not converge after 1 iterations (residual = 83.18357644249869, tolerance = 1.0e-7).
```

Sharper monitors need more iterations, so raise `maxiter` before loosening `tol`.

See the [documentation](https://albert-de-montserrat.github.io/PoissonGrids.jl/dev/)
for the grid equation, the equidistribution principle behind it, and worked
examples for each monitor.

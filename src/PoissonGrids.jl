module PoissonGrids

using StaticArrays
using ForwardDiff: ForwardDiff
import LinearAlgebra: norm

const derivative = Val(true)
const primitive = Val(false)

include("monitors.jl")
export gaussian_monitor, tanh_monitor, window_monitor

include("solver.jl")
export ConvergenceError, solve_grid

end # module PoissonGrids

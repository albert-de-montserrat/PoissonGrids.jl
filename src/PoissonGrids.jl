module PoissonGrids

using ForwardDiff, StaticArrays
import LinearAlgebra: norm

using DifferentiationInterface
using ForwardDiff: ForwardDiff

const derivative = Val(true)
const primitive = Val(false)

include("monitors.jl")
export gaussian_monitor, tanh_monitor

include("solver.jl")
export solve_grid

end # module PoissonGrids

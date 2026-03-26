"""
    gaussian_monitor(α, xc, σ)

Construct a Gaussian monitor function.

The returned callable is

    M(x) = 1 + α * exp(-((x - xc)^2) / σ^2)

and can be passed to [`solve_grid`](@ref) to concentrate grid points around `xc`.

# Arguments
- `α`: amplitude of the refinement bump. Larger values produce stronger clustering.
- `xc`: center of the refined region.
- `σ`: width of the refined region. Smaller values produce narrower clustering.

# Returns
- A scalar function `M(x)` that evaluates the monitor at position `x`.

# Example
```julia
M = gaussian_monitor(5.0, -20.0, 1.0)
u = solve_grid(-50.0, 50.0, M, 127)
```
"""
function gaussian_monitor(α, xc, σ)
    return x -> 1.0 + α * exp(-((x - xc)^2) / σ^2)
end

"""
    tanh_monitor(α, κ, c; direction = :right)

Construct a monotone hyperbolic-tangent monitor function.

This monitor transitions smoothly from `1` to `1 + α` around `c` and is useful
when refinement is desired primarily on one side of an interface.

# Arguments
- `α`: size of the monitor jump.
- `κ`: transition sharpness. Larger values produce a steeper change.
- `c`: center of the transition.

# Keyword Arguments
- `direction = :right`: accepted for API compatibility. In the current
  implementation the returned profile always increases toward larger `x`.

# Returns
- A scalar function `M(x)` that evaluates the monitor at position `x`.

# Example
```julia
M = tanh_monitor(5.0, 1e-1, 0.0; direction = :right)
u = solve_grid(-50.0, 50.0, M, 127)
```
"""
function tanh_monitor(α, κ, c; direction = :right)
    s = direction === :right ? 1 : -1
    x -> (α + α * tanh(κ * (x - c))) / 2 + 1
end

######################## Main ########################
"""
    solve_grid(xmin, xmax, M, nc; verbose = false, tol = 1e-7, maxiter = 20_000)

Solve the 1D grid-generation problem on `[xmin, xmax]` for a prescribed monitor.

The solver starts from a uniform grid in computational space and iteratively moves
the interior vertices so that the physical grid adapts to the monitor function
`M(x)`. Regions where `M` is larger receive more grid points.

# Arguments
- `xmin`, `xmax`: physical domain bounds.
- `M`: scalar monitor function evaluated on the current physical grid.
- `nc`: number of cells in the final grid. The returned vertex array has length
  `nc + 1`.

# Keyword Arguments
- `verbose = false`: print residual information every few nonlinear iterations.
- `tol = 1e-7`: stopping tolerance for the normalized residual.
- `maxiter = 20_000`: maximum number of nonlinear iterations. If convergence is
  not reached within this budget, the solver throws [`ConvergenceError`](@ref).

# Returns
- A vector of grid-vertex coordinates spanning `[xmin, xmax]`.

# Example
```julia
M = gaussian_monitor(10.0, -20.0, 1.0)
u = solve_grid(-50.0, 50.0, M, 127; verbose = true, tol = 1e-8)
```
"""
struct ConvergenceError <: Exception
    maxiter::Int
    residual::Float64
    tolerance::Float64
end

function Base.showerror(io::IO, err::ConvergenceError)
    print(
        io,
        "solve_grid did not converge after $(err.maxiter) iterations ",
        "(residual = $(err.residual), tolerance = $(err.tolerance)).",
    )
end

function solve_grid(xmin, xmax, M, nc; verbose = false, tol = 1.0e-7, maxiter = 20_000)
    validate_solver_inputs(xmin, xmax, nc, tol, maxiter)

    # Resolution
    nv  = nc + 1  # vertices 

    # Computational grid (ξ, η) ∈ [0, 1] × [0, 1]
    ξmin = 0.0
    ξmax = 1.0
    Δξ   = (ξmax - ξmin) / nc

    # initial solution guess (uniform grid)
    u    = collect(LinRange(xmin, xmax, nv))

    validate_monitor_on_grid(M, u)

    # Memory allocations
    r    = zeros(nv)
    ∂u∂τ = zeros(nv)
    r0   = zeros(nv)
    b    = zeros(nv)
    D    = ones(nv)
    G    = ones(nv)
    q    = (;x = zeros(nc), )
    Gq   = (;x = zeros(nc), )

    # Step 1: compute flux Gershgorin first Gq = ∂q∂u
    run_flux!(Gq, u, M, Δξ, derivative)
    # Step 2: compute Poisson Gershgorin first Gq = ∂r∂q*∂q∂u + ∂r∂b*∂b∂u
    Poisson!(G, q, Gq, b, Δξ, derivative)

    # Approximate diagonal PC (does not see BC effect)
    @. D = G / 2

    # Iteration parameters
    CFL    = 0.99
    c_fact = 0.9
    nr0    = 1e0
    n      = √(prod(nc))
    converged = false
    residual = Inf

    λmax = maximum(G ./ D)
    Δτ   = 2 / √(λmax) * CFL
    λmin = 0.0
    c    = 2 * √(λmin) * c_fact
    α    = 2 * Δτ^2 / (2 + c * Δτ)
    β    = (2 - c * Δτ) / (2 + c * Δτ)

    # Iterations
    for iter in 1:maxiter
        copyto!(r0, r)
    
        #-------------------------
        # 1 - flux
        run_flux!(q, u, M, Δξ, primitive)
        # 2 - balance
        Poisson!(r, q, Gq, b, Δξ, primitive)
        #-------------------------
        update_rate!(∂u∂τ, r, D, β)
        update_variable!(u, ∂u∂τ, α)
        validate_monitor_on_grid(M, u)

        if iter == 1 || mod(iter, 5) == 0
            residual = norm(r) / √n
            nr0 = iter == 1 ? residual : nr0
            verbose && @info "iter. $(iter) --- abs. |r| =  $(norm(r)) ---  nr = $(residual / nr0)"
            
            if residual < tol
                converged = true
                break
            end

            # Step 1: compute flux Gershgorin first Gq = ∂q∂u
            run_flux!(Gq, u, M, Δξ, derivative)
            # Step 2: compute Poisson Gershgorin first Gq = ∂r∂q*∂q∂u + ∂r∂b*∂b∂u
            Poisson!(G, q, Gq, b, Δξ, derivative)

            # Approximate diagonal PC (does not see BC effect)
            @. D = G / 2
            
            λmax = maximum(G ./ D)
            Δτ   = 2 / √(maximum(λmax)) * CFL
            λmin = abs((sum(Δτ .* ∂u∂τ .* ((r .- r0) ./ D)))) / sum((Δτ .* ∂u∂τ) .^ 2)
            c    = 2 * √(λmin) * c_fact
            α    = 2 * Δτ^2 / (2 + c * Δτ)
            β    = (2 - c * Δτ) / (2 + c * Δτ)
        end
    end

    if !converged
        residual = isfinite(residual) ? residual : norm(r) / √n
        throw(ConvergenceError(maxiter, residual, float(tol)))
    end

    return u
end

function validate_solver_inputs(xmin, xmax, nc, tol, maxiter)
    xmin isa Real && isfinite(xmin) || throw(ArgumentError("xmin must be a finite real number."))
    xmax isa Real && isfinite(xmax) || throw(ArgumentError("xmax must be a finite real number."))
    xmin < xmax || throw(ArgumentError("xmin must be strictly smaller than xmax."))
    nc isa Integer || throw(ArgumentError("nc must be an integer number of cells."))
    nc > 0 || throw(ArgumentError("nc must be positive."))
    tol isa Real && isfinite(tol) || throw(ArgumentError("tol must be a finite real number."))
    tol > 0 || throw(ArgumentError("tol must be positive."))
    maxiter isa Integer || throw(ArgumentError("maxiter must be an integer."))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive."))
    return nothing
end

function validate_monitor_on_grid(M, u)
    @inbounds for i in 1:length(u)-1
        xmid = (u[i] + u[i + 1]) / 2
        mx = M(xmid)
        mx isa Real && isfinite(mx) || throw(ArgumentError("Monitor values must be finite real numbers."))
        mx > 0 || throw(ArgumentError("Monitor values must stay strictly positive on the grid."))
    end
    return nothing
end


######################## Physics ########################
# Poisson wrapper
function Poisson!(u, q, Gq, b, Δξ, ::Val{AD}) where {AD}
    for I in eachindex(u)[1:end-2]
        # Flux
        𝑞 = get_fluxes(q, I)
        # RHS
        𝑏 = b[I]
        # Call Poisson evaluation
        u[I+1] = if AD
            ∂𝑞∂𝑢 = get_fluxes(Gq, I)
            ∂r∂𝑞 = ForwardDiff.gradient(x -> Poisson_local(x, 𝑏, Δξ), 𝑞)
            sum(abs.(∂r∂𝑞 .* ∂𝑞∂𝑢))
        else
            Poisson_local(𝑞, 𝑏, Δξ)
        end
    end
    return nothing
end

@inline function get_fluxes(q, I)
    qxW = q.x[I]
    qxE = q.x[I + 1]
    return SA[qxW, qxE]
end

# Local Poisson function
# divergence of the flux on each face + source term b
@inline function Poisson_local(q::SVector{N, T}, b, Δξ) where {N, T}
    # Flux components on each face
    qleft  = q[1]
    qright = q[2]
    # Flux divergence
    ∇q = (qright - qleft) / Δξ
    # Residual
    r = -∇q - b
    return r
end

# Flux component

function run_flux!(q, u, M, Δξ, primitive)
    for I in eachindex(q.x)
        @inbounds q.x[I] = qi(u, M, Δξ, primitive, I) 
    end
end

# Flux component wrapper
@inline function qi(u, M, Δξ, ::Val{AD}, I) where {AD}
    # Neighbours
    uW = u[I]
    uE = u[I + 1]
    # Call flux component evaluation
    if AD
        𝑢 = @SVector [uW, uE]
        ∂q∂𝑢 = ForwardDiff.gradient(x -> qi_local(x, M, Δξ), 𝑢)
        return sum(abs.(∂q∂𝑢))
    else
        𝑢 = @SVector [uW, uE]
        return qi_local(𝑢, M, Δξ)
    end
end

@inline function qi_local(u, M, Δξ)
    # Flux component
    uW, uE = u
    Mx = M((uW + uE) / 2) # monitoring function evaluated at flux node
    return -Mx * (uE - uW) / Δξ
end

@inline function update_rate!(∂u∂τ, r, D, β)
    for I in eachindex(r)
        @inbounds ∂u∂τ[I] = r[I] / D[I] + β * ∂u∂τ[I]
    end
end

@inline function update_variable!(u, ∂u∂τ, α)
    for I in eachindex(u)[2:end-1]
        @inbounds u[I] = u[I] + α * ∂u∂τ[I]
    end
end

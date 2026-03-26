# using Metal
using ForwardDiff, StaticArrays
import LinearAlgebra: norm

using DifferentiationInterface
using ForwardDiff: ForwardDiff
using MuladdMacro

using GLMakie

const DAT = Float64
const derivative = Val(true)
const primitive = Val(false)

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

######################## Main ########################
function main(xmin, xmax, M, nc; verbose = false)

    # Resolution
    nv  = nc + 1  # vertices 

    # Computational grid (ξ, η) ∈ [0, 1] × [0, 1]
    ξmin = 0.0
    ξmax = 1.0
    Δξ   = (ξmax - ξmin) / ncx

    # initial solution guess (uniform grid)
    u    = collect(LinRange(xmin, xmax, nv))

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
    tol    = 1.0e-7
    niter  = 20000
    CFL    = 0.99
    c_fact = 0.9
    nr0    = 1e0
    n      = √(prod(nc))

    λmax = maximum(G ./ D)
    Δτ   = 2 / √(λmax) * CFL
    λmin = 0.0
    c    = 2 * √(λmin) * c_fact
    α    = 2 * Δτ^2 / (2 + c * Δτ)
    β    = (2 - c * Δτ) / (2 + c * Δτ)

    # Iterations
    @time for iter in 1:niter
        copyto!(r0, r)
    
        #-------------------------
        # 1 - flux
        run_flux!(q, u, M, Δξ, primitive)
        # 2 - balance
        Poisson!(r, q, Gq, b, Δξ, primitive)
        #-------------------------
        update_rate!(∂u∂τ, r, D, β)
        update_variable!(u, ∂u∂τ, α)

        if iter == 1 || mod(iter, 5) == 0
            nr = norm(r) / √n
            nr0 = iter == 1 ? nr : nr0
            verbose && @info "iter. $(iter) --- abs. |r| =  $(norm(r)) ---  nr = $(nr / nr0)"
            
            nr < tol && break

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

    return u
end

function gaussian_monitor(α, xc, σ)
    return x -> 1.0 + α * exp(-((x - xc)^2) / σ^2)
end

function tanh_monitor(α, κ, c; direction = :right)
    s = direction === :right ? 1 : -1
    x -> (α + α * tanh(κ * (x - c))) / 2 + 1
end

M = tanh_monitor(5, 1e-1 , 0e0)
# scatter(x, M.(x))
# number of centers
nc = 127 

# Physical grid
xmin, xmax = -50.0,  50.0 # kms

# # Monitor function
# xc_ref = -20 # km
# α      = 10e0     # Larger α = stronger refinement
# σ      = 1e0     # Smaller σ = narrower refined region
# M      = gaussian_monitor(α, xc_ref, σ)

u = main(xmin, xmax, M, nc; verbose = true);

fig = Figure(size=(1200, 1200), fontsize=20)
ax1 = Axis(fig[1,1]; xlabel = "x", ylabel = "Δx")
ax2 = Axis(fig[2,1]; xlabel = "x")
scatterlines!(ax1, u[1:end-1], diff(u))
scatterlines!(ax2, u, zero(u))
fig

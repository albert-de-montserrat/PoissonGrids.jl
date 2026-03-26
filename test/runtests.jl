using Test
using PoissonGrids

@testset "Monitor constructors" begin
    @testset "gaussian_monitor" begin
        α, xc, σ = 5.0, -2.0, 0.5
        M = gaussian_monitor(α, xc, σ)

        @test M(xc) == 1.0 + α
        @test M(xc - 0.25) ≈ M(xc + 0.25)
        @test M(xc + 5σ) ≈ 1.0 atol = 1e-10
    end

    @testset "tanh_monitor" begin
        α, κ, c = 4.0, 2.0, 0.0
        M = tanh_monitor(α, κ, c)
        Mleft = tanh_monitor(α, κ, c; direction = :left)

        @test M(c) ≈ 1.0 + α / 2
        @test M(-10.0) ≈ 1.0 atol = 1e-8
        @test M(10.0) ≈ 1.0 + α atol = 1e-8

        # `direction` is currently accepted but does not alter the profile.
        @test Mleft(-1.0) ≈ M(-1.0)
        @test Mleft(1.0) ≈ M(1.0)
    end
end

@testset "solve_grid" begin
    @testset "uniform monitor reproduces a uniform grid" begin
        xmin, xmax, nc = -1.0, 1.0, 8
        u = solve_grid(xmin, xmax, x -> 1.0, nc)

        @test length(u) == nc + 1
        @test first(u) == xmin
        @test last(u) == xmax
        @test u ≈ collect(LinRange(xmin, xmax, nc + 1))
    end

    @testset "adaptive monitor preserves ordering and endpoints" begin
        xmin, xmax, nc = -1.0, 1.0, 32
        M = gaussian_monitor(5.0, 0.0, 0.2)
        u = solve_grid(xmin, xmax, M, nc)
        Δx = diff(u)

        @test length(u) == nc + 1
        @test first(u) == xmin
        @test last(u) == xmax
        @test issorted(u)
        @test all(>(0.0), Δx)
        @test minimum(Δx) < maximum(Δx)
        @test minimum(Δx) < (xmax - xmin) / nc
    end
end

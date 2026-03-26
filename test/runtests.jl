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

        @test Mleft(c) ≈ 1.0 + α / 2
        @test Mleft(-10.0) ≈ 1.0 + α atol = 1e-8
        @test Mleft(10.0) ≈ 1.0 atol = 1e-8
    end

    @testset "window_monitor" begin
        α, κ, b, c = 4.0, 8.0, 0.5, 1.25
        M = window_monitor(α, κ, b, c)

        @test M(c) ≈ 1.0 + α atol = 1e-2
        @test M(-10.0) ≈ 1.0 atol = 1e-8
        @test M(10.0) ≈ 1.0 atol = 1e-8
        @test M(c - 0.25) ≈ M(c + 0.25)
        @test M(c) > M(c + b + 0.5)
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

    @testset "rejects invalid inputs" begin
        M = x -> 1.0

        @test_throws ArgumentError solve_grid(0.0, 0.0, M, 8)
        @test_throws ArgumentError solve_grid(1.0, -1.0, M, 8)
        @test_throws ArgumentError solve_grid(-1.0, 1.0, M, 0)
        @test_throws ArgumentError solve_grid(-1.0, 1.0, M, 8; tol = 0.0)
        @test_throws ArgumentError solve_grid(-1.0, 1.0, M, 8; maxiter = 0)
    end

    @testset "rejects invalid monitor values" begin
        @test_throws ArgumentError solve_grid(-1.0, 1.0, x -> 0.0, 8)
        @test_throws ArgumentError solve_grid(-1.0, 1.0, x -> NaN, 8)
    end

    @testset "throws on non-convergence" begin
        M = gaussian_monitor(5.0, 0.0, 0.2)
        err = try
            solve_grid(-1.0, 1.0, M, 32; maxiter = 1)
            nothing
        catch exc
            exc
        end
        @test err isa ConvergenceError
        @test err.maxiter == 1
        @test err.tolerance == 1.0e-7
        @test isfinite(err.residual)
    end
end

using ExportNURBS
using Test
using Random
using BasicBSpline
using Colors

@testset "ExportNURBS.jl" begin
    Random.seed!(42)

    @testset "luxor" begin
        p = 2 # degree of polynomial
        k = Knots(1:8) # knot vector
        P = FastBSplineSpace(p, k) # B-spline space
        rand_a = [rand(2) for i in 1:dim(P), j in 1:dim(P)]
        a = [[2 * i - 6.5, 2 * j - 6.5] for i in 1:dim(P), j in 1:dim(P)] + rand_a # random generated control points
        M = BSplineSurface([P, P], a) # Define B-spline manifold

        @testset "svg" begin
            save_svg("test.svg",M)
            @test isfile("test.svg")
        end
        @testset "png" begin
            save_png("test.png",M)
            @test isfile("test.png")
        end
        @testset "color-png" begin
            color(u) = rand(RGB)
            save_png("test_color.png", M, color)
            @test isfile("test_color.png")
        end
    end

    @testset "povray" begin
        p = 2 # degree of polynomial
        k = Knots(1:8) # knot vector
        P = FastBSplineSpace(p, k) # B-spline space
        rand_a = [rand(3) for i in 1:dim(P), j in 1:dim(P)]
        a = [[2 * i - 6.5, 2 * j - 6.5,-0.5] for i in 1:dim(P), j in 1:dim(P)] + rand_a # random generated control points
        M = BSplineSurface([P, P], a) # Define B-spline manifold

        @testset "2d3d" begin
            save_pov("test.inc",M)
            @test isfile("test.inc")
        end
    end
end

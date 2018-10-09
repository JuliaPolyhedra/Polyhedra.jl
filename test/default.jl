@testset "Default" begin
    @testset "Library" begin
        function default_library_test(d)
            @test default_library(d(2), Float64) isa Polyhedra.DefaultLibrary
            @test default_library(d(1), Float64) isa IntervalLibrary
            @test default_library(d(1), Int) isa IntervalLibrary
            @test default_library(d(1), BigInt) isa IntervalLibrary
            @test default_library(d(2), AbstractFloat) isa Polyhedra.DefaultLibrary{Float64}
            @test default_library(d(1), AbstractFloat) isa IntervalLibrary{Float64}
        end
        default_library_test(N -> StaticArrays.Size((N,)))
        default_library_test(identity)
    end
    @testset "LP Solver" begin
        @test Polyhedra.default_solver(hrep([HalfSpace((@SVector [1, 2]), 3)])) === nothing
        @test Polyhedra.solver(hrep([HalfSpace((@SVector [1, 2]), 3)]), lp_solver) === lp_solver
    end
end

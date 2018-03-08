@testset "Default" begin
    @testset "Library" begin
        @test default_library(Polyhedra.FullDim{2}(), Float64) isa SimplePolyhedraLibrary
        @test default_library(Polyhedra.FullDim{1}(), Float64) isa IntervalLibrary
        @test default_library(Polyhedra.FullDim{1}(), Float64) isa IntervalLibrary
    end
    @testset "LP Solver" begin
        @test defaultLPsolverfor(hrep([HalfSpace((@SVector [1, 2]), 3)]), lpsolver) === lpsolver
    end
end

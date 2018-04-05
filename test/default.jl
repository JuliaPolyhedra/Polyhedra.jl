@testset "Default" begin
    @testset "Library" begin
        @test default_library(Polyhedra.FullDim{2}(), Float64) isa SimplePolyhedraLibrary
        @test default_library(Polyhedra.FullDim{1}(), Float64) isa IntervalLibrary
        @test default_library(Polyhedra.FullDim{1}(), Float64) isa IntervalLibrary
        @test default_library(Polyhedra.FullDim{1}(), Float64) isa IntervalLibrary
        @test default_library(Polyhedra.FullDim{2}(), AbstractFloat) isa SimplePolyhedraLibrary{Float64}
        @test default_library(Polyhedra.FullDim{1}(), AbstractFloat) isa IntervalLibrary{Float64}
    end
    @testset "LP Solver" begin
        @test Polyhedra.default_solver(hrep([HalfSpace((@SVector [1, 2]), 3)])) isa JuMP.UnsetSolver
        @test Polyhedra.solver(hrep([HalfSpace((@SVector [1, 2]), 3)]), lp_solver) === lp_solver
    end
end

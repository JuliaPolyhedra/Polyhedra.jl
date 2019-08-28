using Test
using JuMP
using Polyhedra

@testset "Model to Polyhedra with default library" begin
    @testset "1 variable" begin
        model = Model()
        @variable(model, x >= -1, Int)
        @constraint(model, x ≤ 1)
        p = polyhedron(model)
        @test p isa Interval{Float64}
        @test nhalfspaces(p) == 2
    end
    @testset "2 variables" begin
        model = Model()
        @variable(model, x[1:2] >= -1, integer = true)
        @constraint(model, sum(x) ≤ 1)
        p = polyhedron(model)
        @test p isa DefaultPolyhedron{Float64}
        @test nhalfspaces(p) == 3
    end
    # TODO add test with binary variables once https://github.com/JuliaOpt/MathOptInterface.jl/issues/704 is done.
end

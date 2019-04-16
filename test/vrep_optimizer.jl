using Test
using Polyhedra
using JuMP
const MOIT = MOI.Test
const MOIB = MOI.Bridges

@testset "Continuous Linear problems with VRepOptimizer" begin
    optimizer = VRepOptimizer{Float64}()
    @testset "SolverName" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "VRep"
    end
    cache = MOIU.UniversalFallback(Polyhedra._MOIModel{Float64}())
    cached = MOIU.CachingOptimizer(cache, optimizer)
    bridged = MOIB.full_bridge_optimizer(cached, Float64)
    config = MOIT.TestConfig(duals=false)
    MOIT.contlineartest(bridged, config,
                        # linear8a and linear12 will be solved by https://github.com/JuliaOpt/MathOptInterface.jl/pull/702
                        ["linear8a", "linear12", "partial_start"])
end

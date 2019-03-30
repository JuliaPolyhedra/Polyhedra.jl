using Test
using Polyhedra
import MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

@testset "Continuous Linear problems with VRepOptimizer" begin
    optimizer = VRepOptimizer{Float64}()
    cache = MOIU.UniversalFallback(Polyhedra._MOIModel{Float64}())
    cached = MOIU.CachingOptimizer(cache, optimizer)
    bridged = MOIB.full_bridge_optimizer(cached, Float64)
    config = MOIT.TestConfig(duals=false)
    MOIT.contlineartest(bridged, config, ["partial_start"])
end

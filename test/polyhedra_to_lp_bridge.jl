using Test
using JuMP
using Polyhedra
@testset "PolyhedraToLPBridge{$T}" for T in [Float64, Int, BigInt, Rational{BigInt}]
    h = HalfSpace(ones(T, 2), one(T)) ∩ HalfSpace(-ones(T, 2), -one(T))
    set = Polyhedra.PolyhedraOptSet(h)
    S = typeof(set)
    Fs = [MOI.VectorOfVariables, MOI.VectorAffineFunction{T}, MOI.VectorQuadraticFunction{T}]
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{T}())
    for F in Fs
        @test !MOI.supports_constraint(mock, F, S)
    end

    @testset "LazyBridgeOptimizer" begin
        bridged = MOI.Bridges.LazyBridgeOptimizer(mock)
        for F in Fs
            @test !MOI.supports_constraint(bridged, F, S)
        end
        MOI.Bridges.add_bridge(bridged, Polyhedra.PolyhedraToLPBridge{T})
        for F in Fs
            @test MOI.supports_constraint(bridged, F, S)
        end
        x = MOI.add_constrained_variables(bridged, set)
     end

    @testset "LazyBridgeOptimizer" begin
        bridged = MOI.Bridges.full_bridge_optimizer(mock, T)
        for F in Fs
            @test !MOI.supports_constraint(bridged, F, S)
        end
        MOI.Bridges.add_bridge(bridged, Polyhedra.PolyhedraToLPBridge{T})
        for F in Fs
            @test MOI.supports_constraint(bridged, F, S)
        end
        x = MOI.add_constrained_variables(bridged, set)
     end
end

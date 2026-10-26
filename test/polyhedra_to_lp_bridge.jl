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

# See https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/306
@testset "ConstraintConflictStatus" begin
    import HiGHS
    p = polyhedron(HalfSpace([1.0, 1.0], 1.0) ∩ HalfSpace([-1.0, 0.0], 0.0) ∩ HalfSpace([0.0, -1.0], 0.0))
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[1:2])
    c = @constraint(model, x in p)
    d = @constraint(model, x[1] >= 2)
    e = @constraint(model, x[2] <= 10)
    optimize!(model)
    @test termination_status(model) == MOI.INFEASIBLE
    compute_conflict!(model)
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), d) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), e) == MOI.NOT_IN_CONFLICT
end

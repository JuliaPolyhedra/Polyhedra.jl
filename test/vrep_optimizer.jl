using Test
using Polyhedra
using JuMP
const MOIB = MOI.Bridges

@testset "Continuous Linear problems with VRepOptimizer" begin
    optimizer = VRepOptimizer{Float64}()
    @testset "SolverName" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "VRep"
    end
    cache = MOIU.UniversalFallback(Polyhedra._MOIModel{Float64}())
    cached = MOIU.CachingOptimizer(cache, optimizer)
    bridged = MOIB.full_bridge_optimizer(cached, Float64)
    config = MOI.Test.Config(
        exclude=Any[
            MOI.ConstraintBasisStatus,
            MOI.VariableBasisStatus,
            MOI.ConstraintDual,
            MOI.DualObjectiveValue,
            MOI.ObjectiveBound,
        ],
    )
    MOI.Test.runtests(
        bridged,
        config,
        exclude = String[
            "test_attribute_RawStatusString",
            "test_attribute_SolveTimeSec",
            "test_attribute_SolverVersion",
            #   MathOptInterface.jl issue #1431
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
        ],
    )
end
@testset "simplex chebyshev center with $T" for T in [Float64, Rational{BigInt}]
    h = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HyperPlane([1, 1], 1)
    lib = DefaultLibrary{T}(VRepOptimizer{T})
    poly = polyhedron(h, lib)
    center, radius = chebyshevcenter(poly)
    @test center ≈ [1/2, 1/2]
    @test radius ≈ 1/2
    @test dim(poly) == 1
end

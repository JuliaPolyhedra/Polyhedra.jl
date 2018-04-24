@testset "Dual Type" begin
    h = LPHRepresentation(spzeros(Int, 2, 2), [1, 2], [3, 4], [4, 5], [6, 7])
    p = @inferred polyhedron(h)
    @test p isa SimplePolyhedron{2, Rational{BigInt}, LPHRepresentation{2, Rational{BigInt}, SparseMatrixCSC{Rational{BigInt},Int}}, Polyhedra.Hull{2, Rational{BigInt}, Vector{Rational{BigInt}}}}
    h = hrep(zeros(2, 2), zeros(2))
    p = @inferred polyhedron(h)
    @test p isa SimplePolyhedron{2, Float64, MixedMatHRep{2, Float64, Matrix{Float64}}, MixedMatVRep{2, Float64, Matrix{Float64}}}
    h = hrep(spzeros(2, 2), zeros(2))
    p = @inferred polyhedron(h)
    @test p isa SimplePolyhedron{2, Float64, MixedMatHRep{2, Float64, SparseMatrixCSC{Float64, Int}}, MixedMatVRep{2, Float64, Matrix{Float64}}}
    v = vrep(zeros(2, 3))
    p = @inferred polyhedron(v)
    @test p isa SimplePolyhedron{3, Float64, MixedMatHRep{3, Float64, Matrix{Float64}}, MixedMatVRep{3, Float64, Matrix{Float64}}}
    v = vrep(spzeros(2, 3))
    p = @inferred polyhedron(v)
    @test p isa SimplePolyhedron{3, Float64, MixedMatHRep{3, Float64, Matrix{Float64}}, MixedMatVRep{3, Float64, SparseMatrixCSC{Float64, Int}}}
end

struct BadPoly{N, T} <: Polyhedra.Polyhedron{N, T}
end

@testset "Unimplemented methods" begin
    p = BadPoly{2, Int}()
    @test_throws ErrorException detectvlinearity!(p)
    @test_throws ErrorException detecthlinearity!(p)
    @test_throws ErrorException eliminate(p, [1], FourierMotzkin())
    @test_throws ErrorException eliminate(p, [1], BlockElimination())
    @test_throws ErrorException hrep(p)
    @test_throws ErrorException vrep(p)
    @test_throws ErrorException loadpolyhedron!(p, "a", "b")
    @test_throws ErrorException loadpolyhedron!(p, "a", "ine")
    @test_throws ErrorException loadpolyhedron!(p, "a", "ext")
end

@testset "SimplePolyhedron constructor with nothing" begin
    vr = convexhull([-1, 0], [0, -1]) + conichull([1, 1], [-1, 1])
    hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1)
    p = SimplePolyhedron{2, Int, typeof(hr), typeof(vr)}(hr, nothing, lp_solver)
    @test hrep(p) === hr
    @test !vrepiscomputed(p)
    @test p.solver === lp_solver
    p = SimplePolyhedron{2, Int, typeof(hr), typeof(vr)}(nothing, vr, lp_solver)
    @test !hrepiscomputed(p)
    @test vrep(p) === vr
    @test p.solver === lp_solver
end

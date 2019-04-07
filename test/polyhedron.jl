@testset "Dual Type" begin
    model = JuMP.Model()
    x = @variable(model, [1:2])
    @constraint(model, sum(x) == 1)
    h = hrep(model)
    p = polyhedron(h)
    @test p isa Polyhedra.DefaultPolyhedron{Float64, LPHRep{Float64}, Polyhedra.Hull{Float64, Vector{Float64}, Int}}
    h = hrep(zeros(Int, 2, 2), zeros(Int, 2))
    p = polyhedron(h)
    @test p isa Polyhedra.DefaultPolyhedron{Rational{BigInt}, MixedMatHRep{Rational{BigInt}, Matrix{Rational{BigInt}}}, MixedMatVRep{Rational{BigInt}, Matrix{Rational{BigInt}}}}
    h = hrep(spzeros(2, 2), zeros(2))
    p = polyhedron(h)
    @test p isa Polyhedra.DefaultPolyhedron{Float64, MixedMatHRep{Float64, SparseMatrixCSC{Float64, Int}}, MixedMatVRep{Float64, Matrix{Float64}}}
    v = vrep(zeros(2, 3))
    p = polyhedron(v)
    @test p isa Polyhedra.DefaultPolyhedron{Float64, MixedMatHRep{Float64, Matrix{Float64}}, MixedMatVRep{Float64, Matrix{Float64}}}
    v = vrep(spzeros(2, 3))
    p = polyhedron(v)
    @test p isa Polyhedra.DefaultPolyhedron{Float64, MixedMatHRep{Float64, Matrix{Float64}}, MixedMatVRep{Float64, SparseMatrixCSC{Float64, Int}}}
end

struct BadPoly{T} <: Polyhedra.Polyhedron{T}
end

@testset "Unimplemented methods" begin
    p = BadPoly{Int}()
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

@testset "Polyhedra.DefaultPolyhedron constructor with nothing" begin
    vr = convexhull([-1, 0], [0, -1]) + conichull([1, 1], [-1, 1])
    hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1)
    p = Polyhedra.DefaultPolyhedron{Int, typeof(hr), typeof(vr)}(hr, nothing, lp_solver)
    @test hrep(p) === hr
    @test !vrepiscomputed(p)
    @test p.solver === lp_solver
    p = Polyhedra.DefaultPolyhedron{Int, typeof(hr), typeof(vr)}(nothing, vr, lp_solver)
    @test !hrepiscomputed(p)
    @test vrep(p) === vr
    @test p.solver === lp_solver
end

using Test

using Polyhedra

include("solvers.jl")

function _triangle(T)
    return hrep([
        HalfSpace(T[-1, 0], zero(T)),
        HalfSpace(T[0, -1], zero(T)),
        HalfSpace(T[1, 1], one(T)),
    ])
end

function _square(T)
    return hrep([
        HalfSpace(T[-1, 0], zero(T)),
        HalfSpace(T[0, -1], zero(T)),
        HalfSpace(T[1, 0], one(T)),
        HalfSpace(T[0, 1], one(T)),
    ])
end

# Edge of the triangle that is in the hyperplane `x + y = 1`
function _edge(T)
    return hrep(
        [HyperPlane(T[1, 1], one(T))],
        [HalfSpace(T[-1, 0], zero(T)), HalfSpace(T[0, -1], zero(T))],
    )
end

@testset "issubset $T" for T in (Float64, Rational{BigInt})
    lib = DefaultLibrary{T}(lp_solver)
    triangle, square, edge = _triangle(T), _square(T), _edge(T)
    vtriangle = convexhull(T[0, 0], T[1, 0], T[0, 1])
    vsquare = convexhull(T[0, 0], T[1, 0], T[0, 1], T[1, 1])
    vedge = convexhull(T[1, 0], T[0, 1])
    @testset "Element" begin
        @test issubset(vtriangle, HalfSpace(T[1, 1], one(T)))
        @test !issubset(vsquare, HalfSpace(T[1, 1], one(T)))
        @test issubset(triangle, HalfSpace(T[1, 1], one(T)), lp_solver)
        @test !issubset(square, HalfSpace(T[1, 1], one(T)), lp_solver)
        @test issubset(edge, HyperPlane(T[1, 1], one(T)), lp_solver)
        @test !issubset(triangle, HyperPlane(T[1, 1], one(T)), lp_solver)
    end
    @testset "V-representation in H-representation" begin
        @test vtriangle ⊆ square
        @test vtriangle ⊆ triangle
        @test !(vsquare ⊆ triangle)
        @test vedge ⊆ edge
        @test vedge ⊆ triangle
        @test !(vtriangle ⊆ edge)
        # A ray is not in a bounded polyhedron
        @test !(vtriangle + conichull(T[1, 0]) ⊆ square)
        orthant =
            hrep([HalfSpace(T[-1, 0], zero(T)), HalfSpace(T[0, -1], zero(T))])
        @test vtriangle + conichull(T[1, 0]) ⊆ orthant
        # A line is in the hyperplane if it is parallel to it
        @test vedge + Line(T[1, -1]) ⊆ hrep([HyperPlane(T[1, 1], one(T))])
        @test !(vedge + Line(T[1, -1]) ⊆ triangle)
    end
    @testset "H-representation in H-representation" begin
        @test issubset(triangle, square, lp_solver)
        @test !issubset(square, triangle, lp_solver)
        @test issubset(edge, triangle, lp_solver)
        @test !issubset(triangle, edge, lp_solver)
    end
    @testset "Polyhedron in polyhedron" begin
        for compute_vp in (false, true), compute_vq in (false, true)
            p = polyhedron(triangle, lib)
            q = polyhedron(square, lib)
            compute_vp && vrep(p)
            compute_vq && vrep(q)
            @test p ⊆ q
            @test !(q ⊆ p)
            @test p ⊆ p
        end
        # The H-representation of `q` is computed if needed
        p = polyhedron(vtriangle, lib)
        q = polyhedron(vsquare, lib)
        @test !hrepiscomputed(q)
        @test p ⊆ q
        @test !(q ⊆ p)
        @test polyhedron(vedge, lib) ⊆ polyhedron(edge, lib)
        @test polyhedron(edge, lib) ⊆ p
        @test !(p ⊆ polyhedron(edge, lib))
        # Polyhedron in representation and representation in polyhedron
        @test p ⊆ square
        @test vtriangle ⊆ q
        @test issubset(triangle, q, lp_solver)
    end
    @testset "Intersection #303" begin
        A = T[1 0; 0 1; -1 0; 0 -1; -1 -1]
        b = T[3, 3, 3, 3, 3]
        p = polyhedron(hrep(A, b), lib)
        @test p ⊆ intersect(p, p)
        @test intersect(p, p) ⊆ p
    end
    @testset "V-representation as superset" begin
        err = ErrorException(
            "`issubset(p, v::VRepresentation)` is not supported as it " *
            "requires the H-representation of `v`. Use a polyhedron for the " *
            "second argument instead, e.g., `issubset(p, polyhedron(v))`.",
        )
        @test_throws err vtriangle ⊆ vsquare
        @test_throws err issubset(triangle, vsquare, lp_solver)
        @test_throws err polyhedron(triangle, lib) ⊆ vsquare
    end
end

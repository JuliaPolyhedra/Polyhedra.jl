using Test
using Polyhedra

include("utils.jl")
include("solvers.jl")

function _test_block_elimination(p, delset)
    q = eliminate(p, delset, BlockElimination())
    @test !vrepiscomputed(q)
    expected = eliminate(p, delset, ProjectGenerators())
    removevredundancy!(expected)
    @test fulldim(q) == fulldim(expected)
    @test isempty(q) == isempty(expected)
    if !isempty(expected)
        # Block elimination may produce redundant halfspaces so we compare
        # the V-representations
        generator_fulltest(q, vrep(expected))
    end
end

@testset "Block elimination with $T" for T in [Rational{BigInt}, Float64]
    lib = DefaultLibrary{T}(lp_solver)
    @testset "Polytope" begin
        # x1 + x2 ≤ 1, x1 - x2 ≤ 0, x1 ≥ 0
        p = polyhedron(hrep(T[1 1; 1 -1; -1 0], T[1, 0, 0]), lib)
        _test_block_elimination(p, [1])
        _test_block_elimination(p, [2])
    end
    @testset "Hyperplanes" begin
        # x1 + x2 + x3 = 1, x ≥ 0
        p = polyhedron(
            hrep(T[1 1 1], T[1], BitSet(1)) ∩
            hrep(-T[1 0 0; 0 1 0; 0 0 1], zeros(T, 3)),
            lib,
        )
        _test_block_elimination(p, [3])
        _test_block_elimination(p, [1, 3])
        q = eliminate(p, [3], BlockElimination())
        @test !hashyperplanes(q)
    end
    @testset "Hyperplane not involving eliminated variables" begin
        # x1 = x2, -1 ≤ x3 ≤ 1
        p = polyhedron(
            hrep(T[1 -1 0; 0 0 1; 0 0 -1], T[0, 1, 1], BitSet(1)),
            lib,
        )
        q = eliminate(p, [3], BlockElimination())
        @test nhyperplanes(q) == 1
        @test !hashalfspaces(q)
        @test [1, 1] in q
        @test !([1, 0] in q)
    end
    @testset "Hyperplane only involving eliminated variables" begin
        # x3 = 0, 0 ≤ x1 ≤ 1, 0 ≤ x2 ≤ 1
        p = polyhedron(
            hrep(
                T[0 0 1; 1 0 0; 0 1 0; -1 0 0; 0 -1 0],
                T[0, 1, 1, 0, 0],
                BitSet(1),
            ),
            lib,
        )
        _test_block_elimination(p, [3])
        @test !hashyperplanes(eliminate(p, [3], BlockElimination()))
    end
    @testset "Only hyperplanes" begin
        # x1 + x2 + x3 = 1, x1 - x2 = 0
        p = polyhedron(hrep(T[1 1 1; 1 -1 0], T[1, 0], BitSet(1:2)), lib)
        q = eliminate(p, [3], BlockElimination())
        @test nhyperplanes(q) == 1
        @test !hashalfspaces(q)
        @test [1, 1] in q
        @test !([1, 0] in q)
    end
    @testset "Eliminate no variable" begin
        p = polyhedron(hrep(T[1 1; 1 -1; -1 0], T[1, 0, 0]), lib)
        _test_block_elimination(p, Int[])
    end
    @testset "Eliminate all variables" begin
        p = polyhedron(hrep(T[1 1; 1 -1; -1 0], T[1, 0, 0]), lib)
        q = eliminate(p, [1, 2], BlockElimination())
        @test fulldim(q) == 0
        @test !hashyperplanes(q)
        @test !hashalfspaces(q)
        # x1 + x2 ≤ -1, x1 + x2 ≥ 1
        p = polyhedron(hrep(T[1 1; -1 -1], T[-1, -1]), lib)
        q = eliminate(p, [1, 2], BlockElimination())
        @test fulldim(q) == 0
        @test nhalfspaces(q) == 1
        @test first(halfspaces(q)).β < 0
    end
    @testset "Unbounded" begin
        # x1 ≥ |x2|, x3 free
        p = polyhedron(hrep(T[-1 1 0; -1 -1 0], T[0, 0]), lib)
        _test_block_elimination(p, [2])
        _test_block_elimination(p, [3])
        _test_block_elimination(p, [2, 3])
    end
    @testset "Empty" begin
        # x1 + x2 ≤ -1, x1 + x2 ≥ 1
        p = polyhedron(hrep(T[1 1; -1 -1], T[-1, -1]), lib)
        q = eliminate(p, [2], BlockElimination())
        @test isempty(q)
    end
    @testset "Full space" begin
        p = polyhedron(hrep(HalfSpace{T,Vector{T}}[]; d = 3), lib)
        q = eliminate(p, [2], BlockElimination())
        @test fulldim(q) == 2
        @test !hashyperplanes(q)
        @test !hashalfspaces(q)
    end
    @testset "Default algorithm" begin
        p = polyhedron(hrep(T[1 1; 1 -1; -1 0], T[1, 0, 0]), lib)
        @test !supportselimination(p, BlockElimination())
    end
end

@testset "LiftedHRepresentation" begin
    # 0 ≤ x1 ≤ 1, 0 ≤ x2 ≤ 1
    h = LiftedHRepresentation([1 -1 0; 1 0 -1; 0 1 0; 0 0 1])
    p = polyhedron(h)
    generator_fulltest(p, convexhull([0, 0], [1, 0], [0, 1], [1, 1]))
    v = Polyhedra.dualfullspace(h, 2, Int)
    @test v isa LiftedVRepresentation{Int}
    @test npoints(v) == 1
    @test nlines(v) == 2
    @test !hasrays(v)
end

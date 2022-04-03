using LinearAlgebra, Test
using StaticArrays
using Polyhedra

function test_cartesian_product()
    p1 = polyhedron(hrep([HalfSpace([1.0], 1.0), HalfSpace([-1.0], 0.0)]))
    p2 = polyhedron(hrep([HalfSpace([1.0], 1.0), HalfSpace([-1.0], 0.0)]))
    p = Polyhedra.hcartesianproduct(p1, p2)
    @test p isa DefaultPolyhedron{Float64,
                                  Polyhedra.Intersection{Float64,StaticArrays.SArray{Tuple{2},Float64,1,2},StaticArrays.Size{(2,)}},
                                  Polyhedra.Hull{Float64,StaticArrays.SArray{Tuple{2},Float64,1,2},StaticArrays.Size{(2,)}}}
    expected = HalfSpace((@SVector [1.0, 0.0]), 1.0) ∩ HalfSpace((@SVector [-1.0, 0.0]), 0.0) ∩ HalfSpace((@SVector [0.0, 1.0]), 1.0) ∩ HalfSpace((@SVector [0.0, -1.0]), 0.0)
    inequality_fulltest(p, expected)
end

function _test_sethrep(fh, fv, args...)
    #             x_1 ≤ 2   x_1 ≥ -1 <=> -x_1 ≤ 1
    h1 = HalfSpace([1], 2) ∩ HalfSpace([-1], 1)
    v1 = convexhull([-1], [2])

    ph = polyhedron(h1)
    @test volume(ph) == 3
    inequality_fulltest(ph, h1)
    generator_fulltest(ph, v1)

    pv = polyhedron(v1)
    @test volume(pv) == 3
    inequality_fulltest(pv, h1)
    generator_fulltest(pv, v1)

    #             x_1 ≤ 1   x_1 ≥ -1 <=> -x_1 ≤ 1
    h2 = HalfSpace([1], 1) ∩ HalfSpace([-1], 1)
    v2 = convexhull([-1], [1])

    # Test that the three fields have been modified
    fh(ph, h2, args...)
    @test volume(ph) == 2
    inequality_fulltest(ph, h2)
    generator_fulltest(ph, v2)
    fv(pv, v2, args...)
    @test volume(pv) == 2
    inequality_fulltest(pv, h2)
    generator_fulltest(pv, v2)
end
function test_sethrep()
    _test_sethrep(Polyhedra.sethrep!, Polyhedra.setvrep!)
    _test_sethrep(Polyhedra.sethrep!, Polyhedra.setvrep!, Polyhedra.UNKNOWN_REDUNDANCY)
    _test_sethrep(Polyhedra.resethrep!, Polyhedra.resetvrep!)
    _test_sethrep(Polyhedra.resethrep!, Polyhedra.resetvrep!, Polyhedra.UNKNOWN_REDUNDANCY)
end

function test_redundancy()
    h = HalfSpace([1], 2) ∩ HalfSpace([-1], 1)
    p = polyhedron(h)
    function f()
        removevredundancy!(p)
        removevredundancy!(p, nothing; strongly = false, planar = false)
        detectvlinearity!(p)
        detectvlinearity!(p, nothing; verbose=1)
        removehredundancy!(p)
        removehredundancy!(p, nothing; strongly = false, planar = false)
        detecthlinearity!(p)
        detecthlinearity!(p, nothing; verbose=1)
    end
    alloc_test(f, 0)
end

@testset "Interval tests" begin
    # Closed interval
    h = HalfSpace([-1], 3.) ∩ HalfSpace([1.], 8)
    v = convexhull([-3.], [8.])

    v0 = [[1.], [8], [-2], [4], [-1], [-3], [5]]
    pp = polyhedron(vrep(v0), Polyhedra.DefaultLibrary{Float64}())
    @test !hrepiscomputed(pp)
    @test vrepiscomputed(pp)
    d = StaticArrays.Size{(1,)}()
    for p in (Interval{Float64, SVector{1, Float64}, typeof(d)}(pp),
              Interval{Float64, SVector{1, Float64}, typeof(d)}(d, Polyhedra.vreps(pp)...))
        @test similar_library(pp, 1) == IntervalLibrary{Float64}()
        @test similar_library(pp, StaticArrays.Size((1,))) == IntervalLibrary{Float64}()
        @test library(p) == IntervalLibrary{Float64}()
        @test similar_library(p, 2) == Polyhedra.DefaultLibrary{Float64}()
        @test similar_library(p, StaticArrays.Size((2,))) == Polyhedra.DefaultLibrary{Float64}()
        @test hrepiscomputed(p)
        @test vrepiscomputed(p)
        @test p isa Interval{Float64}
        @test !isempty(p)
        @test volume(p) == 11
        @test center_of_mass(p) ≈ [2.5]
        @test dim(p) == 1
        inequality_fulltest(p, h)
        generator_fulltest(p, v)
    end

    pp = polyhedron(h, Polyhedra.DefaultLibrary{Float64}())
    @test hrepiscomputed(pp)
    @test !vrepiscomputed(pp)
    for p in (Interval{Float64, SVector{1, Float64}, typeof(d)}(pp),
              Interval{Float64, SVector{1, Float64}, typeof(d)}(d, Polyhedra.hreps(pp)...))
        @test hrepiscomputed(p)
        @test vrepiscomputed(p)
        @test p isa Interval{Float64}
        @test !isempty(p)
        @test iszero(surface(p))
        @test volume(p) == 11
        @test center_of_mass(p) ≈ [2.5]
        @test dim(p) == 1
        generator_fulltest(p, v)
        inequality_fulltest(p, h)
    end

    # Singleton
    h = intersect(HyperPlane([1.], 2))
    v = convexhull([2.])

    v0 = reshape([2., 2., 2.], 3, 1)
    p = polyhedron(vrep(v0))
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test iszero(surface(p))
    @test iszero(volume(p))
    @test dim(p) == 0
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(h)
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test iszero(surface(p))
    @test iszero(volume(p))
    @test dim(p) == 0
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    # Left-open interval
    h = intersect(HalfSpace([1.], 2))
    v = v + conichull([-1.])

    r0 = reshape([-1., -1], 2, 1)
    p = polyhedron(vrep(v0, r0))
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test iszero(surface(p))
    @test dim(p) == 1
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(HalfSpace([1.], 3) ∩ HalfSpace([5.], 10) ∩ HalfSpace([2.], 6))
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test iszero(surface(p))
    @test dim(p) == 1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    # Right-open interval
    h = intersect(HalfSpace([-1.], -2))
    v = convexhull([2.]) + conichull([1.])

    p = polyhedron(convexhull([2.], [3.]) + conichull([1.], [3.]))
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test iszero(surface(p))
    @test dim(p) == 1
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(HalfSpace([-3.], -6) ∩ HalfSpace([-2.], -2))
    # TODO test value of p
    p = polyhedron(h)
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test iszero(surface(p))
    @test dim(p) == 1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    # Empty
    h = HyperPlane([0], 1) ∩ HyperPlane([1], 0)
    v = vrep(Line{Int, StaticArrays.SVector{1, Int}}[])

    p = polyhedron(v)
    @test p isa Interval{Int}
    @test isempty(p)
    @test iszero(surface(p))
    @test iszero(volume(p))
    @test dim(p) == -1
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(HalfSpace([1], 1) ∩ HalfSpace([-1], -2))
    @test p isa Interval{Int}
    @test isempty(p)
    @test iszero(surface(p))
    @test iszero(volume(p))
    @test dim(p) == -1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    p = polyhedron(intersect(HyperPlane([0], 1)))
    @test p isa Interval{Int}
    @test isempty(p)
    @test iszero(surface(p))
    @test iszero(volume(p))
    @test dim(p) == -1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    p = polyhedron(intersect(HalfSpace([0], -1)))
    @test p isa Interval{Int}
    @test isempty(p)
    @test iszero(surface(p))
    @test iszero(volume(p))
    @test dim(p) == -1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    # Line
    h = hrep(HyperPlane{Int, SVector{1, Int}}[])
    v = conichull(Line([1]))

    p = polyhedron(v)
    @test p isa Interval{Int}
    @test !isempty(p)
    @test iszero(surface(p))
    @test dim(p) == 1
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(h)
    @test p isa Interval{Int}
    @test !isempty(p)
    @test iszero(surface(p))
    @test dim(p) == 1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    # Symmetric interval
    h = HalfSpace([1], 1) ∩ HalfSpace([-1], 1)
    v = convexhull([1], [-1])

    p = polyhedron(convexhull([-1], [1], [0]))
    @test p isa Interval{Int}
    @test !isempty(p)
    @test iszero(surface(p))
    @test volume(p) == 2
    @test dim(p) == 1
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(h)
    @test p isa Interval{Int}
    @test !isempty(p)
    @test iszero(surface(p))
    @test volume(p) == 2
    @test dim(p) == 1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    @testset "Cartesian product (#132)" begin
        test_cartesian_product()
    end

    @testset "sethrep! (#267)" begin
        test_sethrep()
    end

    @testset "test_redundancy (#298)" begin
        test_redundancy()
    end
end

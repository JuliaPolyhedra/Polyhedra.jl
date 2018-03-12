@testset "Interval tests" begin
    # Closed interval
    h = HalfSpace([-1], 3.) ∩ HalfSpace([1.], 8)
    v = convexhull([-3.], [8.])

    v0 = [[1.], [8], [-2], [4], [-1], [-3], [5]]
    pp = polyhedron(vrep(v0), SimplePolyhedraLibrary{Float64}())
    p = Interval{Float64, SVector{1, Float64}}(pp)
    @test !hrepiscomputed(pp)
    @test vrepiscomputed(pp)
    @test hrepiscomputed(p)
    @test vrepiscomputed(p)
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test volume(p) == 11
    @test dim(p) == 1
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    pp = polyhedron(h, SimplePolyhedraLibrary{Float64}())
    p = Interval{Float64, SVector{1, Float64}}(pp)
    @test hrepiscomputed(pp)
    @test !vrepiscomputed(pp)
    @test hrepiscomputed(p)
    @test vrepiscomputed(p)
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test iszero(surface(p))
    @test volume(p) == 11
    @test dim(p) == 1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

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

    p = polyhedron(h)
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

    p = polyhedron(h)
    @test p isa Interval{Float64}
    @test !isempty(p)
    @test iszero(surface(p))
    @test dim(p) == 1
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    # Empty
    h = HyperPlane([0], 1) ∩ HyperPlane([1], 0)
    v = vrep(SymPoint{1, Int, Vector{Int}}[])

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
end

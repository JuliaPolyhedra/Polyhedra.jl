@testset "Interval tests" begin
    # Closed interval
    h = HalfSpace([-1], 3.) ∩ HalfSpace([1.], 8)
    v = convexhull([-3.], [8.])

    v0 = [[1.], [8], [-2], [4], [-1], [-3], [5]]
    p = polyhedron(vrep(v0))
    @test p isa Interval{Float64}
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(h)
    @test p isa Interval{Float64}
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    # Singleton
    h = intersect(HyperPlane([1.], 2))
    v = convexhull([2.])

    v0 = reshape([2., 2., 2.], 3, 1)
    p = polyhedron(vrep(v0))
    @test p isa Interval{Float64}
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(h)
    @test p isa Interval{Float64}
    generator_fulltest(p, v)
    inequality_fulltest(p, h)

    # Left-open interval
    h = intersect(HalfSpace([1.], 2))
    v = v + conichull([-1.])

    r0 = reshape([-1., -1], 2, 1)
    p = polyhedron(vrep(v0, r0))
    @test p isa Interval{Float64}
    inequality_fulltest(p, h)
    generator_fulltest(p, v)

    p = polyhedron(h)
    @test p isa Interval{Float64}
    generator_fulltest(p, v)
    inequality_fulltest(p, h)
end

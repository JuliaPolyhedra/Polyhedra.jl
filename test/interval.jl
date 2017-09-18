@testset "Interval tests" begin
    v0 = reshape([1., 8, -2, 4, -1, -3, 5], 7, 1)
    v1 = reshape([-3., 8], 2, 1)
    p = polyhedron(SimpleVRepresentation(v0))
    @test p isa Interval{Float64}
    A = reshape([1., -1], 2, 1)
    b = [8, 3]
    inequality_fulltest(p, A, b, IntSet())
    generator_fulltest(p, v1)

    p = polyhedron(SimpleHRepresentation(A, b, IntSet()))
    @test p isa Interval{Float64}
    generator_fulltest(p, v1)
    inequality_fulltest(p, A, b, IntSet())

    v0 = reshape([2., 2., 2.], 3, 1)
    v1 = reshape([2.], 1, 1)
    p = polyhedron(SimpleVRepresentation(v0))
    @test p isa Interval{Float64}
    A = reshape([1.], 1, 1)
    b = [2]
    inequality_fulltest(p, A, b, IntSet([1]))
    generator_fulltest(p, v1)

    p = polyhedron(SimpleHRepresentation(A, b, IntSet([1])))
    @test p isa Interval{Float64}
    generator_fulltest(p, v1)
    inequality_fulltest(p, A, b, IntSet([1]))

    r0 = reshape([-1., -1], 2, 1)
    r1 = reshape([-1.], 1, 1)
    p = polyhedron(SimpleVRepresentation(v0, r0))
    @test p isa Interval{Float64}
    inequality_fulltest(p, A, b, IntSet())
    generator_fulltest(p, v1, r1)

    p = polyhedron(SimpleHRepresentation(A, b, IntSet()))
    @test p isa Interval{Float64}
    generator_fulltest(p, v1, r1)
    inequality_fulltest(p, A, b, IntSet())
end

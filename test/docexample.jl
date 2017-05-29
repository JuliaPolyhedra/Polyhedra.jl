function doctest(lib::PolyhedraLibrary)
    A = [1 1;1 -1;-1 0]
    b = [1,0,0]
    V = [1//2 1//2; 0 1; 0 0]
    hrep = SimpleHRepresentation(A, b)
    p = polyhedron(hrep, lib)
    solver = defaultLPsolverfor(p, lpsolver)
    @test in(HalfSpace([1, -1], 0), p, solver)
    @test !(in(HalfSpace([1, 1], 0), p, solver))
    @test !(in(HyperPlane([-1, 0], 0), p, solver))
    inequality_fulltest(p, A, b, IntSet())
    generator_fulltest(p, V)
    @test hrepiscomputed(p)
    @test vrepiscomputed(p) # computation triggered in generator_fulltest
    @test [1//4, 1//2] in p
    @test !([1, 1//2] in p)
    @test [0, 1] in p
    @test !(Ray([0, 1]) in p)
    Ap = [1; -1]
    bp = [1, 0]
    inequality_fulltest(eliminate(p, [1]), Ap, bp, IntSet())
    inequality_fulltest(project(p, [2]), Ap, bp, IntSet())
    inequality_fulltest(eliminate(p, [1], Val{:ProjectGenerators}), Ap, bp, IntSet())
    inequality_fulltest(project(p, [2], Val{:ProjectGenerators}), Ap, bp, IntSet())
    if implementseliminationmethod(p, Val{:BlockElimination})
        inequality_fulltest(eliminate(p, [1], Val{:BlockElimination}), Ap, bp, IntSet())
        inequality_fulltest(project(p, [2], Val{:BlockElimination}), Ap, bp, IntSet())
    end
    inequality_fulltest(project(p, reshape([0, 1], 2, 1)), Ap, bp, IntSet())
    inequality_fulltest(fixandeliminate(p, 2, eltype(p) <: AbstractFloat ? 1/4 : 1//4), reshape([1, -1], 2, 1), [1/4, 0], IntSet())
    # tries to create CDDPolyhedron of BigFloat, needs to let CDD decide of the promotion
    #inequality_fulltest(fixandeliminate(p, 2, 1/4, reshape([1, -1], 2, 1), [1/4, 0], IntSet())
end

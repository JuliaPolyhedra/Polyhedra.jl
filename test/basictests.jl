function basicpolyhedrontests(lib::PolyhedraLibrary)
    @testset "Polyhedron eltype conversion tests with $(typeof(lib))" begin
        A = [1 0; 0 1]
        b = [0, 0]
        linset = IntSet([2])
        R = [1 1]
        h = SimpleHRepresentation(A, b, linset)
        v = SimpleVRepresentation(A, R)

        p = copy(polyhedron(h, lib))
        @test hrepiscomputed(p)
        @test !vrepiscomputed(p)
        p = Polyhedron{2, Float64}(p)
        @test hrepiscomputed(p)
        @test !vrepiscomputed(p)
        inequality_fulltest(SimpleHRepresentation(hrep(p)), A, b, linset)
        p = copy(polyhedron(v, lib))
        @test !hrepiscomputed(p)
        @test vrepiscomputed(p)
        p = Polyhedron{2, Float64}(p)
        @test !hrepiscomputed(p)
        @test vrepiscomputed(p)
        generator_fulltest(SimpleVRepresentation(vrep(p)), A, R)
    end
    @testset "Polyhedron cartesian product tests with $(typeof(lib))" begin
        pv1 = polyhedron(SimpleVRepresentation([1 2; 3 4], [5 6]), lib)
        pv2 = polyhedron(SimpleVRepresentation([5 6; 7 8], [9 1; 2 3]), lib)
        qv = pv1 * pv2
        generator_fulltest(SimpleVRepresentation(vrep(qv)), [1 2 5 6; 1 2 7 8; 3 4 5 6; 3 4 7 8], [5 6 0 0; 0 0 9 1; 0 0 2 3])
        ph1 = polyhedron(SimpleHRepresentation([1 2; 3 4], [5, 6]), lib)
        ph2 = polyhedron(SimpleHRepresentation([7 8; 9 1], [2, 3]), lib)
        qh = ph1 * ph2
        inequality_fulltest(SimpleHRepresentation(hrep(qh)), [1 2 0 0; 3 4 0 0; 0 0 7 8; 0 0 9 1], [5, 6, 2, 3], IntSet())
    end
end

function basicpolyhedrontests(lib::PolyhedraLibrary)
    @testset "Polyhedron eltype conversion tests with $(typeof(lib))" begin
        A = [1 0; 0 1]
        b = [0, 0]
        linset = IntSet([2])
        R = [1 1]
        h = SimpleHRepresentation(A, b, linset)
        v = SimpleVRepresentation(A, R)
        p = polyhedron(h, lib)
        @test hrepiscomputed(p)
        @test !vrepiscomputed(p)
        p = Polyhedron{2, Float64}(p)
        @test hrepiscomputed(p)
        @test !vrepiscomputed(p)
        inequality_fulltest(SimpleHRepresentation(hrep(p)), A, b, linset)
        p = polyhedron(v, lib)
        @test !hrepiscomputed(p)
        @test vrepiscomputed(p)
        p = Polyhedron{2, Float64}(p)
        @test !hrepiscomputed(p)
        @test vrepiscomputed(p)
        generator_fulltest(SimpleVRepresentation(vrep(p)), A, R)
    end
end

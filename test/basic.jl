function basictest(lib::Polyhedra.Library)
    @testset "Polyhedron eltype conversion tests with $(typeof(lib))" begin
        A = [1 0; 0 1]
        b = [0, 0]
        linset = BitSet([2])
        R = [1 1]
        h = hrep(A, b, linset)
        v = vrep(A, R)

        p = copy(polyhedron(h, lib))
        @test hrepiscomputed(p)
        @test !vrepiscomputed(p)
        p = Polyhedron{Float64}(p)
        @test hrepiscomputed(p)
        @test !vrepiscomputed(p)
        inequality_fulltest(hrep(p), A, b, linset)
        p = copy(polyhedron(v, lib))
        @test !hrepiscomputed(p)
        @test vrepiscomputed(p)
        p = Polyhedron{Float64}(p)
        @test !hrepiscomputed(p)
        @test vrepiscomputed(p)
        generator_fulltest(vrep(p), A, R)
    end
    @testset "Polyhedron cartesian product tests with $(typeof(lib))" begin
        v1 = vrep([1 2; 3 4], [5 6])
        v2 = vrep([5 6; 7 8], [9 1; 2 3])
        v = v1 * v2
        expv = vrep([1 2 5 6; 1 2 7 8; 3 4 5 6; 3 4 7 8], [5 6 0 0; 0 0 9 1; 0 0 2 3])
        generator_fulltest(v, expv)
        pv1 = polyhedron(v1, lib)
        pv2 = polyhedron(v2, lib)
        qv = pv1 * pv2
        generator_fulltest(vrep(qv), expv)
        h1 = hrep([1 2; 3 4], [5, 6])
        h2 = hrep([7 8; 9 1], [2, 3])
        h = h1 * h2
        exph = hrep([1 2 0 0; 3 4 0 0; 0 0 7 8; 0 0 9 1], [5, 6, 2, 3])
        inequality_fulltest(h, exph)
        ph1 = polyhedron(h1, lib)
        ph2 = polyhedron(h2, lib)
        qh = ph1 * ph2
        inequality_fulltest(hrep(qh), exph)
    end
end

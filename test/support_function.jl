function support_function_test(p::Rep)
    @test support_function([1, 1], p) ≈ 1
    @test support_function([1, 0], p) ≈ 1
    @testset "Scale" begin
        @test support_function([1, 1], 2p) ≈ 2
        @test support_function([1, 0], 2p) ≈ 2
    end
end
function support_function_test(lib::Polyhedra.Library)
    v = convexhull([1, 0], [0, 1], [-1, -1])
    @testset "V-representation" begin
        support_function_test(v)
    end
    p = polyhedron(v, lib)
    @testset "Polyhedron with V-representation" begin
        support_function_test(p)
    end
    @test !hrepiscomputed(p)
    h = hrep(polyhedron(v, lib))
    # We don't test just with the H-rep since it needs a solver
    p = polyhedron(h, lib)
    @testset "Polyhedron with H-representation" begin
        support_function_test(p)
    end
    @test !vrepiscomputed(p)
end

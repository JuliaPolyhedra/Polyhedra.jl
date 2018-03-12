mutable struct TestLib <: PolyhedraLibrary
end
mutable struct TestPoly{N, T} <: Polyhedron{N, T}
end

@testset "Testing error for unimplemented algo" begin
    Polyhedra.polyhedron{N, T}(h::Representation{N, T}, ::TestLib) = TestPoly{N, T}()

    h = hrep([1 0; 0 1], [0, 0])
    v = vrep([1 0; 0 1], [1 1])
    p = polyhedron(h, TestLib())

    @test_throws ErrorException push!(p, h)
    @test_throws ErrorException push!(p, v)
    @test_throws ErrorException hrepiscomputed(p)
    @test_throws ErrorException hrep(p)
    @test_throws ErrorException vrepiscomputed(p)
    @test_throws ErrorException vrep(p)
    @test_throws ErrorException loadpolyhedron!(p, "test", "ine")
    @test_throws ErrorException loadpolyhedron!(p, "test", "ext")
    for algo in (FourierMotzkin(), BlockElimination())
        @test !supportselimination(p, algo)
        @test_throws ErrorException eliminate(p, [1], algo)
    end
end

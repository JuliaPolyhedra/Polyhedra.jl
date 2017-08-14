mutable struct TestLib <: PolyhedraLibrary
end
mutable struct TestPoly{N, T} <: Polyhedron{N, T}
end

@testset "Testing error for unimplemented method" begin
    Polyhedra.polyhedron{N, T}(h::Representation{N, T}, ::TestLib) = TestPoly{N, T}()

    h = SimpleHRepresentation([1 0; 0 1], [0, 0])
    v = SimpleVRepresentation([1 0; 0 1], [1 1])
    p = polyhedron(h, TestLib())

    @test_throws ErrorException push!(p, h)
    @test_throws ErrorException push!(p, v)
    @test_throws ErrorException hrepiscomputed(p)
    @test_throws ErrorException hrep(p)
    @test_throws ErrorException vrepiscomputed(p)
    @test_throws ErrorException vrep(p)
    @test_throws ErrorException loadpolyhedron!(p, "test", "ine")
    @test_throws ErrorException loadpolyhedron!(p, "test", "ext")
    for method in (Val{:FourierMotzkin}, Val{:BlockElimination})
        @test !implementseliminationmethod(p, method)
        @test_throws ErrorException eliminate(p, [1], method)
    end
end

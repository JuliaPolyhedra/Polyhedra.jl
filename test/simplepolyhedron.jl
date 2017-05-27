type TestLib <: PolyhedraLibrary
end
type TestPoly{N, T} <: Polyhedron{N, T}
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

tests = Tuple{String, Function}[]
push!(tests, ("Hypercube in 2 dimensions", lib->hypercubetest(lib, 2)))
push!(tests, ("Simplex in 2 dimensions", lib->simplextest(lib, 2))) # linearity
push!(tests, ("Simplex with the origin in 2 dimensions", lib->simplexorigtest(lib, 2)))
push!(tests, ("Cross Polytope in 2 dimensions", lib->crosspolytopetest(lib, 2)))
push!(tests, ("The ex1 example", ex1test))
push!(tests, ("Infeasible in 2 dimensions", lib->infeasibletest(lib, 2)))
push!(tests, ("Non full-dimensional", nonfulldimensionaltest)) #linearity

@testset "JuMP tests for SimplePolyhedron in float arithmetic" begin
    runtests(SimplePolyhedraLibrary(), tests)
end

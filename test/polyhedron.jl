struct BadPoly{N, T} <: Polyhedra.Polyhedron{N, T}
end

@testset "Unimplemented methods" begin
    p = BadPoly{2, Int}()
    @test_throws ErrorException detectvlinearities!(p)
    @test_throws ErrorException detecthlinearities!(p)
    @test_throws ErrorException eliminate(p, [1], FourierMotzkin())
    @test_throws ErrorException eliminate(p, [1], BlockElimination())
end

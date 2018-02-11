@testset "Duplicate removal" begin
    h = removeduplicates(SimpleHRepresentation([1 1; -1 -1; 1 0; 0 -1], [1, -1, 1, 0]))
    @test nhyperplanes(h) == 1
    @test first(hyperplanes(h)) == HyperPlane([1, 1], 1)
    @test first(hyperplanes(h)) == HyperPlane([-1, -1], -1)
    @test nhalfspaces(h) == 1
    h = removeduplicates(SimpleHRepresentation([1 1; -1 -1], [2, -1]))
    @test !hashyperplanes(h)
    @test nhalfspaces(h) == 2
end

@testset "Duplicate removal" begin
    h = removeduplicates(SimpleHRepresentation([1 1; -1 -1; 1 0; 0 -1], [1, -1, 1, 0]))
    @test neqs(h) == 1
    @test first(eqs(h)) == HyperPlane([1, 1], 1)
    @test first(eqs(h)) == HyperPlane([-1, -1], -1)
    @test nineqs(h) == 1
    h = removeduplicates(SimpleHRepresentation([1 1; -1 -1], [2, -1]))
    @test !haseqs(h)
    @test nineqs(h) == 2
end

@testset "VecRep" begin
    @testset "Intersection" begin
        h = HalfSpace([1, 0], 1) ∩ HyperPlane([0, 1], 1)
        @test nhyperplanes(h) == 1
        @test nhalfspaces(h) == 1
        hs = h.halfspaces
        hp = h.hyperplanes
        _hp = h.hyperplanes.hyperplanes
        intersect!(h, HalfSpace([1, 1], 2))
        @test h.halfspaces === hs
        @test h.hyperplanes === hp
        @test h.hyperplanes.hyperplanes === _hp
        @test nhyperplanes(h) == 1
        @test nhalfspaces(h) == 2
        intersect!(h, HyperPlane([1, 1], 3))
        @test h.halfspaces === hs
        @test h.hyperplanes === hp
        @test h.hyperplanes.hyperplanes === _hp
        @test nhyperplanes(h) == 2
        @test nhalfspaces(h) == 2
    end
end

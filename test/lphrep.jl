@testset "LPHRep" begin
    hs = HalfSpace([1, 0], 1)
    hp = HyperPlane([0, 1], 1)
    h = hs ∩ hp
    model = Model()
    @variable(model, [1:2] in h)
    lph = hrep(model)
    @test nhalfspaces(lph) == 1
    @test first(halfspaces(lph)) == hs
    @test nhyperplanes(lph) == 1
    @test first(hyperplanes(lph)) == hp
end

using Polyhedra
using Test
h = Polyhedra.HyperPlanesIntersection(2, [Polyhedra.HyperPlane{Int}([0, 1], 1)])
let
    i = 0
    for a in hyperplanes(h)
        i += 1
    end
    @test i == 0
end

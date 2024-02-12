using Polyhedra
using Aqua

@testset "aqua" begin
    Aqua.test_ambiguities(Polyhedra, broken=true)
    Aqua.test_all(Polyhedra, ambiguities=false)
end

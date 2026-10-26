using Polyhedra
using Aqua
using Test

@testset "Aqua" begin
    # Scalar and linear-map multiplication methods overlap for Line/Ray pairs,
    # whose product is undefined.
    Aqua.test_all(Polyhedra; ambiguities = (exclude = [*],))
end

using Test
using Polyhedra
using StaticArrays

function cartesian_interval_test(lib::Polyhedra.Library)
    G = [
        -1 0
        0 -1
        1 1
    ]
    F = [0, 0, 1]
    S = polyhedron(hrep(G, F), lib)
    U = polyhedron(convexhull([-2], [2]))
    @test Polyhedra.typed_fulldim(U) == StaticArrays.Size{(1,)}()
    # Mixes `S` whose `FullDim` is `2` and `U`
    # whose `FullDim` has type `StaticArrays.Size`.
    poly = S * U
    @test Polyhedra.typed_fulldim(poly) == 3
    h =
        HalfSpace([-1, 0, 0], 0) ∩ HalfSpace([0, -1, 0], 0) ∩
        HalfSpace([1, 1, 0], 1) ∩ HalfSpace([0, 0, 1], 2) ∩
        HalfSpace([0, 0, -1], 2)
    inequality_fulltest(poly, h)
    v = convexhull(
        [0, 0, 2],
        [0, 1, 2],
        [1, 0, 2],
        [0, 0, -2],
        [0, 1, -2],
        [1, 0, -2],
    )
    return generator_fulltest(poly, v)
end

if abspath(PROGRAM_FILE) == @__FILE__
    include("utils.jl")
    include("solvers.jl")
    @testset "Cartesian interval" begin
        cartesian_interval_test(DefaultLibrary{Float64}(lp_solver))
    end
end

function test_lphrep()
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

# Test that the order of variables is kept, see
# https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/338
# `MOI.Utilities.default_copy_to` will copy `y` first
# as it is a constrained variable. This tests that we do
# our custom `copy` that preserves order.
function test_order_lphrep()
    model = Model()
    @variable(model, x)
    @variable(model, y <= 1)
    h = hrep(model)
    @show h
    @test nhalfspaces(h) == 1
    @test first(halfspaces(h)) == HalfSpace([0, 1], 1)
end

using SparseArrays, Test, JuMP, Polyhedra

function _hypercube(n, base_name)
    model = Model()
    @variable(model, [1:n], lower_bound = 0, upper_bound = 1, base_name = base_name)
    h = hrep(model)
    @test dimension_names(h) == ["$base_name[$i]" for i in 1:n]
    return h
end
function hypercube_name()
    h = _hypercube(1, "x") * _hypercube(2, "y")
    @test dimension_names(h) == ["x[1]", "y[1]", "y[2]"]
    @test dimension_names(spzeros(3, 4) \ h) === nothing
end
@testset "Hypercube name" begin
    hypercube_name()
end

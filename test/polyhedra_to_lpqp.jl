using JuMP

@testset "computeoffsets" begin
    m = Model()
    n = 3
    @variable(m, x[1:n] >= 0)
    @constraint(m, sum(x[i] for i=1:n) == 1)
    lp = hrep(m)
    coloffset, rowoffset = Polyhedra.computeoffsets(lp)
    @test coloffset == [[-2], [-3], [-4]]
    @test rowoffset == [[1]]
end

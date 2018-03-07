@testset "computeoffsets" begin
    m = Model()
    @variable(m, x[1:n] >= 0)
    @constraint(m, sum(x[i] for i=1:n) == 1)
    lp = hrep(m)
    coloffset, rowoffset = computeoffset(lp)
    @show coloffset
    @show rowoffset
end

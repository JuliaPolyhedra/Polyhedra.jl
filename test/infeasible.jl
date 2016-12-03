using JuMP

# Inspired from Joey Huchette's test in ConvexHull.jl
function infeasibletest(lib::PolyhedraLibrary, n)
    m = Model()
    @variable(m, 0 ≤ x[1:n] ≤ 1)
    @constraint(m, x[1] ≥ 2)

    poly = polyhedron(m, lib)

    @test !hasrays(poly)
    @test !haspoints(poly)
    @test !hasvreps(poly)

    @test !haseqs(poly)
    @test hasineqs(poly)
    @test hashreps(poly)

    # haseqs has triggered a detection of linearity but
    # it shouldn't affect the V-representation
    @test !hasrays(poly)
    @test !haspoints(poly)
    @test !hasvreps(poly)
end

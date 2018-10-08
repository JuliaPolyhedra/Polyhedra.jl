using JuMP

# Inspired from Joey Huchette's test in ConvexHull.jl
function infeasibletest(lib::Polyhedra.Library, n)
    m = Model()
    @variable(m, 0 ≤ x[1:n] ≤ 1)
    @constraint(m, x[1] ≥ 2)

    poly = polyhedron(m, lib)

    @test !haspoints(poly)
    @test !haslines(poly)
    @test !hasrays(poly)
    @test !hasallrays(poly)

    @test !hashyperplanes(poly)
    @test hashalfspaces(poly)
    @test hasallhalfspaces(poly)

    # haseqs has triggered a detection of linearity but
    # it shouldn't affect the V-representation
    @test !haspoints(poly)
    @test !haslines(poly)
    @test !hasrays(poly)
    @test !hasallrays(poly)
end

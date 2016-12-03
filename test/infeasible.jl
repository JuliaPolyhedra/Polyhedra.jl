using JuMP

# Inspired from Joey Huchette's test in ConvexHull.jl
function infeasibletest(lib::PolyhedraLibrary, n)
    m = Model()
    @variable(m, 0 ≤ x[1:n] ≤ 1)
    @constraint(m, x[1] ≥ 2)

    lphrep = LPHRepresentation(m)
    poly = polyhedron(lphrep, lib)

    @fact hasrays(poly) --> false
    @fact haspoints(poly) --> false
    @fact hasvreps(poly) --> false

    @fact haseqs(poly) --> false
    @fact hasineqs(poly) --> true
    @fact hashreps(poly) --> true

    # haseqs has triggered a detection of linearity but
    # it shouldn't affect the V-representation
    @fact hasrays(poly) --> false
    @fact haspoints(poly) --> false
    @fact hasvreps(poly) --> false
end

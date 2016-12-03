using JuMP

# Inspired from Joey Huchette's test in ConvexHull.jl
function simplextest(lib::PolyhedraLibrary, n)
    m = Model()
    @variable(m, x[1:n] >= 0)
    @constraint(m, sum(x[i] for i=1:n) == 1)

    poly = polyhedron(m, lib)

    @test npoints(poly) == n
    @test nrays(poly) == 0

    V = zeros(Int, n, n)
    for k in 1:n
        V[k, k] = 1
    end
    generator_fulltest(poly, V)
end

# Inspired from Joey Huchette's test in ConvexHull.jl
function simplexorigtest(lib::PolyhedraLibrary, n)
    m = Model()
    @variable(m, x[1:n] >= 0)
    @constraint(m, sum(x[i] for i=1:n) ≤ 1)

    poly = polyhedron(m, lib)

    @test npoints(poly) == n+1
    @test nrays(poly) == 0

    V = zeros(Int, n+1, n)
    for k in 1:n
        V[k, k] = 1
    end
    generator_fulltest(poly, V)
end

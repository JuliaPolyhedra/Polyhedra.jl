using JuMP

# Inspired from Joey Huchette's test in ConvexHull.jl
function jumpsimplextest(lib::Polyhedra.Library, n)
    m = Model()
    @variable(m, x[1:n] >= 0)
    @constraint(m, sum(x[i] for i=1:n) == 1)

    poly = polyhedron(m, lib)

    V = zeros(Int, n, n)
    for k in 1:n
        V[k, k] = 1
    end
    generator_fulltest(poly, V)
end

# Inspired from Joey Huchette's test in ConvexHull.jl
function simplexorigtest(lib::Polyhedra.Library, n)
    m = Model()
    @variable(m, x[1:n] >= 0)
    @constraint(m, sum(x[i] for i=1:n) ≤ 1)

    poly = polyhedron(m, lib)

    V = zeros(Int, n+1, n)
    for k in 1:n
        V[k, k] = 1
    end
    generator_fulltest(poly, V)
end

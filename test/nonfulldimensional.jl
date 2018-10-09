using JuMP
using Combinatorics

# Inspired from Joey Huchette's test in ConvexHull.jl
function nonfulldimensionaltest(lib::Polyhedra.Library)
    m = Model()
    @variable(m, x[1:3] ≥ 1)
    @constraints(m, begin
        x[1] == 2
        x[2] ≤ 2
    end)

    poly = polyhedron(m, lib)

    V = [2 1 1;
         2 2 1]
    R = [0 0 1]
    generator_fulltest(poly, V, R)
end

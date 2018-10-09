using JuMP

# Inspired from Joey Huchette's test in ConvexHull.jl
function ex1test(lib::Polyhedra.Library)
    m = Model()
    @variable(m, x)
    @variable(m, y)
    @constraints(m, begin
        12 + 2x - y ≥ 0
        -6 - x + 2y ≥ 0
        -3 + x + y ≥ 0
        1 + x ≥ 0
    end)

    poly = polyhedron(m, lib)

    V = [ 0  3;
         -1  4;
         -1 10]

    R = [1 2;
         2 1]

    generator_fulltest(poly, V, R)
end

using JuMP
using Combinatorics

# Inspired from Joey Huchette's test in ConvexHull.jl
function crosspolytopetest(lib::Polyhedra.Library, n)
    model = Model()

    @variable(model, x[1:n])

    for k in 0:n
        for S in combinations(1:n, k)
            Sᶜ = setdiff(1:n, S)
            @constraint(model, sum(x[i] for i in S) - sum(x[i] for i in Sᶜ) ≤ 1)
        end
    end

    poly = polyhedron(model, lib)

    V = zeros(Int, 2n, n)
    for i in 1:n
        V[2i-1, i] = 1
        V[2i, i] = -1
    end
    generator_fulltest(poly, V)
end

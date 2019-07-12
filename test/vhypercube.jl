using LinearAlgebra
using Combinatorics

using Polyhedra

# Hypercube from `α * ones(n)` to `β * ones(n)`.
function vhypercubetest(lib::Polyhedra.Library, n, α, β)
    row = 0
    V = α * ones(Int, 2^n, n)
    for k in 0:n
        for p in combinations(1:n, k)
            row += 1
            V[row, p] = β * ones(Int, length(p))
        end
    end
    v = vrep(V)

    poly = polyhedron(v, lib)

    h = hrep([Matrix(1I, n, n); Matrix(-1I, n, n)], [β * ones(Int, n); -α * ones(Int, n)])

    inequality_fulltest(poly, h)
    generator_fulltest(poly, v)
end

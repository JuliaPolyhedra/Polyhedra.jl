A = [1 1; 0 1] - [1; 1] * [1.17 1.03]

using Test
using Polyhedra
Wv = vrep([[x, y] for x in [-1.0, 1.0] for y in [-1.0, 1.0]])

using GLPK
using JuMP
lib = DefaultLibrary{Float64}(GLPK.Optimizer)
W = polyhedron(Wv, lib)

function Fs(s::Integer, verbose=1)
    @assert s ≥ 1
    F = W
    A_W = W
    for i in 1:(s-1)
        A_W = A * A_W
        F += A_W
        if verbose ≥ 1
            println("Number of points after adding A^$i * W: ", npoints(F))
        end
        removevredundancy!(F)
        if verbose ≥ 1
            println("Number of points after removing redundant ones: ", npoints(F))
        end
    end
    return F
end

@test npoints(Fs(4)) == 16

function αo(s)
    A_W = A^s \ W
    hashyperplanes(A_W) && error("HyperPlanes not supported")
    return maximum([Polyhedra.support_function(h.a, W) / h.β for h in halfspaces(A_W)])
end

α = αo(10)
@test α ≈ 1.91907e-5

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


using Polyhedra

using SparseArrays
using LinearAlgebra
using Test

using StaticArrays

include("utils.jl")

include("solvers.jl")

include("elements.jl")
include("comp.jl")
include("representation.jl")
include("polyhedron.jl")

include("redundancy.jl")
include("doubledescription.jl")

include("interval.jl")

include("polyhedra_to_lpqp.jl")
include("default.jl")

# TODO fix on Julia v1.0
#include("show.jl")

include("polyhedra.jl")
for (arith, T) in (("floating point", Float64), ("exact", Rational{BigInt}))
    @testset "Polyhedra tests in $arith arithmetic" begin
        # TODO readd empty
        polyhedratest(SimplePolyhedraLibrary{T}(lp_solver), ["empty", "board"])
    end
end

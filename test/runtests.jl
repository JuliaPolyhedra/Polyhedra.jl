using Polyhedra
using Base.Test

import MultivariatePolynomials
MP = MultivariatePolynomials

using StaticArrays

include("utils.jl")
include("config.jl")
include("alltests.jl")

include("libraries.jl")

include("elements.jl")
include("comp.jl")
include("representation.jl")
include("polyhedron.jl")

include("redundancy.jl")
include("doubledescription.jl")

include("interval.jl")

include("polyhedra_to_lpqp.jl")
include("default.jl")

include("show.jl")

for arith in (("floating point", Float64), ("exact", Rational{BigInt}))
    @testset "Polyhedra tests in $arith arithmetic" begin
        polyhedratest(SimplePolyhedraLibrary{T}())
    end
end

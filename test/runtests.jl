using Polyhedra
using Base.Test

import MultivariatePolynomials
MP = MultivariatePolynomials

using StaticArrays

include("alltests.jl")

include("redundancy.jl")

include("libraries.jl")
if !isempty(lp_solvers)
    lpsolver = first(lp_solvers)
end

include("default.jl")
include("element.jl")
include("rep.jl")
include("interval.jl")
include("polyhedra_to_lpqp.jl")
include("show.jl")

for lib in libraries
    basicpolyhedrontests(lib)
end

#include("simplepolyhedron.jl")

for (testname, testfun) in alltests
    @testset "$testname tests" begin
        for arith in [("floating point", float_libraries), ("exact", exact_libraries)]
            arithname, arith_libraries = arith
            @testset "$testname tests in $arithname arithmetic with $(typeof(library))" for library in arith_libraries
                testfun(library)
            end
            length(arith_libraries) == 0 && warn("$testname tests in $arithname arithmetics test not run!")
        end
    end
end

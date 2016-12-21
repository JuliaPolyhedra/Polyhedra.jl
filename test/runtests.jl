using Polyhedra
using Base.Test

include("default.jl")
include("badargs.jl")

include("alltests.jl")

include("libraries.jl")

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

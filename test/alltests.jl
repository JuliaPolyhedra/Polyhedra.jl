using Polyhedra

include("jump.jl")
include("misc.jl")

function runtests(lib::PolyhedraLibrary, tests=alltests)
    for (testname, testfun) in tests
        @testset "$testname tests" begin
            testfun(lib)
        end
    end
end

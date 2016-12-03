using Polyhedra
using FactCheck

include("badargs.jl")

include("alltests.jl")

include("libraries.jl")

tests = Tuple{String, Function}[]
push!(tests, ("Hypercube in 2 dimensions", lib->hypercubetest(lib, 2)))
push!(tests, ("Simplex in 2 dimensions", lib->simplextest(lib, 2)))
push!(tests, ("Simplex with the origin in 2 dimensions", lib->simplexorigtest(lib, 2)))
push!(tests, ("Cross Polytope in 2 dimensions", lib->crosspolytopetest(lib, 2)))
push!(tests, ("The ex1 example", ex1test))
push!(tests, ("Infeasible in 2 dimensions", lib->infeasibletest(lib, 2)))
push!(tests, ("Non full-dimensional", nonfulldimensionaltest))
push!(tests, ("Simplex", simplextest))
push!(tests, ("Permutahedron", permutahedrontest))
push!(tests, ("Board", boardtest))

for test in tests
    testname, testfun = test
    for arith in [("floating point", float_libraries), ("exact", exact_libraries)]
        arithname, arith_libraries = arith
        facts("$testname tests in $arithname arithmetic") do
        for library in arith_libraries
        context("With library $(typeof(library))") do
            testfun(library)
        end; end; end
        length(arith_libraries) == 0 && warn("$testname tests in $arithname arithmetics test not run!")
    end
end

FactCheck.exitstatus()

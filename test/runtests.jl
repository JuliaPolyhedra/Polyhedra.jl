using Polyhedra
using FactCheck

facts("Simple Representation with bad arguments") do
    A = [1 1; -1 0; 0 -1]
    b = [1, 0, 0]
    linset = IntSet([1])

    @fact_throws ErrorException SimpleHRepresentation(A, [0, 0], linset)
    @fact_throws ErrorException SimpleHRepresentation(A, [0, 0], IntSet([4]))
    ine = SimpleHRepresentation(A, b, linset)
    @fact fulldim(ine) --> 2
    @fact eltype(ine) --> Int

    V = [0 1; 1 0]
    @fact_throws ErrorException SimpleVRepresentation(V, [1 0 0], IntSet(), IntSet())
    @fact_throws ErrorException SimpleVRepresentation(V, [1 1], IntSet(), IntSet([2]))
    @fact_throws ErrorException SimpleVRepresentation(V, [1 1], IntSet([4]), IntSet())
    @fact_throws ErrorException SimpleVRepresentation(V, IntSet([4]))
    ext = SimpleVRepresentation(V)
    @fact fulldim(ext) --> 2
    @fact eltype(ext) --> Int
end

include("alltests.jl")

include("libraries.jl")

for test in [("Simplex", simplextest), ("Permutahedron", permutahedrontest), ("Board", boardtest)]
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

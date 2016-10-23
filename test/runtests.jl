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

FactCheck.exitstatus()

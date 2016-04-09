using Polyhedra
using Base.Test

A = [1 1; -1 0; 0 -1]
b = [1, 0, 0]
linset = IntSet([1])

@test_throws ErrorException SimpleHRepresentation(A, [0, 0], linset)
@test_throws ErrorException SimpleHRepresentation(A, [0, 0], IntSet([4]))
ine = SimpleHRepresentation(A, b, linset)
@test fulldim(ine) == 2
@test eltype(ine) == Int

V = [0 1; 1 0]
@test_throws ErrorException SimpleVRepresentation(V, [1 0 0], IntSet([]), IntSet([]))
@test_throws ErrorException SimpleVRepresentation(V, [1 1], IntSet([]), IntSet([2]))
@test_throws ErrorException SimpleVRepresentation(V, [1 1], IntSet([4]), IntSet([]))
@test_throws ErrorException SimpleVRepresentation(V, IntSet([4]))
ext = SimpleVRepresentation(V)
@test fulldim(ext) == 2
@test eltype(ext) == Int

using Polyhedra
using BenchmarkTools

let
    A = [1 1;1 -1;-1 0]
    b = [1,0,0]
    h = SimpleHRepresentation{2, Float64}(A, b)

    display(@benchmark(doubledescription($h)))
end

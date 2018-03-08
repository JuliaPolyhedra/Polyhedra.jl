using Polyhedra
using BenchmarkTools

# Polyhedra v0.2.3 with SimpleHRepresentation (now renamed MixedMatHRep)
# BenchmarkTools.Trial:
#   memory estimate:  150.98 KiB
#   allocs estimate:  3388
#   --------------
#   minimum time:     567.867 μs (0.00% GC)
#   median time:      579.217 μs (0.00% GC)
#   mean time:        610.671 μs (3.98% GC)
#   maximum time:     5.010 ms (75.87% GC)
#   --------------
#   samples:          8156
#   evals/sample:     1
#
# Polyhedra v0.3.0 with Intersection

let
    h = hrep([HalfSpace([1., 1.], 1.), HalfSpace([ 1., -1.], 0.), HalfSpace([-1.,  0.], 0.)])
    @show typeof(doubledescription(h))
    display(@benchmark(doubledescription($h)))
end

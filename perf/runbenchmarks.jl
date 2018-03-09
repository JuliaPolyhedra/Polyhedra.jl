using Polyhedra
using BenchmarkTools
using StaticArrays

# Intersection/Vector
# BenchmarkTools.Trial:
#   memory estimate:  133.73 KiB
#   allocs estimate:  2424
#   --------------
#   minimum time:     454.533 μs (0.00% GC)
#   median time:      461.111 μs (0.00% GC)
#   mean time:        488.875 μs (4.38% GC)
#   maximum time:     4.360 ms (87.57% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
# Intersection/SVector # FIXME why does it need more allocs ? It should be the opposite
# BenchmarkTools.Trial:
#   memory estimate:  572.52 KiB
#   allocs estimate:  21094
#   --------------
#   minimum time:     499.253 μs (0.00% GC)
#   median time:      612.736 μs (0.00% GC)
#   mean time:        866.356 μs (29.38% GC)
#   maximum time:     12.673 ms (93.11% GC)
#   --------------
#   samples:          5749
#   evals/sample:     1
# MixedMatHRep
# BenchmarkTools.Trial:
#   memory estimate:  152.77 KiB
#   allocs estimate:  2893
#   --------------
#   minimum time:     541.129 μs (0.00% GC)
#   median time:      550.350 μs (0.00% GC)
#   mean time:        586.288 μs (4.79% GC)
#   maximum time:     5.012 ms (84.70% GC)
#   --------------
#   samples:          8505
#   evals/sample:     1

let
    h = hrep([HalfSpace([ 1.,  1], 1),
              HalfSpace([ 1., -1], 0),
              HalfSpace([-1.,  0], 0)])
    display(@benchmark(doubledescription($h)))
end

let
    h = hrep([HalfSpace((@SVector [ 1,  1]), 1),
              HalfSpace((@SVector [ 1, -1]), 0),
              HalfSpace((@SVector [-1,  0]), 0)])
    display(@benchmark(doubledescription($h)))
end

let
    h = hrep([ 1  1
               1 -1
              -1  0],
             [1., 0, 0])
    display(@benchmark(doubledescription($h)))
end

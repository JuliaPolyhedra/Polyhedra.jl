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
#   memory estimate:  148.64 KiB
#   allocs estimate:  2775
#   --------------
#   minimum time:     479.105 μs (0.00% GC)
#   median time:      486.236 μs (0.00% GC)
#   mean time:        508.813 μs (3.73% GC)
#   maximum time:     3.555 ms (83.56% GC)
#   --------------
#   samples:          9809
#   evals/sample:     1

using BenchmarkTools
using Polyhedra

function vector(T)
    h = hrep([HalfSpace(T[ 1,  1], T(1)),
              HalfSpace(T[ 1, -1], T(0)),
              HalfSpace(T[-1,  0], T(0))])
    display(@benchmark(doubledescription($h)))
end

using StaticArrays
function static(T)
    h = hrep([HalfSpace((@SVector T[ 1,  1]), T(1)),
              HalfSpace((@SVector T[ 1, -1]), T(0)),
              HalfSpace((@SVector T[-1,  0]), T(0))])
    display(@benchmark(doubledescription($h)))
    @profview_allocs for _ in 1:10000
        doubledescription(h)
    end
end

function matrix(T)
    h = hrep(T[
        1  1
        1 -1
       -1  0],
        T[1, 0, 0],
    )
    display(@benchmark(doubledescription($h)))
end

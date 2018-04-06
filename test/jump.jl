include("hypercube.jl")
include("crosspolytope.jl")
include("ex1.jl")
include("infeasible.jl")
include("jumpsimplex.jl")
include("nonfulldimensional.jl")

const jumptests = Dict("hypercube"          => lib->hypercubetest(lib, 2),
                       "jumpsimplex"        => lib->jumpsimplextest(lib, 2),
                       "simplexorig"        => lib->simplexorigtest(lib, 2),
                       "crosspolytope"      => lib->crosspolytopetest(lib, 2),
                       "ex1"                => ex1test,
                       "infeasible"         => lib->infeasibletest(lib, 2),
                       "nonfulldimensional" => nonfulldimensionaltest)

@polytestset jump

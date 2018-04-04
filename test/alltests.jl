using Polyhedra

# Basic
include("basictests.jl")

# JuMP
include("hypercube.jl")
include("crosspolytope.jl")
include("ex1.jl")
include("infeasible.jl")
include("jumpsimplex.jl")
include("nonfulldimensional.jl")

# non-JuMP
include("simplex.jl")
include("permutahedron.jl")
include("board.jl")
include("docexample.jl")
include("decompose.jl")
include("issue48.jl")
include("empty.jl")

alltests = Tuple{String, Function}[]
push!(alltests, ("Hypercube in 2 dimensions", lib->hypercubetest(lib, 2)))
push!(alltests, ("Simplex in 2 dimensions", lib->simplextest(lib, 2)))
push!(alltests, ("Simplex with the origin in 2 dimensions", lib->simplexorigtest(lib, 2)))
push!(alltests, ("Cross Polytope in 2 dimensions", lib->crosspolytopetest(lib, 2)))
push!(alltests, ("The ex1 example", ex1test))
push!(alltests, ("Infeasible in 2 dimensions", lib->infeasibletest(lib, 2)))
push!(alltests, ("Non full-dimensional", nonfulldimensionaltest))
push!(alltests, ("Doc example", doctest))
push!(alltests, ("Simplex", simplextest))
push!(alltests, ("Permutahedron", permutahedrontest))
push!(alltests, ("Board", boardtest))
push!(alltests, ("Decompose", decompose))
push!(alltests, ("Issue 48", issue48test))
push!(alltests, ("Empty", emptytest))

function runtests(lib::PolyhedraLibrary, tests=alltests)
    for (testname, testfun) in tests
        @testset "$testname tests" begin
            testfun(lib)
        end
    end
end

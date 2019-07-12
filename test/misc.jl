include("basic.jl")
include("simplex.jl")
include("permutahedron.jl")
include("board.jl")
include("docexample.jl")
include("issue48.jl")
include("empty.jl")
include("sparse.jl")
include("sparserect.jl")
include("recipe.jl")
include("support_function.jl")
include("vhypercube.jl")

const misctests = Dict("basic" => basictest,
                       "doc" => doctest,
                       "simplex" => simplextest,
                       "permutahedron" => permutahedrontest,
                       "board" => boardtest,
                       "issue48" => issue48test,
                       "empty" => emptytest,
                       "sparse" => sparsetest,
                       "sparserect" => sparserecttest,
                       "recipe" => recipetest,
                       "support_function" => support_function_test,
                       "vhypercubetest1c" => lib -> vhypercubetest(lib, 1, -1, 1),
                       "vhypercubetest1u" => lib -> vhypercubetest(lib, 1,  1, 2),
                       "vhypercubetest2c" => lib -> vhypercubetest(lib, 2, -1, 1),
                       "vhypercubetest2u" => lib -> vhypercubetest(lib, 2,  1, 2),
                       "vhypercubetest3c" => lib -> vhypercubetest(lib, 3, -1, 1),
                       "vhypercubetest3u" => lib -> vhypercubetest(lib, 3,  1, 2))

@polytestset misc

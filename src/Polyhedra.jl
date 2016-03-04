module Polyhedra

abstract PolyhedraLibrary
abstract Polyhedron
export PolyhedraLibrary, Polyhedron

include("description.jl")
include("operations.jl")

end # module

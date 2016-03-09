module Polyhedra

using GeometryTypes

abstract PolyhedraLibrary
abstract Polyhedron <: GeometryPrimitive
export PolyhedraLibrary, Polyhedron

include("description.jl")
include("operations.jl")

end # module

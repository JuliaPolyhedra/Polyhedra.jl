__precompile__()

module Polyhedra

using FixedSizeArrays, GeometryTypes

export PolyhedraLibrary, Polyhedron, getlibrary, getlibraryfor

abstract PolyhedraLibrary
abstract Polyhedron{N,T} <: GeometryPrimitive{N,T}

# Libraries can implement this function to be the default for a certain type T
# getlibraryfor{T<:Real}(::Type{T})
# This function should be implemented by libraries
# or it would be like not implementing similar for AbstractArray
getlibraryfor{T}(p::Polyhedron, ::Type{T}) = getlibraryfor(T)
getlibrary(p::Polyhedron) = getlibraryfor(p, eltype(p))

import Base.intersect, Base.+, Base.*, Base.isempty, Base.copy, Base.push!, Base.length, Base.eltype, Base.start, Base.done, Base.next

# Definitions
include("elements.jl")
include("mycomp.jl")
include("representation.jl")
include("repop.jl")
include("operations.jl")

# Implementations
include("lphrep.jl")
include("simplerep.jl")
include("liftedrep.jl")
include("simplepolyhedron.jl")

# Optimization
include("opt.jl")
include("lpqp_to_polyhedra.jl")
include("polyhedra_to_lpqp.jl")
include("simplevrepsolver.jl")
include("default.jl")

# Visualization
include("show.jl")
include("decompose.jl")

end # module

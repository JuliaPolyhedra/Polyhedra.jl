module Polyhedra

using GeometryTypes

export PolyhedraLibrary, Polyhedron, getlibrary, getlibraryfor

abstract PolyhedraLibrary
abstract Polyhedron{N,T} <: GeometryPrimitive{N,T}

# Libraries can implement this function to be the default for a certain type T
# getlibraryfor{T<:Real}(::Type{T})
# This function should be implemented by libraries
# or it would be like not implementing similar for AbstractArray
getlibraryfor{T}(p::Polyhedron, ::Type{T}) = getlibraryfor(T)
getlibrary(p::Polyhedron) = getlibraryfor(p, eltype(p))

const threshold = 1e-8
myeqzero{T<:Real}(x::T) = x == zero(T)
myeqzero{T<:AbstractFloat}(x::T) = -threshold < x < threshold
myeq{T<:Real}(x::T, y::T) = x == y
myeq{T<:AbstractFloat}(x::T, y::T) = y < x+threshold && x < y+threshold
myeq{T<:Real, S<:Real}(x::T, y::S) = myeq(promote(x, y)...)
myeqzero{T<:Real}(x::Vector{T}) = myeqzero(sum(abs(x)))

include("representation.jl")
include("operations.jl")
include("decompose.jl")

end # module

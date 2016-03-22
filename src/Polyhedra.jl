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
mylt{T<:Real}(x::T, y::T) = x < y
mylt{T<:AbstractFloat}(x::T, y::T) = x+threshold < y
mylt{S<:Real,T<:Real}(x::S, y::T) = mylt(promote(x, y)...)
mygt{S<:Real, T<:Real}(x::S, y::T) = mylt(y, x)
myleq{T<:Real}(x::T, y::T) = x <= y
myleq{T<:AbstractFloat}(x::T, y::T) = x < y+threshold
myleq{S<:Real,T<:Real}(x::S, y::T) = myleq(promote(x, y)...)
mygt{S<:Real,T<:Real}(x::S, y::T) = mylt(y, x)
mygeq{T<:Real}(x::T, y::T) = myleq(y, x)
mypos{T<:Real}(x::T) = mygt(x, zero(T))
myneg{T<:Real}(x::T) = mylt(x, zero(T))
mynonneg{T<:Real}(x::T) = mygeq(x, zero(T))
mynonpos{T<:Real}(x::T) = myleq(x, zero(T))
myeqzero{T<:Real}(x::T) = x == zero(T)
myeqzero{T<:AbstractFloat}(x::T) = -threshold < x < threshold
myeq{T<:Real}(x::T, y::T) = x == y
myeq{T<:AbstractFloat}(x::T, y::T) = y < x+threshold && x < y+threshold
myeq{T<:Real, S<:Real}(x::T, y::S) = myeq(promote(x, y)...)
myeqzero{T<:Real}(x::Vector{T}) = myeqzero(sum(abs(x)))

import Base.intersect, Base.+, Base.*, Base.isempty, Base.copy, Base.push!

include("representation.jl")
include("operations.jl")
include("decompose.jl")

end # module

__precompile__()

module Polyhedra

using GeometryTypes

using MultivariatePolynomials

export PolyhedraLibrary, Polyhedron, getlibrary, getlibraryfor, FullDim

abstract type PolyhedraLibrary end
abstract type Polyhedron{N,T} <: GeometryPrimitive{N,T} end

import Base: intersect, ==, +, *, \, /, isempty, copy, push!, length, eltype, start, done, next

using StaticArrays
using StaticArrays.FixedSizeArrays: FixedVector

# Similar to StaticArrays.Size
struct FullDim{N}
    function FullDim{N}() where N
        new{N::Int}()
    end
end
fulldim(::FullDim{N}) where N = N
+(d1::FullDim{N1}, d2::FullDim{N2}) where {N1, N2} = FullDim{N1+N2}()

FullDim(v::AbstractVector) = FullDim{length(v)}()
FullDim(v::StaticArrays.StaticArray{S, T, 1}) where {S, T} = FullDim{S[1]}()
MultivariatePolynomials.coefficienttype(::Union{AbstractVector{T}, Type{<:AbstractVector{T}}}) where T = T
similar_type(::Type{<:Vector}, ::FullDim, ::Type{T}) where T = Vector{T}
similar_type(::Type{SparseVector{S, IT}}, ::FullDim, ::Type{T}) where {S, IT, T} = SparseVector{T, IT}
# We define the following method to avoid the fallback to call FullDim(Vector{...})
similar_type(::Type{<:Vector}, ::Type{T}) where T = Vector{T}
similar_type(::Type{SparseVector{S, IT}}, ::Type{T}) where {S, IT, T} = SparseVector{T, IT}
function similar_type(::Type{SAT}, ::FullDim{D}, ::Type{Tout}) where {SAT <: StaticArrays.StaticArray, D, Tout}
    StaticArrays..similar_type(SAT, Tout, Size(D))
end

similar_type(::Type{<:Point}, ::FullDim{N}, ::Type{T}) where {N, T} = Point{N,T}
similar_type(::Type{<:Vec}, ::FullDim{N}, ::Type{T}) where {N, T} = Vec{N,T}

# Definitions
include("elements.jl")
include("mycomp.jl")
include("representation.jl")

function similar_type(::Type{ET}, ::Type{Tout}) where {Tout, ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, FullDim(ET), Tout)
end
function similar_type(::Type{ET}, d::FullDim) where {ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, d, coefficienttype(ET))
end

include("repop.jl")
include("center.jl")
include("repelemop.jl")
include("operations.jl")
include("aff.jl")
include("redundancy.jl")
include("projection.jl")

# Implementations
include("lphrep.jl")
include("jump.jl")
include("vecrep.jl")
include("simplerep.jl")
include("liftedrep.jl")
include("doubledescription.jl")
include("interval.jl") # 1D polyhedron
include("simplepolyhedron.jl")

# Optimization
importall MathProgBase.SolverInterface
include("opt.jl")
include("lpqp_to_polyhedra.jl")
include("polyhedra_to_lpqp.jl")
include("simplevrepsolver.jl")
include("default.jl")

# Visualization
include("show.jl")
include("recipe.jl")
include("decompose.jl")

end # module

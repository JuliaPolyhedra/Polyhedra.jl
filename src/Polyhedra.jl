__precompile__()

module Polyhedra

using Compat
using Compat.SparseArrays
using Compat.LinearAlgebra

using MultivariatePolynomials

export PolyhedraLibrary, Polyhedron

abstract type PolyhedraLibrary end
abstract type Polyhedron{N, T} <: GeometryPrimitive{N, T} end

import MathProgBase
const MPB = MathProgBase
const MPBSI = MPB.SolverInterface

import StaticArrays
function similar_type(SAT::Type{<:StaticArrays.SVector},
                      size::StaticArrays.Size, T::Type)
    StaticArrays.similar_type(SAT, T, size)
end

"""
    FullDim(p)::FullDim

Similar to [`fulldim`](@ref) but used for type stability with the vector type.
If the vector type is `StaticArrays.SVector` then it returns a
`StaticArrays.Size`.
"""
function FullDim end

"""
    fulldim(p)::Int

Returns the full dimension of the polyhedron, polyhedron representation,
polyhedron representation element or vector.
"""
function fulldim end

# The dimension is Int for Vector, SparseVector and
# StaticArrays.Size for SVector
# This allows similar_type to be compiled as only one method for Vector
# and SparseVector and be type stable for SVector
const FullDim = Union{Int, StaticArrays.Size}
fulldim(N::Int) = N
fulldim(::StaticArrays.Size{N}) where N = N
FullDim(v::AbstractVector) = length(v)
FullDim(v::Union{StaticArrays.SVector{N}, Type{<:StaticArrays.SVector{N}}}) where N = StaticArrays.Size(v)

MultivariatePolynomials.coefficienttype(::Union{AbstractVector{T}, Type{<:AbstractVector{T}}}) where T = T
similar_type(::Type{<:Vector}, ::Int, ::Type{T}) where T = Vector{T}
similar_type(::Type{SparseVector{S, IT}}, ::Int, ::Type{T}) where {S, IT, T} = SparseVector{T, IT}

emptymatrix(::Type{MT}, m, n) where {MT<:AbstractMatrix} = MT(undef, m, n)
emptymatrix(::Type{SparseMatrixCSC{T, Int}}, m, n) where T = spzeros(T, m, n)
similar_type(::Type{<:Matrix}, ::Type{T}) where T = Matrix{T}
similar_type(::Type{SparseMatrixCSC{S, I}}, ::Type{T}) where {S, I, T} = SparseMatrixCSC{T, I}

# Interface/Definitions
include("elements.jl")
include("comp.jl")
include("representation.jl")
include("indices.jl")
include("incidence.jl")
include("iterators.jl")
include("polyhedron.jl")

# For retro-compatibility with CDD and LRS, remove in v0.3.4
vvectortype(RepT::Type{<:VRep}) = vectortype(RepT)
hvectortype(RepT::Type{<:HRep}) = vectortype(RepT)
# TODO Only define vectortype for the type for Polyhedron v0.4
function vectortype(RepT::Type{<:Polyhedron})
    @assert hvectortype(RepT) == vvectortype(RepT)
    hvectortype(RepT)
end
vectortype(HRepT::Type{<:HRepresentation}) = hvectortype(HRepT)
vectortype(VRepT::Type{<:VRepresentation}) = vvectortype(VRepT)
vectortype(p::Rep) = vectortype(typeof(p))

vectortype(::Type{<:AbstractSparseArray{T}}) where T = SparseVector{T, Int}
vectortype(::Type{<:AbstractMatrix{T}}) where T = Vector{T}

const arraytype = vectortype

hmatrixtype(RepT::Type{<:HRep}, T::Type) = matrixtype(similar_type(hvectortype(RepT), T))
vmatrixtype(RepT::Type{<:VRep}, T::Type) = matrixtype(similar_type(vvectortype(RepT), T))
matrixtype(::Type{<:AbstractVector{T}}) where T = Matrix{T}
matrixtype(::Type{<:AbstractSparseVector{T}}) where T = SparseMatrixCSC{T, Int}

function similar_type(::Type{ET}, ::Type{Tout}) where {Tout, ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, FullDim(ET), Tout)
end
function similar_type(::Type{ET}, d::FullDim) where {ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, d, coefficienttype(ET))
end

# Operations
include("repop.jl")
include("center.jl")
include("repelemop.jl")
include("aff.jl")
include("redundancy.jl")
include("projection.jl")

# Implementations of representations
include("vecrep.jl")
include("mixedrep.jl")
include("lphrep.jl")
include("jump.jl")
include("matrep.jl")
include("liftedrep.jl")
include("doubledescription.jl") # FIXME move it after projection.jl once it stops depending on LiftedRep
include("interval.jl") # 1D polyhedron
include("simplepolyhedron.jl")

# Optimization
include("opt.jl")
include("lpqp_to_polyhedra.jl")
include("polyhedra_to_lpqp.jl")
include("vrepsolver.jl")
include("default.jl")

# Visualization
include("show.jl")
include("recipe.jl")
include("decompose.jl")

end # module

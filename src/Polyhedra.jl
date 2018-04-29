__precompile__()

module Polyhedra

using GeometryTypes

using MultivariatePolynomials

export PolyhedraLibrary, Polyhedron, FullDim

abstract type PolyhedraLibrary end
abstract type Polyhedron{N,T} <: GeometryPrimitive{N,T} end

using StaticArrays
using StaticArrays.FixedSizeArrays: FixedVector

using MathProgBase
const MPB = MathProgBase

# Similar to StaticArrays.Size
struct FullDim{N}
    function FullDim{N}() where N
        new{N::Int}()
    end
end
fulldim(::FullDim{N}) where N = N
Base.:+(d1::FullDim{N1}, d2::FullDim{N2}) where {N1, N2} = FullDim{N1+N2}()

FullDim(v::AbstractVector) = FullDim{length(v)}()
FullDim(::Union{StaticArrays.SVector{N}, Type{<:StaticArrays.SVector{N}}}) where N = FullDim{N}()
MultivariatePolynomials.coefficienttype(::Union{AbstractVector{T}, Type{<:AbstractVector{T}}}) where T = T
similar_type(::Type{<:Vector}, ::FullDim, ::Type{T}) where T = Vector{T}
similar_type(::Type{SparseVector{S, IT}}, ::FullDim, ::Type{T}) where {S, IT, T} = SparseVector{T, IT}
# We define the following method to avoid the fallback to call FullDim(Vector{...})
similar_type(::Type{<:Vector}, ::Type{T}) where T = Vector{T}
similar_type(::Type{SparseVector{S, IT}}, ::Type{T}) where {S, IT, T} = SparseVector{T, IT}
function similar_type(::Type{SAT}, ::FullDim{D}, ::Type{Tout}) where {SAT <: StaticArrays.SVector, D, Tout}
    StaticArrays.similar_type(SAT, Tout, StaticArrays.Size(D))
end

similar_type(::Type{<:Point}, ::FullDim{N}, ::Type{T}) where {N, T} = Point{N,T}
similar_type(::Type{<:Vec}, ::FullDim{N}, ::Type{T}) where {N, T} = Vec{N,T}

emptymatrix(::Type{MT}, m, n) where {MT<:AbstractMatrix} = MT(m, n)
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
importall MathProgBase.SolverInterface
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

__precompile__()

module Polyhedra

using Compat
using Compat.SparseArrays

abstract type Polyhedron{T} end

similar_type(::Type{<:Vector}, ::Int, ::Type{T}) where T = Vector{T}

# Interface/Definitions
#include("dimension.jl")
import StaticArrays
const FullDim = Union{Int, StaticArrays.Size}
FullDim(::Type{<:AbstractVector}) = -1 # Shouldn't hurt as it will not be used
FullDim(v::AbstractVector) = length(v)
FullDim_convert(::Type{Int}, d::Int) = d
#include("elements.jl")
abstract type HRepElement{T, AT} end
struct HalfSpace{T, AT <: AbstractVector{T}} <: HRepElement{T, AT}
    a::AT
    β::T
    function HalfSpace{T, AT}(a::AT, β::T) where {T, AT<:AbstractVector{T}}
        new{T, AT}(a, β)
    end
end

struct HyperPlane{T, AT<:AbstractVector{T}} <: HRepElement{T, AT}
    a::AT
    β::T
    function HyperPlane{T, AT}(a::AT, β::T) where {T, AT<:AbstractVector{T}}
        new{T, AT}(a, β)
    end
end

HyperPlane{T}(a::AT, β::T) where {T, AT <: AbstractVector{T}} = HyperPlane{T, AT}(a, β)
HyperPlane{T}(a::AbstractVector, β) where {T} = HyperPlane{T}(_vec(T, a), T(β))

for ElemT in [:HalfSpace, :HyperPlane]
    @eval begin
        function similar_type(::Type{$ElemT{T, AT}}, dout::FullDim, ::Type{Tout}) where {T, AT, Tout}
            return $ElemT{Tout, similar_type(AT, dout, Tout)}
        end
    end
end


abstract type Representation{T <: Real} end
abstract type HRepresentation{T} <: Representation{T} end
abstract type VRepresentation{T} <: Representation{T} end
const  Rep{T} = Union{ Representation{T}, Polyhedron{T}}
const HRep{T} = Union{HRepresentation{T}, Polyhedron{T}}
const VRep{T} = Union{VRepresentation{T}, Polyhedron{T}}

include("indices.jl")
include("iterators.jl")

vvectortype(rep::Rep) = vvectortype(typeof(rep))
hvectortype(rep::Rep) = hvectortype(typeof(rep))
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

function similar_type(::Type{ET}, ::Type{Tout}) where {Tout, ET<:Union{HRepElement, Rep}}
    similar_type(ET, FullDim(ET), Tout)
end
function similar_type(::Type{ET}, d::FullDim) where {ET<:Union{HRepElement, Rep}}
    similar_type(ET, d, coefficient_type(ET))
end

# Operations
include("aff.jl")
FullDim(rep::Type{<:HyperPlanesIntersection}) = FullDim(hvectortype(rep))
# collect(::Vector) does a copy.
# lazy_collect avoid this copy in case `v` is a Vector
lazy_collect(v::Vector) = v
lazy_collect(v) = collect(v)
function FullDim_rec(it::It, its::Union{Rep, It}...)
    if isempty(it)
        return FullDim_rec(its...)
    else
        return FullDim(first(it))
    end
end

end # module

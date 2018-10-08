export Representation, HRepresentation, VRepresentation, fulldim
export RepIterator
export Rep

"""
    Representation{T<:Real}

Supertype for H-(or V-)representations with coefficient type `T`.
"""
abstract type Representation{T <: Real} end
"""
    HRepresentation{T<:Real}

Supertype for H-representations with coefficient type `T`.
"""
abstract type HRepresentation{T} <: Representation{T} end
"""
    VRepresentation{T<:Real}

Supertype for V-representations coefficient type `T`.
"""
abstract type VRepresentation{T} <: Representation{T} end

export Rep, HRep, VRep

const  Rep{T} = Union{ Representation{T}, Polyhedron{T}}
const HRep{T} = Union{HRepresentation{T}, Polyhedron{T}}
const VRep{T} = Union{VRepresentation{T}, Polyhedron{T}}

Base.broadcastable(x::Union{Rep, Library}) = Ref(x)

"""
    coefficient_type(rep::Rep)

Returns the type of the coefficients used in the representation of `rep`.
"""
coefficient_type(rep::Union{Rep{T}, Type{<:Rep{T}}}) where {T} = T

FullDim(rep::Type{<:VRep}) = FullDim(vvectortype(rep))
FullDim(rep::Type{<:HRep}) = FullDim(hvectortype(rep))
FullDim(rep::Type{<:Polyhedron}) = FullDim(hvectortype(rep))

# Check that it is either empty or it has a point
vconsistencyerror() = error("Non-empty V-representation must contain at least one point. If it is a polyhedral cone, the origin should be added.")
function checkvconsistency(vrep::VRep)
    if !haspoints(vrep) && hasallrays(vrep)
        vconsistencyerror()
    end
end
checkvconsistency(p::Polyhedron) = vrepiscomputed(p) && checkvconsistency(vrep(p))

# This method solves the ambiguity with the following methods and the general method
# Base.convert{T}(::Type{T}, p::T) = p
Base.convert(::Type{T}, p::T) where {T<:HRepresentation} = p
Base.convert(::Type{T}, p::T) where {T<:VRepresentation} = p

Base.convert(RepT::Type{<:HRep}, p::HRepresentation)            = hconvert(RepT, p)
Base.convert(RepT::Type{<:HRepresentation}, p::HRep)            = hconvert(RepT, p)
# avoid ambiguity
Base.convert(RepT::Type{<:HRepresentation}, p::HRepresentation) = hconvert(RepT, p)

Base.copy(rep::HRepresentation) = typeof(rep)(FullDim_hreps(rep)...)
Base.copy(rep::VRepresentation) = typeof(rep)(FullDim_vreps(rep)...)

function Polyhedron{S}(p::Polyhedron{T}) where {S, T}
    RepT = similar_type(typeof(p), S)
    if !hrepiscomputed(p) && vrepiscomputed(p)
        vconvert(RepT, p)
    else
        hconvert(RepT, p)
    end
end

Base.convert(RepT::Type{<:VRep}, p::VRepresentation)            = vconvert(RepT, p)
Base.convert(RepT::Type{<:VRepresentation}, p::VRep)            = vconvert(RepT, p)
# avoid ambiguity
Base.convert(RepT::Type{<:VRepresentation}, p::VRepresentation) = vconvert(RepT, p)

# Used by VRepPolyhedraModel
Base.convert(::Type{VRep}, p::VRepresentation) = p

change_coefficient_type(p::Rep{T}, ::Type{T}) where {T} = p
change_coefficient_type(p::Rep, T::Type) = convert(similar_type(typeof(p), T), p)

VRepresentation{T}(v::VRepresentation) where {T} = convert(similar_type(typeof(v), FullDim(v), T), v)
HRepresentation{T}(h::HRepresentation) where {T} = convert(similar_type(typeof(h), FullDim(h), T), h)

VRep{T}(v::VRepresentation) where {T} = VRepresentation{T}(v)
VRep{T}(p::Polyhedron) where {T} = Polyhedron{T}(p)
HRep{T}(h::HRepresentation) where {T} = HRepresentation{T}(h)
HRep{T}(p::Polyhedron) where {T} = Polyhedron{T}(p)

# FIXME it does not get called. The calls always go throug vconvert and hconvert. Use changecoefficient_type instead
#function Base.convert{T, RepT<:Representation}(::Type{Representation{T}}, rep::RepT)
#  if fulldim(RepT) != N
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficient_type(RepT, T), rep)
#end
#function Base.convert{T, RepT<:HRepresentation}(::Type{HRepresentation{T}}, rep::RepT)
#  if fulldim(RepT) != N
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficient_type(RepT, T), rep)
#end
#function Base.convert{S, RepT<:VRepresentation}(::Type{VRepresentation{S}}, rep::RepT)
#  if fulldim(RepT) != M
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficient_type(RepT, S), rep)
#end

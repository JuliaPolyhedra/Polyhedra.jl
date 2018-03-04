export Representation, HRepresentation, VRepresentation, coefficienttype, fulldim
export RepIterator
export Rep

abstract type Representation{N, T <: Real} end
abstract type HRepresentation{N,T} <: Representation{N,T} end
abstract type VRepresentation{N,T} <: Representation{N,T} end

export Rep, HRep, VRep

const  Rep{N,T} = Union{ Representation{N,T}, Polyhedron{N,T}}
const HRep{N,T} = Union{HRepresentation{N,T}, Polyhedron{N,T}}
const VRep{N,T} = Union{VRepresentation{N,T}, Polyhedron{N,T}}

Base.copy(rep::Rep)            = error("copy not implemented for $(typeof(rep))")

"""
    coefficienttype(rep::Rep)

Returns the type of the coefficients used in the representation of `rep`.
"""
MultivariatePolynomials.coefficienttype(rep::Union{Rep{N,T}, Type{<:Rep{N,T}}}) where {N,T} = T

"""
    fulldim(rep::Rep)

Returns the dimension of the space in which the representation is defined.
That is, a straight line in a 3D space has `fulldim` 3.
"""
fulldim(p) = fulldim(FullDim(p))
FullDim(rep::Union{Rep{N}, Type{<:Rep{N}}}) where N = FullDim{N}()

# Check that it is either empty or it has a point
vconsistencyerror() = error("Non-empty V-representation must contain at least one point. If it is a polyhedral cone, the origin should be added.")
function checkvconsistency(vrep::VRep)
    if !hasallpoints(vrep) && hasallrays(vrep)
        vconsistencyerror()
    end
end
checkvconsistency(p::Polyhedron) = vrepiscomputed(p) && checkvconsistency(p)

# This method solves the ambiguity with the following methods and the general method
# Base.convert{T}(::Type{T}, p::T) = p
Base.convert{T<:HRepresentation}(::Type{T}, p::T) = p
Base.convert{T<:VRepresentation}(::Type{T}, p::T) = p

Base.convert(::Type{RepTout}, p::HRepresentation) where RepTout<:HRep = hconvert(RepTout, p)
Base.convert(::Type{RepTout}, p::HRep) where {RepTout<:HRepresentation} = hconvert(RepTout, p)
# avoid ambiguity
Base.convert(::Type{RepTout}, p::HRepresentation) where {RepTout<:HRepresentation} = hconvert(RepTout, p)

function Polyhedron{N, S}(p::Polyhedron{N, T}) where {N, S, T}
    RepTout = similar_type(typeof(p), S)
    if !hrepiscomputed(p) && vrepiscomputed(p)
        vconvert(RepTout, p)
    else
        hconvert(RepTout, p)
    end
end

#function vconvert{RepTout<:VRep, RepTin<:VRep}(::Type{RepTout}, p::RepTin)
#    Nin  = fulldim(RepTin)
#    Nout = fulldim(RepTout)
#    if Nin != Nout
#        error("Different dimension")
#    end
#    Tin  = eltype(RepTin)
#    Tout = eltype(RepTout)
#    if Tin == Tout
#        f = nothing
#    else
#        f = (i,x) -> similar_type(typeof(x), Tout)(x)
#    end
#    if decomposedvfast(p)
#        RepTout(PointIterator{Nout,Tout,Nin,Tin}([p], f), RayIterator{Nout,Tout,Nin,Tin}([p], f))
#    else
#        RepTout(VRepIterator{Nout,Tout,Nin,Tin}([p], f))
#    end
#end

Base.convert(::Type{RepTout}, p::VRepresentation) where {RepTout<:VRep} = vconvert(RepTout, p)
Base.convert(::Type{RepTout}, p::VRep) where {RepTout<:VRepresentation} = vconvert(RepTout, p)
# avoid ambiguity
Base.convert(::Type{RepTout}, p::VRepresentation) where {RepTout<:VRepresentation} = vconvert(RepTout, p)

MultivariatePolynomials.changecoefficienttype(p::Rep{N,T}, ::Type{T}) where {N,T} = p
function MultivariatePolynomials.changecoefficienttype(p::RepTin, ::Type{Tout}) where {RepTin<:Rep, Tout}
    RepTout = similar_type(RepTin, Tout)
    RepTout(p)
end

VRepresentation{N, T}(v::RepTin) where {N, T, RepTin} = similar_type(RepTin, FullDim{N}(), T)(v)
HRepresentation{N, T}(h::RepTin) where {N, T, RepTin} = similar_type(RepTin, FullDim{N}(), T)(h)

VRep{N, T}(v::VRepresentation) where {N, T} = VRepresentation{N, T}(v)
VRep{N, T}(p::Polyhedron) where {N, T} = Polyhedron{N, T}(p)
HRep{N, T}(h::HRepresentation) where {N, T} = HRepresentation{N, T}(h)
HRep{N, T}(p::Polyhedron) where {N, T} = Polyhedron{N, T}(p)

# FIXME it does not get called. The calls always go throug vconvert and hconvert. Use changecoefficienttype instead
#function Base.convert{N, T, RepT<:Representation}(::Type{Representation{N, T}}, rep::RepT)
#  if fulldim(RepT) != N
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficienttype(RepT, T), rep)
#end
#function Base.convert{N, T, RepT<:HRepresentation}(::Type{HRepresentation{N, T}}, rep::RepT)
#  if fulldim(RepT) != N
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficienttype(RepT, T), rep)
#end
#function Base.convert{M, S, RepT<:VRepresentation}(::Type{VRepresentation{M, S}}, rep::RepT)
#  if fulldim(RepT) != M
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficienttype(RepT, S), rep)
#end

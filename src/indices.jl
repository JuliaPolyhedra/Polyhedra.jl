# Index of representation element of type ElemT
"""
    Index{N, T, ElemT}

Index of an element of type `ElemT` in a `Rep{N, T}`.
"""
struct Index{N, T, ElemT}
    value::Int
end

# Type of the value associated to this index
valuetype(::Union{Index{N, T, ElemT}, Type{Index{N, T, ElemT}}}) where {N, T, ElemT} = ElemT

const HyperPlaneIndex{N, T} = Index{N, T, <:HyperPlane{N, T}}
const HalfSpaceIndex{N, T} = Index{N, T, <:HalfSpace{N, T}}
const HIndex{N, T} = Union{HyperPlaneIndex{N, T}, HalfSpaceIndex{N, T}}

#const SymPointIndex{N, T} = Index{N, T, <:SymPoint{N, T}}
const PointIndex{N, T} = Index{N, T, <:AbstractPoint{N, T}}
#const PIndex{N, T} = Union{SymPointIndex{N, T}, PointIndex{N, T}}
const PIndex{N, T} = PointIndex{N, T}
const LineIndex{N, T} = Index{N, T, <:Line{N, T}}
const RayIndex{N, T} = Index{N, T, <:Ray{N, T}}
const RIndex{N, T} = Union{LineIndex{N, T}, RayIndex{N, T}}
const VIndex{N, T} = Union{PIndex{N, T}, RIndex{N, T}}
islin(::Union{Index{N, T, ElemT}, Type{Index{N, T, ElemT}}}) where {N, T, ElemT} = islin(ElemT)
ispoint(::Union{Index{N, T, ElemT}, Type{Index{N, T, ElemT}}}) where {N, T, ElemT} = ispoint(ElemT)

"""
    Indices{N, T, ElemT, RepT<:Rep{N, T}}

Iterator over the indices of the elements of type `ElemT` of the field `rep`.
"""
struct Indices{N, T, ElemT, RepT<:Rep{N, T}}
    rep::RepT
    function Indices{N, T, ElemT}(rep) where {N, T, ElemT}
        new{N, T, ElemT, typeof(rep)}(rep)
    end
end

Base.eltype(::Indices{N, T, ElemT}) where {N, T, ElemT} = Index{N, T, ElemT}
valuetype(idxs::Indices) = valuetype(eltype(idxs))

const HyperPlaneIndices{N, T, RepT} = Indices{N, T, <:HyperPlane{N, T}, RepT}
const HalfSpaceIndices{N, T, RepT} = Indices{N, T, <:HalfSpace{N, T}, RepT}
const HIndices{N, T, RepT} = Union{HyperPlaneIndices{N, T, RepT}, HalfSpaceIndices{N, T, RepT}}

#const SymPointIndices{N, T, RepT} = Indices{N, T, <:SymPoint{N, T}, RepT}
const PointIndices{N, T, RepT} = Indices{N, T, <:AbstractPoint{N, T}, RepT}
#const PIndices{N, T, RepT} = Union{SymPointIndices{N, T, RepT}, PointIndices{N, T, RepT}}
const PIndices{N, T, RepT} = PointIndices{N, T, RepT}
const LineIndices{N, T, RepT} = Indices{N, T, <:Line{N, T}, RepT}
const RayIndices{N, T, RepT} = Indices{N, T, <:Ray{N, T}, RepT}
const RIndices{N, T, RepT} = Union{LineIndices{N, T, RepT}, RayIndices{N, T, RepT}}
const VIndices{N, T, RepT} = Union{PIndices{N, T, RepT}, RIndices{N, T, RepT}}

function Base.next(idxs::Indices{N, T, ElemT}, idx::Index{N, T, ElemT}) where {N, T, ElemT}
    nextidx = nextindex(idxs.rep, idx)
    idx, nextidx
end

repfor(p, ::Type{<:HRepElement}) = hrep(p)
repfor(p, ::Type{<:VRepElement}) = vrep(p)
Base.length(idxs::Indices{N, T, ElemT, <:Polyhedron{N, T}}) where {N, T, ElemT} = length(Indices{N, T, ElemT}(repfor(idxs.rep, ElemT)))
Base.isempty(idxs::Indices{N, T, ElemT, <:Polyhedron{N, T}}) where {N, T, ElemT} = isempty(Indices{N, T, ElemT}(repfor(idxs.rep, ElemT)))
Base.start(idxs::Indices{N, T, ElemT, <:Polyhedron{N, T}}) where {N, T, ElemT} = start(Indices{N, T, ElemT}(repfor(idxs.rep, ElemT)))
Base.done(idxs::Indices{N, T, ElemT, <:Polyhedron{N, T}}, idx::Index{N, T, ElemT}) where {N, T, ElemT} = done(Indices{N, T, ElemT}(repfor(idxs.rep, ElemT)), idx)
Base.get(p::Polyhedron{N, T}, idx::Index{N, T, ElemT}) where {N, T, ElemT} = get(repfor(p, ElemT), idx)
nextindex(p::Polyhedron{N, T}, idx::Index{N, T, ElemT}) where {N, T, ElemT} = nextindex(repfor(p, ElemT), idx)

"""
The representation `rep` does not contain any `elem`.
"""
macro norepelem(rep, elem)
    idxs = Symbol(string(elem) * "Indices")
    idx = Symbol(string(elem) * "Index")
    quote
        Base.length(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = 0
        Base.isempty(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = true
        Base.start(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = eltype(idxs)(0)
        Base.done(idxs::$idxs{N, T, <:$rep{N, T}}, ::$idx{N, T}) where {N, T} = true
    end
end

abstract type HAffineSpace{N, T} <: HRepresentation{N, T} end
@norepelem HAffineSpace HalfSpace

_promote_reptype(P1::Type{<:HAffineSpace}, ::Type{<:HAffineSpace}) = P1
_promote_reptype(P1::Type{<:HAffineSpace}, ::Type{<:HRep}) = hreptype(P1)

abstract type VPolytope{N, T} <: VRepresentation{N, T} end
@norepelem VPolytope Line
@norepelem VPolytope Ray

_promote_reptype(P1::Type{<:VPolytope}, ::Type{<:VPolytope}) = P1
_promote_reptype(P1::Type{<:VPolytope}, ::Type{<:VRep}) = vreptype(P1)

abstract type VSymPolytope{N, T} <: VPolytope{N, T} end
@norepelem VSymPolytope Point

abstract type VCone{N, T} <: VRepresentation{N, T} end
#@norepelem VCone SymPoint
# See issue #28
Base.length(idxs::PointIndices{N, T, <:VCone{N, T}}) where {N, T} = hasallrays(idxs.rep) ? 1 : 0
Base.isempty(idxs::PointIndices{N, T, <:VCone{N, T}}) where {N, T} = !hasallrays(idxs.rep)
Base.start(idxs::PointIndices{N, T, <:VCone{N, T}}) where {N, T} = eltype(idxs)(hasallrays(idxs.rep) ? 1 : 2)
Base.done(::PointIndices{N, T, <:VCone{N, T}}, idx::PointIndex{N, T}) where {N, T} = idx.value > 1
Base.get(L::VCone{N, T}, ::PointIndex{N, T}) where {N, T} = origin(arraytype(L), FullDim{N}())
nextindex(::VCone{N, T}, idx::PointIndex{N, T}) where {N, T} = typeof(idx)(idx.value + 1)

_promote_reptype(P1::Type{<:VCone}, ::Type{<:VCone}) = P1
_promote_reptype(P1::Type{<:VCone}, ::Type{<:VRep}) = vreptype(P1)

abstract type VLinearSpace{N, T} <: VCone{N, T} end
@norepelem VLinearSpace Ray

_promote_reptype(P1::Type{<:VLinearSpace}, ::Type{<:VLinearSpace}) = P1
_promote_reptype(P1::Type{<:VLinearSpace}, ::Type{<:VCone}) = conetype(P1)
_promote_reptype(P1::Type{<:VLinearSpace}, ::Type{<:VRep}) = vreptype(P1)

"""
The representation `rep` contain the elements `elem` inside a vector in the field `field`.
"""
macro vecrepelem(rep, elem, field)
    idxs = Symbol(string(elem) * "Indices")
    idx = Symbol(string(elem) * "Index")
    esc(quote
        Base.length(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = length(idxs.rep.$field)
        Base.isempty(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = isempty(idxs.rep.$field)
        Base.start(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = eltype(idxs)(1)
        Base.done(idxs::$idxs{N, T, <:$rep{N, T}}, idx::$idx{N, T}) where {N, T} = idx.value > length(idxs)
        Base.get(rep::$rep{N, T}, idx::$idx{N, T}) where {N, T} = rep.$field[idx.value]
        nextindex(::$rep{N, T}, idx::$idx{N, T}) where {N, T} = typeof(idx)(idx.value + 1)
    end)
end

"""
The representation `rep` contain the elements `elem` inside a representation in the field `field`.
"""
macro subrepelem(rep, elem, field)
    idxst = Symbol(string(elem) * "Indices")
    idxs = :(Polyhedra.$idxst)
    idxt = Symbol(string(elem) * "Index")
    idx = :(Polyhedra.$idxt)
    subidxs = :(Polyhedra.Indices{N, T, Polyhedra.valuetype(idxs)}(idxs.rep.$field))
    esc(quote
        Base.length(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = length($subidxs)
        Base.isempty(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = isempty($subidxs)
        Base.start(idxs::$idxs{N, T, <:$rep{N, T}}) where {N, T} = start($subidxs)
        Base.done(idxs::$idxs{N, T, <:$rep{N, T}}, idx::$idx{N, T}) where {N, T} = done($subidxs, idx)
        Base.get(rep::$rep{N, T}, idx::$idx{N, T}) where {N, T} = get(rep.$field, idx)
        Polyhedra.nextindex(rep::$rep{N, T}, idx::$idx{N, T}) where {N, T} = Polyhedra.nextindex(rep.$field, idx)
    end)
end

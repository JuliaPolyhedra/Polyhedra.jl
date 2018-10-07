# Index of representation element of type ElemT
"""
    Index{T, ElemT}

Index of an element of type `ElemT` in a `Rep{T}`.
"""
struct Index{T, ElemT}
    value::Int
end

# Type of the value associated to this index
valuetype(::Union{Index{T, ElemT}, Type{Index{T, ElemT}}}) where {T, ElemT} = ElemT

const HyperPlaneIndex{T} = Index{T, <:HyperPlane{T}}
const HalfSpaceIndex{T} = Index{T, <:HalfSpace{T}}
const HIndex{T} = Union{HyperPlaneIndex{T}, HalfSpaceIndex{T}}

"""
    Indices{T, ElemT, RepT<:Rep{T}}

Iterator over the indices of the elements of type `ElemT` of the field `rep`.
"""
struct Indices{T, ElemT, RepT<:Rep{T}}
    rep::RepT
    function Indices{T, ElemT}(rep) where {T, ElemT}
        new{T, ElemT, typeof(rep)}(rep)
    end
end

Base.eltype(::Indices{T, ElemT}) where {T, ElemT} = Index{T, ElemT}
valuetype(idxs::Indices) = valuetype(eltype(idxs))

const HyperPlaneIndices{T, RepT} = Indices{T, <:HyperPlane{T}, RepT}
const HalfSpaceIndices{T, RepT} = Indices{T, <:HalfSpace{T}, RepT}
const HIndices{T, RepT} = Union{HyperPlaneIndices{T, RepT}, HalfSpaceIndices{T, RepT}}

undouble_it(idx::Nothing) = nothing
double_it(idx::Nothing) = nothing
undouble_it(idx::NTuple{2, Index}) = idx[1]
double_it(idx::Index) = idx, idx
function Base.iterate(idxs::Indices)
    return double_it(startindex(idxs))
end
function Base.iterate(idxs::Indices{T, ElemT},
                      idx::Index{T, ElemT}) where {T, ElemT}
    return double_it(nextindex(idxs.rep, idx))
end

# For polyhedron, redirect to hrep or vrep depending on whether it is an
# HRepElement or VRepElement
repfor(p::Polyhedron, ::Type{<:HRepElement}) = hrep(p)
function startindex(idxs::Indices{T, ElemT,
                                  <:Polyhedron{T}}) where {T, ElemT}
    return startindex(Indices{T, ElemT}(repfor(idxs.rep, ElemT)))
end
function nextindex(p::Polyhedron{T}, idx::Index{T, ElemT}) where {T, ElemT}
    return nextindex(repfor(p, ElemT), idx)
end
function Base.length(idxs::Indices{T, ElemT, <:Polyhedron{T}}) where {T, ElemT}
    return length(Indices{T, ElemT}(repfor(idxs.rep, ElemT)))
end
function Base.isempty(idxs::Indices{T, ElemT, <:Polyhedron{T}}) where {T, ElemT}
    return isempty(Indices{T, ElemT}(repfor(idxs.rep, ElemT)))
end
function Base.get(p::Polyhedron{T}, idx::Index{T, ElemT}) where {T, ElemT}
    return get(repfor(p, ElemT), idx)
end

"""
The representation `rep` does not contain any `elem`.
"""
macro norepelem(rep, elem)
    idxs = Symbol(string(elem) * "Indices")
    idx = Symbol(string(elem) * "Index")
    quote
        Base.length(idxs::$idxs{T, <:$rep{T}}) where {T} = 0
        Base.isempty(idxs::$idxs{T, <:$rep{T}}) where {T} = true
        Base.iterate(idxs::$idxs{T, <:$rep{T}}) where {T} = nothing
    end
end

abstract type HAffineSpace{T} <: HRepresentation{T} end
@norepelem HAffineSpace HalfSpace

_promote_reptype(P1::Type{<:HAffineSpace}, ::Type{<:HAffineSpace}) = P1
_promote_reptype(P1::Type{<:HAffineSpace}, ::Type{<:HRep}) = hreptype(P1)

"""
The representation `rep` contain the elements `elem` inside a vector in the field `field`.
"""
macro vecrepelem(rep, elem, field)
    idxs = Symbol(string(elem) * "Indices")
    idx = Symbol(string(elem) * "Index")
    esc(quote
        Base.length(idxs::$idxs{T, <:$rep{T}}) where {T} = length(idxs.rep.$field)
        Base.isempty(idxs::$idxs{T, <:$rep{T}}) where {T} = isempty(idxs.rep.$field)
        function Polyhedra.startindex(idxs::$idxs{T, <:$rep{T}}) where {T}
            if isempty(idxs.rep.$field)
                return nothing
            else
                return eltype(idxs)(1)
            end
        end
        Base.get(rep::$rep{T}, idx::$idx{T}) where {T} = rep.$field[idx.value]
        function Polyhedra.nextindex(rep::$rep{T}, idx::$idx{T}) where {T}
            if idx.value >= length(rep.$field)
                return nothing
            else
                return typeof(idx)(idx.value + 1)
            end
        end
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
    subidxs = :(Polyhedra.Indices{T, Polyhedra.valuetype(idxs)}(idxs.rep.$field))
    esc(quote
        Base.length(idxs::$idxs{T, <:$rep{T}}) where {T} = length($subidxs)
        Base.isempty(idxs::$idxs{T, <:$rep{T}}) where {T} = isempty($subidxs)
        Polyhedra.startindex(idxs::$idxs{T, <:$rep{T}}) where {T} = Polyhedra.startindex($subidxs)
        Base.get(rep::$rep{T}, idx::$idx{T}) where {T} = get(rep.$field, idx)
        Polyhedra.nextindex(rep::$rep{T}, idx::$idx{T}) where {T} = Polyhedra.nextindex(rep.$field, idx)
    end)
end

# Index of representation element of type ElemT
"""
    RepElementIndex{N, T, ElemT}

Index of an element of type `ElemT` in a `Rep{N, T}`.
"""
struct RepElementIndex{N, T, ElemT}
    value::Int
end

# Type of the value associated to this index
valuetype(::RepElementIndex{N, T, ElemT}) where {N, T, ElemT} = ElemT

const HyperPlaneIndex{N, T} = RepElementIndex{N, T, <:HyperPlane{N, T}}
const HalfSpaceIndex{N, T} = RepElementIndex{N, T, <:HalfSpace{N, T}}
const HIndex{N, T} = Union{HyperPlaneIndex{N, T}, HalfSpaceIndex{N, T}}

const SymPointIndex{N, T} = RepElementIndex{N, T, <:SymPoint{N, T}}
const PointIndex{N, T} = RepElementIndex{N, T, <:MyPoint{N, T}}
const PIndex{N, T} = Union{SymPointIndex{N, T}, PointIndex{N, T}}
const LineIndex{N, T} = RepElementIndex{N, T, <:Line{N, T}}
const RayIndex{N, T} = RepElementIndex{N, T, <:Ray{N, T}}
const RIndex{N, T} = Union{LineIndex{N, T}, RayIndex{N, T}}
const VIndex{N, T} = Union{PIndex{N, T}, RIndex{N, T}}
islin(::Union{RepElementIndex{N, T, ElemT}, Type{RepElementIndex{N, T, ElemT}}}) where {N, T, ElemT} = islin(ElemT)
ispoint(::Union{RepElementIndex{N, T, ElemT}, Type{RepElementIndex{N, T, ElemT}}}) where {N, T, ElemT} = ispoint(ElemT)

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

Base.eltype(::Indices{N, T, ElemT}) where {N, T, ElemT} = RepElementIndex{N, T, ElemT}

const HyperPlaneIndices{N, T, RepT} = Indices{N, T, <:HyperPlane{N, T}, RepT}
const HalfSpaceIndices{N, T, RepT} = Indices{N, T, <:HalfSpace{N, T}, RepT}
const HIndices{N, T, RepT} = Union{HyperPlaneIndices{N, T, RepT}, HalfSpaceIndices{N, T, RepT}}

const SymPointIndices{N, T, RepT} = Indices{N, T, <:SymPoint{N, T}, RepT}
const PointIndices{N, T, RepT} = Indices{N, T, <:MyPoint{N, T}, RepT}
const PIndices{N, T, RepT} = Union{SymPointIndices{N, T, RepT}, PointIndices{N, T, RepT}}
const LineIndices{N, T, RepT} = Indices{N, T, <:Line{N, T}, RepT}
const RayIndices{N, T, RepT} = Indices{N, T, <:Ray{N, T}, RepT}
const RIndices{N, T, RepT} = Union{LineIndices{N, T, RepT}, RayIndices{N, T, RepT}}
const VIndices{N, T, RepT} = Union{PIndices{N, T, RepT}, RIndices{N, T, RepT}}

function Base.next(idx::Indices{N, T, ElemT}, state::RepElementIndex{N, T, ElemT}) where {N, T, ElemT}
    nextidx = nextindex(idx.rep, state)
    nextidx, nextidx
end

repfor(p, ::Type{<:HRepElement}) = hrep(p)
repfor(p, ::Type{<:VRepElement}) = vrep(p)
Base.length(idxs::Indices{N, T, ElemT, <:Polyhedron{N, T}}) where {N, T, ElemT} = length(Indices{N, T, ElemT}(repfor(idxs.rep, ElemT)))
Base.isempty(idxs::Indices{N, T, ElemT, <:Polyhedron{N, T}}) where {N, T, ElemT} = isempty(Indices{N, T, ElemT}(repfor(idxs.rep, ElemT)))
Base.start(idxs::Indices{N, T, ElemT, <:Polyhedron{N, T}}) where {N, T, ElemT} = start(Indices{N, T, ElemT}(repfor(idxs.rep, ElemT)))
Base.done(idxs::Indices{N, T, ElemT, <:Polyhedron{N, T}}, idx::RepElementIndex{N, T, ElemT}) where {N, T, ElemT} = done(Indices{N, T, ElemT}(repfor(idxs.rep, ElemT)), idx)
Base.get(p::Polyhedron{N, T}, idx::RepElementIndex{N, T, ElemT}) where {N, T, ElemT} = get(repfor(p, ElemT), idx)
nextindex(p::Polyhedron{N, T}, idx::RepElementIndex{N, T, ElemT}) where {N, T, ElemT} = nextindex(repfor(p, ElemT), idx)

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
        Base.get(L::$rep{N, T}, idx::$idx{N, T}) where {N, T} = L.$field[idx.value]
        nextindex(L::$rep{N, T}, idx::$idx{N, T}) where {N, T} = typeof(idx)(idx.value + 1)
    end)
end

isincident(v::VRepElement, h::HRepElement; tol) = in(v, hyperplane(h); tol)
isincident(h::HRepElement, v::VRepElement; tol) = in(v, hyperplane(h); tol)
isincident(p::Polyhedron, v::Index, h::Index; tol) = isincident(get(p, v), get(p, h); tol)

abstract type Incident{T, ElemT<:RepElement{T}, PT<:Polyhedron{T}, IdxT<:Index{T}} end

function Base.get(p::Polyhedron{T}, inc::Incident{T, ElemT}; tol) where {T, ElemT}
    el = get(p, inc.idx)
    incs = _inctype(inc)[]
    for idx in Indices{T, ElemT}(p)
        eli = get(p, idx)
        if isincident(el, eli; tol)
            push!(incs, _incel(inc, idx, eli))
        end
    end
    return incs
end

struct IncidentElements{T, ElemT, PT, IdxT} <: Incident{T, ElemT, PT, IdxT}
    p::PT
    idx::IdxT
end
IncidentElements{T, ElemT}(p, idx) where {T, ElemT<:RepElement{T}} = IncidentElements{T, ElemT, typeof(p), typeof(idx)}(p, idx)
_inctype(inc::IncidentElements{T, ElemT}) where {T, ElemT} = ElemT
_incel(inc::IncidentElements, idx, eli) = eli

struct IncidentIndices{T, ElemT, PT, IdxT} <: Incident{T, ElemT, PT, IdxT}
    p::PT
    idx::IdxT
end
IncidentIndices{T, ElemT}(p, idx) where {T, ElemT<:RepElement{T}} = IncidentIndices{T, ElemT, typeof(p), typeof(idx)}(p, idx)
_inctype(inc::IncidentIndices{T, ElemT}) where {T, ElemT} = Index{T, ElemT}
_incel(inc::IncidentIndices, idx, eli) = idx

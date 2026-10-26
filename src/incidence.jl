isincident(v::VRepElement, h::HRepElement; tol) = in(v, hyperplane(h); tol)
isincident(h::HRepElement, v::VRepElement; tol) = in(v, hyperplane(h); tol)
function isincident(p::Polyhedron, v::Index, h::Index; tol)
    return isincident(get(p, v), get(p, h); tol)
end

abstract type Incident{T,ElemT<:RepElement{T},PT<:Polyhedron{T},IdxT<:Index{T}} end

function Base.get(
    p::Polyhedron{T},
    inc::Incident{T,ElemT};
    tol = _default_tol(T),
) where {T,ElemT}
    el = get(p, inc.idx)
    incs = _inctype(inc)[]
    for idx in Indices{T,ElemT}(p)
        eli = get(p, idx)
        if isincident(el, eli; tol)
            push!(incs, _incel(inc, idx, eli))
        end
    end
    return incs
end

struct IncidentElements{T,ElemT,PT,IdxT} <: Incident{T,ElemT,PT,IdxT}
    p::PT
    idx::IdxT
end
function IncidentElements{T,ElemT}(p, idx) where {T,ElemT<:RepElement{T}}
    return IncidentElements{T,ElemT,typeof(p),typeof(idx)}(p, idx)
end
_inctype(inc::IncidentElements{T,ElemT}) where {T,ElemT} = ElemT
_incel(inc::IncidentElements, idx, eli) = eli

struct IncidentIndices{T,ElemT,PT,IdxT} <: Incident{T,ElemT,PT,IdxT}
    p::PT
    idx::IdxT
end
function IncidentIndices{T,ElemT}(p, idx) where {T,ElemT<:RepElement{T}}
    return IncidentIndices{T,ElemT,typeof(p),typeof(idx)}(p, idx)
end
_inctype(inc::IncidentIndices{T,ElemT}) where {T,ElemT} = Index{T,ElemT}
_incel(inc::IncidentIndices, idx, eli) = idx

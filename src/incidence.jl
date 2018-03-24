isincident(v::VRepElement, h::HRepElement) = v in hyperplane(h)
isincident(h::HRepElement, v::VRepElement) = v in hyperplane(h)

abstract type Incident{N, T, ElemT<:Element{N, T}, PT<:Polyhedron{N, T}, IdxT<:Index{N, T}} end

function Base.get(p::Polyhedron{N, T}, inc::Incident{N, T, ElemT}) where {N, T, ElemT}
    el = get(p, inc.idx)
    incs = _inctype(inc)[]
    for idx in Indices{N, T, ElemT}(p)
        eli = get(p, idx)
        if isincident(el, eli)
            push!(incs, _incel(inc, idx, eli))
        end
    end
    incs
end

struct IncidentElements{N, T, ElemT, PT, IdxT} <: Incident{N, T, ElemT, PT, IdxT}
    p::PT
    idx::IdxT
end
IncidentElements{N, T, ElemT}(p, idx) where {N, T, ElemT<:Element{N, T}} = IncidentElements{N, T, ElemT, typeof(p), typeof(idx)}(p, idx)
_inctype(inc::IncidentElements{N, T, ElemT}) where {N, T, ElemT} = ElemT
_incel(inc::IncidentElements, idx, eli) = eli

struct IncidentIndices{N, T, ElemT, PT, IdxT} <: Incident{N, T, ElemT, PT, IdxT}
    p::PT
    idx::IdxT
end
IncidentIndices{N, T, ElemT}(p, idx) where {N, T, ElemT<:Element{N, T}} = IncidentIndices{N, T, ElemT, typeof(p), typeof(idx)}(p, idx)
_inctype(inc::IncidentIndices{N, T, ElemT}) where {N, T, ElemT} = Index{N, T, ElemT}
_incel(inc::IncidentIndices, idx, eli) = idx

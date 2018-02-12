# Useful for CDD, LRS and matrix representations

# HyperPlanes and HalfSpaces mixed
abstract type MixedHRep{N, T} <: HRepresentation{N, T} end
# SymPoints and Points mixed and/or Lines and Rays mixed
abstract type MixedVRep{N, T} <: VRepresentation{N, T} end

const MixedRep{N, T} = Union{MixedHRep{N, T}, MixedVRep{N, T}}

# This could be defined for all types Indices but it is not recommended for reps other than Mixed to use this implementation of length as it is inefficient
# a MethodError is more helpful than a hidden inefficiency
function Base.length(idxs::Indices{N, T, ElemT, <:MixedRep{N, T}}) where {N, T, ElemT}
    count = 0
    for idx in idxs
        count += 1
    end
    count
end

function mixednext(rep::MixedRep{N, T}, idx::IdxT) where {N, T, ElemT, IdxT<:Index{N, T, ElemT}}
    idx = IdxT(idx.value+1)
    while !done(Indices{N, T, ElemT}(rep), idx) && !isvalid(rep, idx)
        idx = IdxT(idx.value+1)
    end
    idx
end
Base.start(idx::Indices{N, T, ElemT, <:MixedRep{N, T}}) where {N, T, ElemT} = mixednext(idx.rep, Index{N, T, ElemT}(0))
nextindex(rep::MixedRep{N, T}, idx::Index{N, T}) where {N, T} = mixednext(rep, idx)

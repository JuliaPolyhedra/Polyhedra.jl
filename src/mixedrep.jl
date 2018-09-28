# Useful for CDD, LRS and matrix representations

# HyperPlanes and HalfSpaces mixed
abstract type MixedHRep{T} <: HRepresentation{T} end
# SymPoints and Points mixed and/or Lines and Rays mixed
abstract type MixedVRep{T} <: VRepresentation{T} end

const MixedRep{T} = Union{MixedHRep{T}, MixedVRep{T}}

# This could be defined for all types Indices but it is not recommended for reps other than Mixed to use this implementation of length as it is inefficient
# a MethodError is more helpful than a hidden inefficiency
function mixedlength(idxs::Indices)
    count = 0
    for idx in idxs
        count += 1
    end
    count
end
Base.length(idxs::Indices{T, ElemT, <:MixedRep{T}}) where {T, ElemT} = mixedlength(idxs)

function mixednext(rep::MixedRep{T}, idx::IdxT) where {T, ElemT, IdxT<:Index{T, ElemT}}
    idx = IdxT(idx.value+1)
    while !done(Indices{T, ElemT}(rep), idx) && !isvalid(rep, idx)
        idx = IdxT(idx.value+1)
    end
    if done(Indices{T, ElemT}(rep), idx)
        return nothing
    else
        return idx
    end
end
startindex(idx::Indices{T, ElemT, <:MixedRep{T}}) where {T, ElemT} = mixednext(idx.rep, Index{T, ElemT}(0))
nextindex(rep::MixedRep{T}, idx::Index{T}) where {T} = mixednext(rep, idx)

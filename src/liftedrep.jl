export LiftedHRepresentation, LiftedVRepresentation

# H-Represenation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LiftedHRepresentation{N, T} <: MixedHRep{N, T}
    # Ax >= 0, it is [b -A] * [z; x] where z = 1
    A::AbstractMatrix{T}
    linset::IntSet

    function LiftedHRepresentation{N, T}(A::AbstractMatrix, linset::IntSet=IntSet()) where {N, T}
        if !isempty(linset) && last(linset) > size(A, 1)
            error("The elements of linset should be between 1 and the number of rows of A")
        end
        if size(A, 2) != N+1
            error("dimension does not match")
        end
        new{N, T}(A, linset)
    end
end

similar_type{N,T}(::Type{<:LiftedHRepresentation}, ::FullDim{N}, ::Type{T}) = LiftedHRepresentation{N,T}
arraytype(p::Union{LiftedHRepresentation{N, T}, Type{LiftedHRepresentation{N, T}}}) where {N, T} = Vector{T}

linset(rep::LiftedHRepresentation) = copy(rep.linset)

LiftedHRepresentation(A::AbstractMatrix{T}, linset::IntSet=IntSet()) where {T <: Real} = LiftedHRepresentation{size(A,2)-1,T}(A, linset)
LiftedHRepresentation(h::HRepresentation{N,T}) where {N,T} = LiftedHRepresentation{N,T}(h)

function LiftedHRepresentation{N, T}(hyperplanes::ElemIt{<:HyperPlane{N, T}}, halfspaces::ElemIt{<:HalfSpace{N, T}}) where {N, T}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = Matrix{T}(nhrep, N+1)
    linset = IntSet(1:nhyperplane)
    for (i, h) in enumerate(hyperplanes)
        A[i,2:end] = -h.a
        A[i,1] = h.β
    end
    for (i, h) in enumerate(halfspaces)
        A[nhyperplane+i,2:end] = -h.a
        A[nhyperplane+i,1] = h.β
    end
    LiftedHRepresentation{N, T}(A, linset)
end

Base.copy(ine::LiftedHRepresentation{N,T}) where {N,T} = LiftedHRepresentation{N,T}(copy(ine.A), copy(ine.linset))

Base.isvalid(hrep::LiftedHRepresentation{N, T}, idx::HIndex{N, T}) where {N, T} = 0 < idx.value <= size(hrep.A, 1) && (idx.value in hrep.linset) == islin(idx)
Base.done(idxs::HIndices{N, T, <:LiftedHRepresentation{N, T}}, idx::HIndex{N, T}) where {N, T} = idx.value > size(idxs.rep.A, 1)
Base.get(hrep::LiftedHRepresentation{N, T}, idx::HIndex{N, T}) where {N, T} = valuetype(idx)(-hrep.A[idx.value,2:end], hrep.A[idx.value,1])

Base.getindex(h::LiftedHRepresentation, I::AbstractArray) = LiftedHRepresentation(h.A[I, :], filterintset(h.linset, I))

# V-Represenation

mutable struct LiftedVRepresentation{N,T} <: MixedVRep{N,T}
    R::AbstractMatrix{T} # each row is a vertex if the first element is 1 and a ray otherwise
    linset::IntSet

    function LiftedVRepresentation{N, T}(R::AbstractMatrix, linset::IntSet=IntSet([])) where {N, T}
        if length(R) > 0 && size(R, 2) != N+1
            error("dimension does not match")
        end
        if !isempty(linset) && last(linset) > size(R, 1)
            error("The elements of linset should be between 1 and the number of rows of R")
        end
        new{N, T}(R, linset)
    end
end

similar_type{N,T}(::Type{<:LiftedVRepresentation}, ::FullDim{N}, ::Type{T}) = LiftedVRepresentation{N,T}
arraytype(p::Union{LiftedVRepresentation{N, T}, Type{LiftedVRepresentation{N, T}}}) where {N, T} = Vector{T}

linset(rep::LiftedVRepresentation) = copy(rep.linset)

LiftedVRepresentation(R::AbstractMatrix{T}, linset::IntSet=IntSet()) where {T <: Real} = LiftedVRepresentation{size(R,2)-1,T}(R, linset)
LiftedVRepresentation(v::VRepresentation{N,T}) where {N,T} = LiftedVRepresentation{N,T}(v)

function LiftedVRepresentation{N, T}(vits::VIt{N, T}...) where {N, T}
    sympoints, points, lines, rays = fillvits(vits...)
    nsympoint = length(sympoints)
    npoint = length(points)
    nline = length(lines)
    nray = length(rays)
    nvrep = nsympoint + npoint + nline + nray
    R = Matrix{T}(nvrep, N+1)
    linset = IntSet()
    function _fill(offset, z, ps)
        for (i, p) in enumerate(ps)
            R[offset + i,2:end] = coord(p)
            R[offset + i,1] = z
            if islin(p)
                push!(linset, offset + i)
            end
        end
    end
    _fill(0, one(T), sympoints)
    _fill(nsympoint, one(T), points)
    _fill(nsympoint+npoint, zero(T), lines)
    _fill(nsympoint+npoint+nline, zero(T), rays)
    LiftedVRepresentation{N, T}(R, linset)
end

Base.copy(ext::LiftedVRepresentation{N,T}) where {N,T} = LiftedVRepresentation{N,T}(copy(ext.R), copy(ext.linset))

nvreps(ext::LiftedVRepresentation) = size(ext.R, 1)

function isrowpoint(ext::LiftedVRepresentation{N,T}, i) where {N,T}
    ispoint = ext.R[i,1]
    @assert ispoint == zero(T) || ispoint == one(T)
    ispoint == one(T)
end

function Base.isvalid(vrep::LiftedVRepresentation{N, T}, idx::VIndex{N, T}) where {N, T}
    0 < idx.value <= nvreps(vrep) && isrowpoint(vrep, idx.value) == ispoint(idx) && (idx.value in vrep.linset) == islin(idx)
end
Base.done(idxs::VIndices{N, T, <:LiftedVRepresentation{N, T}}, idx::VIndex{N, T}) where {N, T} = idx.value > size(idxs.rep.R, 1)
Base.get(vrep::LiftedVRepresentation{N, T}, idx::VIndex{N, T}) where {N, T} = valuetype(idx)(vrep.R[idx.value,2:end])

Base.getindex(v::LiftedVRepresentation, I::AbstractArray) = LiftedVRepresentation(v.R[I, :], filterintset(v.linset, I))

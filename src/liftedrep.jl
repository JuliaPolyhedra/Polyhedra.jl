export LiftedHRepresentation, LiftedVRepresentation

# H-Represenation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LiftedHRepresentation{N, T, MT<:AbstractMatrix{T}} <: MixedHRep{N, T}
    # Ax >= 0, it is [b -A] * [z; x] where z = 1
    A::MT
    linset::IntSet

    function LiftedHRepresentation{N, T, MT}(A::MT, linset::IntSet=IntSet()) where {N, T, MT}
        if !isempty(linset) && last(linset) > size(A, 1)
            error("The elements of linset should be between 1 and the number of rows of A")
        end
        if size(A, 2) != N+1
            error("dimension does not match")
        end
        new{N, T, MT}(A, linset)
    end
end

similar_type(::Type{LiftedHRepresentation{M, S, MT}}, ::FullDim{N}, ::Type{T}) where {M, S, N, T, MT} = LiftedHRepresentation{N, T, similar_type(MT, T)}
hvectortype(p::Type{LiftedHRepresentation{N, T, MT}}) where {N, T, MT} = vectortype(MT)

LiftedHRepresentation{N, T}(A::AbstractMatrix{T}, linset::IntSet=IntSet()) where {N, T} = LiftedHRepresentation{N, T, typeof(A)}(A, linset)
LiftedHRepresentation{N, T}(A::AbstractMatrix, linset::IntSet=IntSet()) where {N, T} = LiftedHRepresentation{N, T}(AbstractMatrix{T}(A), linset)
LiftedHRepresentation(A::AbstractMatrix{T}, linset::IntSet=IntSet()) where T = LiftedHRepresentation{size(A, 2) - 1, T}(A, linset)
function LiftedHRepresentation(h::HRepresentation{N, T}) where {N, T}
    LiftedHRepresentation{N, T, hvectortype(typeof(h)) <: AbstractSparseVector ? SparseMatrixCSC{T, Int} : Matrix{T}}(h)
end

function LiftedHRepresentation{N, T, MT}(hyperplanes::ElemIt{<:HyperPlane{N, T}}, halfspaces::ElemIt{<:HalfSpace{N, T}}) where {N, T, MT}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = emptymatrix(MT, nhrep, N+1)
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

# V-Representation

mutable struct LiftedVRepresentation{N, T, MT<:AbstractMatrix{T}} <: MixedVRep{N, T}
    R::MT # each row is a vertex if the first element is 1 and a ray otherwise
    linset::IntSet

    function LiftedVRepresentation{N, T, MT}(R::MT, linset::IntSet=IntSet([])) where {N, T, MT}
        if length(R) > 0 && size(R, 2) != N+1
            error("dimension does not match")
        end
        if !isempty(linset) && last(linset) > size(R, 1)
            error("The elements of linset should be between 1 and the number of rows of R")
        end
        new{N, T, MT}(R, linset)
    end
end

similar_type(::Type{LiftedVRepresentation{M, S, MT}}, ::FullDim{N}, ::Type{T}) where {M, S, N, T, MT} = LiftedVRepresentation{N, T, similar_type(MT, T)}
vvectortype(p::Type{LiftedVRepresentation{N, T, MT}}) where {N, T, MT} = vectortype(MT)

LiftedVRepresentation{N, T}(R::AbstractMatrix{T}, linset::IntSet=IntSet()) where {N, T} = LiftedVRepresentation{N, T, typeof(R)}(R, linset)
LiftedVRepresentation{N, T}(R::AbstractMatrix, linset::IntSet=IntSet()) where {N, T} = LiftedVRepresentation{N, T}(AbstractMatrix{T}(R), linset)
LiftedVRepresentation(R::AbstractMatrix{T}, linset::IntSet=IntSet()) where T = LiftedVRepresentation{size(R, 2) - 1, T}(R, linset)
LiftedVRepresentation(v::VRepresentation{N,T}) where {N, T} = LiftedVRepresentation{N, T, Matrix{T}}(v)

function LiftedVRepresentation{N, T, MT}(vits::VIt{N, T}...) where {N, T, MT}
    points, lines, rays = fillvits(FullDim{N}(), vits...)
    npoint = length(points)
    nline = length(lines)
    nray = length(rays)
    nvrep = npoint + nline + nray
    R = emptymatrix(MT, nvrep, N+1)
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
    _fill(0, one(T), points)
    _fill(npoint, zero(T), lines)
    _fill(npoint+nline, zero(T), rays)
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
    isp = isrowpoint(vrep, idx.value)
    isl = (idx.value in vrep.linset)
    @assert !isp || !isl # if isp && isl, it is a symmetric point but it is not allowed to mix symmetric points and points
    0 < idx.value <= nvreps(vrep) && isp == ispoint(idx) && isl == islin(idx)
end
Base.done(idxs::VIndices{N, T, <:LiftedVRepresentation{N, T}}, idx::VIndex{N, T}) where {N, T} = idx.value > size(idxs.rep.R, 1)
Base.get(vrep::LiftedVRepresentation{N, T}, idx::VIndex{N, T}) where {N, T} = valuetype(idx)(vrep.R[idx.value,2:end])

export LiftedHRepresentation, LiftedVRepresentation

# H-Represenation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LiftedHRepresentation{T, MT<:AbstractMatrix{T}} <: MixedHRep{T}
    # Ax >= 0, it is [b -A] * [z; x] where z = 1
    A::MT
    linset::BitSet

    function LiftedHRepresentation{T, MT}(A::MT, linset::BitSet=BitSet()) where {T, MT}
        if !isempty(linset) && last(linset) > size(A, 1)
            error("The elements of linset should be between 1 and the number of rows of A")
        end
        new{T, MT}(A, linset)
    end
end
FullDim(rep::LiftedHRepresentation) = size(rep.A, 2) - 1

similar_type(::Type{LiftedHRepresentation{S, MT}}, ::FullDim, ::Type{T}) where {S, T, MT} = LiftedHRepresentation{T, similar_type(MT, T)}
hvectortype(p::Type{LiftedHRepresentation{T, MT}}) where {T, MT} = vectortype(MT)

LiftedHRepresentation{T}(A::AbstractMatrix{T}, linset::BitSet=BitSet()) where {T} = LiftedHRepresentation{T, typeof(A)}(A, linset)
LiftedHRepresentation{T}(A::AbstractMatrix, linset::BitSet=BitSet()) where {T} = LiftedHRepresentation{T}(AbstractMatrix{T}(A), linset)
LiftedHRepresentation(A::AbstractMatrix{T}, linset::BitSet=BitSet()) where T = LiftedHRepresentation{T}(A, linset)

LiftedHRepresentation(h::HRepresentation{T}) where {T} = LiftedHRepresentation{T}(h)
LiftedHRepresentation{T}(h::HRepresentation) where {T} = convert(LiftedHRepresentation{T, hmatrixtype(typeof(h), T)}, h)

function LiftedHRepresentation{T, MT}(d::FullDim, hyperplanes::HyperPlaneIt{T},
                                      halfspaces::HalfSpaceIt{T}) where {T, MT}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = emptymatrix(MT, nhrep, fulldim(d)+1)
    linset = BitSet(1:nhyperplane)
    for (i, h) in enumerate(hyperplanes)
        A[i,2:end] = -h.a
        A[i,1] = h.β
    end
    for (i, h) in enumerate(halfspaces)
        A[nhyperplane+i,2:end] = -h.a
        A[nhyperplane+i,1] = h.β
    end
    LiftedHRepresentation{T}(A, linset)
end

Base.copy(ine::LiftedHRepresentation{T}) where {T} = LiftedHRepresentation{T}(copy(ine.A), copy(ine.linset))

Base.isvalid(hrep::LiftedHRepresentation{T}, idx::HIndex{T}) where {T} = 0 < idx.value <= size(hrep.A, 1) && (idx.value in hrep.linset) == islin(idx)
done(idxs::HIndices{T, <:LiftedHRepresentation{T}}, idx::HIndex{T}) where {T} = idx.value > size(idxs.rep.A, 1)
Base.get(hrep::LiftedHRepresentation{T}, idx::HIndex{T}) where {T} = valuetype(idx)(-hrep.A[idx.value,2:end], hrep.A[idx.value,1])

# V-Representation

mutable struct LiftedVRepresentation{T, MT<:AbstractMatrix{T}} <: MixedVRep{T}
    R::MT # each row is a vertex if the first element is 1 and a ray otherwise
    linset::BitSet

    function LiftedVRepresentation{T, MT}(R::MT, linset::BitSet=BitSet([])) where {T, MT}
        if !isempty(linset) && last(linset) > size(R, 1)
            error("The elements of linset should be between 1 and the number of rows of R")
        end
        new{T, MT}(R, linset)
    end
end
FullDim(rep::LiftedVRepresentation) = size(rep.R, 2) - 1

similar_type(::Type{LiftedVRepresentation{S, MT}}, ::FullDim, ::Type{T}) where {S, T, MT} = LiftedVRepresentation{T, similar_type(MT, T)}
vvectortype(::Type{LiftedVRepresentation{T, MT}}) where {T, MT} = vectortype(MT)

LiftedVRepresentation{T}(R::AbstractMatrix{T}, linset::BitSet=BitSet()) where {T} = LiftedVRepresentation{T, typeof(R)}(R, linset)
LiftedVRepresentation{T}(R::AbstractMatrix, linset::BitSet=BitSet()) where {T} = LiftedVRepresentation{T}(AbstractMatrix{T}(R), linset)
LiftedVRepresentation(R::AbstractMatrix{T}, linset::BitSet=BitSet()) where T = LiftedVRepresentation{T}(R, linset)

LiftedVRepresentation(v::VRepresentation{T}) where {T} = LiftedVRepresentation{T}(v)
LiftedVRepresentation{T}(v::VRepresentation) where {T} = convert(LiftedVRepresentation{T, vmatrixtype(typeof(v), T)}, v)

function LiftedVRepresentation{T, MT}(d::FullDim,
                                      vits::VIt{T}...) where {T, MT}
    points, lines, rays = fillvits(d, vits...)
    npoint = length(points)
    nline = length(lines)
    nray = length(rays)
    nvrep = npoint + nline + nray
    R = emptymatrix(MT, nvrep, fulldim(d) + 1)
    linset = BitSet()
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
    LiftedVRepresentation{T}(R, linset)
end

Base.copy(ext::LiftedVRepresentation{T}) where {T} = LiftedVRepresentation{T}(copy(ext.R), copy(ext.linset))

nvreps(ext::LiftedVRepresentation) = size(ext.R, 1)

function isrowpoint(ext::LiftedVRepresentation{T}, i) where {T}
    ispoint = ext.R[i,1]
    @assert ispoint == zero(T) || ispoint == one(T)
    ispoint == one(T)
end

function Base.isvalid(vrep::LiftedVRepresentation{T}, idx::VIndex{T}) where {T}
    isp = isrowpoint(vrep, idx.value)
    isl = (idx.value in vrep.linset)
    @assert !isp || !isl # if isp && isl, it is a symmetric point but it is not allowed to mix symmetric points and points
    0 < idx.value <= nvreps(vrep) && isp == ispoint(idx) && isl == islin(idx)
end
done(idxs::VIndices{T, <:LiftedVRepresentation{T}}, idx::VIndex{T}) where {T} = idx.value > size(idxs.rep.R, 1)
Base.get(vrep::LiftedVRepresentation{T}, idx::VIndex{T}) where {T} = valuetype(idx)(vrep.R[idx.value,2:end])

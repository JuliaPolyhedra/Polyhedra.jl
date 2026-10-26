export LiftedHRepresentation, LiftedVRepresentation

# H-Represenation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LiftedHRepresentation{T,MT<:AbstractMatrix{T}} <: MixedHRep{T}
    # Ax >= 0, it is [b -A] * [z; x] where z = 1
    A::MT
    linset::BitSet

    function LiftedHRepresentation{T,MT}(
        A::MT,
        linset::BitSet = BitSet(),
    ) where {T,MT}
        if !isempty(linset) && last(linset) > size(A, 1)
            error(
                "The elements of linset should be between 1 and the number of rows of A",
            )
        end
        return new{T,MT}(A, linset)
    end
end
FullDim(rep::LiftedHRepresentation) = size(rep.A, 2) - 1

function similar_type(
    ::Type{LiftedHRepresentation{S,MT}},
    ::FullDim,
    ::Type{T},
) where {S,T,MT}
    return LiftedHRepresentation{T,similar_type(MT, T)}
end
hvectortype(p::Type{LiftedHRepresentation{T,MT}}) where {T,MT} = vectortype(MT)

function LiftedHRepresentation{T}(
    A::AbstractMatrix{T},
    linset::BitSet = BitSet(),
) where {T}
    return LiftedHRepresentation{T,typeof(A)}(A, linset)
end
function LiftedHRepresentation{T}(
    A::AbstractMatrix,
    linset::BitSet = BitSet(),
) where {T}
    return LiftedHRepresentation{T}(AbstractMatrix{T}(A), linset)
end
function LiftedHRepresentation(
    A::AbstractMatrix{T},
    linset::BitSet = BitSet(),
) where {T}
    return LiftedHRepresentation{T}(A, linset)
end

function LiftedHRepresentation(h::HRepresentation{T}) where {T}
    return LiftedHRepresentation{T}(h)
end
function LiftedHRepresentation{T}(h::HRepresentation) where {T}
    return convert(LiftedHRepresentation{T,hmatrixtype(typeof(h), T)}, h)
end

function LiftedHRepresentation{T,MT}(
    d::FullDim,
    hyperplanes::HyperPlaneIt{T},
    halfspaces::HalfSpaceIt{T},
) where {T,MT}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = emptymatrix(MT, nhrep, fulldim(d)+1)
    linset = BitSet(1:nhyperplane)
    for (i, h) in enumerate(hyperplanes)
        A[i, 2:end] = -h.a
        A[i, 1] = h.β
    end
    for (i, h) in enumerate(halfspaces)
        A[nhyperplane+i, 2:end] = -h.a
        A[nhyperplane+i, 1] = h.β
    end
    return LiftedHRepresentation{T}(A, linset)
end

function Base.copy(ine::LiftedHRepresentation{T}) where {T}
    return LiftedHRepresentation{T}(copy(ine.A), copy(ine.linset))
end

function Base.isvalid(hrep::LiftedHRepresentation{T}, idx::HIndex{T}) where {T}
    return 0 < idx.value <= size(hrep.A, 1) &&
           (idx.value in hrep.linset) == islin(idx)
end
function done(
    idxs::HIndices{T,<:LiftedHRepresentation{T}},
    idx::HIndex{T},
) where {T}
    return idx.value > size(idxs.rep.A, 1)
end
function Base.get(hrep::LiftedHRepresentation{T}, idx::HIndex{T}) where {T}
    return valuetype(idx)(-hrep.A[idx.value, 2:end], hrep.A[idx.value, 1])
end

# V-Representation

mutable struct LiftedVRepresentation{T,MT<:AbstractMatrix{T}} <: MixedVRep{T}
    R::MT # each row is a vertex if the first element is 1 and a ray otherwise
    linset::BitSet

    function LiftedVRepresentation{T,MT}(
        R::MT,
        linset::BitSet = BitSet([]),
    ) where {T,MT}
        if !isempty(linset) && last(linset) > size(R, 1)
            error(
                "The elements of linset should be between 1 and the number of rows of R",
            )
        end
        return new{T,MT}(R, linset)
    end
end
FullDim(rep::LiftedVRepresentation) = size(rep.R, 2) - 1

function similar_type(
    ::Type{LiftedVRepresentation{S,MT}},
    ::FullDim,
    ::Type{T},
) where {S,T,MT}
    return LiftedVRepresentation{T,similar_type(MT, T)}
end
vvectortype(::Type{LiftedVRepresentation{T,MT}}) where {T,MT} = vectortype(MT)

function LiftedVRepresentation{T}(
    R::AbstractMatrix{T},
    linset::BitSet = BitSet(),
) where {T}
    return LiftedVRepresentation{T,typeof(R)}(R, linset)
end
function LiftedVRepresentation{T}(
    R::AbstractMatrix,
    linset::BitSet = BitSet(),
) where {T}
    return LiftedVRepresentation{T}(AbstractMatrix{T}(R), linset)
end
function LiftedVRepresentation(
    R::AbstractMatrix{T},
    linset::BitSet = BitSet(),
) where {T}
    return LiftedVRepresentation{T}(R, linset)
end

function LiftedVRepresentation(v::VRepresentation{T}) where {T}
    return LiftedVRepresentation{T}(v)
end
function LiftedVRepresentation{T}(v::VRepresentation) where {T}
    return convert(LiftedVRepresentation{T,vmatrixtype(typeof(v), T)}, v)
end

function LiftedVRepresentation{T,MT}(d::FullDim, vits::VIt{T}...) where {T,MT}
    points, lines, rays = fillvits(d, vits...)
    npoint = length(points)
    nline = length(lines)
    nray = length(rays)
    nvrep = npoint + nline + nray
    R = emptymatrix(MT, nvrep, fulldim(d) + 1)
    linset = BitSet()
    function _fill(offset, z, ps)
        for (i, p) in enumerate(ps)
            R[offset+i, 2:end] = coord(p)
            R[offset+i, 1] = z
            if islin(p)
                push!(linset, offset + i)
            end
        end
    end
    _fill(0, one(T), points)
    _fill(npoint, zero(T), lines)
    _fill(npoint+nline, zero(T), rays)
    return LiftedVRepresentation{T}(R, linset)
end

function Base.copy(ext::LiftedVRepresentation{T}) where {T}
    return LiftedVRepresentation{T}(copy(ext.R), copy(ext.linset))
end

nvreps(ext::LiftedVRepresentation) = size(ext.R, 1)

function isrowpoint(ext::LiftedVRepresentation{T}, i) where {T}
    ispoint = ext.R[i, 1]
    @assert ispoint == zero(T) || ispoint == one(T)
    return ispoint == one(T)
end

function Base.isvalid(vrep::LiftedVRepresentation{T}, idx::VIndex{T}) where {T}
    isp = isrowpoint(vrep, idx.value)
    isl = (idx.value in vrep.linset)
    @assert !isp || !isl # if isp && isl, it is a symmetric point but it is not allowed to mix symmetric points and points
    return 0 < idx.value <= nvreps(vrep) &&
           isp == ispoint(idx) &&
           isl == islin(idx)
end
function done(
    idxs::VIndices{T,<:LiftedVRepresentation{T}},
    idx::VIndex{T},
) where {T}
    return idx.value > size(idxs.rep.R, 1)
end
function Base.get(vrep::LiftedVRepresentation{T}, idx::VIndex{T}) where {T}
    return valuetype(idx)(vrep.R[idx.value, 2:end])
end

export MixedMatHRep, MixedMatVRep

# H-Representation
"""
    hrep(A::AbstractMatrix, b::AbstractVector, linset::BitSet=BitSet())

Creates an H-representation for the polyhedron defined by the inequalities ``\\langle A_i, x \\rangle = b_i`` if `i in linset`
and ``\\langle A_i, x \\rangle \\le b_i`` otherwise where ``A_i`` is the ``i``th row of `A`, i.e. `A[i,:]` and ``b_i`` is `b[i]`.
"""
hrep(A::AbstractMatrix, b::AbstractVector, linset::BitSet=BitSet()) = MixedMatHRep(A, b, linset)

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct MixedMatHRep{T, MT<:AbstractMatrix{T}} <: MixedHRep{T}
    # Ax <= b
    A::MT
    b::Vector{T}
    linset::BitSet

    function MixedMatHRep{T, MT}(A::MT, b::AbstractVector, linset::BitSet) where {T, MT}
        if size(A, 1) != length(b)
            error("The length of b must be equal to the number of rows of A")
        end
        if !isempty(linset) && last(linset) > length(b)
            error("The elements of linset should be between 1 and the number of rows of A/length of b")
        end
        if size(A, 2) != N
            error("dimension does not match")
        end
        new{T, typeof(A)}(A, b, linset)
    end
end

similar_type(::Type{<:MixedMatHRep{M, S, MT}}, ::FullDim{N}, ::Type{T}) where {M, S, T, MT} = MixedMatHRep{T, similar_type(MT, T)}
hvectortype(p::Type{MixedMatHRep{T, MT}}) where {T, MT} = vectortype(MT)
fulltype(::Type{MixedMatHRep{T, MT}}) where {T, MT} = MixedMatHRep{T, MT}

MixedMatHRep{T}(A::AbstractMatrix{T}, b::AbstractVector{T}, linset::BitSet) where {T} = MixedMatHRep{T, typeof(A)}(A, b, linset)
MixedMatHRep(A::AbstractMatrix{T}, b::AbstractVector{T}, linset::BitSet) where T = MixedMatHRep{size(A, 2), T, typeof(A)}(A, b, linset)
MixedMatHRep{U}(A::AbstractMatrix{S}, b::AbstractVector{T}, linset::BitSet) where {S, T, U} = MixedMatHRep{U}(AbstractMatrix{U}(A), AbstractVector{U}(b), linset)
function MixedMatHRep(A::AbstractMatrix{S}, b::AbstractVector{T}, linset::BitSet) where {S, T}
    U = promote_type(S, T)
    MixedMatHRep(AbstractMatrix{U}(A), AbstractVector{U}(b), linset)
end
MixedMatHRep(A::AbstractMatrix{T}, b::AbstractVector{T}, linset::BitSet) where T <: Real = MixedMatHRep{size(A,2),T}(A, b, linset)

MixedMatHRep(h::HRep{T}) where {T} = MixedMatHRep{T}(h)
MixedMatHRep{T}(h::HRep{N}) where {T} = MixedMatHRep{T, hmatrixtype(typeof(h), T)}(h)

MixedMatHRep{T}(hits::HIt{T}...) where {T} = MixedMatHRep{T, Matrix{T}}(hits...) # FIXME required by CDD and LRS tests
function MixedMatHRep{T, MT}(hyperplanes::HyperPlaneIt{T}, halfspaces::HalfSpaceIt{T}) where {T, MT}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = emptymatrix(MT, nhrep, N)
    b = Vector{T}(undef, nhrep)
    linset = BitSet(1:nhyperplane)
    for (i, h) in enumerate(hyperplanes)
        A[i,:] = h.a
        b[i] = h.β
    end
    for (i, h) in enumerate(halfspaces)
        A[nhyperplane+i,:] = h.a
        b[nhyperplane+i] = h.β
    end
    MixedMatHRep{T, MT}(A, b, linset)
end

Base.copy(ine::MixedMatHRep{T, MT}) where {T, MT} = MixedMatHRep{T, MT}(copy(ine.A), copy(ine.b), copy(ine.linset))

Base.isvalid(hrep::MixedMatHRep{T}, idx::HIndex{T}) where {T} = 0 < idx.value <= size(hrep.A, 1) && (idx.value in hrep.linset) == islin(idx)
Base.done(idxs::HIndices{T, <:MixedMatHRep{T}}, idx::HIndex{T}) where {T} = idx.value > size(idxs.rep.A, 1)
Base.get(hrep::MixedMatHRep{T}, idx::HIndex{T}) where {T} = valuetype(idx)(hrep.A[idx.value,:], hrep.b[idx.value])

# V-Representation
"""
    vrep(V::AbstractMatrix, R::AbstractMatrix, Rlinset::BitSet=BitSet())

Creates a V-representation for the polyhedron defined by the points ``V_i``,
lines ``R_i`` if `i in Rlinset` and rays ``R_i`` otherwise
where ``V_i`` (resp. ``R_i``) is the ``i``th row of `V` (resp. `R`), i.e. `V[i,:]` (resp. `R[i,:]`).
"""
vrep(V::AbstractMatrix, R::AbstractMatrix, Rlinset::BitSet=BitSet()) = MixedMatVRep(V, R, Rlinset)
vrep(V::AbstractMatrix) = vrep(V, similar(V, 0, size(V, 2)))

mutable struct MixedMatVRep{T, MT<:AbstractMatrix{T}} <: MixedVRep{T}
    V::MT # each row is a vertex
    R::MT # each row is a ray
    Rlinset::BitSet

    function MixedMatVRep{T, MT}(V::MT, R::MT, Rlinset::BitSet) where {T, MT}
        if iszero(size(V, 1)) && !iszero(size(R, 1))
            vconsistencyerror()
        end
        if (length(R) > 0 && size(R, 2) != N) || (length(V) > 0 && size(V, 2) != N)
            error("dimension does not match")
        end
        if !isempty(Rlinset) && last(Rlinset) > size(R, 1)
            error("The elements of Rlinset should be between 1 and the number of rows of R")
        end
        new{T, MT}(V, R, Rlinset)
    end
end

similar_type(::Type{<:MixedMatVRep{M, S, MT}}, ::FullDim{N}, ::Type{T}) where {M, S, T, MT} = MixedMatVRep{T, similar_type(MT, T)}
vvectortype(p::Type{MixedMatVRep{T, MT}}) where {T, MT} = vectortype(MT)
fulltype(::Type{MixedMatVRep{T, MT}}) where {T, MT} = MixedMatVRep{T, MT}

MixedMatVRep{T}(V::MT, R::MT, Rlinset::BitSet) where {T, MT<:AbstractMatrix{T}} = MixedMatVRep{T, MT}(V, R, Rlinset)
MixedMatVRep(V::MT, R::MT, Rlinset::BitSet) where {T, MT<:AbstractMatrix{T}} = MixedMatVRep{size(V, 2), T, MT}(V, R, Rlinset)
MixedMatVRep{U}(V::AbstractMatrix{S}, R::AbstractMatrix{T}, Rlinset::BitSet) where {S, T, U} = MixedMatVRep{U}(AbstractMatrix{U}(V), AbstractMatrix{U}(R), Rlinset)
function MixedMatVRep(V::AbstractMatrix{S}, R::AbstractMatrix{T}, Rlinset::BitSet) where {S, T}
    U = promote_type(S, T)
    MixedMatVRep(AbstractMatrix{U}(V), AbstractMatrix{U}(R), Rlinset)
end

MixedMatVRep(v::VRep{T}) where {T} = MixedMatVRep{T}(v)
MixedMatVRep{T}(v::VRep{N}) where {T} = MixedMatVRep{T, vmatrixtype(typeof(v), T)}(v)

MixedMatVRep{T}(vits::VIt{T}...) where {T} = MixedMatVRep{T, Matrix{T}}(vits...) # FIXME required by CDD and LRS tests
function MixedMatVRep{T, MT}(vits::VIt{T}...) where {T, MT}
    points, lines, rays = fillvits(FullDim{N}(), vits...)
    npoint = length(points)
    nline = length(lines)
    nray = length(rays)
    V = emptymatrix(MT, npoint, N)
    R = emptymatrix(MT, nline + nray, N)
    Rlinset = BitSet()
    function _fill!(M, linset, offset, ps)
        for (i, p) in enumerate(ps)
            M[offset + i, :] = coord(p)
            if islin(p)
                push!(linset, offset + i)
            end
        end
    end
    _fill!(V, nothing, 0, points)
    _fill!(R, Rlinset, 0, lines)
    _fill!(R, Rlinset, nline, rays)
    MixedMatVRep{T, MT}(V, R, Rlinset)
end

Base.copy(ext::MixedMatVRep{T, MT}) where {T, MT} = MixedMatVRep{T, MT}(copy(ext.V), copy(ext.R), copy(ext.Rlinset))

_mat(hrep::MixedMatVRep, ::PIndex) = hrep.V
_mat(hrep::MixedMatVRep, ::RIndex) = hrep.R
_islin(hrep::MixedMatVRep, ::PIndex) = false
_islin(hrep::MixedMatVRep, idx::RIndex) = idx.value in hrep.Rlinset

Base.isvalid(vrep::MixedMatVRep{T}, idx::VIndex{T}) where {T} = 0 < idx.value <= size(_mat(vrep, idx), 1) && _islin(vrep, idx) == islin(idx)
Base.done(idxs::VIndices{T, <:MixedMatVRep{T}}, idx::VIndex{T}) where {T} = idx.value > size(_mat(idxs.rep, idx), 1)
Base.get(vrep::MixedMatVRep{T}, idx::VIndex{T}) where {T} = valuetype(idx)(_mat(vrep, idx)[idx.value, :])

# sparse halfspaces does not mean sparse points
dualtype(::Type{MixedMatHRep{T, MT}}) where {T, MT} = MixedMatVRep{T, Matrix{T}}
dualtype(::Type{MixedMatVRep{T, MT}}) where {T, MT} = MixedMatHRep{T, Matrix{T}}
function dualfullspace(h::MixedMatHRep, ::FullDim{N}, ::Type{T}) where {T}
    MixedMatVRep{T}(zeros(T, 1, N), eye(T, N), BitSet(1:N))
end

export SimpleHRepresentation, SimpleVRepresentation

function SimpleHRepresentation(args...)
    Base.depwarn("`SimpleHRepresentation` is deprecated, it has been renamed to `MixedMatHRep`. It is recommended to use `hrep` to build a `MixedMatHRep`, e.g. `hrep([1 2; 3 4], [5, 6])` for `x + 2y ≤ 5` and `3x + 4y ≤ 6`.", :SimpleHRepresentation)
    hrep(args...)
end
function SimpleVRepresentation(args...)
    Base.depwarn("`SimpleVRepresentation` is deprecated, it has been renamed to `MixedMatVRep` It is recommended to use `vrep` to build a `MixedMatVRep`, e.g. `vrep([1 2; 3 4])` for the convex hull of `(1, 2)` and `(3, 4)`.", :SimpleVRepresentation)
    vrep(args...)
end

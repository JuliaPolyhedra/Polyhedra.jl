export MixedMatHRep, MixedMatVRep

# H-Representation
"""
    hrep(A::AbstractMatrix, b::AbstractVector, linset::IntSet=IntSet())

Creates an H-representation for the polyhedron defined by the inequalities ``\\langle A_i, x \\rangle = b_i`` if `i in linset`
and ``\\langle A_i, x \\rangle \\le b_i`` otherwise where ``A_i`` is the ``i``th row of `A`, i.e. `A[i,:]` and ``b_i`` is `b[i]`.
"""
hrep(A::AbstractMatrix, b::AbstractVector, linset::IntSet=IntSet()) = MixedMatHRep(A, b, linset)

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct MixedMatHRep{N, T} <: MixedHRep{N, T}
    # Ax <= b
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    linset::IntSet

    function MixedMatHRep{N, T}(A::AbstractMatrix, b::AbstractVector, linset::IntSet) where {N, T}
        if size(A, 1) != length(b)
            error("The length of b must be equal to the number of rows of A")
        end
        if !isempty(linset) && last(linset) > length(b)
            error("The elements of linset should be between 1 and the number of rows of A/length of b")
        end
        if size(A, 2) != N
            error("dimension does not match")
        end
        new{N, T}(A, b, linset)
    end
end

similar_type{N,T}(::Type{<:MixedMatHRep}, ::FullDim{N}, ::Type{T}) = MixedMatHRep{N,T}
arraytype(p::Union{MixedMatHRep{N, T}, Type{MixedMatHRep{N, T}}}) where {N, T} = Vector{T}
fulltype(::Type{MixedMatHRep{N, T}}) where {N, T} = MixedMatHRep{N, T}

function MixedMatHRep(A::AbstractMatrix{S}, b::AbstractVector{T}, linset::IntSet) where {S <: Real, T <: Real}
    U = promote_type(S, T)
    MixedMatHRep{size(A,2),U}(AbstractMatrix{U}(A), AbstractVector{U}(b), linset)
end
MixedMatHRep(A::AbstractMatrix{T}, b::AbstractVector{T}, linset::IntSet) where T <: Real = MixedMatHRep{size(A,2),T}(A, b, linset)

MixedMatHRep(h::HRep{N,T}) where {N,T} = MixedMatHRep{N,T}(h)

function MixedMatHRep{N, T}(hyperplanes::HyperPlaneIt{N, T}, halfspaces::HalfSpaceIt{N, T}) where {N, T}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = Matrix{T}(nhrep, N)
    b = Vector{T}(nhrep)
    linset = IntSet(1:nhyperplane)
    for (i, h) in enumerate(hyperplanes)
        A[i,:] = h.a
        b[i] = h.β
    end
    for (i, h) in enumerate(halfspaces)
        A[nhyperplane+i,:] = h.a
        b[nhyperplane+i] = h.β
    end
    MixedMatHRep{N, T}(A, b, linset)
end

Base.copy(ine::MixedMatHRep{N,T}) where {N,T} = MixedMatHRep{N,T}(copy(ine.A), copy(ine.b), copy(ine.linset))

Base.isvalid(hrep::MixedMatHRep{N, T}, idx::HIndex{N, T}) where {N, T} = 0 < idx.value <= size(hrep.A, 1) && (idx.value in hrep.linset) == islin(idx)
Base.done(idxs::HIndices{N, T, <:MixedMatHRep{N, T}}, idx::HIndex{N, T}) where {N, T} = idx.value > size(idxs.rep.A, 1)
Base.get(hrep::MixedMatHRep{N, T}, idx::HIndex{N, T}) where {N, T} = valuetype(idx)(hrep.A[idx.value,:], hrep.b[idx.value])

function filterintset(J::IntSet, I)
    K = IntSet()
    for (idx, i) in enumerate(I)
        if i in J
            push!(K, idx)
        end
    end
    K
end
Base.getindex(h::MixedMatHRep, I::AbstractArray) = MixedMatHRep(h.A[I, :], h.b[I], filterintset(h.linset, I))

# V-Representation
"""
    vrep(V::AbstractMatrix, R::AbstractMatrix, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet())

Creates a V-representation for the polyhedron defined by the symmetric points ``V_i`` if `i in Vlinset`,
points ``V_i`` otherwise, lines ``R_i`` if `i in Rlinset` and rays ``R_i`` otherwise
where ``V_i`` (resp. ``R_i``) is the ``i``th row of `V` (resp. `R`), i.e. `V[i,:]` (resp. `R[i,:]`).
"""
vrep(V::AbstractMatrix, R::AbstractMatrix, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet()) = MixedMatVRep(V, R, Vlinset, Rlinset)
vrep(V::AbstractMatrix, linset::IntSet=IntSet()) = vrep(V, similar(V, 0, size(V, 2)), linset)

mutable struct MixedMatVRep{N,T} <: MixedVRep{N,T}
    V::AbstractMatrix{T} # each row is a vertex
    R::AbstractMatrix{T} # each row is a ray
    Vlinset::IntSet
    Rlinset::IntSet

    function MixedMatVRep{N, T}(V::AbstractMatrix, R::AbstractMatrix, Vlinset::IntSet, Rlinset::IntSet) where {N, T}
        if iszero(size(V, 1)) && !iszero(size(R, 1))
            vconsistencyerror()
        end
        if (length(R) > 0 && size(R, 2) != N) || (length(V) > 0 && size(V, 2) != N)
            error("dimension does not match")
        end
        if !isempty(Vlinset) && last(Vlinset) > size(V, 1)
            error("The elements of Vlinset should be between 1 and the number of rows of V")
        end
        if !isempty(Rlinset) && last(Rlinset) > size(R, 1)
            error("The elements of Rlinset should be between 1 and the number of rows of R")
        end
        new{N, T}(V, R, Vlinset, Rlinset)
    end
end

similar_type{N,T}(::Type{<:MixedMatVRep}, ::FullDim{N}, ::Type{T}) = MixedMatVRep{N,T}
arraytype(p::Union{MixedMatVRep{N, T}, Type{MixedMatVRep{N, T}}}) where {N, T} = Vector{T}
fulltype(::Type{MixedMatVRep{N, T}}) where {N, T} = MixedMatVRep{N, T}

function MixedMatVRep(V::AbstractMatrix{S}, R::AbstractMatrix{T}, Vlinset::IntSet, Rlinset::IntSet) where {S <: Real, T <: Real}
    U = promote_type(S, T)
    MixedMatVRep{size(V,2),U}(AbstractMatrix{U}(V), AbstractMatrix{U}(R), Vlinset, Rlinset)
end

MixedMatVRep(v::VRep{N,T}) where {N,T} = MixedMatVRep{N,T}(v)

function MixedMatVRep{N, T}(sympoints::SymPointIt{N, T}, points::PointIt{N, T}, lines::LineIt{N, T}, rays::RayIt{N, T}) where {N, T}
    nsympoint = length(sympoints)
    npoint = length(points)
    nline = length(lines)
    nray = length(rays)
    V = Matrix{T}(nsympoint + npoint, N)
    R = Matrix{T}(nline + nray, N)
    Vlinset = IntSet()
    Rlinset = IntSet()
    function _fill!(M, linset, offset, ps)
        for (i, p) in enumerate(ps)
            M[offset + i, :] = coord(p)
            if islin(p)
                push!(linset, offset + i)
            end
        end
    end
    _fill!(V, Vlinset, 0, sympoints)
    _fill!(V, Vlinset, nsympoint, points)
    _fill!(R, Rlinset, 0, lines)
    _fill!(R, Rlinset, nline, rays)
    MixedMatVRep{N, T}(V, R, Vlinset, Rlinset)
end

Base.copy(ext::MixedMatVRep{N,T}) where {N,T} = MixedMatVRep{N,T}(copy(ext.V), copy(ext.R), copy(ext.Vlinset), copy(ext.Rlinset))

_mat(hrep::MixedMatVRep, ::PIndex) = hrep.V
_mat(hrep::MixedMatVRep, ::RIndex) = hrep.R
_linset(hrep::MixedMatVRep, ::PIndex) = hrep.Vlinset
_linset(hrep::MixedMatVRep, ::RIndex) = hrep.Rlinset

Base.isvalid(vrep::MixedMatVRep{N, T}, idx::VIndex{N, T}) where {N, T} = 0 < idx.value <= size(_mat(vrep, idx), 1) && (idx.value in _linset(vrep, idx)) == islin(idx)
Base.done(idxs::VIndices{N, T, <:MixedMatVRep{N, T}}, idx::VIndex{N, T}) where {N, T} = idx.value > size(_mat(idxs.rep, idx), 1)
Base.get(vrep::MixedMatVRep{N, T}, idx::VIndex{N, T}) where {N, T} = valuetype(idx)(_mat(vrep, idx)[idx.value, :])

function Base.getindex(h::MixedMatVRep, I::AbstractArray)
    Ir = filter(i -> i <= nrays(h), I)
    Ip = filter(i -> i > nrays(h), I) - nrays(h)
    MixedMatVRep(h.V[Ip,:], h.R[Ir,:], filterintset(h.Vlinset, I), filterintset(h.Rlinset, I))
end

function dualfullspace(h::MixedMatHRep, ::FullDim{N}, ::Type{T}) where {N, T}
    MixedMatVRep{N, T}(zeros(T, 1, N), eye(T, N), IntSet(), IntSet(1:N))
end

export SimpleHRepresentation, SimpleVRepresentation

function SimpleHRepresentation(args...)
    warn("`SimpleHRepresentation` is deprecated, it has been renamed to `MixedMatHRep`")
    MixedMatHRep(args...)
end
function SimpleVRepresentation(args...)
    warn("`SimpleVRepresentation` is deprecated, it has been renamed to `MixedMatVRep`")
    MixedMatVRep(args...)
end

export SimpleHRepresentation, SimpleVRepresentation

# H-Representation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct SimpleHRepresentation{N, T} <: MixedHRep{N, T}
    # Ax <= b
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    linset::IntSet

    function SimpleHRepresentation{N, T}(A::AbstractMatrix, b::AbstractVector, linset::IntSet=IntSet()) where {N, T}
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

similar_type{N,T}(::Type{<:SimpleHRepresentation}, ::FullDim{N}, ::Type{T}) = SimpleHRepresentation{N,T}
arraytype(p::Union{SimpleHRepresentation{N, T}, Type{SimpleHRepresentation{N, T}}}) where {N, T} = Vector{T}

function SimpleHRepresentation(A::AbstractMatrix{S}, b::AbstractVector{T}, linset::IntSet=IntSet()) where {S <: Real, T <: Real}
    U = promote_type(S, T)
    SimpleHRepresentation{size(A,2),U}(AbstractMatrix{U}(A), AbstractVector{U}(b), linset)
end
SimpleHRepresentation(A::AbstractMatrix{T}, b::AbstractVector{T}, linset::IntSet=IntSet()) where {T <: Real} = SimpleHRepresentation{size(A,2),T}(A, b, linset)

SimpleHRepresentation(h::HRep{N,T}) where {N,T} = SimpleHRepresentation{N,T}(h)

function SimpleHRepresentation{N, T}(hyperplanes::ElemIt{<:HyperPlane{N, T}}, halfspaces::ElemIt{<:HalfSpace{N, T}}) where {N, T}
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
    SimpleHRepresentation{N, T}(A, b, linset)
end

Base.copy(ine::SimpleHRepresentation{N,T}) where {N,T} = SimpleHRepresentation{N,T}(copy(ine.A), copy(ine.b), copy(ine.linset))

Base.isvalid(hrep::SimpleHRepresentation{N, T}, idx::HIndex{N, T}) where {N, T} = 0 < idx.value <= size(hrep.A, 1) && (idx.value in hrep.linset) == islin(idx)
Base.done(idxs::HIndices{N, T, <:SimpleHRepresentation{N, T}}, idx::HIndex{N, T}) where {N, T} = idx.value > size(idxs.rep.A, 1)
Base.get(hrep::SimpleHRepresentation{N, T}, idx::HIndex{N, T}) where {N, T} = valuetype(idx)(hrep.A[idx.value,:], hrep.b[idx.value])

function filterintset(J::IntSet, I)
    K = IntSet()
    for (idx, i) in enumerate(I)
        if i in J
            push!(K, idx)
        end
    end
    K
end
Base.getindex(h::SimpleHRepresentation, I::AbstractArray) = SimpleHRepresentation(h.A[I, :], h.b[I], filterintset(h.linset, I))

# V-Representation

mutable struct SimpleVRepresentation{N,T} <: MixedVRep{N,T}
    V::AbstractMatrix{T} # each row is a vertex
    R::AbstractMatrix{T} # each row is a ray
    Vlinset::IntSet
    Rlinset::IntSet

    function SimpleVRepresentation{N, T}(V::AbstractMatrix, R::AbstractMatrix, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet()) where {N, T}
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

similar_type{N,T}(::Type{<:SimpleVRepresentation}, ::FullDim{N}, ::Type{T}) = SimpleVRepresentation{N,T}
arraytype(p::Union{SimpleVRepresentation{N, T}, Type{SimpleVRepresentation{N, T}}}) where {N, T} = Vector{T}

function SimpleVRepresentation(V::AbstractMatrix{S}, R::AbstractMatrix{T}, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet()) where {S <: Real, T <: Real}
    U = promote_type(S, T)
    SimpleVRepresentation{size(V,2),U}(AbstractMatrix{U}(V), AbstractMatrix{U}(R), Vlinset, Rlinset)
end
SimpleVRepresentation(V::AbstractMatrix{T}, linset::IntSet=IntSet()) where {T <: Real} = SimpleVRepresentation{size(V, 2),T}(V, similar(V, 0, size(V, 2)), linset, IntSet())

SimpleVRepresentation(v::VRep{N,T}) where {N,T} = SimpleVRepresentation{N,T}(v)

function SimpleVRepresentation{N, T}(sympoints::ElemIt{<:SymPoint{N, T}}, points::ElemIt{<:MyPoint{N, T}}, lines::ElemIt{<:Line{N, T}}, rays::ElemIt{<:Ray{N, T}}) where {N, T}
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
    SimpleVRepresentation{N, T}(V, R, Vlinset, Rlinset)
end

Base.copy(ext::SimpleVRepresentation{N,T}) where {N,T} = SimpleVRepresentation{N,T}(copy(ext.V), copy(ext.R), copy(ext.Vlinset), copy(ext.Rlinset))

_mat(hrep::SimpleVRepresentation, ::PIndex) = hrep.V
_mat(hrep::SimpleVRepresentation, ::RIndex) = hrep.R
_linset(hrep::SimpleVRepresentation, ::PIndex) = hrep.Vlinset
_linset(hrep::SimpleVRepresentation, ::RIndex) = hrep.Rlinset

Base.isvalid(vrep::SimpleVRepresentation{N, T}, idx::VIndex{N, T}) where {N, T} = 0 < idx.value <= size(_mat(vrep, idx), 1) && (idx.value in _linset(vrep, idx)) == islin(idx)
Base.done(idxs::VIndices{N, T, <:SimpleVRepresentation{N, T}}, idx::VIndex{N, T}) where {N, T} = idx.value > size(_mat(idxs.rep, idx), 1)
Base.get(vrep::SimpleVRepresentation{N, T}, idx::VIndex{N, T}) where {N, T} = valuetype(idx)(_mat(vrep, idx)[idx.value, :])

function Base.getindex(h::SimpleVRepresentation, I::AbstractArray)
    Ir = filter(i -> i <= nrays(h), I)
    Ip = filter(i -> i > nrays(h), I) - nrays(h)
    SimpleVRepresentation(h.V[Ip,:], h.R[Ir,:], filterintset(h.Vlinset, I), filterintset(h.Rlinset, I))
end

export LPHRepresentation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LPHRepresentation{N, T, MT<:AbstractMatrix{T}} <: MixedHRep{N, T}
    # lb <= Ax <= ub
    # l <= x <= u
    A::MT
    l::AbstractVector{T}
    u::AbstractVector{T}
    colleqs::IntSet
    colgeqs::IntSet
    coleqs::IntSet
    lb::AbstractVector{T}
    ub::AbstractVector{T}
    rowleqs::IntSet
    rowgeqs::IntSet
    roweqs::IntSet

    function LPHRepresentation{N, T, MT}(A::MT, l::AbstractVector{T}, u::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T}) where {N, T, MT<:AbstractMatrix{T}}
        if length(l) != length(u) || size(A, 2) != length(l)
            error("The length of l and u must be equal to the number of rows of A")
        end
        if length(lb) != length(ub) || size(A, 1) != length(lb)
            error("The length of lb and ub must be equal to the number of columns of A")
        end
        if size(A, 2) != N
            error("Type dimension does not match the number of rows of A")
        end
        colleqs = IntSet()
        colgeqs = IntSet()
        coleqs = IntSet()
        for i in 1:N
            leq = u[i] < typemax(T)
            geq = l[i] > typemin(T)
            if leq && geq && myeq(l[i], u[i])
                push!(coleqs, i)
            else
                if leq
                    push!(colleqs, i)
                end
                if geq
                    push!(colgeqs, i)
                end
            end
        end
        rowleqs = IntSet()
        rowgeqs = IntSet()
        roweqs = IntSet()
        for i in 1:size(A, 1)
            leq = ub[i] < typemax(T)
            geq = lb[i] > typemin(T)
            if leq && geq && myeq(lb[i], ub[i])
                push!(roweqs, i)
            else
                if leq
                    push!(rowleqs, i)
                end
                if geq
                    push!(rowgeqs, i)
                end
            end
        end
        new{N, T, typeof(A)}(A, l, u, colleqs, colgeqs, coleqs, lb, ub, rowleqs, rowgeqs, roweqs)
    end
end

LPHRepresentation(A::AbstractMatrix{T}, l::AbstractVector{T}, u::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T}) where {T <: Real} = LPHRepresentation{size(A,2), T, typeof(A)}(A, l, u, lb, ub)
function LPHRepresentation(A::AbstractMatrix, l::AbstractVector, u::AbstractVector, lb::AbstractVector, ub::AbstractVector)
    T = promote_type(eltype(A), eltype(l), eltype(u), eltype(lb), eltype(ub))
    AT = AbstractMatrix{T}(A)
    LPHRepresentation{size(A,2), T, typeof(AT)}(AT, AbstractVector{T}(l), AbstractVector{T}(u), AbstractVector{T}(lb), AbstractVector{T}(ub))
end

similar_type(::Type{<:Matrix}, ::Type{T}) where T = Matrix{T}
similar_type(::Type{SparseMatrixCSC{S, I}}, ::Type{T}) where {S, I, T} = SparseMatrixCSC{T, I}
arraytype(::Union{LPHRepresentation{N, T, MT}, Type{LPHRepresentation{N, T, MT}}}) where {N, T, MT} = MT <: AbstractSparseArray ? SparseVector{T, Int} : Vector{T}
similar_type(::Type{LPHRepresentation{M, S, MT}}, ::FullDim{N}, ::Type{T}) where {M, S, N, T, MT} = LPHRepresentation{N, T, similar_type(MT, T)}

LPHRepresentation(rep::LPHRepresentation) = rep
_mattype(::Type{<:AbstractVector{T}}) where T = Matrix{T}
_mattype(::Type{<:AbstractSparseVector{T}}) where T = SparseMatrixCSC{T, Int}
LPHRepresentation(rep::HRep{N, T}) where {N, T} = LPHRepresentation{N, T, _mattype(arraytype(rep))}(rep)

#function LPHRepresentation{N, T}(it::HRepIterator{N, T}) where {N,T}
#    A = Matrix{T}(length(it), N)
#    lb = Vector{T}(length(it))
#    ub = Vector{T}(length(it))
#    MathProgBase.HighLevelInterface.warn_no_inf(T)
#    l = fill(typemin(T), N)
#    u = fill(typemax(T), N)
#    for (i, h) in enumerate(it)
#        A[i,:] = h.a
#        ub[i] = h.β
#        if islin(h)
#            lb[i] = ub[i]
#        else
#            lb[i] = typemin(T)
#        end
#    end
#    LPHRepresentation{N, T}(A, l, u, lb, ub)
#end
function LPHRepresentation{N, T, MT}(hyperplanes::ElemIt{<:HyperPlane{N, T}}, halfspaces::ElemIt{<:HalfSpace{N, T}}) where {N, T, MT}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = MT <: AbstractSparseArray ? spzeros(eltype(MT), nhrep, N) : MT(nhrep, N)
    lb = Vector{T}(nhrep)
    ub = Vector{T}(nhrep)
    MathProgBase.HighLevelInterface.warn_no_inf(T)
    l = fill(typemin(T), N)
    u = fill(typemax(T), N)
    for (i, h) in enumerate(hyperplanes)
        A[i,:] = h.a
        lb[i] = h.β
        ub[i] = h.β
    end
    for (i, h) in enumerate(halfspaces)
        A[nhyperplane+i,:] = h.a
        lb[nhyperplane+i] = typemin(T)
        ub[nhyperplane+i] = h.β
    end
    LPHRepresentation{N, T, MT}(A, l, u, lb, ub)
end

Base.copy(lp::LPHRepresentation{N, T, MT}) where {N, T, MT} = LPHRepresentation{N, T, MT}(copy(lp.A), copy(lp.l), copy(lp.u), copy(lp.colleqs), copy(lp.colgeqs), copy(lp.coleqs), copy(lp.lb), copy(lp.ub), copy(lp.rowleqs), copy(lp.rowgeqs), copy(lp.roweqs))

Base.length(idxs::Indices{N, T, <:HyperPlane{N, T}, LPHRepresentation{N, T}}) where {N, T} = length(idxs.rep.coleqs) + length(idxs.rep.roweqs)
Base.length(idxs::Indices{N, T, <:HalfSpace{N, T}, LPHRepresentation{N, T}}) where {N, T} = length(idxs.rep.colleqs) + length(idxs.rep.colgeqs) + length(idxs.rep.rowleqs) + length(idxs.rep.rowgeqs)

# m rows (constraints) and n columns (variables)
# state : colrow, i, lgeq
# * colrow :
#   - 1 : columns (variables)
#   - 2 : rows (constraints)
#   - 3 : done
# * i : index of row or column
# * lgeq :
#   - 1 : ⟨a, x⟩ ≤ β
#   - 2 : ⟨a, x⟩ ≥ β
#   - 3 : ⟨a, x⟩ = β

# HyperPlane index j :
#   1 <= j <= n   : column j
# n+1 <= j <= n+m : row    j-n
function _index2state(lp, idx::HyperPlaneIndex)
    m = size(lp.A, 1)
    n = size(lp.A, 2)
    j = idx.value
    @assert j >= 1
    if j <= n
        1, j, 3
    elseif j <= n+m
        2, j-n, 3
    else
        3, 0, 0
    end
end

# HalfSpace index k :
# first bit :
# - 0 : lgeq = 1
# - 1 : lgeq = 2
# remaining bits : j
#   1 <= j+1 <= n   : column j
# n+1 <= j+1 <= n+m : row    j-n
function _index2state(lp, idx::HalfSpaceIndex)
    m = size(lp.A, 1)
    n = size(lp.A, 2)
    k = idx.value-1
    lgeq = (k & 1) + 1
    j = (k >> 1) + 1 # drop first bit
    @assert j >= 1
    if j <= n
        1, j, lgeq
    elseif j <= n+m
        2, j-n, lgeq
    else
        3, 0, 0
    end
end

function Base.isvalid(lp::LPHRepresentation{N, T}, idx::HIndex{N, T}) where {N, T}
    colrow, i, lgeq = _index2state(lp, idx)
    if colrow == 1
        lgeqs = (lp.colleqs, lp.colgeqs, lp.coleqs)
    else
        lgeqs = (lp.rowleqs, lp.rowgeqs, lp.roweqs)
    end
    1 <= colrow <= 2 && i in lgeqs[lgeq]
end
Base.done(idxs::HIndices{N, T, <:LPHRepresentation{N, T}}, idx::HIndex{N, T}) where {N, T} = _index2state(idxs.rep, idx)[1] == 3

function getaβ(lp::LPHRepresentation{N, T}, idx::HIndex{N, T}) where {N, T}
    colrow, i, lgeq = _index2state(lp, idx)
    if colrow == 1
        a = origin(arraytype(lp), FullDim{N}())
        a[i] = lgeq == 2 ? -one(T) : one(T)
        β = lgeq == 2 ? -lp.l[i] : lp.u[i]
    else
        @assert colrow == 2
        a = lgeq == 2 ? -lp.A[i,:] : lp.A[i,:]
        β = lgeq == 2 ? -lp.lb[i] : lp.ub[i]
    end
    a, β
end
Base.get(lp::LPHRepresentation{N, T}, idx::HIndex{N, T}) where {N, T} = valuetype(idx)(getaβ(lp, idx)...)

dualfullspace(h::LPHRepresentation, d::FullDim{N}, ::Type{T}, ::Type{AT}) where {N, T, AT} = dualfullspace(Intersection{N, T, AT}, d, T, AT)
dualfullspace(h::LPHRepresentation, d::FullDim, ::Type{T}) where T = dualfullspace(h, d, T, Vector{T})

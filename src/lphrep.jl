export LPHRepresentation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LPHRepresentation{T, MT<:AbstractMatrix{T}} <: MixedHRep{T}
    # lb <= Ax <= ub
    # l <= x <= u
    A::MT
    l::AbstractVector{T}
    u::AbstractVector{T}
    colleqs::BitSet
    colgeqs::BitSet
    coleqs::BitSet
    lb::AbstractVector{T}
    ub::AbstractVector{T}
    rowleqs::BitSet
    rowgeqs::BitSet
    roweqs::BitSet

    function LPHRepresentation{T, MT}(A::MT, l::AbstractVector{T}, u::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T}) where {T, MT<:AbstractMatrix{T}}
        if length(l) != length(u) || size(A, 2) != length(l)
            throw(DimensionMismatch("The length of l and u must be equal to the number of rows of A"))
        end
        if length(lb) != length(ub) || size(A, 1) != length(lb)
            throw(DimensionMismatch("The length of lb and ub must be equal to the number of columns of A"))
        end
        colleqs = BitSet()
        colgeqs = BitSet()
        coleqs = BitSet()
        for i in 1:size(A, 2)
            leq = u[i] < typemax(T)
            geq = l[i] > typemin(T)
            if leq && geq && _isapprox(l[i], u[i])
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
        rowleqs = BitSet()
        rowgeqs = BitSet()
        roweqs = BitSet()
        for i in 1:size(A, 1)
            leq = ub[i] < typemax(T)
            geq = lb[i] > typemin(T)
            if leq && geq && _isapprox(lb[i], ub[i])
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
        new{T, typeof(A)}(A, l, u, colleqs, colgeqs, coleqs, lb, ub, rowleqs, rowgeqs, roweqs)
    end
end
FullDim(rep::LPHRepresentation) = size(rep.A, 2)

LPHRepresentation(A::AbstractMatrix{T}, l::AbstractVector{T}, u::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T}) where {T <: Real} = LPHRepresentation{T, typeof(A)}(A, l, u, lb, ub)
function LPHRepresentation(A::AbstractMatrix, l::AbstractVector, u::AbstractVector, lb::AbstractVector, ub::AbstractVector)
    T = promote_type(eltype(A), eltype(l), eltype(u), eltype(lb), eltype(ub))
    AT = AbstractMatrix{T}(A)
    LPHRepresentation{T, typeof(AT)}(AT, AbstractVector{T}(l), AbstractVector{T}(u), AbstractVector{T}(lb), AbstractVector{T}(ub))
end

hvectortype(::Type{LPHRepresentation{T, MT}}) where {T, MT} = vectortype(MT)
similar_type(::Type{LPHRepresentation{S, MT}}, ::FullDim, ::Type{T}) where {S, T, MT} = LPHRepresentation{T, similar_type(MT, T)}
fulltype(::Type{LPHRepresentation{T, MT}}) where {T, MT} = LPHRepresentation{T, MT}

LPHRepresentation(h::HRep{T}) where {T} = LPHRepresentation{T}(h)
LPHRepresentation{T}(h::HRep) where {T} = convert(LPHRepresentation{T, hmatrixtype(typeof(h), T)}, h)

#function LPHRepresentation{T}(it::HRepIterator{T}) where {T}
#    A = Matrix{T}(length(it), N)
#    lb = Vector{T}(length(it))
#    ub = Vector{T}(length(it))
#    MPB.HighLevelInterface.warn_no_inf(T)
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
#    LPHRepresentation{T}(A, l, u, lb, ub)
#end
function LPHRepresentation{T, MT}(d::FullDim,
                                  hyperplanes::ElemIt{<:HyperPlane{T}},
                                  halfspaces::ElemIt{<:HalfSpace{T}}) where {T, MT}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    N = fulldim(d)
    A = emptymatrix(MT, nhrep, N)
    lb = Vector{T}(undef, nhrep)
    ub = Vector{T}(undef, nhrep)
    MPB.HighLevelInterface.warn_no_inf(T)
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
    LPHRepresentation{T, MT}(A, l, u, lb, ub)
end

Base.copy(lp::LPHRepresentation{T, MT}) where {T, MT} = LPHRepresentation{T, MT}(copy(lp.A), copy(lp.l), copy(lp.u), copy(lp.colleqs), copy(lp.colgeqs), copy(lp.coleqs), copy(lp.lb), copy(lp.ub), copy(lp.rowleqs), copy(lp.rowgeqs), copy(lp.roweqs))

Base.length(idxs::Indices{T, <:HyperPlane{T}, LPHRepresentation{T}}) where {T} = length(idxs.rep.coleqs) + length(idxs.rep.roweqs)
Base.length(idxs::Indices{T, <:HalfSpace{T}, LPHRepresentation{T}}) where {T} = length(idxs.rep.colleqs) + length(idxs.rep.colgeqs) + length(idxs.rep.rowleqs) + length(idxs.rep.rowgeqs)

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

function Base.isvalid(lp::LPHRepresentation{T}, idx::HIndex{T}) where {T}
    colrow, i, lgeq = _index2state(lp, idx)
    if colrow == 1
        lgeqs = (lp.colleqs, lp.colgeqs, lp.coleqs)
    else
        lgeqs = (lp.rowleqs, lp.rowgeqs, lp.roweqs)
    end
    1 <= colrow <= 2 && i in lgeqs[lgeq]
end

function getaβ(lp::LPHRepresentation{T}, idx::HIndex{T}) where {T}
    colrow, i, lgeq = _index2state(lp, idx)
    if colrow == 1
        a = origin(hvectortype(typeof(lp)), fulldim(lp))
        a[i] = lgeq == 2 ? -one(T) : one(T)
        β = lgeq == 2 ? -lp.l[i] : lp.u[i]
    else
        @assert colrow == 2
        a = lgeq == 2 ? -lp.A[i,:] : lp.A[i,:]
        β = lgeq == 2 ? -lp.lb[i] : lp.ub[i]
    end
    a, β
end
done(idxs::HIndices{T, <:LPHRepresentation{T}}, idx::HIndex{T}) where {T} = _index2state(idxs.rep, idx)[1] == 3
Base.get(lp::LPHRepresentation{T}, idx::HIndex{T}) where {T} = valuetype(idx)(getaβ(lp, idx)...)

dualtype(::Type{<:LPHRepresentation{T}}, ::Type{AT}) where {T, AT} = dualtype(Intersection{T, AT, Int}, AT)
dualfullspace(h::LPHRepresentation, d::FullDim, ::Type{T}, ::Type{AT}) where {T, AT} = dualfullspace(Intersection{T, AT, Int}, d, T, AT)

struct HSumSpace{T} <: HRepresentation{T}
    d::Int
end
hvectortype(::Type{HSumSpace{T}}) where T = SparseVector{T, Int}
Base.length(idxs::HyperPlaneIndices{T, HSumSpace{T}}) where T = idxs.rep.d
Base.isempty(idxs::HyperPlaneIndices{T, HSumSpace{T}}) where T = idxs.rep.d <= 0
startindex(idxs::HyperPlaneIndices{T, HSumSpace{T}}) where T = isempty(idxs) ? nothing : eltype(idxs)(1)
function Base.get(h::HSumSpace{T}, idx::HyperPlaneIndex) where T
    i = idx.value
    a = sparsevec(i:h.d:i+2h.d, StaticArrays.SVector(one(T), -one(T), -one(T)), 3h.d)
    return HyperPlane(a, zero(T))
end
function nextindex(h::HSumSpace, idx::HyperPlaneIndex)
    if idx.value < h.d
        return typeof(idx)(idx.value + 1)
    else
        return nothing
    end
end
@norepelem HSumSpace HalfSpace

struct Projection{S, I}
    set::S
    dimensions::I
end
fulldim(p::Projection) = length(p.dimensions)
dimension_names(p::Projection) = dimension_names(p.set)[p.dimensions]

function __pad_map(i, d, n, start)
    if i <= start
        return i
    elseif i <= start + n
        return nothing
    else
        return i - n
    end
end
function __pad_map(i, d, n)
    if n < 0
        if i <= n
            return nothing
        else
            return i - n
        end
    else
        if i <= d
            return i
        else
            return nothing
        end
    end
end
function _pad_map(i, d, padding...)
    j = __pad_map(i, d, padding...)
    return j === nothing ? nothing : (1, j)
end

# TODO should be cartesian product with FullSpace
function zeropad(p::HRep{T}, padding...) where T
    d = map_fulldim(N -> N + abs(padding[1]), FullDim(p))
    f = (i, el) -> zeropad(el, padding...)
    return similar(p, d, T, hmap(f, d, T, p)...;
                   dimension_map = i -> _pad_map(i, fulldim(p), padding...))
end

function Base.intersect(p1::Projection{S1, <:Base.OneTo}, p2::Projection{S2, <:Base.OneTo}) where {S1, S2}
    d = length(p1.dimensions)
    d == length(p2.dimensions) || throw(DimensionMismatch())
    p = intersect(zeropad(p1.set, fulldim(p2.set) - d),
                  zeropad(p2.set, fulldim(p1.set) - d, d))
    return Projection(p, Base.OneTo(d))
end

"""
    convexhull(p1::HRepresentation, p2::HRepresentation)

Returns the Balas [Theorem 3.3, B85] extended H-representation of the convex hull
of `p1` and `p2`.

[B85] Balas, E., 1985.
*Disjunctive programming and a hierarchy of relaxations for discrete optimization problems*.
SIAM Journal on Algebraic Discrete Methods, 6(3), pp.466-486.
"""
function convexhull(p1::HRepresentation, p2::HRepresentation)
    d = FullDim(p1)
    T = promote_coefficient_type((p1, p2))
    d_2 = map_fulldim(N -> 2N, d)
    d_3 = map_fulldim(N -> 3N, d)
    function f(i, el)
        if i == 3
            return zeropad(el, 1)
        elseif i == 4
            return zeropad(el, -d_3)
        else
            if i == 1
                _padded = zeropad(el, neg_fulldim(d))
                padded = zeropad(_padded, d)
                return lift(padded, false)
            else
                padded = zeropad(el, neg_fulldim(d_2))
                return complement_lift(padded, false)
            end
        end
    end
    ei = basis(hvectortype(p1), map_fulldim(N -> 1, d), 1)
    λ = HalfSpace(-ei, zero(T)) ∩ HalfSpace(ei, one(T))
    lifted = similar((p1, p2), map_fulldim(N -> 3N + 1, d), T,
        hmap(f, d, T, p1, p2, HSumSpace{T}(fulldim(p1)), λ)...;
        dimension_map = i -> (1 <= i <= fulldim(d)) ? (1, i) : nothing)
    return Projection(lifted, Base.OneTo(fulldim(p1)))
end

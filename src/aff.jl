# It is easy to go from H-rep of affine space to V-rep of affine space by computing the kernel of a matrix using RowEchelon
# However, it is really worth it since Base.in for an HRepElement in HyperPlanesIntersection and Base.in for an VRepElement in LinesHull are false already.
function remproj(x::HRepElement{<:Integer}, l::HRepElement{<:Integer})
    # generates large numbers but keeps the type integer
    x * dot(coord(l), coord(l)) - l * dot(coord(x), coord(l))
end
function remproj(x, l)
    ll = dot(coord(l), coord(l))
    if isapproxzero(ll)
        x
    else
        simplify(x - l * ((dot(coord(x), coord(l)) / ll)))
    end
end

# H-representation

"""
    hrep(hyperplanes::HyperPlaneIt; d::FullDim)

Creates an affine space of full dimension `d` from the list of hyperplanes
`hyperplanes`.

### Examples
```julia
hrep([HyperPlane([0, 1, 0], 1), HyperPlane([0, 0, 1], 0)])
```
creates the 1-dimensional affine subspace containing all the points ``(x_1, 0, 0)``, i.e. the ``x_1``-axis.

```julia
hrep([HyperPlane([1, 1], 1), HyperPlane([1, 0], 0)])
```
creates the 0-dimensional affine subspace only containing the point ``(0, 1)``.
"""
function hrep(hyperplanes::HyperPlaneIt; d::FullDim = FullDim_rec(hyperplanes))
    return HyperPlanesIntersection(d, hyperplanes)
end

# Representation of an affine space as the intersection of hyperplanes.
# Also called affine set, affine manifold, affine variety, linear variety or flat.
# An affine space L satisfies:
# λx + (1-λ)y ∈ L, ∀x, y ∈ L, ∀ λ ∈ R
# Note that λ is not rhyperplaneuired to be between 0 and 1 as in convex sets.
struct HyperPlanesIntersection{T, AT, D<:FullDim} <: HAffineSpace{T}
    d::D
    # HyperPlanes whose intersection is the affine space
    hyperplanes::Vector{HyperPlane{T, AT}}
    function HyperPlanesIntersection{T, AT, D}(d::FullDim,
                                               hps::HyperPlaneIt{T}) where {T, AT, D}
        new{T, AT, D}(FullDim_convert(D, d), lazy_collect(hps))
    end
end
function HyperPlanesIntersection(d::FullDim,
                                 it::ElemIt{HyperPlane{T, AT}}) where {T, AT}
    return HyperPlanesIntersection{T, AT, typeof(d)}(d, it)
end
HyperPlanesIntersection{T, AT}(d::FullDim) where {T, AT} = HyperPlanesIntersection(d, HyperPlane{T, AT}[])

FullDim(h::HyperPlanesIntersection) = h.d
hvectortype(L::Type{<:HyperPlanesIntersection{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:HyperPlanesIntersection}, d::FullDim, ::Type{T}) where T = HyperPlanesIntersection{T, similar_type(hvectortype(PT), d, T), typeof(d)}

Base.intersect!(L::HyperPlanesIntersection, h::HyperPlane) = push!(L.hyperplanes, h)

@vecrepelem HyperPlanesIntersection HyperPlane hyperplanes

hreptype(::Type{<:HyperPlanesIntersection{T, AT}}) where {T, AT} = Intersection{T, AT}

# Returns an HyperPlanesIntersection representing the affine hull of p.
# The affine hull is defined as
# {λx + (1-λ)y | x, y ∈ p, λ ∈ R}
function affinehull(h::HRep, current=false)
    if !current
        detecthlinearity!(h)
    end
    HyperPlanesIntersection(FullDim(h), hyperplanes(h))
end

function remproj(h::HRepElement, L::HyperPlanesIntersection)
    for hp in hyperplanes(L)
        h = remproj(h, hp)
    end
    h
end
function Base.in(h::HRepElement, L::HyperPlanesIntersection)
    h = remproj(h, L)
    isapproxzero(h)
end

function removeduplicates(L::HyperPlanesIntersection{T, AT}) where {T, AT}
    H = HyperPlanesIntersection{T, AT}(FullDim(L))
    for h in hyperplanes(L)
        if !(h in H)
            intersect!(H, h)
        end
    end
    H
end

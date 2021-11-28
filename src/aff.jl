export HyperPlanesIntersection, LinesHull, affinehull, linespace

# It is easy to go from H-rep of affine space to V-rep of affine space by computing the kernel of a matrix using RowEchelon
# However, it is really worth it since Base.in for an HRepElement in HyperPlanesIntersection and Base.in for an VRepElement in LinesHull are false already.
function remproj(x::RepElement{<:Integer}, l::RepElement{<:Integer})
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
mutable struct HyperPlanesIntersection{T, AT, D<:FullDim} <: HAffineSpace{T}
    d::D
    # HyperPlanes whose intersection is the affine space
    hyperplanes::Vector{HyperPlane{T, AT}}
    orthogonalized::Bool
    function HyperPlanesIntersection{T, AT, D}(d::FullDim,
                                               hps::HyperPlaneIt{T}) where {T, AT, D}
        new{T, AT, D}(FullDim_convert(D, d), lazy_collect(hps), false)
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

function Base.empty!(L::HyperPlanesIntersection{T, AT}) where {T, AT}
    # We want the dimension to be `-1` and the hyperplanes to make
    # an infeasible system so we build the system `x = 0 & 0 = 1`.
    empty!(L.hyperplanes)
    for i in 1:fulldim(L)
        push!(L.hyperplanes, HyperPlane(basis(AT, L.d, i), zero(T)))
    end
    push!(L.hyperplanes, HyperPlane(origin(AT, fulldim(L)), one(T)))
    return L
end
function Base.intersect!(L::HyperPlanesIntersection, h::HyperPlane)
    if L.orthogonalized
        h_p = remproj(h, L)
        if isapproxzero(h_p.a)
            if !isapproxzero(h_p.β)
                # `L` is empty
                empty!(L)
            end
        else
            push!(L.hyperplanes, h_p)
        end
    else
        push!(L.hyperplanes, h)
    end
    return L
end

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

function detecthlinearity!(L::HyperPlanesIntersection, solver = nothing)
    if !L.orthogonalized
        V = detecthlinearity(L)
        L.hyperplanes = V.hyperplanes
        L.orthogonalized = true
    end
    return L
end
function detecthlinearity(L::HyperPlanesIntersection{T, AT}, solver=nothing) where {T, AT}
    H = HyperPlanesIntersection{T, AT}(FullDim(L))
    for h in hyperplanes(L)
        hp = remproj(h, H)
        if !(h in H)
            intersect!(H, h)
        end
    end
    H
end

# V-representation
struct VEmptySpace{T, AT, D <: FullDim} <: VLinearSpace{T}
    d::D
end
FullDim(v::VEmptySpace) = v.d
function emptyspace(v::VRep{T}) where {T}
    d = FullDim(v)
    return VEmptySpace{T, vvectortype(typeof(v)), typeof(d)}(d)
end

Base.in(v::VRepElement, L::VEmptySpace) = isapproxzero(v)

"""
    vrep(lines::LineIt; d::FullDim)

Creates an affine space of full dimension `d` from the list of lines `lines`.

### Examples
```julia
vrep([Line([1, 0, 0]), Line([0, 1, 0])])
```
creates the 2-dimensional affine subspace containing all the points ``(x_1, x_2, 0)``, i.e. the ``x_1````x_2``-plane.
"""
vrep(lines::LineIt; d::FullDim = FullDim_rec(lines)) = LinesHull(d, lines)

# Representation of an affine space containing the origin by the minkowsky sum of lines
mutable struct LinesHull{T, AT, D<:FullDim} <: VLinearSpace{T}
    d::D
    lines::Vector{Line{T, AT}}
    orthogonalized::Bool
    function LinesHull{T, AT, D}(d::FullDim, lines::LineIt{T}) where {T, AT, D}
        new{T, AT, D}(FullDim_convert(D, d), lazy_collect(lines), false)
    end
end
function LinesHull(d::FullDim, it::ElemIt{Line{T, AT}}) where {T, AT}
    return LinesHull{T, AT, typeof(d)}(d, it)
end
LinesHull{T, AT}(d::FullDim) where {T, AT} = LinesHull(d, Line{T, AT}[])

FullDim(v::LinesHull) = v.d
vvectortype(::Type{<:LinesHull{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:LinesHull}, d::FullDim, ::Type{T}) where T = LinesHull{T, similar_type(vvectortype(PT), d, T), typeof(d)}

function convexhull!(L::LinesHull, l::Line)
    if L.orthogonalized
        l_p = remproj(l, L)
        if !isapproxzero(l_p)
            push!(L.lines, l_p)
        end
    else
        push!(L.lines, l)
    end
    return L
end

@vecrepelem LinesHull Line lines

conetype(::Type{LinesHull{T, AT, D}}) where {T, AT, D} = RaysHull{T, AT, D}
vreptype(::Type{LinesHull{T, AT, D}}) where {T, AT, D} = Hull{T, AT, D}

# Returns a LinesHull representing the following set (TODO does it have a name?, does someone has a reference talking about it ?)
# {x | y in p => x + y in p}
function linespace(v::VRep, current=false)
    if !current
        detectvlinearity!(v)
    end
    LinesHull(FullDim(v), lines(v))
end

function detectvlinearity!(L::LinesHull, solver = nothing)
    if !L.orthogonalized
        V = detectvlinearity(L)
        L.lines = V.lines
        L.orthogonalized = true
    end
    return L
end
function detectvlinearity(L::LinesHull{T, AT}, solver = nothing) where {T, AT}
    V = LinesHull{T, AT}(FullDim(L))
    V.orthogonalized = true
    for l in lines(L)
        convexhull!(V, l)
    end
    return V
end

detectlinearity!(h::HRepresentation) = detecthlinearity!(h)
detectlinearity!(v::VRepresentation) = detectvlinearity!(v)
_linearity(h::HRepresentation) = hyperplanes(h)
_linearity(v::VRepresentation) = lines(v)
function remproj(el::RepElement, L::Union{HyperPlanesIntersection, LinesHull})
    detectlinearity!(L)
    # Modified Gram-Schmidt process, it is more stable than
    # the classical Gram-Schmidt process.
    # see https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Numerical_stability
    for lin in _linearity(L)
        el = remproj(el, lin)
    end
    return el
end
function Base.in(el::HRepElement, L::HyperPlanesIntersection)
    return isapproxzero(remproj(el, L))
end
function Base.in(el::VRepElement, L::LinesHull)
    return isapproxzero(remproj(el, L))
end

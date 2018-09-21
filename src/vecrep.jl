# collect(::Vector) does a copy.
# lazy_collect avoid this copy in case `v` is a Vector
lazy_collect(v::Vector) = v
lazy_collect(v) = collect(v)

# H-representation

"""
    hrep(hyperplanes::HyperPlaneIt, halfspaces::HalfSpaceIt)

Creates an H-representation for the polyhedron equal to the intersection of the hyperplanes `hyperplanes` and halfspaces `halfspaces`.

### Examples
For instance, the simplex
```math
\\begin{align*}
  x_1 + x_2 &= 1 \\\\
  x_1 &\\geq 0 \\\\
  x_2 &\\geq 0
\\end{align*}
```
can be created as follows:
```julia
hrep([HalfSpace([-1, 0], 0)], [HyperPlane([1, 1], 1), HalfSpace([0, -1], 0)])
```
"""
hrep(hyperplanes::HyperPlaneIt, halfspaces::HalfSpaceIt) = Intersection(hyperplanes, halfspaces)

"""
    hrep(halfspaces::HalfSpaceIt)

Creates an H-representation for the polyhedron equal to the intersection of the halfspaces `halfspaces`.

### Examples
For instance, the polytope
```math
\\begin{align*}
  x_1 + x_2 &\\leq 1 \\\\
  x_1 - x_2 &\\leq 0 \\\\
  x_1 & \\geq 0.
\\end{align*}
```
can be created as follows:
```julia
hrep([HalfSpace([1, 1], 1), HalfSpace([1, -1], 0), HalfSpace([-1, 0], 0)])
```
"""
hrep(halfspaces::ElemIt{HalfSpace{T, AT}}) where {T, AT} = hrep(HyperPlane{T, AT}[], halfspaces)

mutable struct Intersection{T, AT} <: HRepresentation{T}
    hyperplanes::HyperPlanesIntersection{T, AT}
    halfspaces::Vector{HalfSpace{T, AT}}
    function Intersection{T, AT}(hyperplanes::HyperPlaneIt{T}, halfspaces::HalfSpaceIt{T}) where {T, AT}
        new{T, AT}(HyperPlanesIntersection{T, AT}(hyperplanes), lazy_collect(halfspaces))
    end
end
Intersection(hyperplanes::ElemIt{HyperPlane{T, AT}}, halfspaces::ElemIt{HalfSpace{T, AT}}) where {T, AT} = Intersection{T, AT}(hyperplanes, halfspaces)
FullDim(h::Intersection) = FullDim_rec(h.hyperplanes, h.halfspaces)
hvectortype(::Type{Intersection{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:Intersection}, d::FullDim, ::Type{T}) where {T} = Intersection{T, similar_type(hvectortype(PT), d, T)}

Intersection(h::HRepresentation{T}) where {T} = Intersection{T}(h)
Intersection{T}(h::HRepresentation) where {T} = Intersection{T, similar_type(hvectortype(typeof(h)), T)}(h)

@subrepelem Intersection HyperPlane hyperplanes
@vecrepelem Intersection HalfSpace halfspaces

fulltype(::Type{<:Union{Intersection{T, AT}, HyperPlanesIntersection{T, AT}}}) where {T, AT} = Intersection{T, AT}

# V-representation

#"""
#    vrep(sympoints::SymPointIt)
#
#Creates a V-representation for the symmetric polytope equal to the convex hull of the symmetric points `sympoints`.
#
#### Examples
#The following creates a square
#```julia
#vrep([SymPoint([1, 1])], [SymPoint([1, -1])])
#```
#"""
#vrep(sympoints::SymPointIt) = SymPointsHull(sympoints)
#
#mutable struct SymPointsHull{T, AT} <: VSymPolytope{T}
#    sympoints::Vector{SymPoint{T, AT}}
#    function SymPointsHull{T, AT}(sympoints::SymPointIt{T}) where {T, AT}
#        new{T, AT}(lazy_collect(sympoints))
#    end
#end
#SymPointsHull(ps::ElemIt{SymPoint{T, AT}}) where {T, AT<:AbstractVector{T}} = SymPointsHull{T, AT}(collect(ps))
#vectortype(::Union{SymPointsHull{T, AT}, Type{SymPointsHull{T, AT}}}) where {T, AT} = AT
#similar_type(PT::Type{<:SymPointsHull}, d::FullDim, ::Type{T}) where {T} = SymPointsHull{T, similar_type(vectortype(PT), d, T)}
#
#SymPointsHull{T, AT}(sympoints::SymPointIt, points::PointIt, lines::LineIt, rays::RayIt) where {T, AT} = Hull{T, AT}(sympoints, points, lines, rays)
#SymPointsHull{T, AT}(sympoints::SymPointIt, lines::LineIt, rays::RayIt) where {T, AT} = Hull{T, AT}(sympoints, AT[], lines, rays)
#
#@vecrepelem SymPointsHull SymPoint sympoints
#
## SymPoint's can be split
#removevredundancy(vrep::SymPointsHull, hrep::HRep; kws...) = removevredundancy(PointsHull(sympoints(vrep), points(vrep)), hrep; kws...)

#"""
#    vrep(sympoints::SymPointIt, points::PointIt)
#
#Creates a V-representation for the polytope equal to the convex hull of the symmetric points `sympoints` and points `points`.
#
#### Examples
#The convex hull of ``(0, -1)``, ``(0, 1)`` and ``(1/2, 1/2)`` can be created as follows:
#```julia
#vrep([SymPoint([0, 1])], [[1/2, 1/2]])
#```
#"""
#vrep(sympoints::SymPointIt, points::PointIt) = PointsHull(sympoints, points)

"""
    vrep(points::PointIt)

Creates a V-representation for the polytope equal to the convex hull of the points `points`.

### Examples
The convex hull of ``(0, 0)``, ``(0, 1)`` and ``(1/2, 1/2)`` can be created as follows using exact arithmetic
```julia
vrep([[0, 0], [0, 1], [1//2, 1//2]])
```
or as follows using floating point arithmetic
```julia
vrep([[0, 0], [0, 1], [1/2, 1/2]])
```
"""
vrep(points::PointIt) = PointsHull(points)

mutable struct PointsHull{T, AT} <: VPolytope{T}
    points::Vector{AT}
    function PointsHull{T, AT}(points::PointIt) where {T, AT}
        new{T, AT}(lazy_collect(points))
    end
end
PointsHull(points::ElemIt{StaticArrays.SVector{T}}) where {T} = PointsHull{T, StaticArrays.SVector{T}}(points)
function PointsHull(points::PointIt)
    return PointsHull{coefficienttype(eltype(points)), eltype(points)}(points)
end
FullDim(v::PointsHull) = FullDim_rec(v.points)
vvectortype(::Type{PointsHull{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:PointsHull}, d::FullDim, ::Type{T}) where {T} = PointsHull{T, similar_type(vvectortype(PT), d, T)}

vreptype(::Type{PointsHull{T, AT}}) where {T, AT} = Hull{T, AT}

@vecrepelem PointsHull Point points

"""
    vrep(lines::LineIt, rays::RayIt)

Creates a V-representation for the polyhedral cone equal to the conic hull of the lines `lines` and rays `rays`.

### Examples
```julia
vrep([Line([0, 1])], [Ray([1, 0])])
```
creates a V-representation for the halfspace ``x_1 \\ge 0``.
"""
vrep(lines::LineIt, rays::RayIt) = RaysHull(lines, rays)

"""
    vrep(rays::RayIt)

Creates a V-representation for the polyhedral cone equal to the conic hull of the rays `rays`.

### Examples
```julia
vrep([Ray([1, 0]), Ray([0, 1])])
```
creates a V-representation for positive orthant.
"""
vrep(rays::ElemIt{Ray{T, AT}}) where {T, AT} = vrep(Line{T, AT}[], rays)

mutable struct RaysHull{T, AT} <: VCone{T}
    lines::LinesHull{T, AT}
    rays::Vector{Ray{T, AT}}
    function RaysHull{T, AT}(ls::LineIt{T}, rs::RayIt{T}) where {T, AT}
        new{T, AT}(LinesHull{T, AT}(ls), lazy_collect(rs))
    end
end
function RaysHull(ls::ElemIt{Line{T, AT}}, rs::ElemIt{Ray{T, AT}}) where {T, AT}
    RaysHull{T, AT}(ls, rs)
end
FullDim(v::RaysHull) = FullDim_rec(v.lines, v.rays)
vvectortype(::Type{RaysHull{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:RaysHull}, d::FullDim, ::Type{T}) where {T} = RaysHull{T, similar_type(vvectortype(PT), d, T)}

@vecrepelem RaysHull Ray rays
@subrepelem RaysHull Line lines

vreptype(::Type{RaysHull{T, AT}}) where {T, AT} = Hull{T, AT}

"""
    vrep(points::PointIt, lines::LineIt, rays::RayIt)

Creates a V-representation for the polyhedron equal to the minkowski sum of the convex hull of `points` with the conic hull of `lines` and `rays`.
"""
vrep(points::PointIt, lines::LineIt, rays::RayIt) = Hull(points, lines, rays)

vrep(points::ElemIt{AT}, lines::ElemIt{Line{T, AT}}) where {T, AT} = Hull(points, lines, Ray{T, AT}[])

mutable struct Hull{T, AT} <: VRepresentation{T}
    points::PointsHull{T, AT}
    rays::RaysHull{T, AT}
    function Hull{T, AT}(vits::VIt{T}...) where {T, AT}
        N, points, lines, rays = fillvits(vits...)
        # If points is empty and its eltype is Vector, by doing PointsHull(points), we loose the dimension information
        # If it is non-empty, we still have something type unstable
        new{T, AT}(PointsHull{T, AT}(points), RaysHull(lines, rays))
    end
end
function Hull(points::ElemIt{AT}, lines::ElemIt{Line{T, AT}}, rays::ElemIt{Ray{T, AT}}) where {T, AT}
    Hull{T, AT}(points, lines, rays)
end
FullDim(v::Hull) = FullDim_rec(v.points, v.rays)
vvectortype(::Type{Hull{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:Hull}, d::FullDim, ::Type{T}) where {T} = Hull{T, similar_type(vvectortype(PT), d, T)}

Hull(v::VRepresentation{T}) where {T} = Hull{T}(v)
Hull{T}(v::VRepresentation) where {T} = Hull{T, similar_type(vvectortype(typeof(v)), T)}(v)

@subrepelem Hull Point points
@subrepelem Hull Line rays
@subrepelem Hull Ray rays

fulltype(::Type{<:Union{Hull{T, AT}, PointsHull{T, AT}, LinesHull{T, AT}, RaysHull{T, AT}}}) where {T, AT} = Hull{T, AT}

dualtype(::Type{<:Intersection{T}}, ::Type{AT}) where {T, AT} = Hull{T, AT}
dualtype(::Type{<:Hull{T}}, ::Type{AT}) where {T, AT} = Intersection{T, AT}
const AnyIntersection{T, AT} = Union{Intersection{T, AT}, HyperPlanesIntersection{T, AT}}
function dualfullspace(h::Union{AnyIntersection, Type{<:AnyIntersection}},
                       d::FullDim, ::Type{T}, ::Type{AT}) where {T, AT}
    Hull{T, AT}([origin(AT, fulldim(d))],
                Line{T, AT}.(basis.(AT, d, 1:fulldim(d))),
                Ray{T, AT}[])
end

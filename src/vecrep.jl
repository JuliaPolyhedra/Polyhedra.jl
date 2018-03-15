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
  x_1 + x_2 &= 1 \\
  x_1 &\\geq 0 \\
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
\begin{align*}
  x_1 + x_2 &\leq 1 \\
  x_1 - x_2 &\leq 0 \\
  x_1 & \geq 0.
\end{align*}
```
can be created as follows:
```julia
hrep([HalfSpace([1, 1], 1), HalfSpace([1, -1], 0), HalfSpace([-1, 0], 0)])
```
"""
hrep(halfspaces::ElemIt{HalfSpace{N, T, AT}}) where {N, T, AT} = hrep(HyperPlane{N, T, AT}[], halfspaces)

mutable struct Intersection{N, T, AT} <: HRepresentation{N, T}
    hyperplanes::HAffineSpace{N, T, AT}
    halfspaces::Vector{HalfSpace{N, T, AT}}
    function Intersection{N, T, AT}(hyperplanes::ElemIt{HyperPlane{N, T, AT}}, halfspaces::ElemIt{HalfSpace{N, T, AT}}) where {N, T, AT}
        new{N, T, AT}(HAffineSpace(hyperplanes), lazy_collect(halfspaces))
    end
end
Intersection(hyperplanes::ElemIt{HyperPlane{N, T, AT}}, halfspaces::ElemIt{HalfSpace{N, T, AT}}) where {N, T, AT} = Intersection{N, T, AT}(hyperplanes, halfspaces)
arraytype(::Union{Intersection{N, T, AT}, Type{Intersection{N, T, AT}}}) where {N, T, AT} = AT
similar_type(PT::Type{<:Intersection}, d::FullDim{N}, ::Type{T}) where {N, T} = Intersection{N, T, similar_type(arraytype(PT), d, T)}

@subrepelem Intersection HyperPlane hyperplanes
@vecrepelem Intersection HalfSpace halfspaces

fulltype(::Type{<:Union{Intersection{N, T, AT}, HAffineSpace{N, T, AT}}}) where {N, T, AT} = Intersection{N, T, AT}

# V-representation

"""
    vrep(sympoints::SymPointIt)

Creates a V-representation for the symmetric polytope equal to the convex hull of the symmetric points `sympoints`.

### Examples
The following creates a square
```julia
vrep([SymPoint([1, 1])], [SymPoint([1, -1])])
```
"""
vrep(sympoints::SymPointIt) = SymPointsHull(sympoints)

mutable struct SymPointsHull{N, T, AT} <: VSymPolytope{N, T, AT}
    sympoints::Vector{SymPoint{N, T, AT}}
end
SymPointsHull(ps::ElemIt{SymPoint{N, T, AT}}) where {N, T, AT<:AbstractPoint{N, T}} = SymPointsHull{N, T, AT}(collect(ps))
arraytype(::Union{SymPointsHull{N, T, AT}, Type{SymPointsHull{N, T, AT}}}) where {N, T, AT} = AT
similar_type(PT::Type{<:SymPointsHull}, d::FullDim{N}, ::Type{T}) where {N, T} = SymPointsHull{N, T, similar_type(arraytype(PT), d, T)}

SymPointsHull{N, T, AT}(sympoints::SymPointIt, points::PointIt, lines::LineIt, rays::RayIt) where {N, T, AT} = Hull{N, T, AT}(sympoints, points, lines, rays)
SymPointsHull{N, T, AT}(sympoints::SymPointIt, lines::LineIt, rays::RayIt) where {N, T, AT} = Hull{N, T, AT}(sympoints, AT[], lines, rays)

@vecrepelem SymPointsHull SymPoint sympoints

# SymPoint's can be split
removevredundancy(vrep::SymPointsHull, hrep::HRep; kws...) = removevredundancy(PointsHull(sympoints(vrep), points(vrep)), hrep; kws...)

"""
    vrep(sympoints::SymPointIt, points::PointIt)

Creates a V-representation for the polytope equal to the convex hull of the symmetric points `sympoints` and points `points`.

### Examples
The convex hull of ``(0, -1)``, ``(0, 1)`` and ``(1/2, 1/2)`` can be created as follows:
```julia
vrep([SymPoint([0, 1])], [[1/2, 1/2]])
```
"""
vrep(sympoints::SymPointIt, points::PointIt) = PointsHull(sympoints, points)

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
vrep(points::PointIt) = vrep(sympointtype(points)[], points)
sympointtype(points::ElemIt{StaticArrays.SVector{N, T}}) where {N, T} = SymPoint{N, T, StaticArrays.SVector{N, T}}
function sympointtype(points::PointIt)
    isempty(points) && throw(ArgumentError("Cannot create a V-representation from an empty collection of points represented by $(eltype(points)) as the dimension cannot be computed. Use StaticArrays.SVector to represent points instead"))
    SymPoint{length(first(points)), coefficienttype(eltype(points)), eltype(points)}
end

mutable struct PointsHull{N, T, AT} <: VPolytope{N, T, AT}
    sympoints::SymPointsHull{N, T, AT}
    points::Vector{AT}
    function PointsHull{N, T, AT}(sympoints::ElemIt{SymPoint{N, T, AT}}, points::ElemIt{AT}) where {N, T, AT}
        new{N, T, AT}(SymPointsHull(sympoints), lazy_collect(points))
    end
end
PointsHull(sympoints::ElemIt{SymPoint{N, T, AT}}, points::ElemIt{AT}) where {N, T, AT} = PointsHull{N, T, AT}(sympoints, points)
arraytype(::Union{PointsHull{N, T, AT}, Type{PointsHull{N, T, AT}}}) where {N, T, AT} = AT
similar_type(PT::Type{<:PointsHull}, d::FullDim{N}, ::Type{T}) where {N, T} = PointsHull{N, T, similar_type(arraytype(PT), d, T)}

PointsHull{N, T, AT}(sympoints::SymPointIt, points::PointIt, lines::LineIt, rays::RayIt) where {N, T, AT} = Hull{N, T, AT}(sympoints, points, lines, rays)

@vecrepelem PointsHull Point points
@subrepelem PointsHull SymPoint sympoints

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
vrep(rays::ElemIt{Ray{N, T, AT}}) where {N, T, AT} = vrep(Line{N, T, AT}[], rays)

mutable struct RaysHull{N, T, AT} <: VCone{N, T, AT}
    lines::LinesHull{N, T, AT}
    rays::Vector{Ray{N, T, AT}}
    function RaysHull{N, T, AT}(ls::ElemIt{Line{N, T, AT}}, rs::ElemIt{Ray{N, T, AT}}) where {N, T, AT}
        new{N, T, AT}(LinesHull(ls), lazy_collect(rs))
    end
end
function RaysHull(ls::ElemIt{Line{N, T, AT}}, rs::ElemIt{Ray{N, T, AT}}) where {N, T, AT}
    RaysHull{N, T, AT}(ls, rs)
end
arraytype(::Union{RaysHull{N, T, AT}, Type{RaysHull{N, T, AT}}}) where {N, T, AT} = AT
similar_type(PT::Type{<:RaysHull}, d::FullDim{N}, ::Type{T}) where {N, T} = RaysHull{N, T, similar_type(arraytype(PT), d, T)}

@vecrepelem RaysHull Ray rays
@subrepelem RaysHull Line lines

"""
    vrep(sympoints::SymPointIt, points::PointIt, lines::LineIt, rays::RayIt)

Creates a V-representation for the polyhedron equal to the minkowski sum of the convex hull of `sympoints` and `points` with the conic hull of `lines` and `rays`.
"""
vrep(sympoints::SymPointIt, points::PointIt, lines::LineIt, rays::RayIt) = Hull(sympoints, points, lines, rays)

mutable struct Hull{N, T, AT} <: VRepresentation{N, T}
    points::PointsHull{N, T, AT}
    rays::RaysHull{N, T, AT}
    function Hull{N, T, AT}(vits::VIt{N, T}...) where {N, T, AT}
        sympoints, points, lines, rays = fillvits(vits...)
        new{N, T, AT}(PointsHull(sympoints, points), RaysHull(lines, rays))
    end
end
function Hull(sympoints::ElemIt{SymPoint{N, T, AT}}, points::ElemIt{AT}, lines::ElemIt{Line{N, T, AT}}, rays::ElemIt{Ray{N, T, AT}}) where {N, T, AT}
    Hull{N, T, AT}(sympoints, points, lines, rays)
end
arraytype(::Union{Hull{N, T, AT}, Type{Hull{N, T, AT}}}) where {N, T, AT} = AT
similar_type(PT::Type{<:Hull}, d::FullDim{N}, ::Type{T}) where {N, T} = Hull{N, T, similar_type(arraytype(PT), d, T)}

@subrepelem Hull SymPoint points
@subrepelem Hull Point points
@subrepelem Hull Line rays
@subrepelem Hull Ray rays

fulltype(::Type{<:Union{Hull{N, T, AT}, SymPointsHull{N, T, AT}, PointsHull{N, T, AT}, LinesHull{N, T, AT}, RaysHull{N, T, AT}}}) where {N, T, AT} = Hull{N, T, AT}

dualtype(::Type{<:Intersection{N, T}}, ::Type{AT}) where {N, T, AT} = Hull{N, T, AT}
dualtype(::Type{<:Hull{N, T}}, ::Type{AT}) where {N, T, AT} = Intersection{N, T, AT}
function dualfullspace(h::Union{Intersection, Type{<:Intersection}}, d::FullDim{N}, ::Type{T}, ::Type{AT}) where {N, T, AT}
    Hull{N, T, AT}(SymPoint{N, T, AT}[],
                   [origin(AT, d)],
                   Line{N, T, AT}.(basis.(AT, d, 1:N)),
                   Ray{N, T, AT}[])
end

# collect(::Vector) does a copy.
# lazy_collect avoid this copy in case `v` is a Vector
lazy_collect(v::Vector) = v
lazy_collect(v) = collect(v)

# H-representation

"""
    hrep(hyperplanes::HyperPlaneIt, halfspaces::HalfSpaceIt; d::FullDim)

Creates an H-representation for the polyhedron of full dimension `d` equal to
the intersection of the hyperplanes `hyperplanes` and halfspaces `halfspaces`.

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
```jldoctest
julia> hrep([HyperPlane([1, 1], 1)], [HalfSpace([0, -1], 0), HalfSpace([-1, 0], 0)])
H-representation Polyhedra.Intersection{Int64, Vector{Int64}, Int64}:
1-element iterator of HyperPlane{Int64, Vector{Int64}}:
 HyperPlane([1, 1], 1),
2-element iterator of HalfSpace{Int64, Vector{Int64}}:
 HalfSpace([0, -1], 0)
 HalfSpace([-1, 0], 0)
```
"""
function hrep(hyperplanes::HyperPlaneIt, halfspaces::HalfSpaceIt;
              d::FullDim=FullDim_rec(hyperplanes, halfspaces))
    return Intersection(d, hyperplanes, halfspaces)
end

"""
    hrep(halfspaces::HalfSpaceIt; d::FullDim)

Creates an H-representation for the polyhedron of full dimension `d` equal to
the intersection of the halfspaces `halfspaces`.

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
function hrep(halfspaces::ElemIt{HalfSpace{T, AT}}; kws...) where {T, AT}
    return hrep(HyperPlane{T, AT}[], halfspaces; kws...)
end

mutable struct Intersection{T, AT, D<:FullDim} <: HRepresentation{T}
    hyperplanes::HyperPlanesIntersection{T, AT, D}
    halfspaces::Vector{HalfSpace{T, AT}}
    function Intersection{T, AT, D}(d::FullDim, hyperplanes::HyperPlaneIt{T},
                                    halfspaces::HalfSpaceIt{T}) where {T, AT, D}
        new{T, AT, D}(HyperPlanesIntersection{T, AT, D}(d, hyperplanes),
                      lazy_collect(halfspaces))
    end
end
function Intersection(d::FullDim, hyperplanes::ElemIt{HyperPlane{T, AT}},
                      halfspaces::ElemIt{HalfSpace{T, AT}}) where {T, AT}
    return Intersection{T, AT, typeof(d)}(d, hyperplanes, halfspaces)
end

FullDim(h::Intersection) = FullDim(h.hyperplanes)
hvectortype(::Type{<:Intersection{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:Intersection}, d::FullDim, ::Type{T}) where {T} = Intersection{T, similar_type(hvectortype(PT), d, T), typeof(d)}

Intersection(h::HRepresentation{T}) where {T} = Intersection{T}(h)
Intersection{T}(h::HRepresentation) where {T} = convert(Intersection{T, similar_type(hvectortype(typeof(h)), T), typeof(FullDim(h))}, h)

function Base.intersect!(h::Intersection, hp::HyperPlane)
    intersect!(h.hyperplanes, hp)
    return h
end
function Base.intersect!(h::Intersection, hs::HalfSpace)
    push!(h.halfspaces, hs)
    return h
end

@subrepelem Intersection HyperPlane hyperplanes
@vecrepelem Intersection HalfSpace halfspaces

fulltype(::Type{<:Union{Intersection{T, AT, D}, HyperPlanesIntersection{T, AT, D}}}) where {T, AT, D} = Intersection{T, AT, D}

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
    vrep(points::PointIt; d::FullDim)

Creates a V-representation for the polytope of full dimension `d` equal to the
convex hull of the points `points`.

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
vrep(points::PointIt; d::FullDim = FullDim_rec(points)) = PointsHull(d, points)

mutable struct PointsHull{T, AT, D<:FullDim} <: VPolytope{T}
    d::D
    points::Vector{AT}
    function PointsHull{T, AT, D}(d::FullDim, points::PointIt) where {T, AT, D}
        new{T, AT, D}(FullDim_convert(D, d), lazy_collect(points))
    end
end
function PointsHull(d::FullDim, points::PointIt)
    return PointsHull{coefficient_type(eltype(points)), eltype(points),
                      typeof(d)}(d, points)
end
FullDim(v::PointsHull) = v.d
vvectortype(::Type{<:PointsHull{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:PointsHull}, d::FullDim, ::Type{T}) where {T} = PointsHull{T, similar_type(vvectortype(PT), d, T), typeof(d)}

vreptype(::Type{PointsHull{T, AT, D}}) where {T, AT, D} = Hull{T, AT, D}

convexhull!(v::PointsHull, p::AbstractVector) = push!(v.points, p)

@vecrepelem PointsHull Point points

MA.mutability(::Type{<:PointsHull}) = MA.IsMutable()
function mutable_operate_elements!(
    op::F,
    v::PointsHull,
    args::Vararg{Any,N},
) where {F<:Function,N}
    for i in eachindex(v.points)
        v.points[i] = MA.operate!(op, v.points[i], args...)
    end
    return v
end

"""
    vrep(lines::LineIt, rays::RayIt; d::FullDim)

Creates a V-representation for the polyhedral cone of full dimension `d` equal
to the conic hull of the lines `lines` and rays `rays`.

### Examples
```julia
vrep([Line([0, 1])], [Ray([1, 0])])
```
creates a V-representation for the halfspace ``x_1 \\ge 0``.
"""
function vrep(lines::LineIt, rays::RayIt; d = FullDim_rec(lines, rays))
    return RaysHull(d, lines, rays)
end

"""
    vrep(rays::RayIt)

Creates a V-representation for the polyhedral cone of full dimension `d` equal
to the conic hull of the rays `rays`.

### Examples
```julia
vrep([Ray([1, 0]), Ray([0, 1])])
```
creates a V-representation for positive orthant.
"""
function vrep(rays::ElemIt{Ray{T, AT}}; kws...) where {T, AT}
    return vrep(Line{T, AT}[], rays; kws...)
end

mutable struct RaysHull{T, AT, D<:FullDim} <: VCone{T}
    lines::LinesHull{T, AT, D}
    rays::Vector{Ray{T, AT}}
    function RaysHull{T, AT, D}(d::FullDim, ls::LineIt{T},
                                rs::RayIt{T}) where {T, AT, D}
        new{T, AT, D}(LinesHull{T, AT, D}(d, ls), lazy_collect(rs))
    end
end
function RaysHull(d::FullDim, ls::ElemIt{Line{T, AT}}, rs::ElemIt{Ray{T, AT}}) where {T, AT}
    RaysHull{T, AT, typeof(d)}(d, ls, rs)
end
FullDim(v::RaysHull) = FullDim(v.lines)
vvectortype(::Type{<:RaysHull{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:RaysHull}, d::FullDim, ::Type{T}) where {T} = RaysHull{T, similar_type(vvectortype(PT), d, T), typeof(d)}

convexhull!(v::RaysHull, l::Line) = convexhull!(v.lines, l)
convexhull!(v::RaysHull, r::Ray) = push!(v.rays, r)

@vecrepelem RaysHull Ray rays
@subrepelem RaysHull Line lines

MA.mutability(::Type{<:RaysHull}) = MA.IsMutable()
function mutable_operate_elements!(
    ::typeof(translate),
    ::RaysHull,
    ::AbstractVector,
)
end
function mutable_operate_elements!(
    op::F,
    v::RaysHull,
    args::Vararg{Any,N},
) where {F<:Function,N}
    mutable_operate_elements!(op, v.lines, args...)
    for i in eachindex(v.rays)
        v.rays[i] = MA.operate!(op, v.rays[i], args...)
    end
    return v
end

vreptype(::Type{RaysHull{T, AT, D}}) where {T, AT, D} = Hull{T, AT, D}

"""
    vrep(points::PointIt, lines::LineIt, rays::RayIt;
         d = Polyhedra.FullDim_rec(points, lines, rays))

Creates a V-representation for the polyhedron of full dimension `d` equal to the
minkowski sum of the convex hull of `points` with the conic hull of `lines` and
`rays`.
"""
function vrep(points::PointIt, lines::LineIt, rays::RayIt;
              d::FullDim = FullDim_rec(points, lines, rays))
    return Hull(d, points, lines, rays)
end

function vrep(points::ElemIt{AT}, lines::ElemIt{Line{T, AT}};
              kws...) where {T, AT}
    return vrep(points, lines, Ray{T, AT}[]; kws...)
end

mutable struct Hull{T, AT, D<:FullDim} <: VRepresentation{T}
    points::PointsHull{T, AT, D}
    rays::RaysHull{T, AT, D}
    function Hull{T, AT, D}(d::FullDim, vits::VIt{T}...) where {T, AT, D}
        points, lines, rays = fillvits(d, vits...)
        new{T, AT, D}(PointsHull(d, points), RaysHull(d, lines, rays))
    end
end
function Hull(d::FullDim, points::ElemIt{AT}, lines::ElemIt{Line{T, AT}},
              rays::ElemIt{Ray{T, AT}}) where {T, AT}
    Hull{T, AT, typeof(d)}(d, points, lines, rays)
end
FullDim(v::Hull) = FullDim(v.points)
vvectortype(::Type{<:Hull{T, AT}}) where {T, AT} = AT
similar_type(PT::Type{<:Hull}, d::FullDim, ::Type{T}) where {T} = Hull{T, similar_type(vvectortype(PT), d, T), typeof(d)}

Hull(v::VRepresentation{T}) where {T} = Hull{T}(v)
Hull{T}(v::VRepresentation) where {T} = convert(Hull{T, similar_type(vvectortype(typeof(v)), T), typeof(FullDim(v))}, v)

convexhull!(v::Hull, p::AbstractVector) = convexhull!(v.points, p)
convexhull!(v::Hull, r::VStruct) = convexhull!(v.rays, r)

@subrepelem Hull Point points
@subrepelem Hull Line rays
@subrepelem Hull Ray rays

MA.mutability(::Type{<:Hull}) = MA.IsMutable()
function mutable_operate_elements!(
    op::F,
    v::Hull,
    args::Vararg{Any,N},
) where {F<:Function,N}
    mutable_operate_elements!(op, v.points, args...)
    mutable_operate_elements!(op, v.rays, args...)
    return v
end

fulltype(::Type{<:Union{Hull{T, AT, D}, PointsHull{T, AT, D}, LinesHull{T, AT, D}, RaysHull{T, AT, D}}}) where {T, AT, D} = Hull{T, AT, D}

dualtype(::Type{<:Intersection{T, BT, D}}, ::Type{AT}) where {T, AT, BT, D} = Hull{T, AT, D}
dualtype(::Type{<:Hull{T, BT, D}}, ::Type{AT}) where {T, AT, BT, D} = Intersection{T, AT, D}
const AnyIntersection{T, AT, D} = Union{Intersection{T, AT, D}, HyperPlanesIntersection{T, AT, D}}
function dualfullspace(h::Union{AnyIntersection, Type{<:AnyIntersection}},
                       d::FullDim, ::Type{T}, ::Type{AT}) where {T, AT}
    Hull{T, AT, typeof(d)}(d, [origin(AT, fulldim(d))],
                           Line{T, AT}.(basis.(AT, Ref(d), 1:fulldim(d))),
                           Ray{T, AT}[])
end

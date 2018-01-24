abstract type VPolytope{N, T, AT} <: VRepresentation{N, T} end

@norepelem VPolytope Line
@norepelem VPolytope Ray

mutable struct PointsHull{N, T, PT<:MyPoint{N, T}} <: VRepresentation{N, T}
    points::Vector{PT}
end
PointsHull{N, T}(ps::ElemIt{PT}) where {N, T, PT<:MyPoint{N, T}} = PointsHull{N, T, PT}(collect(ps))
arraytype(::PointsHull{N, T, PT}) where {N, T, PT} = PT

@norepelem PointsHull SymPoint
@vecrepelem PointsHull Point points

mutable struct RaysHull{N, T, AT} <: VCone{N, T, AT}
    lines::VAffineSpace{N, T, AT}
    rays::Vector{Ray{N, T, AT}}
    function RaysHull{N, T, AT}(ls::ElemIt{Line{N, T, AT}}, rs::Vector{Ray{N, T, AT}}) where {N, T, AT}
        new{N, T, AT}(VAffineSpace{N, T, AT}(ls), rs)
    end
end
function RaysHull{N, T, AT}(ls::ElemIt{Line{N, T, AT}}, rs::ElemIt{Ray{N, T, AT}}) where {N, T, AT}
    RaysHull{N, T, AT}(ls, collect(rs))
end
function RaysHull{N, T}(ls::ElemIt{Line{N, T, AT}}, rs::ElemIt{Ray{N, T, AT}}) where {N, T, AT}
    RaysHull{N, T, AT}(ls, rs)
end
arraytype(::RaysHull{N, T, AT}) where {N, T, AT} = AT

@vecrepelem RaysHull Ray rays

for op in (:nlines, :startline, :doneline, :nextline)
    @eval begin
        $op(p::RaysHull, args...) = $op(p.ls, args...)
    end
end
nrays(p::RaysHull) = length(p.rs)
startray(p::RaysHull) = start(p.rs)
nextray(p::RaysHull, state) = next(p.rs, state)
doneray(p::RaysHull, state) = done(p.rs, state)

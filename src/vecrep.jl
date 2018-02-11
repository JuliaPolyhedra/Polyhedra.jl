mutable struct PointsHull{N, T, PT<:MyPoint{N, T}} <: VRepresentation{N, T}
    ps::Vector{PT}
end
PointsHull{N, T}(ps::ElemIt{PT}) where {N, T, PT<:MyPoint{N, T}} = PointsHull{N, T, PT}(collect(ps))
arraytype(::PointsHull{N, T, PT}) where {N, T, PT} = PT

nsympoints(::PointsHull) = 0
startsympoint(::PointsHull) = 0
donesympoint(::PointsHull, state) = true
npoints(p::PointsHull) = length(p.ps)
startpoint(p::PointsHull) = start(p.ps)
nextpoint(p::PointsHull, state) = next(p.ps, state)
donepoint(p::PointsHull, state) = done(p.ps, state)
nrays(::PointsHull) = 0
startray(::PointsHull) = 0
doneray(::PointsHull, state) = true
mutable struct RaysHull{N, T, AT} <: VRepresentation{N, T}
    ls::VAffineSpace{N, T, AT}
    rs::Vector{Ray{N, T, AT}}
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

for op in (:nlines, :startline, :doneline, :nextline)
    @eval begin
        $op(p::RaysHull, args...) = $op(p.ls, args...)
    end
end
nrays(p::RaysHull) = length(p.rs)
startray(p::RaysHull) = start(p.rs)
nextray(p::RaysHull, state) = next(p.rs, state)
doneray(p::RaysHull, state) = done(p.rs, state)
nsympoints(::RaysHull) = 0
startsympoint(::RaysHull) = 0
donesympoint(::RaysHull, state) = true
npoints(::RaysHull) = 0
startpoint(::RaysHull) = 0
donepoint(::RaysHull, state) = true

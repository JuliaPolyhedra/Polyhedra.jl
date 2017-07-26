mutable struct PointsHull{N, T, PT<:AbstractPoint{N, T}} <: VRepresentation{N, T}
    ps::Vector{PT}
end
PointsHull{N, T}(ps::Vector{PT}) where {N, T, PT} = PointsHull{N, T, PT}(ps)
npoints(p::PointsHull) = length(p.ps)
startpoint(p::PointsHull) = start(p.ps)
nextpoint(p::PointsHull, state) = next(p.ps, state)
donepoint(p::PointsHull, state) = done(p.ps, state)
nrays(::PointsHull) = 0
startray(::PointsHull) = 0
doneray(::PointsHull, state) = true
mutable struct RaysHull{N, T} <: VRepresentation{N, T}
    rs::Vector{AbstractRay{N, T}} # FIXME this should be Ray only
end
nrays(p::RaysHull) = length(p.rs)
startray(p::RaysHull) = start(p.rs)
nextray(p::RaysHull, state) = next(p.rs, state)
doneray(p::RaysHull, state) = done(p.rs, state)
npoints(::RaysHull) = 0
startpoint(::RaysHull) = 0
donepoint(::RaysHull, state) = true

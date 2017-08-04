export SimplePolyhedraLibrary, SimplePolyhedron

struct SimplePolyhedraLibrary{T} <: PolyhedraLibrary
end

mutable struct SimplePolyhedron{N, T} <: Polyhedron{N, T}
    hrep::Nullable{HRepresentation{N, T}}
    vrep::Nullable{VRepresentation{N, T}}
end

changefulldim{N, T}(::Type{SimplePolyhedron{N, T}}, n) = SimplePolyhedron{n, T}

getlibraryfor(p::SimplePolyhedron, N::Int, ::Type{T}) where {T} = SimplePolyhedraLibrary{T}()

function SimplePolyhedron{N, T}(rep::HRepresentation{N, T}) where {N, T}
    SimplePolyhedron{N, T}(rep, nothing)
end
function SimplePolyhedron{N, T}(rep::VRepresentation{N, T}) where {N, T}
    SimplePolyhedron{N, T}(nothing, rep)
end

function SimplePolyhedron{N, T}(rep::HRepIterator) where {N, T}
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(rep))
end
function SimplePolyhedron{N, T}(eqs::EqIterator, ineqs::IneqIterator) where {N, T}
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(eqs, ineqs))
end
function SimplePolyhedron{N, T}(eqs::EqIterator) where {N, T}
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(eqs))
end
function SimplePolyhedron{N, T}(ineqs::IneqIterator) where {N, T}
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(ineqs))
end

function SimplePolyhedron{N, T}(rep::VRepIterator) where {N, T}
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(rep))
end
function SimplePolyhedron{N, T}(points::PointIterator, rays::RayIterator) where {N, T}
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(points, rays))
end
function SimplePolyhedron{N, T}(rays::RayIterator) where {N, T}
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(rays))
end
function SimplePolyhedron{N, T}(points::PointIterator) where {N, T}
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(points))
end

function polyhedron(rep::Representation{N}, ::SimplePolyhedraLibrary{T}) where {N, T}
    SimplePolyhedron{N, T}(rep)
end
function polyhedron(repit::Union{HRepIterator{N}, VRepIterator{N}}, lib::SimplePolyhedraLibrary{T}) where {N, T}
    SimplePolyhedron{N, T}(repit)
end
function polyhedron(hps::EqIterator{N}, hss::IneqIterator{N}, ::SimplePolyhedraLibrary{T}) where {N, T}
    SimplePolyhedron{N, T}(hps, hss)
end
function polyhedron(ps::PointIterator{N}, rs::RayIterator{N}, ::SimplePolyhedraLibrary{T}) where {N, T}
    SimplePolyhedron{N, T}(ps, rs)
end

function Base.copy(p::SimplePolyhedron{N, T}) where {N, T}
    if !isnull(p.hrep)
        SimplePolyhedron{N, T}(get(p.hrep))
    else
        SimplePolyhedron{N, T}(get(p.vrep))
    end
end

function Base.push!(p::SimplePolyhedron{N}, ine::HRepresentation{N}) where N
    p.hrep = hrep(p) ∩ ine
    p.vrep = nothing
end
function Base.push!(p::SimplePolyhedron{N}, ext::VRepresentation{N}) where N
    p.vrep = convexhull(vrep(p), ext)
    p.hrep = nothing
end

hrepiscomputed(p::SimplePolyhedron) = !isnull(p.hrep)
function computehrep!(p::SimplePolyhedron)
    # vrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.hrep = doubledescription(get(p.vrep))
end
function hrep(p::SimplePolyhedron)
    if !hrepiscomputed(p)
        computehrep!(p)
    end
    get(p.hrep)
end
vrepiscomputed(p::SimplePolyhedron) = !isnull(p.vrep)
function computevrep!(p::SimplePolyhedron)
    # hrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.vrep = doubledescription(get(p.hrep))
end
function vrep(p::SimplePolyhedron)
    if !vrepiscomputed(p)
        computevrep!(p)
    end
    get(p.vrep)
end
function decomposedhfast(p::SimplePolyhedron)
    if isnull(p.hrep)
        false # what should be done here ?
    else
        decomposedhfast(get(p.hrep))
    end
end
function decomposedvfast(p::SimplePolyhedron)
    if isnull(p.vrep)
        false # what should be done here ?
    else
        decomposedvfast(get(p.vrep))
    end
end

function detecthlinearities!(p::SimplePolyhedron)
    p.hrep = removeduplicates(hrep(p))
end
function detectvlinearities!(p::SimplePolyhedron)
    p.vrep = removeduplicates(vrep(p))
end
function removehredundancy!(p::SimplePolyhedron)
    detectvlinearities!(p)
    detecthlinearities!(p)
    p.hrep = removehredundancy(hrep(p), vrep(p))
end
function removevredundancy!(p::SimplePolyhedron)
    detecthlinearities!(p)
    detectvlinearities!(p)
    p.vrep = removevredundancy(vrep(p), hrep(p))
end

for op in [:nhreps, :starthrep, :neqs, :starteq, :nineqs, :startineq]
    @eval begin
        $op(p::Union{Interval, SimplePolyhedron}) = $op(hrep(p))
    end
end
for op in [:donehrep, :nexthrep, :doneeq, :nexteq, :doneineq, :nextineq]
    @eval begin
        $op(p::Union{Interval, SimplePolyhedron}, state) = $op(hrep(p), state)
    end
end

for op in [:nvreps, :startvrep, :npoints, :startpoint, :nrays, :startray]
    @eval begin
        $op(p::Union{Interval, SimplePolyhedron}) = $op(vrep(p))
    end
end
for op in [:donevrep, :nextvrep, :donepoint, :nextpoint, :doneray, :nextray]
    @eval begin
        $op(p::Union{Interval, SimplePolyhedron}, state) = $op(vrep(p), state)
    end
end

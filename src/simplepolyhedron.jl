export SimplePolyhedraLibrary, SimplePolyhedron

struct SimplePolyhedraLibrary{T} <: PolyhedraLibrary
end

mutable struct SimplePolyhedron{N, T} <: Polyhedron{N, T}
    hrep::Nullable{HRepresentation{N, T}}
    vrep::Nullable{VRepresentation{N, T}}
end

changefulldim{N, T}(::Type{SimplePolyhedron{N, T}}, n) = SimplePolyhedron{n, T}

getlibraryfor{T}(p::SimplePolyhedron, N::Int, ::Type{T}) = SimplePolyhedraLibrary{T}()

function (::Type{SimplePolyhedron{N, T}}){N, T}(rep::HRepresentation{N, T})
    SimplePolyhedron{N, T}(rep, nothing)
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(rep::VRepresentation{N, T})
    SimplePolyhedron{N, T}(nothing, rep)
end

function (::Type{SimplePolyhedron{N, T}}){N, T}(rep::HRepIterator)
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(rep))
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(eqs::EqIterator, ineqs::IneqIterator)
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(eqs, ineqs))
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(eqs::EqIterator)
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(eqs))
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(ineqs::IneqIterator)
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(ineqs))
end

function (::Type{SimplePolyhedron{N, T}}){N, T}(rep::VRepIterator)
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(rep))
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(points::PointIterator, rays::RayIterator)
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(points, rays))
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(rays::RayIterator)
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(rays))
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(points::PointIterator)
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(points))
end

function polyhedron{N, T}(rep::Representation{N}, ::SimplePolyhedraLibrary{T})
    SimplePolyhedron{N, T}(rep)
end
function polyhedron{N, T}(repit::Union{HRepIterator{N}, VRepIterator{N}}, lib::SimplePolyhedraLibrary{T})
    SimplePolyhedron{N, T}(repit)
end
function polyhedron{N, T}(hps::EqIterator{N}, hss::IneqIterator{N}, ::SimplePolyhedraLibrary{T})
    SimplePolyhedron{N, T}(hps, hss)
end
function polyhedron{N, T}(ps::PointIterator{N}, rs::RayIterator{N}, ::SimplePolyhedraLibrary{T})
    SimplePolyhedron{N, T}(ps, rs)
end

function Base.copy{N, T}(p::SimplePolyhedron{N, T})
    if !isnull(p.hrep)
        SimplePolyhedron{N, T}(get(p.hrep))
    else
        SimplePolyhedron{N, T}(get(p.vrep))
    end
end

function Base.push!{N}(p::SimplePolyhedron{N}, ine::HRepresentation{N})
    p.hrep = hrep(p) ∩ ine
    p.vrep = nothing
end
function Base.push!{N}(p::SimplePolyhedron{N}, ext::VRepresentation{N})
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
        $op(p::SimplePolyhedron) = $op(hrep(p))
    end
end
for op in [:donehrep, :nexthrep, :doneeq, :nexteq, :doneineq, :nextineq]
    @eval begin
        $op(p::SimplePolyhedron, state) = $op(hrep(p), state)
    end
end

for op in [:nvreps, :startvrep, :npoints, :startpoint, :nrays, :startray]
    @eval begin
        $op(p::SimplePolyhedron) = $op(vrep(p))
    end
end
for op in [:donevrep, :nextvrep, :donepoint, :nextpoint, :doneray, :nextray]
    @eval begin
        $op(p::SimplePolyhedron, state) = $op(vrep(p), state)
    end
end

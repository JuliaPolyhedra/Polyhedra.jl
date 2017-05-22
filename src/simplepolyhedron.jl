export SimplePolyhedraLibrary, SimplePolyhedron

type SimplePolyhedraLibrary <: PolyhedraLibrary
end

type SimplePolyhedron{N, T} <: Polyhedron{N, T}
    hrep::Nullable{HRepresentation{N, T}}
    vrep::Nullable{VRepresentation{N, T}}
end

getlibraryfor{T}(p::SimplePolyhedron, N::Int, ::Type{T}) = SimplePolyhedraLibrary()

function (::Type{SimplePolyhedron{N, T}}){N, T}(rep::HRepresentation{N, T})
    SimplePolyhedron{N, T}(rep, nothing)
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(rep::VRepresentation{N, T})
    SimplePolyhedron{N, T}(nothing, rep)
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(rep::HRepIterator{N, T})
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(rep))
end
function (::Type{SimplePolyhedron{N, T}}){N, T}(rep::VRepIterator{N, T})
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(rep))
end


function polyhedron{N, T}(rep::Representation{N, T}, ::SimplePolyhedraLibrary)
    SimplePolyhedron{N, T}(rep)
end
function polyhedron{N, T}(repit::Union{HRepIterator{N, T}, VRepIterator{N, T}}, lib::SimplePolyhedraLibrary)
    SimplePolyhedron{N, T}(repit)
end
function polyhedron(lib::SimplePolyhedraLibrary; eqs=nothing, ineqs=nothing, points=nothing, rays=nothing)
    its = [eqs, ineqs, points, rays]
    i = findfirst(x -> !(x === nothing), its)
    if i == 0
        error("polyhedron should be given at least one iterator")
    end
    N = fulldim(its[i])
    T = typeof(its[i]).parameters[2]
    noth = eqs === nothing && ineqs === nothing
    notv = points === nothing && rays === nothing
    if !noth && !notv
        error("SimplePolyhedron constructed with a combination of eqs/ineqs with points/rays")
    end
    if notv
        rep = SimpleHRepresentation{N, T}(eqs=eqs, ineqs=ineqs)
    else
        rep = SimpleVRepresentation{N, T}(points=points, rays=rays)
    end
    polyhedron(rep, lib)
end

function Base.copy{N, T}(p::SimplePolyhedron{N, T})
    if !isnull(p.hrep)
        SimplePolyhedron{N, T}(get(p.hrep))
    else
        SimplePolyhedron{N, T}(get(p.vrep))
    end
end

function Base.push!{N}(p::SimplePolyhedron{N}, ine::HRepresentation{N})
    p.hrep = get(p.hrep) ∩ ine
end

hrepiscomputed(p::SimplePolyhedron) = !isnull(p.hrep)
hrep(p::SimplePolyhedron) = get(p.hrep) # TODO copy
vrepiscomputed(p::SimplePolyhedron) = !isnull(p.vrep)
vrep(p::SimplePolyhedron) = get(p.vrep) # TODO copy
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

for op in [:nhreps, :starthrep, :neqs, :starteq, :nineqs, :startineq]
    @eval begin
        $op(p::SimplePolyhedron) = $op(get(p.hrep))
    end
end
for op in [:donehrep, :nexthrep, :doneeq, :nexteq, :doneineq, :nextineq]
    @eval begin
        $op(p::SimplePolyhedron, state) = $op(get(p.hrep), state)
    end
end

for op in [:nvreps, :startvrep, :npoints, :startpoint, :nrays, :startray]
    @eval begin
        $op(p::SimplePolyhedron) = $op(get(p.vrep))
    end
end
for op in [:donevrep, :nextvrep, :donepoint, :nextpoint, :doneray, :nextray]
    @eval begin
        $op(p::SimplePolyhedron, state) = $op(get(p.vrep), state)
    end
end

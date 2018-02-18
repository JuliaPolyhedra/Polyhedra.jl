export SimplePolyhedraLibrary, SimplePolyhedron

struct SimplePolyhedraLibrary{T} <: PolyhedraLibrary
end

mutable struct SimplePolyhedron{N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}} <: Polyhedron{N, T}
    hrep::Nullable{HRepT}
    vrep::Nullable{VRepT}
    function SimplePolyhedron{N, T, HRepT, VRepT}(hrep::HRepT) where {N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}}
        new{N, T, HRepT, VRepT}(hrep, nothing)
    end
    function SimplePolyhedron{N, T, HRepT, VRepT}(vrep::VRepT) where {N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}}
        new{N, T, HRepT, VRepT}(nothing, vrep)
    end
end

function arraytype(p::SimplePolyhedron{N, T, HRepT, VRepT}) where{N, T, HRepT, VRepT}
    @assert arraytype(HRepT) == arraytype(VRepT)
    arraytype(HRepT)
end

similar_type(::Type{<:SimplePolyhedron{M, S, HRepT, VRepT}}, d::FullDim{N}, ::Type{T}) where {M, S, HRepT, VRepT, N, T} = SimplePolyhedron{N, T, similar_type(HRepT, d, T), similar_type(VRepT, d, T)}

getlibraryfor(p::SimplePolyhedron, N::Int, ::Type{T}) where {T} = SimplePolyhedraLibrary{T}()

function SimplePolyhedron{N, T, HRepT, VRepT}(hits::HIt...) where {N, T, HRepT, VRepT}
    SimplePolyhedron{N, T, HRepT, VRepT}(HRepT(hits...))
end
function SimplePolyhedron{N, T, HRepT, VRepT}(vits::VIt...) where {N, T, HRepT, VRepT}
    SimplePolyhedron{N, T, HRepT, VRepT}(VRepT(vits...))
end

function SimplePolyhedron{N, T}(rep::HRepresentation{N, T}) where {N, T}
    SimplePolyhedron{N, T, typeof(rep), SimpleVRepresentation{N, T}}(rep)
end
function SimplePolyhedron{N, T}(rep::VRepresentation{N, T}) where {N, T}
    SimplePolyhedron{N, T, SimpleHRepresentation{N, T}, typeof(rep)}(rep)
end

#function SimplePolyhedron{N, T}(rep::HRepIterator) where {N, T}
#    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(rep))
#end
function SimplePolyhedron{N, T}(hyperplanes::ElemIt{<:HyperPlane{N, T}}, halfspaces::ElemIt{<:HalfSpace{N, T}}) where {N, T}
    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(hyperplanes, halfspaces))
end
#function SimplePolyhedron{N, T}(hyperplanes::ElemIt{<:HyperPlane{N, T}}) where {N, T}
#    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(hyperplanes))
#end
#function SimplePolyhedron{N, T}(halfspaces::ElemIt{<:HalfSpace{N, T}}) where {N, T}
#    SimplePolyhedron{N, T}(SimpleHRepresentation{N, T}(halfspaces))
#end

#function SimplePolyhedron{N, T}(rep::VRepIterator) where {N, T}
#    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(rep))
#end
function SimplePolyhedron{N, T}(sympoints::ElemIt{<:SymPoint{N, T}}, points::ElemIt{<:MyPoint{N, T}}, lines::ElemIt{<:Line{N, T}}, rays::ElemIt{<:Ray{N, T}}) where {N, T}
    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(sympoints, points, lines, rays))
end
#function SimplePolyhedron{N, T}(rays::RayIterator) where {N, T}
#    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(rays))
#end
#function SimplePolyhedron{N, T}(points::PointIterator) where {N, T}
#    SimplePolyhedron{N, T}(SimpleVRepresentation{N, T}(points))
#end

function polyhedron(rep::Representation{N}, ::SimplePolyhedraLibrary{T}) where {N, T}
    SimplePolyhedron{N, T}(rep)
end
#function polyhedron(repit::Union{HRepIterator{N}, VRepIterator{N}}, lib::SimplePolyhedraLibrary{T}) where {N, T}
#    SimplePolyhedron{N, T}(repit)
#end
function polyhedron(hyperplanes::ElemIt{<:HyperPlane{N, T}}, halfspaces::ElemIt{<:HalfSpace{N, T}}, ::SimplePolyhedraLibrary{T}) where {N, T}
    SimplePolyhedron{N, T}(hps, hss)
end
function polyhedron(sympoints::ElemIt{<:SymPoint{N, T}}, points::ElemIt{<:MyPoint{N, T}}, lines::ElemIt{<:Line{N, T}}, rays::ElemIt{<:Ray{N, T}}, ::SimplePolyhedraLibrary{T}) where {N, T}
    SimplePolyhedron{N, T}(sympoints, points, lines, rays)
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

getlibraryfor{T<:Real}(::Type{T}) = SimplePolyhedraLibrary()

type SimplePolyhedraLibrary <: PolyhedraLibrary
end

type SimplePolyhedron{N, T} <: Polyhedron{N, T}
  hrep::Nullable{HRepresentation{N, T}}
  vrep::Nullable{HRepresentation{N, T}}
end

function polyhedron{N, T}(hrep::HRepresentation{N, T}, ::SimplePolyhedraLibrary)
  SimplePolyhedron{N, T}(hrep, nothing)
end
function polyhedron{N, T}(vrep::VRepresentation{N, T}, ::SimplePolyhedraLibrary)
  SimplePolyhedron{N, T}(nothing, vrep)
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

inequalitiesarecomputed(p::SimplePolyhedron) = !isnull(p.hrep)
getinequalities(p::SimplePolyhedron) = get(p.hrep) # TODO copy
generatorsarecomputed(p::SimplePolyhedron) = !isnull(p.vrep)
getgenerators(p::SimplePolyhedron) = get(p.vrep) # TODO copy

import Base.round, Base.eltype

# TODO affinehull and sparse

export Representation, HRepresentation, VRepresentation, fulldim

abstract Representation{N, T <: Real}
abstract HRepresentation{N,T} <: Representation{N,T}
abstract VRepresentation{N,T} <: Representation{N,T}

typealias  Rep{N,T} Union{ Representation{N,T}, Polyhedron{N,T}}
typealias HRep{N,T} Union{HRepresentation{N,T}, Polyhedron{N,T}}
typealias VRep{N,T} Union{VRepresentation{N,T}, Polyhedron{N,T}}

Base.copy(rep::Representation)    = error("not implemented")

decomposedhfast(rep::Polyhedron)          = error("not implemented")
decomposedvfast(rep::Polyhedron)          = error("not implemented")
decomposedfast(rep::HRepresentation)      = error("not implemented")
decomposedfast(rep::VRepresentation)      = error("not implemented")
decomposedhfast(rep::HRepresentation)     = decomposedfast(rep)
decomposedvfast(rep::VRepresentation)     = decomposedfast(rep)

removevredundancy!(p::VRep)               = error("not implemented")
isvredundant(p::VRep, i::Integer; strict=false, cert=false, solver = defaultLPsolverfor(p))         = error("not implemented")
detectvlinearies!(p::VRep)                = error("not implemented")

removehredundancy!(p::HRep)               = error("not implemented")
ishredundant(p::HRep, i::Integer; strict=false, cert=false, solver = defaultLPsolverfor(p))         = error("not implemented")
detecthlinearities!(p::HRep)              = error("not implemented")

Base.eltype{N,T}(rep::Rep{N,T}) = T
fulldim{N}(rep::Rep{N}) = N
Base.eltype{N,T}(::Type{Rep{N,T}}) = T
fulldim{N}(::Type{Rep{N}}) = N

# type EmptyIterator{N,T}
# end
# length(it::VRepIterator) = 0
# isempty(it::VRepIterator) = true
# eltype{N,T}(it::VRepIterator{N,T}) = T
# start(it::VRepIterator) = 0
# done(it::VRepIterator, state) = true
# next(it::VRepIterator, state) = error("This iterator is empty")

function checknext(it, i, state, donep, startp)
  while i <= length(it.ps) && (i == 0 || donep(it.pas[i], state))
    i += 1
    if i <= length(it.ps)
      state = startp(it.ps[i])
    end
  end
  i > length(it.ps) ? i : (i, state)
end

# convention: ax <= β
for (rep, HorVRep, low) in [(true, :VRep, "vrep"), (false, :VRep, "point"), (false, :VRep, "ray"), (true, :HRep, "hrep"), (false, :HRep, "ineq"), (false, :HRep, "eq")]
  if rep
    up = uppercase(low[1:2]) * low[3:end]
    shortcut = Symbol(low)
    isemp = Symbol("has" * low)
  else
    up = uppercase(low[1:1]) * low[2:end]
    shortcut = Symbol(low * "s")
    isemp = Symbol("has" * low * "s")
  end
  typename = Symbol(up * "Iterator")
  donep = Symbol("done" * low)
  startp = Symbol("start" * low)
  nextp = Symbol("next" * low)
  lenp = Symbol("n" * low * "s")
  @show typename

  @eval begin
    if !$rep
      $lenp(p::$HorVRep)   = error("$lenp not implemented")
      $startp(p::$HorVRep) = error("$startp not implemented")
      $donep(p::$HorVRep)  = error("$donep not implemented")
      $nextp(p::$HorVRep)  = error("$nextp not implemented")
    end

    type $typename{Nin,Nout,T}
      ps::Vector{$HorVRep{Nin,T}}
      f::Nullable{Function}
    end
    $typename{N,T}(ps::Vector{$HorVRep{N,T}}, f::Function=nothing) = $typename{N,N,T}(ps, f)
    $shortcut{N,T}(p::$HorVRep{N,T}, f::Function=nothing) = $typename([p], f)

    Base.length(it::$typename) = sum([$lenp(p) for p in it.ps])
    Base.isempty(it::$typename) = reduce(&, true, [$isemp(p) for p in it.ps])
    Base.eltype{N,T}(it::$typename{N,T}) = T

    Base.start(it::$typename) = checknext(it, 0, nothing, $donep, $startp)
    Base.done(it::$typename, state) = state[1] > length(it.ps)
    function Base.next(it::$typename, state)
      item, newsubstate = $nextp(it.ps[state[1]], state[2])
      newstate = checknext(it, state[1], newsubstate, $donep, $startp)
      (isnull(f) ? item : f(state[1], item), newstate)
    end

  end
end

# Default implementation for hrep and vrep
function checknext(rep::Rep, i, state, donep, startp)
  while i <= 2 && (i == 0 || donep[i](rep, state))
    i += 1
    if i <= 2
      state = startp[i](rep)
    end
  end
  i > 2 ? i : (i, state)
end

#HRep
Base.length(hrep::HRepresentation)  = nhreps(hrep)
Base.isempty(hrep::HRepresentation) = hashreps(hrep)

nhreps(hrep::HRep) = neqs(hrep) + nineqs(hrep)

haseqs(hrep::HRep)   = neqs(hrep) > 0
hasineqs(hrep::HRep) = nhrep(hrep) > 0
hashreps(hrep::HRep) = nhrep(hrep) > 0

starthrep(hrep::HRep) = checknext(hrep, 0, nothing, [doneeq, doneineq], [starteq, startineq])
donehrep(hrep::HRep, state) = state[1] > 2
function nexthrep(hrep::HRep, state)
  nextp = [nexteq, nextineq]
  item, newsubstate = nextp[state[1]](hrep, state[2])
  newstate = checknext(hrep, state[1], newsubstate, [doneeq, doneineq], [starteq, startineq])
  (item, newstate)
end

#VRep
Base.length(vrep::VRepresentation)  = nvreps(vrep)
Base.isempty(vrep::VRepresentation) = hasvreps(vrep)

nvreps(vrep::VRep) = nrays(vrep) + npoints(vrep)

hasrays(vrep::VRep)   = nrays(vrep) > 0
haspoints(vrep::VRep) = npoints(vrep) > 0
hasvreps(vrep::VRep)  = nvreps(vrep) > 0

startvrep(vrep::VRep) = checknext(vrep, 0, nothing, [doneray, donepoint], [startray, startpoint])
donevrep(vrep::VRep, state) = state[1] > 2
function nextvrep(vrep::VRep, state)
  nextp = [nextray, nextpoint]
  item, newsubstate = nextp[state[1]](vrep, state[2])
  newstate = checknext(vrep, state[1], newsubstate, [doneray, donepoint], [startray, startpoint])
  (item, newstate)
end

# Modify type
changeeltype{RepT<:Rep,T}(::Type{RepT}, ::Type{T})  = error("not implemented")
changeeltype{RepT<:Rep}(::Type{RepT}, N)            = error("not implemented")
changeboth{RepT<:Rep,T}(::Type{RepT}, N, ::Type{T}) = error("not implemented")

# Conversion
function Base.convert{Tout<:HRep, Tin<:HRepresentation}(::Type{Tout}, p::Tin)
  if fulldim(Tout) != fulldim(Tin)
    error("Different dimension")
  end
  if eltype(Tout) == eltype(Tin)
    f = nothing
  else
    S = eltype(Tout)
    f2 = (i,x) -> (AbstractVector{S}(x[1]), S(x[2]))
    f3 = (i,x) -> (AbstractVector{S}(x[1]), S(x[2]), x[3])
  end
  if decomposedfast(p)
    Tout(eqs=EqIterator(p, f2), ineqs=IneqIterator(p, f2))
  else
    Tout(HRepIterator(p, f3))
  end
end
function Base.convert{Tout<:VRep, Tin<:VRepresentation}(::Type{Tout}, p::Tin)
  if fulldim(Tout) != fulldim(Tin)
    error("Different dimension")
  end
  if eltype(Tout) == eltype(Tin)
    f = nothing
  else
    f2 = (i,x) -> (AbstractVector{eltype(Tout)}(x[1]), x[2])
    f3 = (i,x) -> (AbstractVector{eltype(Tout)}(x[1]), x[2], x[3])
  end
  if decomposedfast(p)
    Tout(points=PointIterator(p, f2), rays=RayIterator(p, f2))
  else
    Tout(VRepIterator(p, f3))
  end
end

function Base.convert{N, T, RepT<:Representation}(::Type{Representation{N, T}}, rep::RepT)
  if fulldim(RepT) != N
    error("Cannot convert representations of the same dimension")
  end
  Base.convert(changeeltype(RepT, T), rep)
end
function Base.convert{N, T, RepT<:HRepresentation}(::Type{HRepresentation{N, T}}, rep::RepT)
  if fulldim(RepT) != N
    error("Cannot convert representations of the same dimension")
  end
  Base.convert(changeeltype(RepT, T), rep)
end
function Base.convert{N, T, RepT<:VRepresentation}(::Type{VRepresentation{N, T}}, rep::RepT)
  if fulldim(RepT) != N
    error("Cannot convert representations of the same dimension")
  end
  Base.convert(changeeltype(RepT, T), rep)
end

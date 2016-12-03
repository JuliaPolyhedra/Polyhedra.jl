import Base.round, Base.eltype

# TODO affinehull and sparse

export Representation, HRepresentation, VRepresentation, fulldim
export HRepIterator, EqIterator, IneqIterator, VRepIterator, RayIterator, PointIterator
export Rep
export changeeltype, changefulldim, changeboth

abstract Representation{N, T <: Real}
abstract HRepresentation{N,T} <: Representation{N,T}
abstract VRepresentation{N,T} <: Representation{N,T}

export Rep, HRep, VRep

typealias  Rep{N,T} Union{ Representation{N,T}, Polyhedron{N,T}}
typealias HRep{N,T} Union{HRepresentation{N,T}, Polyhedron{N,T}}
typealias VRep{N,T} Union{VRepresentation{N,T}, Polyhedron{N,T}}

Base.copy(rep::Rep)            = error("copy not implemented for $(typeof(rep))")

export decomposedhfast, decomposedvfast, decomposedfast

decomposedhfast(p::Polyhedron)          = error("decomposedhfast not implemented for $(typeof(p))")
decomposedvfast(p::Polyhedron)          = error("decomposedvfast not implemented for $(typeof(p))")
decomposedfast(rep::HRepresentation)      = error("decomposedfast not implemented for $(typeof(rep))")
decomposedfast(rep::VRepresentation)      = error("decomposedfast not implemented for $(typeof(rep))")
decomposedhfast(rep::HRepresentation)     = decomposedfast(rep)
decomposedvfast(rep::VRepresentation)     = decomposedfast(rep)

export removevredundancy!, removehredundancy!, detecthlinearities!, detectvlinearities!, isvredundant, ishredundant

removevredundancy!(p::VRep)               = error("removevredundancy! not implemented for $(typeof(p))")
isvredundant(p::VRep, i::Integer; strongly=false, cert=false, solver = defaultLPsolverfor(p))         = error("not implemented for $(typeof(p))")
detectvlinearities!(p::VRep)                = error("detectvlinearities! not implemented for $(typeof(p))")

removehredundancy!(p::HRep)               = error("removehredundancy! not implemented for $(typeof(p))")
ishredundant(p::HRep, i::Integer; strongly=false, cert=false, solver = defaultLPsolverfor(p))         = error("ishredundant not implemented for $(typeof(p))")
detecthlinearities!(p::HRep)              = error("detecthlinearities! not implemented for $(typeof(p))")

Base.eltype{N,T}(rep::Rep{N,T}) = T
fulldim{N}(rep::Rep{N}) = N
Base.eltype{RepT<:Rep}(::Type{RepT}) = RepT.parameters[2]
fulldim{RepT<:Rep}(::Type{RepT}) = RepT.parameters[1]

# type EmptyIterator{N,T}
# end
# length(it::VRepIterator) = 0
# isempty(it::VRepIterator) = true
# eltype{N,T}(it::VRepIterator{N,T}) = T
# start(it::VRepIterator) = 0
# done(it::VRepIterator, state) = true
# next(it::VRepIterator, state) = error("This iterator is empty")

function checknext(it, i, state, donep, startp)
  while i <= length(it.ps) && (i == 0 || donep(it.ps[i], state))
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
  else
    up = uppercase(low[1:1]) * low[2:end]
  end
  typename = Symbol(up * "Iterator")
  donep = Symbol("done" * low)
  startp = Symbol("start" * low)
  nextp = Symbol("next" * low)
  shortcuts = low * "s"
  shortcut = Symbol(shortcuts)
  lenp = Symbol("n" * shortcuts)
  isemp = Symbol("has" * shortcuts)

  @eval begin
    export $shortcut, $lenp, $startp, $donep, $nextp, $isemp
    if !$rep
      $lenp(p::$HorVRep)   = error("$($lenp) not implemented for $(typeof(p))")
      $startp(p::$HorVRep) = error("$($startp) not implemented for $(typeof(p))")
      $donep(p::$HorVRep)  = error("$($donep) not implemented for $(typeof(p))")
      $nextp(p::$HorVRep)  = error("$($nextp) not implemented for $(typeof(p))")
    end

    type $typename{Nout, Tout, Nin, Tin}
      ps::Vector
      f::Nullable{Function}
      function $typename{RepT<:$HorVRep}(ps::Vector{RepT}, f)
        new(ps, f)
      end
    end
    $typename{RepT<:$HorVRep}(ps::Vector{RepT}, f=nothing) = $typename{fulldim(RepT),eltype(RepT),fulldim(RepT),eltype(RepT)}(ps, f)
    $shortcut{N,T}(p::$HorVRep{N,T}, f=nothing) = $typename([p], f)

    Base.length(it::$typename) = sum([$lenp(p) for p in it.ps])
    Base.isempty(it::$typename) = reduce(&, true, [$isemp(p) for p in it.ps])
    fulldim{N}(it::$typename{N}) = N
    Base.eltype{N,T}(it::$typename{N,T}) = T

    Base.start(it::$typename) = checknext(it, 0, nothing, $donep, $startp)
    Base.done(it::$typename, state) = state[1] > length(it.ps)
    function Base.next(it::$typename, state)
      item, newsubstate = $nextp(it.ps[state[1]], state[2])
      newstate = checknext(it, state[1], newsubstate, $donep, $startp)
      (isnull(it.f) ? item : get(it.f)(state[1], item), newstate)
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
hasineqs(hrep::HRep) = nhreps(hrep) > 0
hashreps(hrep::HRep) = nhreps(hrep) > 0

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

# Linearity Set
export linset
function linset(rep::HRepresentation)
  s = IntSet()
  for (i,h) in enumerate(hreps(rep))
    if islin(h)
      push!(s, i)
    end
  end
  s
end
function linset(rep::VRepresentation)
  s = IntSet()
  for (i,v) in enumerate(vreps(rep))
    if islin(v)
      push!(s, i)
    end
  end
  s
end

# Modify type
changeeltype{RepT<:Rep,T}(::Type{RepT}, ::Type{T})  = error("changeeltype not implemented for $(RepT)")
changefulldim{RepT<:Rep}(::Type{RepT}, N)           = error("changefulldim not implemented for $(RepT)")
changeboth{RepT<:Rep,T}(::Type{RepT}, N, ::Type{T}) = error("changeboth not implemented for $(RepT)")
lazychangeeltype{RepT<:Rep,T}(::Type{RepT}, ::Type{T}) = eltype(RepT) == T ? RepT : changeleltype(RepT, T)
lazychangefulldim{RepT<:Rep}(::Type{RepT}, N)          = fulldim(RepT) == N ? RepT : changelfulldim(RepT, N)
function lazychangeboth{RepT<:Rep,T}(::Type{RepT}, N, ::Type{T})
    if eltype(RepT) == T
        lazychangefulldim(RepT, N)
    elseif fulldim(RepT) == N
        changeeltype(RepT, T)
    else
        changeboth(RepT, N, T)
    end
end

# Conversion
changeeltype{S,N}(::Type{S}, x::HalfSpace{N}) = HalfSpace{N,S}(changeeltype(S, x.a), S(β))
changeeltype{S,N}(::Type{S}, x::HyperPlane{N}) = HyperPlane{N,S}(changeeltype(S, x.a), S(β))

changeeltype{S,N}(::Type{S}, x::Vec{N}) = Vec{N,S}(x)
changeeltype{S,N}(::Type{S}, x::Ray{N}) = Ray{N,S}(changeeltype(S, x.r))
changeeltype{S,N}(::Type{S}, x::Line{N}) = Line{N,S}(changeeltype(S, x.r))
changeeltype{S,N}(::Type{S}, x::Point{N}) = Point{N,S}(x)
changeeltype{S,N}(::Type{S}, x::SymPoint{N}) = SymPoint{N,S}(changeeltype(S, x))
changeeltype{S}(::Type{S}, x::AbstractVector) = AbstractVector{S}(x)

# This method solves the ambiguity with the following methods and the general method
# Base.convert{T}(::Type{T}, p::T) = p
Base.convert{T<:HRepresentation}(::Type{T}, p::T) = p
Base.convert{T<:VRepresentation}(::Type{T}, p::T) = p

function hconvert{RepTout<:HRep, RepTin<:HRep}(::Type{RepTout}, p::RepTin)
  Nin  = fulldim(RepTin)
  Nout = fulldim(RepTout)
  if Nin != Nout
    error("Different dimension")
  end
  Tin  = eltype(RepTin)
  Tout = eltype(RepTout)
  if Tin == Tout
    f = nothing
  else
    f = (i,x) -> changeeltype(typeof(x), Tout)(x)
  end
  if decomposedhfast(p)
    RepTout(eqs=EqIterator{Nout,Tout,Nin,Tin}([p], f), ineqs=IneqIterator{Nout,Tout,Nin,Tin}([p], f))
  else
    RepTout(HRepIterator{Nout,Tout,Nin,Tin}([p], f))
  end
end

Base.convert{RepTout<:HRep, RepTin<:HRepresentation}(::Type{RepTout}, p::RepTin) = hconvert(RepTout, p)
Base.convert{RepTout<:HRepresentation, RepTin<:HRep}(::Type{RepTout}, p::RepTin) = hconvert(RepTout, p)
# avoid ambiguity
Base.convert{RepTout<:HRepresentation, RepTin<:HRepresentation}(::Type{RepTout}, p::RepTin) = hconvert(RepTout, p)

function vconvert{RepTout<:VRep, RepTin<:VRep}(::Type{RepTout}, p::RepTin)
  Nin  = fulldim(RepTin)
  Nout = fulldim(RepTout)
  if Nin != Nout
    error("Different dimension")
  end
  Tin  = eltype(RepTin)
  Tout = eltype(RepTout)
  if Tin == Tout
    f = nothing
  else
    f = (i,x) -> changeeltype(typeof(x), Tout)(x)
  end
  if decomposedvfast(p)
    RepTout(points=PointIterator{Nout,Tout,Nin,Tin}([p], f), rays=RayIterator{Nout,Tout,Nin,Tin}([p], f))
  else
    RepTout(VRepIterator{Nout,Tout,Nin,Tin}([p], f))
  end
end

Base.convert{RepTout<:VRep, RepTin<:VRepresentation}(::Type{RepTout}, p::RepTin) = vconvert(RepTout, p)
Base.convert{RepTout<:VRepresentation, RepTin<:VRep}(::Type{RepTout}, p::RepTin) = vconvert(RepTout, p)
# avoid ambiguity
Base.convert{RepTout<:VRepresentation, RepTin<:VRepresentation}(::Type{RepTout}, p::RepTin) = vconvert(RepTout, p)

changeeltype{N,T}(p::Rep{N,T}, ::Type{T}) = p
function changeeltype{RepTin<:Rep, Tout}(p::RepTin, ::Type{Tout})
    RepTout = changeeltype(RepTin, Tout)
    RepTout(p)
end

# FIXME it does not get called. The calls always go throug vconvert and hconvert. Use changeeltype instead
#function Base.convert{N, T, RepT<:Representation}(::Type{Representation{N, T}}, rep::RepT)
#  if fulldim(RepT) != N
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changeeltype(RepT, T), rep)
#end
#function Base.convert{N, T, RepT<:HRepresentation}(::Type{HRepresentation{N, T}}, rep::RepT)
#  if fulldim(RepT) != N
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changeeltype(RepT, T), rep)
#end
#function Base.convert{M, S, RepT<:VRepresentation}(::Type{VRepresentation{M, S}}, rep::RepT)
#  if fulldim(RepT) != M
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changeeltype(RepT, S), rep)
#end

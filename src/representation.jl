import Base.eltype

export Representation, HRepresentation, VRepresentation, fulldim
export AbstractRepIterator, AbstractHRepIterator, HRepIterator, EqIterator, IneqIterator, AbstractVRepIterator, VRepIterator, LineIterator, RayIterator, PointIterator
export Rep
export changeeltype, changefulldim, changeboth

abstract type Representation{N, T <: Real} end
abstract type HRepresentation{N,T} <: Representation{N,T} end
abstract type VRepresentation{N,T} <: Representation{N,T} end

export Rep, HRep, VRep

const  Rep{N,T} = Union{ Representation{N,T}, Polyhedron{N,T}}
const HRep{N,T} = Union{HRepresentation{N,T}, Polyhedron{N,T}}
const VRep{N,T} = Union{VRepresentation{N,T}, Polyhedron{N,T}}

Base.copy(rep::Rep)            = error("copy not implemented for $(typeof(rep))")

export decomposedhfast, decomposedvfast, decomposedfast

decomposedhfast(p::Polyhedron)        = error("decomposedhfast not implemented for $(typeof(p))")
decomposedvfast(p::Polyhedron)        = error("decomposedvfast not implemented for $(typeof(p))")
decomposedfast(rep::HRepresentation)  = error("decomposedfast not implemented for $(typeof(rep))")
decomposedfast(rep::VRepresentation)  = error("decomposedfast not implemented for $(typeof(rep))")
decomposedhfast(rep::HRepresentation) = decomposedfast(rep)
decomposedvfast(rep::VRepresentation) = decomposedfast(rep)

Base.eltype(rep::Rep{N,T}) where {N,T} = T
"""
    fulldim(rep::Rep)

Returns the dimension of the space in which the representation is defined.
That is, a straight line in a 3D space has `fulldim` 3.
"""
fulldim(rep::Rep{N}) where {N} = N
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
abstract type AbstractRepIterator{N, T} end
abstract type AbstractHRepIterator{N, T} <: AbstractRepIterator{N, T} end
abstract type AbstractVRepIterator{N, T} <: AbstractRepIterator{N, T} end

for (rep, HorVRep, elt, low) in [(true, :VRep, :VRepElement, "vrep"), (false, :VRep, :AbstractPoint, "point"), (false, :VRep, :Line, "line"), (false, :VRep, :AbstractRay, "ray"), (true, :HRep, :HRepElement, "hrep"), (false, :HRep, :HalfSpace, "ineq"), (false, :HRep, :HyperPlane, "eq")]
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
    abstractit = Symbol("Abstract" * string(HorVRep) * "Iterator")

    @eval begin
        export $shortcut, $lenp, $startp, $donep, $nextp, $isemp
        if !$rep
            $lenp(p::$HorVRep)   = error("$($lenp) not implemented for $(typeof(p))")
            $startp(p::$HorVRep) = error("$($startp) not implemented for $(typeof(p))")
            $donep(p::$HorVRep)  = error("$($donep) not implemented for $(typeof(p))")
            $nextp(p::$HorVRep)  = error("$($nextp) not implemented for $(typeof(p))")
        end

        struct $typename{Nout, Tout, Nin, Tin} <: $abstractit{Nout, Tout}
            ps::Vector
            f::Nullable{Function}
            function $typename{Nout, Tout, Nin, Tin}(ps::Vector, f=nothing) where {Nout, Tout, Nin, Tin}
                new{Nout, Tout, Nin, Tin}(ps, f)
            end
        end
        function $typename(ps::Vector{RepT}, f=nothing) where {RepT<:$HorVRep}
            $typename{fulldim(RepT),eltype(RepT),fulldim(RepT),eltype(RepT)}(ps, f)
        end
        function $shortcut{N,T}(p::$HorVRep{N,T}, f=nothing)
            $typename([p], f)
        end

        Base.length(it::$typename) = sum([$lenp(p) for p in it.ps])
        Base.isempty(it::$typename) = !reduce(|, false, $isemp.(it.ps))
        fulldim{N}(it::$typename{N}) = N
        Base.eltype{N, T}(it::$typename{N, T}) = $elt{N, T}

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

"""
    hreps(hr::HRep)

Returns an iterator over the elements of the H-representation.

### Note

This is type unstable as the iterator returns both halfspaces and hyperplanes.
It is therefore more efficient to call [`eqs`](@ref) and [`ineqs`](@ref) separately.
"""
function hreps end

"""
    eqs(hr::HRep)

Returns an iterator over the hyperplanes of the H-representation.
"""
function eqs end

"""
    ineqs(hr::HRep)

Returns an iterator over the halfspaces of the H-representation.
"""
function ineqs end

"""
    nhreps(hr::HRep)

Returns the number of halfspaces and hyperplanes of the H-representation.

### Note

Note that it does not do redundancy removal so it is not the minimal number of halfspace and hyperplanes needed to represent the polyhedron, it is simply the number that are currently used.
"""
nhreps(hrep::HRep) = neqs(hrep) + nineqs(hrep)

"""
    neqs(hr::HRep)

Returns the number of hyperplanes of the H-representation.
"""
function neqs end

"""
    nineqs(hr::HRep)

Returns the number of halfspaces of the H-representation.
"""
function nineqs end

"""
    haseqs(hr::HRep)

Returns whether the H-representation contain any hyperplane.
"""
haseqs(hrep::HRep)   = neqs(hrep) > 0
"""
    hasineqs(hr::HRep)

Returns whether the H-representation contain any halfspace.
"""
hasineqs(hrep::HRep) = nineqs(hrep) > 0
"""
    hashreps(hr::HRep)

Returns whether the H-representation contain any halfspace or hyperplane.
"""
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
#nlines(vrep::VRep) = sum(map(islin, rays(vrep))) # TODO: call detectvlinearity! before

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
lazychangefulldim{RepT<:Rep}(::Type{RepT}, N)          = fulldim(RepT) == N ? RepT : changefulldim(RepT, N)
# TODO in Julia v0.6, we can do {N1, T, RepT<:Rep{N2, T}} and {N, T1, RepT<:Rep{N, T2}} and remove lazychangeboth
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
        RepTout(EqIterator{Nout,Tout,Nin,Tin}([p], f), IneqIterator{Nout,Tout,Nin,Tin}([p], f))
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
        RepTout(PointIterator{Nout,Tout,Nin,Tin}([p], f), RayIterator{Nout,Tout,Nin,Tin}([p], f))
    else
        RepTout(VRepIterator{Nout,Tout,Nin,Tin}([p], f))
    end
end

Base.convert{RepTout<:VRep, RepTin<:VRepresentation}(::Type{RepTout}, p::RepTin) = vconvert(RepTout, p)
Base.convert{RepTout<:VRepresentation, RepTin<:VRep}(::Type{RepTout}, p::RepTin) = vconvert(RepTout, p)
# avoid ambiguity
Base.convert{RepTout<:VRepresentation, RepTin<:VRepresentation}(::Type{RepTout}, p::RepTin) = vconvert(RepTout, p)

changeeltype(p::Rep{N,T}, ::Type{T}) where {N,T} = p
function changeeltype(p::RepTin, ::Type{Tout}) where {RepTin<:Rep, Tout}
    RepTout = changeeltype(RepTin, Tout)
    RepTout(p)
end

VRepresentation{N, T}(v::RepTin) where {N, T, RepTin} = lazychangeboth(RepTin, N, T)(v)
HRepresentation{N, T}(h::RepTin) where {N, T, RepTin} = lazychangeboth(RepTin, N, T)(h)

VRep{N, T}(v::VRepresentation) where {N, T} = VRepresentation{N, T}(v)
VRep{N, T}(p::Polyhedron) where {N, T} = Polyhedron{N, T}(p)
HRep{N, T}(h::HRepresentation) where {N, T} = HRepresentation{N, T}(h)
HRep{N, T}(p::Polyhedron) where {N, T} = Polyhedron{N, T}(p)

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

export Representation, HRepresentation, VRepresentation, fulldim
export RepIterator
export Rep

abstract type Representation{N, T <: Real} end
abstract type HRepresentation{N,T} <: Representation{N,T} end
abstract type VRepresentation{N,T} <: Representation{N,T} end

export Rep, HRep, VRep

const  Rep{N,T} = Union{ Representation{N,T}, Polyhedron{N,T}}
const HRep{N,T} = Union{HRepresentation{N,T}, Polyhedron{N,T}}
const VRep{N,T} = Union{VRepresentation{N,T}, Polyhedron{N,T}}

Base.copy(rep::Rep)            = error("copy not implemented for $(typeof(rep))")

MultivariatePolynomials.coefficienttype(rep::Union{Rep{N,T}, Type{<:Rep{N,T}}}) where {N,T} = T
"""
    fulldim(rep::Rep)

Returns the dimension of the space in which the representation is defined.
That is, a straight line in a 3D space has `fulldim` 3.
"""
fulldim(p) = fulldim(FullDim(p))
FullDim(rep::Union{Rep{N}, Type{<:Rep{N}}}) where N = FullDim{N}()

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

#abstract type AbstractRepIterator{N, T} end
#abstract type AbstractHRepIterator{N, T} <: AbstractRepIterator{N, T} end
#abstract type AbstractVRepIterator{N, T} <: AbstractRepIterator{N, T} end

# Subtyping AbstractVector{ElemT} make Base think that RepIterator implements indexing e.g. for copy!
struct RepIterator{N, T, ElemT, PT}
    ps::PT
    f::Nullable{Function}
    function RepIterator{N, T, ElemT}(ps::PT, f=nothing) where {N, T, ElemT, PT<:Tuple}
        new{N, T, ElemT, PT}(ps, f)
    end
end

function Base.first(it::RepIterator)
    next(it, start(it))[1]
end
FullDim(it::RepIterator{N}) where N = FullDim{N}()
Base.eltype(it::RepIterator{N, T, ElemT}) where {N, T, ElemT} = ElemT
function typed_map(f, d::FullDim{N}, ::Type{T}, it::RepIterator{Nin, Tin, ElemT}) where {N, T, Nin, Tin, ElemT}
    @assert isnull(it.f)
    RepIterator{N, T, similar_type(ElemT, d, T)}(it.ps, f)
end

function RepIterator{N, T}(it::RepIterator) where {N, T}
    typed_map((i,x) -> similar_type(typeof(x), T)(x), FullDim{N}(), T, it)
end

# FIXME the variables need to be defined outside of the local scope of for
#       for Julia to know them inside the """ ... """ of the docstrings
HorV = HorVRep = horvrep = singular = singularlin = plural = plurallin = lenp = isnotemptyp = repexem = listexem = :_

for (isVrep, elt, singular) in [(true, :SymPoint, :sympoint), (true, :MyPoint, :point),
                                (true, :Line, :line), (true, :Ray, :ray),
                                (false, :HyperPlane, :hyperplane), (false, :HalfSpace, :halfspace)]
    if isVrep
        HorV = :V
        HorVRep = :VRep
        horvrep = :vrep
    else
        HorV = :H
        HorVRep = :HRep
        horvrep = :hrep
    end
    typename = :(RepIterator{N, T, <:$elt})
    singularstr = string(singular)
    elemtype = Symbol(singularstr * "type")
    donep = Symbol("done" * singularstr)
    startp = Symbol("start" * singularstr)
    nextp = Symbol("next" * singularstr)
    pluralstr = singularstr * "s"
    plural = Symbol(pluralstr)
    lenp = Symbol("n" * pluralstr)
    isnotemptyp = Symbol("has" * pluralstr)
    mapit = Symbol("map" * pluralstr)

    @eval begin
        export $plural, $lenp, $isnotemptyp, $startp, $donep, $nextp

        """
            $plural($horvrep::$HorVRep)

        Returns an iterator over the $plural of the $HorV-representation `$horvrep`.
        """
        function $plural end

        """
            $lenp($horvrep::$HorVRep)

        Returns the number of $plural of the $HorV-representation `$horvrep`.
        """
        function $lenp end

        """
            $isnotemptyp($horvrep::$HorVRep)

        Returns whether the $HorV-representation `$horvrep` has any $singular.
        """
        $isnotemptyp($horvrep::$HorVRep) = !iszero($lenp($horvrep))

        function $startp end
        function $donep end
        function $nextp end

        if $singularstr == "point"
            $elemtype(p::$HorVRep) = arraytype(p)
        else
            $elemtype(p::$HorVRep{N, T}) where {N, T} = $elt{N, T, arraytype(p)}
        end

        function $plural(p::$HorVRep{N, T}...) where {N, T}
            ElemT = promote_type($elemtype.(p)...)
            RepIterator{N, T, ElemT}(p)
        end

        function $mapit(f::Function, d::FullDim{N}, ::Type{T}, p::$HorVRep...) where {N, T}
            ElemT = promote_type(similar_type.($elemtype.(p), d, T)...)
            RepIterator{N, T, ElemT}(p, f)
        end

        Base.length(it::$typename) where {N, T} = sum($lenp, it.ps)
        Base.isempty(it::$typename) where {N, T} = !any($isnotemptyp.(it.ps))

        Base.start(it::$typename) where {N, T} = checknext(it, 0, nothing, $donep, $startp)
        Base.done(it::$typename, state) where {N, T} = state[1] > length(it.ps)
        function Base.next(it::$typename, state) where {N, T}
            item, newsubstate = $nextp(it.ps[state[1]], state[2])
            newstate = checknext(it, state[1], newsubstate, $donep, $startp)
            (isnull(it.f) ? item : get(it.f)(state[1], item), newstate)
        end
    end
end

# Combines an element type with its linear version.
# e.g. combines points with the sympoints by splitting sympoints in two points.
struct AllRepIterator{N, T, ElemT, LinElemT, PT}
    itlin::RepIterator{N, T, LinElemT, PT}
    it::RepIterator{N, T, ElemT, PT}
end

Base.eltype(it::AllRepIterator{N, T, ElemT}) where {N, T, ElemT} = ElemT
Base.length(it::AllRepIterator) = 2length(it.itlin) + length(it.it)
Base.isempty(it::AllRepIterator) = isempty(it.itlin) && isempty(it.it)

function checknext(it::AllRepIterator, i, state)
    while i <= 3 && ((i <= 2 && done(it.itlin, state)) || (i == 3 && done(it.it, state)))
        i += 1
        if i <= 2
            state = start(it.itlin)
        elseif i == 3
            state = start(it.it)
        end # Otherwise we leave state as it is to be type stable
    end
    i, state
end

splitlin(h::HyperPlane, i) = (i == 1 ? HalfSpace(h.a, h.β) : HalfSpace(-h.a, -h.β))
splitlin(s::SymPoint, i) = (i == 1 ? coord(s) : -coord(s))
splitlin(l::Line, i) = (i == 1 ? Ray(coord(l)) : Ray(-coord(l)))

Base.start(it::AllRepIterator) = checknext(it, 1, start(it.itlin))
Base.done(::AllRepIterator, state) = state[1] > 3
function Base.next(it::AllRepIterator, istate)
    i, state = istate
    if i <= 2
        @assert i >= 1
        itemlin, newstate = next(it.itlin, state)
        item = splitlin(itemlin, i)
    else
        @assert i == 3
        item, newstate = next(it.it, state)
    end
    newistate = checknext(it, i, newstate)
    (item, newistate)
end

for (isVrep, singularlin, singular, repexem, listexem) in [(true, :sympoint, :point, "convexhull(SymPoint([1, 0]), [0, 1])", "[1, 0], [-1, 0], [0, 1]"),
                                                           (true, :line, :ray, "Line([1, 0]) + Ray([0, 1])", "Ray([1, 0]), Ray([-1, 0]), Ray([0, 1])"),
                                                           (false, :hyperplane, :halfspace, "HyperPlane([1, 0], 1) ∩ HalfSpace([0, 1], 1)", "HalfSpace([1, 0]), HalfSpace([-1, 0]), HalfSpace([0, 1])")]
    if isVrep
        HorV = :V
        HorVRep = :VRep
        horvrep = :vrep
    else
        HorV = :H
        HorVRep = :HRep
        horvrep = :hrep
    end
    pluralstrlin = string(singularlin) * "s"
    plurallin = Symbol(pluralstrlin)
    lenplin = Symbol("n" * pluralstrlin)
    isnotemptyplin = Symbol("has" * pluralstrlin)
    pluralstr = string(singular) * "s"
    plural = Symbol(pluralstr)
    lenp = Symbol("n" * pluralstr)
    isnotemptyp = Symbol("has" * pluralstr)
    allpluralstr = "all" * pluralstr
    allplural = Symbol(allpluralstr)
    alllenp = Symbol("n" * allpluralstr)
    allisnotemptyp = Symbol("has" * allpluralstr)
    @eval begin
        export $allplural, $alllenp, $allisnotemptyp

        """
            all$plural($horvrep::$HorVRep)

        Returns an iterator over the $plural and $plurallin in the $HorV-representation `$horvrep` splitting $plurallin in two $plural.

        ### Examples

        ```julia
        $horvrep = $repexem
        collect(all$plural($horvrep)) # Returns [$listexem]
        ```
        """
        $allplural(p::$HorVRep...) = AllRepIterator($plurallin(p...), $plural(p...))

        """
            nall$plural($horvrep::$HorVRep)

        Returns the number of $plural plus twice the number of $plurallin in the $HorV-representation `$horvrep`, i.e. `length(all$plural($horvrep))`
        """
        $alllenp($horvrep::$HorVRep) = 2 * $lenplin($horvrep) + $lenp($horvrep)

        """
            hasall$plural($horvrep::$HorVRep)

        Returns whether the $HorV-representation `$horvrep` contains any $singular or $singularlin.
        """
        $allisnotemptyp($horvrep::$HorVRep) = $isnotemptyplin($horvrep) || $isnotemptyp($horvrep)
    end
end

const ElemIt{ElemT} = Union{AllRepIterator{<:Any, <:Any, ElemT}, RepIterator{<:Any, <:Any, ElemT}, AbstractVector{ElemT}}

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

function hreps(p::HRep{N, T}...) where {N, T}
    hyperplanes(p...), halfspaces(p...)
end

function hmap(f, d::FullDim, ::Type{T}, p::HRep...) where T
    maphyperplanes(f, d, T, p...), maphalfspaces(f, d, T, p...)
end

function hconvert(::Type{RepTout}, p::HRep{N, T}) where {N, T, RepTout<:HRep{N, T}}
    RepTout(hreps(p)...)
end
function hconvert(::Type{RepTout}, p::HRep{N}) where {N, T, RepTout<:HRep{N, T}}
    RepTout(RepIterator{N, T}.(hreps(p))...)
end
#function hconvert(::Type{RepTout}, p::HRep{N, Tin}){N, Tin, Tout, RepTout<:HRep{N, Tout}}
#    if Tin == Tout
#        f = nothing
#    else
#        f = (i,x) -> similar_type(typeof(x), Tout)(x)
#    end
#    if decomposedhfast(p)
#        RepTout(EqIterator{Nout,Tout}((p,), f), IneqIterator{Nout,Tout}((p,), f))
#    else
#        RepTout(HRepIterator{Nout,Tout}((p,), f))
#    end
#end

# This method solves the ambiguity with the following methods and the general method
# Base.convert{T}(::Type{T}, p::T) = p
Base.convert{T<:HRepresentation}(::Type{T}, p::T) = p
Base.convert{T<:VRepresentation}(::Type{T}, p::T) = p

Base.convert(::Type{RepTout}, p::HRepresentation) where RepTout<:HRep = hconvert(RepTout, p)
Base.convert(::Type{RepTout}, p::HRep) where {RepTout<:HRepresentation} = hconvert(RepTout, p)
# avoid ambiguity
Base.convert(::Type{RepTout}, p::HRepresentation) where {RepTout<:HRepresentation} = hconvert(RepTout, p)

function Polyhedron{N, S}(p::Polyhedron{N, T}) where {N, S, T}
    RepTout = similar_type(typeof(p), S)
    if !hrepiscomputed(p) && vrepiscomputed(p)
        vconvert(RepTout, p)
    else
        hconvert(RepTout, p)
    end
end

function vreps(p...)
    preps(p...)..., rreps(p...)...
end
function preps(p::VRep{N, T}...) where {N, T}
    sympoints(p...), points(p...)
end
function rreps(p::VRep{N, T}...) where {N, T}
    lines(p...), rays(p...)
end

function vmap(f, d::FullDim, ::Type{T}, p::VRep...) where T
    mapsympoints(f, d, T, p...), mappoints(f, d, T, p...), maplines(f, d, T, p...), maprays(f, d, T, p...)
end

function vconvert(::Type{RepTout}, p::VRep{N, T}) where {N, T, RepTout<:VRep{N, T}}
    RepTout(vreps(p)...)
end
function vconvert(::Type{RepTout}, p::VRep{N}) where {N, T, RepTout<:VRep{N, T}}
    RepTout(RepIterator{N, T}.(vreps(p))...)
end

#function vconvert{RepTout<:VRep, RepTin<:VRep}(::Type{RepTout}, p::RepTin)
#    Nin  = fulldim(RepTin)
#    Nout = fulldim(RepTout)
#    if Nin != Nout
#        error("Different dimension")
#    end
#    Tin  = eltype(RepTin)
#    Tout = eltype(RepTout)
#    if Tin == Tout
#        f = nothing
#    else
#        f = (i,x) -> similar_type(typeof(x), Tout)(x)
#    end
#    if decomposedvfast(p)
#        RepTout(PointIterator{Nout,Tout,Nin,Tin}([p], f), RayIterator{Nout,Tout,Nin,Tin}([p], f))
#    else
#        RepTout(VRepIterator{Nout,Tout,Nin,Tin}([p], f))
#    end
#end

Base.convert(::Type{RepTout}, p::VRepresentation) where {RepTout<:VRep} = vconvert(RepTout, p)
Base.convert(::Type{RepTout}, p::VRep) where {RepTout<:VRepresentation} = vconvert(RepTout, p)
# avoid ambiguity
Base.convert(::Type{RepTout}, p::VRepresentation) where {RepTout<:VRepresentation} = vconvert(RepTout, p)

MultivariatePolynomials.changecoefficienttype(p::Rep{N,T}, ::Type{T}) where {N,T} = p
function MultivariatePolynomials.changecoefficienttype(p::RepTin, ::Type{Tout}) where {RepTin<:Rep, Tout}
    RepTout = similar_type(RepTin, Tout)
    RepTout(p)
end

VRepresentation{N, T}(v::RepTin) where {N, T, RepTin} = similar_type(RepTin, FullDim{N}(), T)(v)
HRepresentation{N, T}(h::RepTin) where {N, T, RepTin} = similar_type(RepTin, FullDim{N}(), T)(h)

VRep{N, T}(v::VRepresentation) where {N, T} = VRepresentation{N, T}(v)
VRep{N, T}(p::Polyhedron) where {N, T} = Polyhedron{N, T}(p)
HRep{N, T}(h::HRepresentation) where {N, T} = HRepresentation{N, T}(h)
HRep{N, T}(p::Polyhedron) where {N, T} = Polyhedron{N, T}(p)

# FIXME it does not get called. The calls always go throug vconvert and hconvert. Use changecoefficienttype instead
#function Base.convert{N, T, RepT<:Representation}(::Type{Representation{N, T}}, rep::RepT)
#  if fulldim(RepT) != N
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficienttype(RepT, T), rep)
#end
#function Base.convert{N, T, RepT<:HRepresentation}(::Type{HRepresentation{N, T}}, rep::RepT)
#  if fulldim(RepT) != N
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficienttype(RepT, T), rep)
#end
#function Base.convert{M, S, RepT<:VRepresentation}(::Type{VRepresentation{M, S}}, rep::RepT)
#  if fulldim(RepT) != M
#    error("Cannot convert representations of the same dimension")
#  end
#  Base.convert(changecoefficienttype(RepT, S), rep)
#end

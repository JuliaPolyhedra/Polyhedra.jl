# Subtyping AbstractVector{ElemT} make Base think that RepIterator implements indexing e.g. for copy!
abstract type AbstractRepIterator{N, T, ElemT, PT} end

# For a RepIterator with only one representation, we fallback to the indexing interface, the index being the iterator state
# A representation can overwrite this if it can do something more efficient or if it simply does not support indexing
const SingleRepIterator{N, T, ElemT, RepT} = AbstractRepIterator{N, T, ElemT, Tuple{RepT}}
# With MapRepIterator, we could have RepT<:Rep{N', T'} with N' != N or T' != T, with eachindex, we need to replace everything with N', T'
Base.eachindex(it::SingleRepIterator{<:Any, <:Any, ElemT, RepT}) where {N, T, ElemT, RepT<:Rep{N, T}} = Indices{N, T, similar_type(ElemT, FullDim{N}(), T)}(it.ps[1])
Base.start(it::SingleRepIterator{<:Any, <:Any, ElemT, RepT})     where {N, T, ElemT, RepT<:Rep{N, T}} = start(eachindex(it))::Index{N, T, similar_type(ElemT, FullDim{N}(), T)}
Base.done(it::SingleRepIterator, idx::Index) = done(eachindex(it), idx)
Base.next(it::SingleRepIterator, idx::Index) = mapitem(it, 1, get(it.ps[1], idx)), nextindex(it.ps[1], idx)

# If there are multiple representations, we need to iterate.
# Builds a SingleRepIterator{ElemT} from p
function checknext(it::AbstractRepIterator, i::Int, state)
    while i <= length(it.ps) && (i == 0 || done(repit(it, i), state))
        i += 1
        if i <= length(it.ps)
            state = start(repit(it, i))
        end
    end
    i > length(it.ps) ? (i, nothing) : (i, state)
end
Base.start(it::AbstractRepIterator) = checknext(it, 0, nothing)
Base.done(it::AbstractRepIterator, state) = state[1] > length(it.ps)

function Base.next(it::AbstractRepIterator, state)
    i, substate = state
    item, newsubstate = next(repit(it, i), substate)
    newstate = checknext(it, i, newsubstate)
    item, newstate
end

# RepIterator
struct RepIterator{N, T, ElemT, PT<:Tuple{Vararg{Rep{N, T}}}} <: AbstractRepIterator{N, T, ElemT, PT}
    ps::PT
    function RepIterator{N, T, ElemT}(ps::PT) where {N, T, ElemT, PT<:Tuple{Vararg{Rep{N, T}}}}
        new{N, T, ElemT, PT}(ps)
    end
end
repit(it::RepIterator{N, T, ElemT}, i::Int) where {N, T, ElemT} = RepIterator{N, T, ElemT}((it.ps[i],))
mapitem(it::RepIterator, i, item) = item

# MapRepIterator
struct MapRepIterator{N, T, ElemT, PT} <: AbstractRepIterator{N, T, ElemT, PT}
    ps::PT
    f::Function
    function MapRepIterator{N, T, ElemT}(ps::PT, f::Function) where {N, T, ElemT, PT<:Tuple}
        new{N, T, ElemT, PT}(ps, f)
    end
end
# N, T are the output dim and coeftype, it can be different that the ones of is.ps[i]
repit(it::MapRepIterator{N, T, ElemT}, i::Int) where {N, T, ElemT} = MapRepIterator{N, T, ElemT}((it.ps[i],), (j, item) -> it.f(i, item)) # j should be 1
mapitem(it::MapRepIterator, i, item) = it.f(i, item)

function Base.first(it::AbstractRepIterator)
    next(it, start(it))[1]
end
FullDim(it::AbstractRepIterator{N}) where N = FullDim{N}()
Base.eltype(it::AbstractRepIterator{N, T, ElemT}) where {N, T, ElemT} = ElemT
function typed_map(f, d::FullDim{N}, ::Type{T}, it::RepIterator{Nin, Tin, ElemT}) where {N, T, Nin, Tin, ElemT}
    MapRepIterator{N, T, similar_type(ElemT, d, T)}(it.ps, f)
end

function RepIterator{N, T}(it::RepIterator) where {N, T}
    typed_map((i,x) -> similar_type(typeof(x), T)(x), FullDim{N}(), T, it)
end

# FIXME the variables need to be defined outside of the local scope of for
#       for Julia to know them inside the """ ... """ of the docstrings
HorV = HorVRep = horvrep = singular = singularlin = plural = plurallin = lenp = isnotemptyp = repexem = listexem = :_

for (isVrep, elt, singular) in [(true, :SymPoint, :sympoint), (true, :AbstractPoint, :point),
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
    typename = :(AbstractRepIterator{N, T, <:$elt})
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
        export $plural, $lenp, $isnotemptyp, $startp, $donep, $nextp, $elemtype

        """
            $plural($horvrep::$HorVRep)

        Returns an iterator over the $plural of the $HorV-representation `$horvrep`.
        """
        function $plural end

        """
            $lenp($horvrep::$HorVRep)

        Returns the number of $plural of the $HorV-representation `$horvrep`.
        """
        $lenp($horvrep::$HorVRep{N, T}) where {N, T} = length(Indices{N, T, $elemtype($horvrep)}($horvrep))

        """
            $isnotemptyp($horvrep::$HorVRep)

        Returns whether the $HorV-representation `$horvrep` has any $singular.
        """
        $isnotemptyp($horvrep::$HorVRep{N, T}) where {N, T} = !isempty(Indices{N, T, $elemtype($horvrep)}($horvrep))

        $elemtype(p::Polyhedron) = $elemtype($horvrep(p))
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
            MapRepIterator{N, T, ElemT}(p, f)
        end

        Base.length(it::$typename) where {N, T} = sum($lenp, it.ps)
        Base.isempty(it::$typename) where {N, T} = !any($isnotemptyp.(it.ps))
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

const ElemIt{ElemT} = Union{AllRepIterator{<:Any, <:Any, ElemT}, AbstractRepIterator{<:Any, <:Any, ElemT}, AbstractVector{ElemT}}
const HyperPlaneIt{N, T} = ElemIt{<:HyperPlane{N, T}}
const HalfSpaceIt{N, T} = ElemIt{<:HalfSpace{N, T}}
const HIt{N, T} = Union{HyperPlaneIt{N, T}, HalfSpaceIt{N, T}}

const SymPointIt{N, T} = ElemIt{<:SymPoint{N, T}}
const PointIt{N, T} = ElemIt{<:AbstractPoint{N, T}}
const PIt{N, T} = Union{SymPointIt{N, T}, PointIt{N, T}}
const LineIt{N, T} = ElemIt{<:Line{N, T}}
const RayIt{N, T} = ElemIt{<:Ray{N, T}}
const RIt{N, T} = Union{LineIt{N, T}, RayIt{N, T}}
const VIt{N, T} = Union{PIt{N, T}, RIt{N, T}}

function fillvits(sympoints::ElemIt{SymPoint{N, T, AT}}, points::ElemIt{AT}=AT[], lines::ElemIt{Line{N, T, AT}}=Line{N, T, AT}[], rays::ElemIt{Ray{N, T, AT}}=Ray{N, T, AT}[]) where {N, T, AT}
    if isempty(sympoints) && isempty(points) && !(isempty(lines) && isempty(rays))
        vconsistencyerror()
    end
    sympoints, points, lines, rays
end
function fillvits(lines::ElemIt{Line{N, T, AT}}, rays::ElemIt{Ray{N, T, AT}}=Ray{N, T, AT}[]) where {N, T, AT}
    sps = SymPoint{N, T, AT}[]
    if isempty(lines) && isempty(rays)
        ps = AT[]
    else
        ps = [origin(AT, FullDim{N}())]
    end
    sps, ps, lines, rays
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

vreps(p::VPolytope...) = preps(p...)
vreps(p::VCone...) = rreps(p...)
vreps(p...) = preps(p...)..., rreps(p...)...
preps(p::VSymPolytope...) = tuple(sympoints(p...))
preps(p::VRep...) = sympoints(p...), points(p...)
rreps(p::VPolytope...) = tuple()
rreps(p::VRep...) = lines(p...), rays(p...)

function vmap(f, d::FullDim, ::Type{T}, p::VPolytope...) where T
    mapsympoints(f, d, T, p...), mappoints(f, d, T, p...)
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

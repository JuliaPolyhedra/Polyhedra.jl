# Subtyping AbstractVector{ElemT} make Base think that RepIterator implements indexing e.g. for copy!
abstract type AbstractRepIterator{T, ElemT, PT} end

# For a RepIterator with only one representation, we fallback to the indexing interface, the index being the iterator state
# A representation can overwrite this if it can do something more efficient or if it simply does not support indexing
const SingleRepIterator{T, ElemT, RepT} = AbstractRepIterator{T, ElemT, Tuple{RepT}}
# With MapRepIterator, we could have fulldim(RepT) != N or T' != T, with eachindex, we need to replace everything with fulldim(RepT), T'
Base.eachindex(it::SingleRepIterator{<:Any, ElemT, RepT}) where {T, ElemT, RepT<:Rep{T}} = Indices{T, similar_type(ElemT, FullDim(RepT), T)}(it.ps[1])
Base.start(it::SingleRepIterator{<:Any, ElemT, RepT})     where {T, ElemT, RepT<:Rep{T}} = start(eachindex(it))::Index{T, similar_type(ElemT, FullDim(RepT), T)}
Base.done(it::SingleRepIterator, idx::Index) = done(eachindex(it), idx)
function Base.next(it::SingleRepIterator, idx::Index)
    #mapitem(it, 1, get(it.ps[1], idx)), nextindex(it.ps[1], idx)
    x = get(it.ps[1], idx)
    a = mapitem(it, 1, x)
    b = nextindex(it.ps[1], idx)
    a, b
end

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
struct RepIterator{T, ElemT, PT<:Tuple{Vararg{Rep{T}}}} <: AbstractRepIterator{T, ElemT, PT}
    ps::PT
    function RepIterator{T, ElemT}(ps::PT) where {T, ElemT, PT<:Tuple{Vararg{Rep{T}}}}
        new{T, ElemT, PT}(ps)
    end
end
repit(it::RepIterator{T, ElemT}, i::Int) where {T, ElemT} = RepIterator{T, ElemT}((it.ps[i],))
mapitem(it::RepIterator, i, item) = item

# MapRepIterator
struct MapRepIterator{T, ElemT, PT} <: AbstractRepIterator{T, ElemT, PT}
    ps::PT
    f::Function
    function MapRepIterator{T, ElemT}(ps::PT, f::Function) where {T, ElemT, PT<:Tuple}
        new{T, ElemT, PT}(ps, f)
    end
end
# T is the coeftype, it can be different that the ones of is.ps[i]
repit(it::MapRepIterator{T, ElemT}, i::Int) where {T, ElemT} = MapRepIterator{T, ElemT}((it.ps[i],), (j, item) -> it.f(i, item)) # j should be 1
mapitem(it::MapRepIterator, i, item) = it.f(i, item)

function Base.first(it::AbstractRepIterator)
    next(it, start(it))[1]
end
FullDim(it::AbstractRepIterator{T, ElemT}) where {T, ElemT} = FullDim(ElemT)
Base.eltype(it::AbstractRepIterator{T, ElemT}) where {T, ElemT} = ElemT
function typed_map(f, d::FullDim, ::Type{T}, it::RepIterator{Tin, ElemT}) where {T, Tin, ElemT}
    MapRepIterator{T, similar_type(ElemT, d, T)}(it.ps, f)
end

function RepIterator{T}(it::RepIterator) where {T}
    typed_map((i,x) -> similar_type(typeof(x), T)(x), FullDim(it), T, it)
end

# FIXME the variables need to be defined outside of the local scope of for
#       for Julia to know them inside the """ ... """ of the docstrings
singular = HorV = HorVRep = horvrep = singularlin = plural = plurallin = lenp = isnotemptyp = repexem = listexem = :_

for (isVrep, elt, loop_singular) in [(true, :AbstractVector, :point),
                                     (true, :Line, :line), (true, :Ray, :ray),
                                     (false, :HyperPlane, :hyperplane), (false, :HalfSpace, :halfspace)]
    global singular = loop_singular
    if isVrep
        vectortype = :vvectortype
        global HorV = :V
        global HorVRep = :VRep
        global horvrep = :vrep
    else
        vectortype = :hvectortype
        global HorV = :H
        global HorVRep = :HRep
        global horvrep = :hrep
    end
    typename = :(AbstractRepIterator{T, <:$elt})
    singularstr = string(singular)
    elemtype = Symbol(singularstr * "type")
    donep = Symbol("done" * singularstr)
    startp = Symbol("start" * singularstr)
    nextp = Symbol("next" * singularstr)
    pluralstr = singularstr * "s"
    global plural = Symbol(pluralstr)
    global lenp = Symbol("n" * pluralstr)
    global isnotemptyp = Symbol("has" * pluralstr)
    mapit = Symbol("map" * pluralstr)
    inc = Symbol("incident" * pluralstr)
    incidx = Symbol("incident" * singularstr * "indices")

    @eval begin
        export $plural, $lenp, $isnotemptyp, $startp, $donep, $nextp, $elemtype
        export $inc, $incidx

        """
            $plural($horvrep::$HorVRep)

        Returns an iterator over the $plural of the $HorV-representation `$horvrep`.
        """
        function $plural end

        """
            incident$plural(p::Polyhedron, idx)

        Returns the list of $plural incident to idx for the polyhedron `p`.
        """
        $inc(p::Polyhedron{T}, idx) where {T} = get(p, IncidentElements{T, $elemtype(p)}(p, idx))

        """
            incident$(singular)indices(p::Polyhedron, idx)

        Returns the list of the indices of $plural incident to idx for the polyhedron `p`.
        """
        $incidx(p::Polyhedron{T}, idx) where {T} = get(p, IncidentIndices{T, $elemtype(p)}(p, idx))

        """
            $lenp($horvrep::$HorVRep)

        Returns the number of $plural of the $HorV-representation `$horvrep`.
        """
        $lenp($horvrep::$HorVRep{T}) where {T} = length(Indices{T, $elemtype($horvrep)}($horvrep))

        """
            $isnotemptyp($horvrep::$HorVRep)

        Returns whether the $HorV-representation `$horvrep` has any $singular.
        """
        $isnotemptyp($horvrep::$HorVRep{T}) where {T} = !isempty(Indices{T, $elemtype($horvrep)}($horvrep))

        $elemtype(p::Polyhedron) = $elemtype($horvrep(p))
        if $singularstr == "point"
            $elemtype(p::$HorVRep) = $vectortype(typeof(p))
        else
            $elemtype(p::$HorVRep{T}) where {T} = $elt{T, $vectortype(typeof(p))}
        end

        function $plural(p::$HorVRep{T}...) where {T}
            ElemT = promote_type($elemtype.(p)...)
            RepIterator{T, ElemT}(p)
        end

        function $mapit(f::Function, d::FullDim, ::Type{T}, p::$HorVRep...) where {T}
            ElemT = promote_type(similar_type.($elemtype.(p), d, T)...)
            MapRepIterator{T, ElemT}(p, f)
        end

        Base.length(it::$typename) where {T} = sum($lenp, it.ps)
        Base.isempty(it::$typename) where {T} = !any($isnotemptyp.(it.ps))
    end
end

# Combines an element type with its linear version.
# e.g. combines rays with the lines by splitting lines in two points.
struct AllRepIterator{T, ElemT, LinElemT, PT}
    itlin::RepIterator{T, LinElemT, PT}
    it::RepIterator{T, ElemT, PT}
end

Base.eltype(it::AllRepIterator{T, ElemT}) where {T, ElemT} = ElemT
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
#splitlin(s::SymPoint, i) = (i == 1 ? coord(s) : -coord(s))
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

for (isVrep, loop_singularlin,
     loop_singular, loop_repexem,
     loop_listexem) in [#(true, :sympoint, :point, "convexhull(SymPoint([1, 0]), [0, 1])", "[1, 0], [-1, 0], [0, 1]"),
                   (true, :line, :ray, "Line([1, 0]) + Ray([0, 1])",
                    "Ray([1, 0]), Ray([-1, 0]), Ray([0, 1])"),
                   (false, :hyperplane, :halfspace,
                    "HyperPlane([1, 0], 1) ∩ HalfSpace([0, 1], 1)",
                    "HalfSpace([1, 0]), HalfSpace([-1, 0]), HalfSpace([0, 1])")]
    global singularlin = loop_singularlin
    global singular = loop_singular
    global repexem = loop_repexem
    global listexem = loop_listexem
    if isVrep
        global HorV = :V
        global HorVRep = :VRep
        global horvrep = :vrep
    else
        global HorV = :H
        global HorVRep = :HRep
        global horvrep = :hrep
    end
    pluralstrlin = string(singularlin) * "s"
    global plurallin = Symbol(pluralstrlin)
    lenplin = Symbol("n" * pluralstrlin)
    isnotemptyplin = Symbol("has" * pluralstrlin)
    pluralstr = string(singular) * "s"
    global plural = Symbol(pluralstr)
    global lenp = Symbol("n" * pluralstr)
    global isnotemptyp = Symbol("has" * pluralstr)
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

const ElemIt{ElemT} = Union{AllRepIterator{<:Any, ElemT}, AbstractRepIterator{<:Any, ElemT}, AbstractVector{ElemT}}
const HyperPlaneIt{T} = ElemIt{<:HyperPlane{T}}
const HalfSpaceIt{T} = ElemIt{<:HalfSpace{T}}
const HIt{T} = Union{HyperPlaneIt{T}, HalfSpaceIt{T}}

const PointIt{T} = ElemIt{<:AbstractVector{T}}
const PIt{T} = PointIt{T}
const LineIt{T} = ElemIt{<:Line{T}}
const RayIt{T} = ElemIt{<:Ray{T}}
const RIt{T} = Union{LineIt{T}, RayIt{T}}
const VIt{T} = Union{PIt{T}, RIt{T}}

const It{T} = Union{HIt{T}, VIt{T}}

function fillvits(points::ElemIt{AT}, lines::ElemIt{Line{T, AT}}=Line{T, AT}[], rays::ElemIt{Ray{T, AT}}=Ray{T, AT}[]) where {T, AT<:AbstractVector{T}}
    if isempty(points)
        if isempty(lines) && isempty(rays)
            N = 0
        else
            vconsistencyerror()
        end
    else
        N = fulldim(first(points))
    end
    return N, points, lines, rays
end
function fillvits(lines::ElemIt{Line{T, AT}}, rays::ElemIt{Ray{T, AT}}=Ray{T, AT}[]) where {T, AT}
    if isempty(lines) && isempty(rays)
        N = 0
        points = AT[]
    else
        if isempty(lines)
            N = fulldim(first(rays))
        else
            N = fulldim(first(lines))
        end
        points = [origin(AT, N)]
    end
    return N, points, lines, rays
end

hreps(p::HRep{T}...) where {T} = hyperplanes(p...), halfspaces(p...)
hreps(p::HAffineSpace{T}...) where {T} = tuple(hyperplanes(p...))

hmap(f, d::FullDim, ::Type{T}, p::HRep...) where T = maphyperplanes(f, d, T, p...), maphalfspaces(f, d, T, p...)
hmap(f, d::FullDim, ::Type{T}, p::HAffineSpace...) where T = tuple(maphyperplanes(f, d, T, p...))

hconvert(RepT::Type{<:HRep{T}}, p::HRep{T}) where {T} = constructpolyhedron(RepT, (p,), hreps(p)...)
hconvert(RepT::Type{<:HRep{T}}, p::HRep)    where {T} = constructpolyhedron(RepT, (p,), RepIterator{T}.(hreps(p))...)

vreps(p...) = preps(p...)..., rreps(p...)...
preps(p::VRep...) = tuple(points(p...))
preps(p::VCone...) = tuple()
rreps(p::VRep...) = lines(p...), rays(p...)
rreps(p::VLinearSpace...) = tuple(lines(p...))
rreps(p::VPolytope...) = tuple()

vmap(f, d::FullDim, ::Type{T}, p::VRep...) where T = pmap(f, d, T, p...)..., rmap(f, d, T, p...)...
pmap(f, d::FullDim, ::Type{T}, p::VRep...) where T = tuple(mappoints(f, d, T, p...))
pmap(f, d::FullDim, ::Type, p::VCone...) = tuple()
rmap(f, d::FullDim, ::Type{T}, p::VRep...) where T = maplines(f, d, T, p...), maprays(f, d, T, p...)
rmap(f, d::FullDim, ::Type{T}, p::VLinearSpace...) where T = tuple(maplines(f, d, T, p...))
rmap(f, d::FullDim, ::Type, p::VPolytope...) = tuple()

vconvert(RepT::Type{<:VRep{T}}, p::VRep{T}) where {T} = constructpolyhedron(RepT, (p,), vreps(p)...)
vconvert(RepT::Type{<:VRep{T}}, p::VRep)    where {T} = constructpolyhedron(RepT, (p,), RepIterator{T}.(vreps(p))...)

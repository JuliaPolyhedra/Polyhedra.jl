abstract type AbstractRepIterator{T, ElemT} end

# Subtyping AbstractVector{ElemT} make Base think that RepIterator implements indexing e.g. for copy!
abstract type AbstractSingleRepIterator{T, ElemT, PT} <: AbstractRepIterator{T, ElemT} end

# SingleRepIterator
struct SingleRepIterator{T, ElemT, PT<:Rep{T}} <: AbstractSingleRepIterator{T, ElemT, PT}
    p::PT
    function SingleRepIterator{T, ElemT}(p::PT) where {T, ElemT, PT<:Rep{T}}
        new{T, ElemT, PT}(p)
    end
end
mapitem(it::SingleRepIterator, item) = item
iterator(T::Type, ElemT::Type, p::Rep) = SingleRepIterator{T, ElemT}(p)

# SingleMapRepIterator
struct SingleMapRepIterator{T, ElemT, PT<:Rep} <: AbstractSingleRepIterator{T, ElemT, PT}
    p::PT
    f::Function
    function SingleMapRepIterator{T, ElemT}(p::PT, f::Function) where {T, ElemT, PT<:Rep}
        new{T, ElemT, PT}(p, f)
    end
end
mapitem(it::SingleMapRepIterator, item) = it.f(1, item)
iterator(T::Type, ElemT::Type, f::Function, p::Rep) = SingleMapRepIterator{T, ElemT}(p, f)

# For a RepIterator with only one representation, we fallback to the indexing interface, the index being the iterator state
# A representation can overwrite this if it can do something more efficient or if it simply does not support indexing
function Base.eachindex(it::AbstractSingleRepIterator{<:Any, ElemT,
                                              RepT}) where {T, ElemT, RepT<:Rep{T}}
    return Indices{T, similar_type(ElemT, FullDim(RepT), T)}(it.p)
end
element_and_index(it, idx::Nothing) = nothing
function element_and_index(it::AbstractSingleRepIterator, idx::Index)
    return mapitem(it, get(it.p, idx)), idx
end
function Base.iterate(it::AbstractSingleRepIterator)
    idx = undouble_it(iterate(eachindex(it)))
    return element_and_index(it, idx)
end
function Base.iterate(it::AbstractSingleRepIterator, idx::Index)
    idx = undouble_it(iterate(eachindex(it), idx))::Union{Nothing, typeof(idx)}
    return element_and_index(it, idx)
end

abstract type AbstractMultiRepIterator{T, ElemT, PT} <: AbstractRepIterator{T, ElemT} end

# If there are multiple representations, we need to iterate.
# Builds a SingleRepIterator{ElemT} from p
function checknext(it::AbstractMultiRepIterator, i::Int, item_state)
    while i <= length(it.ps) && item_state === nothing
        i += 1
        if i <= length(it.ps)
            item_state = iterate(repit(it, i))
        end
    end
    if item_state === nothing
        return nothing
    else
        return item_state[1], (i, item_state[2])
    end
end
Base.iterate(it::AbstractMultiRepIterator) = checknext(it, 0, nothing)
function Base.iterate(it::AbstractMultiRepIterator, state)
    return checknext(it, state[1], iterate(repit(it, state[1]), state[2]))
end

# RepIterator
struct RepIterator{T, ElemT, PT<:Tuple{Vararg{Rep{T}}}} <: AbstractMultiRepIterator{T, ElemT, PT}
    ps::PT
    function RepIterator{T, ElemT}(ps::PT) where {T, ElemT, PT<:Tuple{Vararg{Rep{T}}}}
        new{T, ElemT, PT}(ps)
    end
end
repit(it::RepIterator{T, ElemT}, i::Int) where {T, ElemT} = SingleRepIterator{T, ElemT}(it.ps[i])

iterator(T::Type, ElemT::Type, p::Rep...) = RepIterator{T, ElemT}(p)

# MapRepIterator
struct MapRepIterator{T, ElemT, PT} <: AbstractMultiRepIterator{T, ElemT, PT}
    ps::PT
    f::Function
    function MapRepIterator{T, ElemT}(ps::PT, f::Function) where {T, ElemT, PT<:Tuple}
        new{T, ElemT, PT}(ps, f)
    end
end
# T is the coeftype, it can be different that the ones of is.ps[i]
repit(it::MapRepIterator{T, ElemT}, i::Int) where {T, ElemT} = SingleMapRepIterator{T, ElemT}(it.ps[i], (j, item) -> it.f(i, item)) # j should be 1

iterator(T::Type, ElemT::Type, f::Function, p::Rep...) = MapRepIterator{T, ElemT}(p, f)

function typed_map(f::Function, ::Type{T}, it::SingleRepIterator{Tin, ElemT}) where {T, Tin, ElemT}
    SingleMapRepIterator{T, similar_type(ElemT, T)}(it.p, f)
end
function typed_map(f::Function, ::Type{T}, it::RepIterator{Tin, ElemT}) where {T, Tin, ElemT}
    MapRepIterator{T, similar_type(ElemT, T)}(it.ps, f)
end

Base.eltype(it::AbstractRepIterator{T, ElemT}) where {T, ElemT} = ElemT
function change_coefficient_type(it::AbstractRepIterator, T::Type)
    typed_map((i, x) -> convert(similar_type(typeof(x), T), x), T, it)
end

# FIXME the variables need to be defined outside of the local scope of for
#       for Julia to know them inside the """ ... """ of the docstrings
singular = HorV = HorVRep = horvrep = singularlin = plural = plurallin = lenp = isnotemptyp = repexem = listexem = :_

for (isVrep, elt, loop_singular) in [(true, :AbstractVector, :point),
                                     (true, :Line, :line), (true, :Ray, :ray),
                                     (false, :HyperPlane, :hyperplane), (false, :HalfSpace, :halfspace)]
    global singular = loop_singular
    if isVrep
        vectortype_fun = :vvectortype
        global HorV = :V
        global HorVRep = :VRep
        global horvrep = :vrep
    else
        vectortype_fun = :hvectortype
        global HorV = :H
        global HorVRep = :HRep
        global horvrep = :hrep
    end
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
            $elemtype(p::$HorVRep) = $vectortype_fun(typeof(p))
        else
            $elemtype(p::$HorVRep{T}) where {T} = $elt{T, $vectortype_fun(typeof(p))}
        end

        function $plural(p::$HorVRep{T}...) where {T}
            ElemT = promote_type($elemtype.(p)...)
            iterator(T, ElemT, p...)
        end

        function $mapit(f::Function, d::FullDim, ::Type{T}, p::$HorVRep...) where {T}
            ElemT = promote_type(similar_type.($elemtype.(p), Ref(d), T)...)
            iterator(T, ElemT, f, p...)
        end

        Base.length(it::AbstractSingleRepIterator{T, <:$elt}) where {T} = $lenp(it.p)
        Base.length(it::AbstractMultiRepIterator{T, <:$elt}) where {T} = sum($lenp, it.ps)
        Base.isempty(it::AbstractSingleRepIterator{T, <:$elt}) where {T} = !any($isnotemptyp(it.p))
        Base.isempty(it::AbstractMultiRepIterator{T, <:$elt}) where {T}  = !any($isnotemptyp.(it.ps))
    end
end

# Combines an element type with its linear version.
# e.g. combines rays with the lines by splitting lines in two points.
struct AllRepIterator{T, ElemT, LinElemT, LRT<:Union{AbstractRepIterator{T, LinElemT}}, RT<:Union{AbstractRepIterator{T, ElemT}}}
    itlin::LRT
    it::RT
end

Base.eltype(it::AllRepIterator{T, ElemT}) where {T, ElemT} = ElemT
Base.length(it::AllRepIterator) = 2length(it.itlin) + length(it.it)
Base.isempty(it::AllRepIterator) = isempty(it.itlin) && isempty(it.it)

function checknext(it::AllRepIterator, i, item_state)
    while i <= 3 && item_state === nothing
        i += 1
        if i <= 2
            @assert i >= 1
            item_state = iterate(it.itlin)
        elseif i == 3
            item_state = iterate(it.it)
        end
    end
    if item_state === nothing
        return nothing
    else
        if i <= 2
            @assert i >= 1
            item = splitlin(item_state[1], i)
        else
            @assert i == 3
            item = item_state[1]
        end
        return item, (i, item_state[2])
    end
end

splitlin(h::HyperPlane, i) = (i == 1 ? HalfSpace(h.a, h.β) : HalfSpace(-h.a, -h.β))
#splitlin(s::SymPoint, i) = (i == 1 ? coord(s) : -coord(s))
splitlin(l::Line, i) = (i == 1 ? Ray(coord(l)) : Ray(-coord(l)))

Base.iterate(it::AllRepIterator) = checknext(it, 0, nothing)
function Base.iterate(it::AllRepIterator, istate)
    i, state = istate
    if i <= 2
        @assert i >= 1
        item_state = iterate(it.itlin, state)
    else
        @assert i == 3
        item_state = iterate(it.it, state)
    end
    return checknext(it, i, item_state)
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

function fillvits(d::FullDim, points::ElemIt{AT},
                  lines::ElemIt{Line{T, AT}}=Line{T, AT}[],
                  rays::ElemIt{Ray{T, AT}}=Ray{T, AT}[]) where {T, AT<:AbstractVector{T}}
    if isempty(points) && !(isempty(lines) && isempty(rays))
        vconsistencyerror()
    end
    return points, lines, rays
end
function fillvits(d::FullDim, lines::ElemIt{Line{T, AT}},
                  rays::ElemIt{Ray{T, AT}}=Ray{T, AT}[]) where {T, AT}
    d = FullDim_rec(lines, rays)
    N = fulldim(d)
    if isempty(lines) && isempty(rays)
        points = AT[]
    else
        points = [origin(AT, N)]
    end
    return points, lines, rays
end

FullDim_hreps(p...) = FullDim(p[1]), hreps(p...)...
FullDim_vreps(p...) = FullDim(p[1]), vreps(p...)...

hreps(p::HRep{T}...) where {T} = hyperplanes(p...), halfspaces(p...)
hreps(p::HAffineSpace{T}...) where {T} = tuple(hyperplanes(p...))

hmap(f, d::FullDim, ::Type{T}, p::HRep...) where T = maphyperplanes(f, d, T, p...), maphalfspaces(f, d, T, p...)
hmap(f, d::FullDim, ::Type{T}, p::HAffineSpace...) where T = tuple(maphyperplanes(f, d, T, p...))

hconvert(RepT::Type{<:HRep{T}}, p::HRep{T}) where {T} = constructpolyhedron(RepT, FullDim(p), (p,), hreps(p)...)
hconvert(RepT::Type{<:HRep{T}}, p::HRep)    where {T} = constructpolyhedron(RepT, FullDim(p), (p,), change_coefficient_type.(hreps(p), T)...)

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

vconvert(RepT::Type{<:VRep{T}}, p::VRep{T}) where {T} = constructpolyhedron(RepT, FullDim(p), (p,), vreps(p)...)
function vconvert(RepT::Type{<:VRep{T}}, p::VRep)    where {T}
    constructpolyhedron(RepT, FullDim(p), (p,), change_coefficient_type.(vreps(p), T)...)
end

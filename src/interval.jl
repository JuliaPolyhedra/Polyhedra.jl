export IntervalLibrary, Interval

"""
    IntervalLibrary{T}

Default library for polyhedra of dimension 1. Many aspect of polyhedral computation become trivial in one dimension. This library exploits this fact.
The library is also used as a fallback for libraries that do not support 1-dimensional polyhedra (e.g. qhull). That is projecting a polyhedron using such library produces a polyhedron using `IntervalLibrary`.
"""
struct IntervalLibrary{T} <: Library
end

similar_library(lib::IntervalLibrary, d::FullDim, ::Type{T}) where T = default_library(d, T) # default_library allows to fallback to DefaultLibrary if d is not FullDim(1)

mutable struct Interval{T, AT, D} <: Polyhedron{T}
    hrep::Intersection{T, AT, D}
    vrep::Hull{T, AT, D}
    length::T
end

Interval{T, AT, D}(d::FullDim, it::HIt...) where {T, AT, D} = _hinterval(Intersection{T, AT, D}(d, it...), AT, D)
Interval{T, AT, D}(d::FullDim, it::VIt...) where {T, AT, D} = _vinterval(Hull{T, AT, D}(d, it...), AT, D)

# If AT is an SVector, it will be StaticArrays.Size{(1,)}
# otherwise, it will be 1
FullDim(p::Interval) = FullDim(p.hrep)

library(::Union{Interval{T}, Type{<:Interval{T}}}) where T = IntervalLibrary{T}()

hvectortype(::Type{<:Interval{T, AT}}) where {T, AT} = AT
vvectortype(::Type{<:Interval{T, AT}}) where {T, AT} = AT

similar_type(::Type{<:Interval}, d::FullDim, ::Type{T}) where T = default_type(d, T)

surface(::Interval{T}) where {T} = zero(T)
volume(p::Interval) = p.length
Base.isempty(p::Interval) = isempty(p.vrep)

function Interval{T, AT, D}(haslb::Bool, lb::T, hasub::Bool, ub::T, isempty::Bool) where {T, AT, D}
    if haslb && hasub && _gt(lb, ub)
        isempty = true
    end
    hps = HyperPlane{T, AT}[]
    hss = HalfSpace{T, AT}[]
    ps = AT[]
    ls = Line{T, AT}[]
    rs = Ray{T, AT}[]
    if !isempty
        if hasub
            push!(ps, StaticArrays.SVector(ub))
            if haslb && _isapprox(lb, ub)
                push!(hps, HyperPlane(StaticArrays.SVector(one(T)), ub))
            else
                push!(hss, HalfSpace(StaticArrays.SVector(one(T)), ub))
            end
            if !haslb
                push!(rs, Ray(StaticArrays.SVector(-one(T))))
            end
        else
            if haslb
                push!(rs, Ray(StaticArrays.SVector(one(T))))
            else
                push!(ps, origin(AT, 1))
                push!(ls, Line(StaticArrays.SVector(one(T))))
            end
        end
        if haslb
            if !_isapprox(lb, ub)
                push!(hss, HalfSpace(StaticArrays.SVector(-one(T)), -lb))
                push!(ps, StaticArrays.SVector(lb))
            end
        end
    else
        # The dimension should be -1 so 1 - nhyperplanes == -1 so nhyperplanes == 2
        push!(hps, HyperPlane(StaticArrays.SVector(one(T)), zero(T)))
        push!(hps, HyperPlane(StaticArrays.SVector(zero(T)), one(T)))
    end
    h = hrep(hps, hss)
    v = vrep(ps, ls, rs)
    volume = isempty ? zero(T) : (haslb && hasub ? max(zero(T), ub - lb) : -one(T))
    Interval{T, AT, D}(h, v, volume)
end

function _hinterval(rep::HRep{T}, ::Type{AT}, D) where {T, AT}
    haslb = false
    lb = zero(T)
    hasub = false
    ub = zero(T)
    empty = false
    function _setlb(newlb)
        if !haslb
            haslb = true
            lb = T(newlb)
        else
            lb = T(max(lb, newlb))
        end
    end
    function _setub(newub)
        if !hasub
            hasub = true
            ub = T(newub)
        else
            ub = T(min(ub, newub))
        end
    end
    for hp in hyperplanes(rep)
        α = hp.a[1]
        if isapproxzero(α)
            if !isapproxzero(hp.β)
                empty = true
            end
        else
            _setlb(hp.β / α)
            _setub(hp.β / α)
        end
    end
    for hs in halfspaces(rep)
        α = hs.a[1]
        if isapproxzero(α)
            if hs.β < 0
                empty = true
            end
        elseif α < 0
            _setlb(hs.β / α)
        else
            _setub(hs.β / α)
        end
    end
    Interval{T, AT, D}(haslb, lb, hasub, ub, empty)
end

function _vinterval(v::VRep{T}, ::Type{AT}, D) where {T, AT}
    haslb = true
    lb = zero(T)
    hasub = true
    ub = zero(T)
    isempty = true
    for p in points(v)
        x = coord(p)[1]
        if isempty
            isempty = false
            lb = x
            ub = x
        else
            lb = min(lb, x)
            ub = max(ub, x)
        end
    end
    for l in lines(v)
        if !isapproxzero(l)
            isempty = false
            haslb = false
            hasub = false
        end
    end
    for r in rays(v)
        x = coord(r)[1]
        if !isapproxzero(x)
            isempty = false
            if x > 0
                hasub = false
            else
                haslb = false
            end
        end
    end
    Interval{T, AT, D}(haslb, lb, hasub, ub, isempty)
end

Interval{T, AT, D}(p::HRepresentation{T}) where {T, AT, D} = _hinterval(p, AT, D)
Interval{T, AT, D}(p::VRepresentation{T}) where {T, AT, D} = _vinterval(p, AT, D)
function Interval{T, AT, D}(p::Polyhedron{T}) where {T, AT, D}
    if hrepiscomputed(p)
        _hinterval(p, AT, D)
    else
        _vinterval(p, AT, D)
    end
end

function polyhedron(rep::Rep{T}, ::IntervalLibrary{T}) where T
    Interval{T, StaticArrays.SVector{1, T}, StaticArrays.Size{(1,)}}(rep)
end

hrep(p::Interval) = p.hrep
vrep(p::Interval) = p.vrep

hrepiscomputed(::Interval) = true
vrepiscomputed(::Interval) = true

# Nothing to do
function detecthlinearity!(::Interval) end
function detectvlinearity!(::Interval) end
function removehredundancy!(::Interval) end
function removevredundancy!(::Interval) end

export IntervalLibrary, Interval

struct IntervalLibrary{T} <: PolyhedraLibrary
end

mutable struct Interval{T, AT} <: Polyhedron{1, T}
    hrep::Intersection{1, T, AT}
    vrep::Hull{1, T, AT}
    length::T
end

arraytype(p::Interval{T, AT}) where {T, AT} = AT

surface{T}(::Interval{T}) = zero(T)
volume(p::Interval) = p.length
Base.isempty(p::Interval) = isempty(p.vrep)

function Interval{T, AT}(haslb::Bool, lb::T, hasub::Bool, ub::T, isempty::Bool) where {T, AT}
    if haslb && hasub && _gt(lb, ub)
        isempty = true
    end
    hps = HyperPlane{1, T, AT}[]
    hss = HalfSpace{1, T, AT}[]
    sps = SymPoint{1, T, AT}[]
    ps = AT[]
    ls = Line{1, T, AT}[]
    rs = Ray{1, T, AT}[]
    if !isempty
        if hasub
            if haslb && _isapprox(lb, -ub)
                push!(sps, SymPoint(SVector(ub)))
            else
                push!(ps, SVector(ub))
            end
            if haslb && _isapprox(lb, ub)
                push!(hps, HyperPlane(SVector(one(T)), ub))
            else
                push!(hss, HalfSpace(SVector(one(T)), ub))
            end
            if !haslb
                push!(rs, Ray(SVector(-one(T))))
            end
        else
            if haslb
                push!(rs, Ray(SVector(one(T))))
            else
                push!(ls, Line(SVector(one(T))))
            end
        end
        if haslb
            if !_isapprox(lb, ub)
                push!(hss, HalfSpace(SVector(-one(T)), -lb))
                if !_isapprox(lb, -ub)
                    push!(ps, SVector(lb))
                end
            end
        end
    else
        # The dimension should be -1 so 1 - nhyperplanes == -1 so nhyperplanes == 2
        push!(hps, HyperPlane(SVector(one(T)), zero(T)))
        push!(hps, HyperPlane(SVector(zero(T)), one(T)))
    end
    h = hrep(hps, hss)
    v = vrep(sps, ps, ls, rs)
    volume = haslb && hasub ? max(zero(T), ub - lb) : -one(T)
    Interval{T, AT}(h, v, volume)
end

function _hinterval(rep::HRep{1, T}, ::Type{AT}) where {T, AT}
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
                lb = T(Inf)
                ub = T(-Inf)
            end
        elseif α < 0
            _setlb(hs.β / α)
        else
            _setub(hs.β / α)
        end
    end
    Interval{T, AT}(haslb, lb, hasub, ub, empty)
end

function _vinterval(v::VRep{1, T}, ::Type{AT}) where {T, AT}
    haslb = true
    lb = zero(T)
    hasub = true
    ub = zero(T)
    isempty = true
    for p in allpoints(v)
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
    Interval{T, AT}(haslb, lb, hasub, ub, isempty)
end

Interval{T, AT}(p::HRepresentation{1, T}) where {T, AT} = _hinterval(p, AT)
Interval{T, AT}(p::VRepresentation{1, T}) where {T, AT} = _vinterval(p, AT)
function Interval{T, AT}(p::Polyhedron{1, T}) where {T, AT}
    if hrepiscomputed(p)
        _hinterval(p, AT)
    else
        _vinterval(p, AT)
    end
end

function polyhedron(rep::Rep{1, T}, ::IntervalLibrary{T}) where T
    Interval{T, SVector{1, T}}(rep)
end

hrep(p::Interval) = p.hrep
vrep(p::Interval) = p.vrep

hrepiscomputed(::Interval) = true
vrepiscomputed(::Interval) = true

# Nothing to do
function detecthlinearities!(::Interval) end
function detectvlinearities!(::Interval) end
function removehredundancy!(::Interval) end
function removevredundancy!(::Interval) end

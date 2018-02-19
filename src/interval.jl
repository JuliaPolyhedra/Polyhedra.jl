export IntervalLibrary, Interval

struct IntervalLibrary{T} <: PolyhedraLibrary
end

# TODO use vecrep instead of simplerep
mutable struct Interval{T} <: Polyhedron{1, T}
    hrep::SimpleHRepresentation{1, T}
    vrep::SimpleVRepresentation{1, T}
    length::T
end

arraytype(p::Interval) = arraytype(p.hrep)

area{T}(::Interval{T}) = zero(T)
volume(p::Interval) = p.length
Base.isempty(p::Interval) = isempty(p.vrep)

function Interval{T}(haslb::Bool, lb::T, hasub::Bool, ub::T, isempty::Bool) where T
    if haslb && hasub && mygt(lb, ub)
        isempty = true
    end
    A = Matrix{Float64}(0, 1)
    b = Float64[]
    linset = IntSet()
    V = Matrix{Float64}(0, 1)
    Vlinset = IntSet()
    R = Matrix{Float64}(0, 1)
    Rlinset = IntSet()
    if !isempty
        if hasub
            A = [A; 1]
            push!(b, ub)
            V = [V; ub]
        else
            R = [R; 1]
        end
        if haslb
            if myeq(lb, ub)
                push!(linset, 1)
            else
                A = [A; -1]
                push!(b, -lb)
                if myeq(lb, -ub)
                    push!(Vlinset, 1)
                else
                    V = [V; lb]
                end
            end
        else
            if hasub
                R = [R; -1]
            else
                push!(Rlinset, 1)
            end
        end
    else
        A = [A; 0]
        push!(b, -1)
        push!(linset, 1)
    end
    h = SimpleHRepresentation{1, Float64}(A, b, linset)
    v = SimpleVRepresentation{1, Float64}(V, R, Vlinset, Rlinset)
    volume = haslb && hasub ? max(zero(T), ub - lb) : -one(T)
    Interval{T}(h, v, volume)
end

function _hinterval(rep::HRep{1, T}) where T
    haslb = false
    lb = zero(T)
    hasub = false
    ub = zero(T)
    empty = false
    function _setlb(newlb)
        if !haslb
            haslb = true
            lb = newlb
        else
            lb = max(lb, newlb)
        end
    end
    function _setub(newub)
        if !hasub
            hasub = true
            ub = newub
        else
            ub = min(ub, newub)
        end
    end
    for hp in hyperplanes(rep)
        α = hp.a[1]
        if myeqzero(α)
            if !myeqzero(hp.β)
                empty = true
            end
        else
            _setlb(hp.β / α)
            _setub(hp.β / α)
        end
    end
    for hs in halfspaces(rep)
        α = hs.a[1]
        if myeqzero(α)
            if hs.β < 0
                lb = Inf
                ub = -Inf
            end
        elseif α < 0
            _setlb(hs.β / α)
        else
            _setub(hs.β / α)
        end
    end
    Interval{T}(haslb, lb, hasub, ub, empty)
end

function _vinterval(v::VRep{1, T}) where T
    haslb = true
    lb = zero(T)
    hasub = true
    ub = zero(T)
    isempty = true
    for r in rays(v)
        x = coord(r)[1]
        if !iszero(x)
            isempty = false
            if islin(r)
                haslb = false
                hasub = false
            else
                if x > 0
                    hasub = false
                else
                    haslb = false
                end
            end
        end
    end
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
        if islin(p)
          lb = min(lb, -x)
          ub = max(ub, -x)
        end
    end
    Interval{T}(haslb, lb, hasub, ub, isempty)
end

Interval{T}(p::HRepresentation{1, T}) where T = _hinterval(p)
Interval{T}(p::VRepresentation{1, T}) where T = _vinterval(p)
function Interval{T}(p::Polyhedron{1, T}) where T
    if hrepiscomputed(p)
        _hinterval(p)
    else
        _vinterval(p)
    end
end

function polyhedron{T}(rep::Rep{1, T}, ::IntervalLibrary{T})
    Interval{T}(rep)
end

hrep(p::Interval) = p.hrep
vrep(p::Interval) = p.vrep

hrepiscomputed(::Interval) = true
vrepiscomputed(::Interval) = true

function detecthlinearities!(::Interval) end
function detectvlinearities!(::Interval) end
function removehredundancy!(::Interval) end
function removevredundancy!(::Interval) end

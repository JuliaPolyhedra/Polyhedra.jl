using LinearAlgebra, Test
using Polyhedra

_isapprox(x::Real, y::Real) = _isapprox(promote(x, y)...)
_isapprox(x::T, y::T) where {T<:Real} = x == y
function _isapprox(x::T, y::T) where {T<:AbstractFloat}
    return y < x+1024*eps(T) && x < y+1024*eps(T)
end
function _isapprox(
    x::AbstractVector{S},
    y::AbstractVector{T},
) where {S<:Real,T<:Real}
    return _isapprox(promote(x, y)...)
end
_isapprox(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Real} = x == y
function _isapprox(
    x::AbstractVector{T},
    y::AbstractVector{T},
) where {T<:AbstractFloat}
    return _isapproxzero(norm(x - y))
end
_isapproxzero(x::T) where {T<:Real} = _isapprox(x, zero(T))

function _parallel(x, y)
    @assert length(x) == length(y)
    cstchosen = false
    cst = zero(Base.promote_op(/, eltype(y), eltype(x)))
    for i in 1:length(x)
        if _isapproxzero(x[i])
            _isapproxzero(y[i]) || return false
        elseif _isapproxzero(y[i])
            _isapproxzero(x[i]) || return false
        else
            c = y[i] / x[i]
            if cstchosen
                _isapprox(cst, c) || return false
            else
                cst = c
            end
        end
    end
    return true
end

function inaffspace(r, s, L, par::Bool = true)
    x = Polyhedra.coord(r)
    y = Polyhedra.coord(s)
    for l in L
        z = Polyhedra.coord(l)
        # remove component
        x = x * dot(z, z) - z * dot(z, x)
        y = y * dot(z, z) - z * dot(z, y)
    end
    return par ? _parallel(x, y) : _isapprox(r, s)
end

function inequality_fulltest(hr::HRepresentation, exp::HRepresentation)
    @test nhyperplanes(hr) == nhyperplanes(exp)
    @test nhalfspaces(hr) == nhalfspaces(exp)

    _lift(h::Polyhedra.HRepElement) = [h.β; h.a]
    aff = map(_lift, hyperplanes(hr))
    inaff(x, y) = inaffspace(x, y, aff)

    for h in hyperplanes(exp)
        hl = _lift(h)
        @test any(hp -> inaff(hl, _lift(hp)), hyperplanes(hr))
    end
    for h in halfspaces(exp)
        hl = _lift(h)
        @test any(hs -> inaff(hl, _lift(hs)), halfspaces(hr))
    end
end
function inequality_fulltest(h::HRepresentation, hrepargs...)
    return inequality_fulltest(h, hrep(hrepargs...))
end
function inequality_fulltest(p::Polyhedron, hrepargs...)
    detecthlinearity!(p)
    removehredundancy!(p)
    return inequality_fulltest(hrep(p), hrepargs...)
end

function generator_fulltest(vr::VRepresentation, exp::VRepresentation)
    @test npoints(vr) == npoints(exp)
    @test nlines(vr) == nlines(exp)
    @test nrays(vr) == nrays(exp)
    linspace = collect(lines(vr))
    inaff(x, y, args...) = inaffspace(x, y, linspace, args...)
    for l in lines(exp)
        @test any(ll -> inaff(l, ll), lines(vr))
    end
    for r in rays(exp)
        @test any(s -> inaff(r, s), rays(vr))
    end
    for p in points(exp)
        @test any(q -> inaff(p, q, false), points(vr))
    end
end
function generator_fulltest(v::VRepresentation, vrepargs...)
    return generator_fulltest(v, vrep(vrepargs...))
end
function generator_fulltest(p::Polyhedron, vrepargs...)
    detectvlinearity!(p)
    removevredundancy!(p)
    return generator_fulltest(vrep(p), vrepargs...)
end

function alloc_test(f, n)
    f() # compile
    @test n == @allocated f()
end

using Test
using Polyhedra

# Allocating size for allocating a `BigInt`.
# Half size on 32-bit.
const BIGINT_ALLOC = Sys.WORD_SIZE == 64 ? 48 : 24

function alloc_test(f, n)
    f() # compile
    @test n == @allocated f()
end

function rec_identical(vr::Polyhedra.LinesHull, vr2::Polyhedra.LinesHull)
    @test vr === vr2
    @test vr.lines === vr2.lines
end
function rec_identical(vr::Polyhedra.RaysHull, vr2::Polyhedra.RaysHull)
    @test vr === vr2
    @test vr2.rays === vr.rays
    rec_identical(vr.lines, vr2.lines)
end
function rec_identical(vr::Polyhedra.PointsHull, vr2::Polyhedra.PointsHull)
    @test vr === vr2
    @test vr2.points === vr.points
end

function rec_identical(vr::Polyhedra.Hull, vr2::Polyhedra.Hull)
    @test vr2 === vr
    rec_identical(vr.points, vr2.points)
    rec_identical(vr.rays, vr2.rays)
end

_isapprox(x::Real, y::Real) = _isapprox(promote(x, y)...)
_isapprox(x::T, y::T) where {T<:Real} = x == y
_isapprox(x::T, y::T) where {T<:AbstractFloat} = y < x+1024*eps(T) && x < y+1024*eps(T)
_isapprox(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real,T<:Real} = _isapprox(promote(x, y)...)
_isapprox(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Real} = x == y
_isapprox(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:AbstractFloat} = _isapproxzero(norm(x - y))
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

function inaffspace(r, s, L, par::Bool=true)
    x = Polyhedra.coord(r)
    y = Polyhedra.coord(s)
    for l in L
        z = Polyhedra.coord(l)
        # remove component
        x = x * dot(z, z) - z * dot(z, x)
        y = y * dot(z, z) - z * dot(z, y)
    end
    par ? _parallel(x, y) : _isapprox(r, s)
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
    inequality_fulltest(h, hrep(hrepargs...))
end
function inequality_fulltest(p::Polyhedron, hrepargs...)
    detecthlinearity!(p)
    removehredundancy!(p)
    inequality_fulltest(hrep(p), hrepargs...)
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
    generator_fulltest(v, vrep(vrepargs...))
end
function generator_fulltest(p::Polyhedron, vrepargs...)
    detectvlinearity!(p)
    removevredundancy!(p)
    generator_fulltest(vrep(p), vrepargs...)
end

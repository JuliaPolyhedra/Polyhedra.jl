isapproxzero(x::T; kws...) where {T<:Real} = x == zero(T)
isapproxzero(x::T; ztol=Base.rtoldefault(Float64)) where {T<:AbstractFloat} = abs(x) < ztol
isapproxzero(x::AbstractVector{T}; kws...) where {T<:Real} = isapproxzero(maximum(abs.(x)); kws...)
#isapproxzero(x::Union{SymPoint, Ray, Line}; kws...) = isapproxzero(coord(x); kws...)
isapproxzero(x::Union{Ray, Line}; kws...) = isapproxzero(coord(x); kws...)
isapproxzero(h::HRepElement; kws...) = isapproxzero(h.a; kws...) && isapproxzero(h.β; kws...)

_isapprox(x::Union{T, AbstractArray{T}}, y::Union{T, AbstractArray{T}}) where {T<:Union{Integer, Rational}} = x == y
# I check with zero because isapprox(0, 1e-100) is false...
# but isapprox(1e-100, 2e-100) should be false
_isapprox(x, y) = (isapproxzero(x) ? isapproxzero(y) : (isapproxzero(y) ? isapproxzero(x) : isapprox(x, y)))

#Base.isapprox(r::SymPoint, s::SymPoint) = _isapprox(coord(r), coord(s)) || _isapprox(coord(r), -coord(s))
function _scaleray(r::Union{Line, Ray}, s::Union{Line, Ray})
    cr = coord(r)
    cs = coord(s)
    cr * sum(abs.(cs)), cs * sum(abs.(cr))
end
function Base.isapprox(r::Line, s::Line)
    rs, ss = _scaleray(r, s)
    _isapprox(rs, ss) || _isapprox(rs, -ss)
end
function Base.isapprox(r::Ray, s::Ray)
    _isapprox(_scaleray(r, s)...)
end
# TODO check that Vec in GeometryTypes also does that
function _scalehp(h1, h2)
    s1 = sum(abs.(h1.a)) + abs(h1.β)
    s2 = sum(abs.(h2.a)) + abs(h2.β)
    (h1.a*s2, h1.β*s2), (h2.a*s1, h2.β*s1)
end
function Base.:(==)(h1::HyperPlane, h2::HyperPlane)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    (a1 == a2 && β1 == β2) || (a1 == -a2 && β1 == -β2)
end
function Base.isapprox(h1::HyperPlane, h2::HyperPlane)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    (_isapprox(a1, a2) && _isapprox(β1, β2)) || (_isapprox(a1, -a2) && _isapprox(β1, -β2))
end
function Base.:(==)(h1::HalfSpace, h2::HalfSpace)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    a1 == a2 && β1 == β2
end
function Base.isapprox(h1::HalfSpace, h2::HalfSpace)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    _isapprox(a1, a2) && _isapprox(β1, β2)
end

_lt(x::T, y::T) where {T<:Real} = x < y
_lt(x::T, y::T) where {T<:AbstractFloat} = x < y && !_isapprox(x, y)
_lt(x::S, y::T) where {S<:Real,T<:Real} = _lt(promote(x, y)...)
_gt(x::S, y::T) where {S<:Real, T<:Real} = _lt(y, x)
_leq(x::T, y::T) where {T<:Real} = x <= y
_leq(x::T, y::T) where {T<:AbstractFloat} = x <= y || _isapprox(x, y)
_leq(x::S, y::T) where {S<:Real,T<:Real} = _leq(promote(x, y)...)
_geq(x::T, y::T) where {T<:Real} = _leq(y, x)
_pos(x::T) where {T<:Real} = _gt(x, zero(T))
_neg(x::T) where {T<:Real} = _lt(x, zero(T))
_nonneg(x::T) where {T<:Real} = _geq(x, zero(T))
_nonpos(x::T) where {T<:Real} = _leq(x, zero(T))

function _promote_type(::Type{T1}, ::Type{T2}) where {T1,T2}
    if T1 <: AbstractFloat || T2 <: AbstractFloat
        Float64
    else
        Rational{BigInt}
    end
end

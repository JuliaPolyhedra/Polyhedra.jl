# We use `Zero` so that it is allocation-free, while `zero(Rational{BigInt})` would be costly
_default_tol(::Type{<:Union{Integer, Rational}}) = MA.Zero()
# `rtoldefault(BigFloat)` is quite slow and allocates a lot so we hope that
# `tol` will be called at some upper level function and passed in here instead
# of created everytime
_default_tol(::Type{T}) where {T<:AbstractFloat} = Base.rtoldefault(T)

isapproxzero(x::T; kws...) where {T<:Real} = x == zero(T)
isapproxzero(x::T; tol=Base.rtoldefault(T)) where {T<:AbstractFloat} = abs(x) < tol
isapproxzero(x::AbstractVector{T}; kws...) where {T<:Real} = isapproxzero(maximum(abs, x); kws...)
#isapproxzero(x::Union{SymPoint, Ray, Line}; kws...) = isapproxzero(coord(x); kws...)
isapproxzero(x::Union{Ray, Line}; kws...) = isapproxzero(coord(x); kws...)
isapproxzero(h::HRepElement; kws...) = isapproxzero(h.a; kws...) && isapproxzero(h.β; kws...)

_isapprox(x::Union{T, AbstractArray{T}}, y::Union{T, AbstractArray{T}}; tol) where {T<:Union{Integer, Rational}} = x == y
# I check with zero because isapprox(0, 1e-100) is false...
# but isapprox(1e-100, 2e-100) should be false
_isapprox(x, y; tol) = (isapproxzero(x; tol) ? isapproxzero(y; tol) : (isapproxzero(y; tol) ? isapproxzero(x; tol) : isapprox(x, y; rtol = tol)))

#Base.isapprox(r::SymPoint, s::SymPoint) = _isapprox(coord(r), coord(s)) || _isapprox(coord(r), -coord(s))
function _scaleray(r::Union{Line, Ray}, s::Union{Line, Ray})
    cr = coord(r)
    cs = coord(s)
    cr * sum(abs.(cs)), cs * sum(abs.(cr))
end
function _isapprox(r::Line, s::Line; tol)
    rs, ss = _scaleray(r, s)
    _isapprox(rs, ss; tol) || _isapprox(rs, -ss; tol)
end
function _isapprox(r::Ray, s::Ray; tol)
    _isapprox(_scaleray(r, s)...; tol)
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
function _isapprox(h1::HyperPlane, h2::HyperPlane; tol)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    (_isapprox(a1, a2; tol) && _isapprox(β1, β2; tol)) || (_isapprox(a1, -a2; tol) && _isapprox(β1, -β2; tol))
end
function Base.:(==)(h1::HalfSpace, h2::HalfSpace)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    a1 == a2 && β1 == β2
end
function _isapprox(h1::HalfSpace, h2::HalfSpace; tol)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    _isapprox(a1, a2; tol) && _isapprox(β1, β2; tol)
end

_lt(x::T, y::T; tol) where {T<:Real} = x < y
_lt(x::T, y::T; tol) where {T<:AbstractFloat} = x < y && !_isapprox(x, y; tol)
_lt(x::S, y::T; tol) where {S<:Real,T<:Real} = _lt(promote(x, y)...; tol)
_gt(x::S, y::T; tol) where {S<:Real, T<:Real} = _lt(y, x; tol)
_leq(x::T, y::T; tol) where {T<:Real} = x <= y
_leq(x::T, y::T; tol) where {T<:AbstractFloat} = x <= y || _isapprox(x, y; tol)
_leq(x::S, y::T; tol) where {S<:Real,T<:Real} = _leq(promote(x, y)...; tol)
_geq(x::T, y::T; tol) where {T<:Real} = _leq(y, x; tol)
_pos(x::T; tol) where {T<:Real} = _gt(x, zero(T); tol)
_neg(x::T; tol) where {T<:Real} = _lt(x, zero(T); tol)
_nonneg(x::T; tol) where {T<:Real} = _geq(x, zero(T); tol)
_nonpos(x::T; tol) where {T<:Real} = _leq(x, zero(T); tol)

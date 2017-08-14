const threshold = 1e-8

myeqzero(x::T) where {T<:Real} = x == zero(T)
myeqzero(x::T) where {T<:AbstractFloat} = -threshold < x < threshold
myeqzero(x::AbstractVector{T}) where {T<:Real} = myeqzero(sum(abs, x))
myeqzero(x::Union{SymPoint, Ray, Line}) = myeqzero(coord(x))
myeqzero(h::HRepElement) = myeqzero(h.a) && myeqzero(h.β)

myeq(x::Union{T, AbstractArray{T}}, y::Union{T, AbstractArray{T}}) where {T<:Union{Integer, Rational}} = x == y
# I check with zero because isapprox(0, 1e-100) is false...
# but isapprox(1e-100, 2e-100) should be false
myeq(x, y) = (myeqzero(x) ? myeqzero(y) : (myeqzero(y) ? myeqzero(x) : isapprox(x, y)))

Base.isapprox(r::SymPoint, s::SymPoint) = myeq(coord(r), coord(s)) || myeq(coord(r), -coord(s))
function _scaleray(r::Union{Line, Ray}, s::Union{Line, Ray})
    cr = coord(r)
    cs = coord(s)
    cr * sum(abs.(cs)), cs * sum(abs.(cr))
end
function Base.isapprox(r::Line, s::Line)
    rs, ss = _scaleray(r, s)
    myeq(rs, ss) || myeq(rs, -ss)
end
function Base.isapprox(r::Ray, s::Ray)
    myeq(_scaleray(r, s)...)
end
# TODO check that Vec in GeometryTypes also does that
function _scalehp(h1, h2)
    s1 = sum(abs.(h1.a)) + abs(h1.β)
    s2 = sum(abs.(h2.a)) + abs(h2.β)
    (h1.a*s2, h1.β*s2), (h2.a*s1, h2.β*s1)
end
function (==)(h1::HyperPlane, h2::HyperPlane)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    (a1 == a2 && β1 == β2) || (a1 == -a2 && β1 == -β2)
end
function Base.isapprox(h1::HyperPlane, h2::HyperPlane)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    (myeq(a1, a2) && myeq(β1, β2)) || (myeq(a1, -a2) && myeq(β1, -β2))
end
function (==)(h1::HalfSpace, h2::HalfSpace)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    a1 == a2 && β1 == β2
end
function Base.isapprox(h1::HalfSpace, h2::HalfSpace)
    (a1, β1), (a2, β2) = _scalehp(h1, h2)
    myeq(a1, a2) && myeq(β1, β2)
end

mylt(x::T, y::T) where {T<:Real} = x < y
mylt(x::T, y::T) where {T<:AbstractFloat} = x < y && !myeq(x, y)
mylt(x::S, y::T) where {S<:Real,T<:Real} = mylt(promote(x, y)...)
mygt(x::S, y::T) where {S<:Real, T<:Real} = mylt(y, x)
myleq(x::T, y::T) where {T<:Real} = x <= y
myleq(x::T, y::T) where {T<:AbstractFloat} = x <= y || myeq(x, y)
myleq(x::S, y::T) where {S<:Real,T<:Real} = myleq(promote(x, y)...)
mygeq(x::T, y::T) where {T<:Real} = myleq(y, x)
mypos(x::T) where {T<:Real} = mygt(x, zero(T))
myneg(x::T) where {T<:Real} = mylt(x, zero(T))
mynonneg(x::T) where {T<:Real} = mygeq(x, zero(T))
mynonpos(x::T) where {T<:Real} = myleq(x, zero(T))

function mypromote_type{T1,T2}(::Type{T1}, ::Type{T2})
    if T1 <: AbstractFloat || T2 <: AbstractFloat
        Float64
    else
        Rational{BigInt}
    end
end

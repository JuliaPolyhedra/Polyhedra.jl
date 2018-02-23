import Base: getindex, vec, dot, cross, -, +, /
import GeometryTypes.Point
export Point
export HRepElement, HalfSpace, HyperPlane
export VRepElement, AbstractPoint, SymPoint, AbstractRay, Ray, Line
export islin, isray, ispoint, ispoint, coord, lift, simplify

mypoint{T}(::Type{T}, a::AbstractVector) = AbstractArray{T}(a)
mypoint{T}(::Type{T}, a::AbstractVector{T}) = a
mypoint{N,T}(::Type{T}, a::Point{N}) = Point{N,T}(a)
mypoint{N,T}(::Type{T}, a::Point{N,T}) = a
const MyVec{N,T} = Union{Vec{N,T},AbstractVector{T}}
myvec{T}(::Type{T}, a::AbstractVector) = AbstractArray{T}(a)
myvec{T}(::Type{T}, a::AbstractVector{T}) = a
myvec{N,T}(::Type{T}, a::Vec{N}) = Vec{N,T}(a)
myvec{N,T}(::Type{T}, a::Vec{N,T}) = a

mydot(a, b) = sum(a .* b)
mydot(a::AbstractVector, b::AbstractVector) = dot(a, b)
mydot(f::FixedVector, v::FixedVector) = dot(f, v)
#mydot{T<:Union{Ray,Line}}(a::T, b::T) = mydot(a.r, b.r)
#mydot{T<:Union{Ray,Line}}(a, b::T) = mydot(a, b.r)
#mydot{T<:Union{Ray,Line}}(a::T, b) = mydot(a.r, b)

vecconv{T}(::Type{T}, a::AbstractVector) = AbstractVector{T}(a)
vecconv{T}(::Type{T}, a::FixedVector) = FixedVector{T}(a)

abstract type HRepElement{N,T} end

"""
    struct HalfSpace{N, T, AT} <: HRepElement{N, T}
        a::AT
        β::T
    end

An halfspace defined by the set of points ``x`` such that ``\\langle a, x \\rangle \\le \\beta``.
"""
struct HalfSpace{N, T, AT<:MyVec{N, T}} <: HRepElement{N,T}
    a::AT
    β::T
    function HalfSpace{N, T, AT}(a::AT, β::T) where {N, T, AT<:MyVec{N, T}}
        new{N, T, AT}(a, β)
    end
end

HalfSpace{N, T}(a::AT, β::T) where {N, T, AT <: MyVec{N, T}} = HalfSpace{N, T, AT}(a, β)
HalfSpace{N, T}(a::MyVec, β) where {N, T} = HalfSpace{N, T}(myvec(T, a), T(β))

"""
    struct HyperPlane{N, T, AT} <: HRepElement{N, T}
        a::AT
        β::T
    end

An hyperplane defined by the set of points ``x`` such that ``\\langle a, x \\rangle = \\beta``.
"""
struct HyperPlane{N, T, AT<:MyVec{N, T}} <: HRepElement{N,T}
    a::AT
    β::T
    function HyperPlane{N, T, AT}(a::AT, β::T) where {N, T, AT<:MyVec{N, T}}
        new{N, T, AT}(a, β)
    end
end

HyperPlane{N, T}(a::AT, β::T) where {N, T, AT <: MyVec{N, T}} = HyperPlane{N, T, AT}(a, β)
HyperPlane{N, T}(a::MyVec, β) where {N, T} = HyperPlane{N, T}(myvec(T, a), T(β))

_HalfSpace(a::MyVec{N, T}, β::T, d::FullDim{N}) where {N, T} = HalfSpace{N, T}(a, β)
_HalfSpace(a::MyVec{N, S}, β::T, d::FullDim{N}) where {N, S, T} = HalfSpace{N, promote_type(S, T)}(a, β)
HalfSpace(a, β) = _HalfSpace(a, β, FullDim(a))
_HyperPlane(a::MyVec{N, T}, β::T, d::FullDim{N}) where {N, T} = HyperPlane{N, T}(a, β)
_HyperPlane(a::MyVec{N, S}, β::T, d::FullDim{N}) where {N, S, T} = HyperPlane{N, promote_type(S, T)}(a, β)
HyperPlane(a, β) = _HyperPlane(a, β, FullDim(a))

Base.convert(::Type{HalfSpace{N, T, AT}}, h::HalfSpace{N}) where {N, T, AT} = HalfSpace{N, T, AT}(AT(h.a), T(h.β))
Base.convert(::Type{HyperPlane{N, T, AT}}, h::HyperPlane{N}) where {N, T, AT} = HyperPlane{N, T, AT}(AT(h.a), T(h.β))
Base.convert(::Type{HalfSpace{N, T, AT}}, h::HalfSpace{N, T, AT}) where {N, T, AT} = h
Base.convert(::Type{HyperPlane{N, T, AT}}, h::HyperPlane{N, T, AT}) where {N, T, AT} = h

islin(::Union{HalfSpace, Type{<:HalfSpace}}) = false
islin(::Union{HyperPlane, Type{<:HyperPlane}}) = true

(-)(h1::HRepElement, h2::HRepElement) = HalfSpace(h1.a - h2.a, h1.β - h2.β)
(-)(h1::HyperPlane, h2::HyperPlane) = HyperPlane(h1.a - h2.a, h1.β - h2.β)

(*)(h::HyperPlane, α::Real) = HyperPlane(h.a * α, h.β * α)
(*)(α::Real, h::HyperPlane) = HyperPlane(α * h.a, α * h.β)
(*)(h::HalfSpace, α::Real) = HalfSpace(h.a * α, h.β * α)
(*)(α::Real, h::HalfSpace) = HalfSpace(α * h.a, α * h.β)

function (/)(h::ElemT, P::Matrix) where {N, T, ElemT<:HRepElement{N, T}}
    Tout = mypromote_type(T, eltype(P))
    ElemTout = similar_type(ElemT, FullDim{size(P, 2)}(), Tout)
    ElemTout(Matrix{Tout}(P) * myvec(Tout, h.a), Tout(h.β))
end
function zeropad(h::ElemT, n::Integer) where {N, T, ElemT<:HRepElement{N, T}}
    if n == 0
        h
    else
        ElemTout = similar_type(ElemT, FullDim{N+abs(n)}())
        if n < 0
            aout = [zeros(T, -n); h.a]
        else
            aout = [h.a; zeros(T, n)]
        end
        ElemTout(aout, h.β)
    end
end

# Point: -> A same Rep should always return the same of the two types so that when points and sympoints will have different accessors it will be type stable
# Point{N, T} or AbstractVector{T}
# Linear Point:
# SymPoint{N,T}
# Ray:
# Ray{N, T}
# Linear Ray:
# Line{N, T}

"""
    const AbstractPoint{N, T} = Union{Point{N, T}, AbstractVector{T}}

A point in dimension `N` and of coefficient type `T`.
"""
const AbstractPoint{N, T} = Union{Point{N, T}, AbstractVector{T}}

"""
    struct SymPoint{N, T, AT <: AbstractPoint{N, T}}
        a::AT
    end

The convex hull of `a` and `-a`.
"""
struct SymPoint{N, T, AT <: AbstractPoint{N, T}}
    a::AT
    function SymPoint{N, T}(a::AT) where {N, T, AT<:AbstractPoint{N, T}}
        new{N, T, AT}(a)
    end
end
SymPoint{N, T, AT}(sympoint::SymPoint) where {N, T, AT} = SymPoint{N, T, AT}(AT(sympoint.a))
SymPoint{N, T}(a::AbstractPoint) where {N,T} = SymPoint{N,T}(mypoint(T, a))

"""
    struct Ray{N, T, AT <: MyVec{N, T}}
        a::AT
    end

The conic hull of `a`, i.e. the set of points `λa` where `λ` is any nonnegative real number.
"""
struct Ray{N, T, AT <: MyVec{N, T}}
    a::AT
    function Ray{N, T, AT}(a::AT) where {N, T, AT<:MyVec{N, T}}
        new{N, T, AT}(a)
    end
end
Ray{N, T, AT}(ray::Ray) where {N, T, AT} = Ray{N, T, AT}(AT(ray.a))
Ray{N, T}(a::AT) where {N, T, AT<:MyVec{N, T}} = Ray{N, T, AT}(a)
Ray{N, T}(a::MyVec) where {N,T} = Ray{N,T}(myvec(T, a))

"""
    struct Line{N, T, AT <: MyVec{N, T}}
        a::AT
    end

The conic hull of `a` and `-a`, i.e. the set of points `λa` where `λ` is any real number.
"""
struct Line{N, T, AT<:MyVec{N, T}}
    a::AT
    function Line{N, T, AT}(a::AT) where {N, T, AT<:MyVec{N, T}}
        new{N, T, AT}(a)
    end
end
Line{N, T, AT}(line::Line) where {N, T, AT} = Line{N, T, AT}(AT(line.a))
Line{N, T}(a::AT) where {N, T, AT<:MyVec{N, T}} = Line{N, T, AT}(a)

getindex(x::Union{SymPoint,Ray,Line}, i) = x.a[i]
vec(x::Union{SymPoint,Ray,Line}) = vec(x.a)

(-)(h::ElemT) where {ElemT<:Union{HyperPlane, HalfSpace}} = ElemT(-h.a, -h.β)
(-)(elem::ElemT) where {ElemT<:Union{SymPoint,Ray,Line}} = ElemT(-coord(elem))
# Used in remproj
(-)(p::AbstractVector, l::Line) = p - coord(l)
# Ray - Line is done in remproj
(-)(r::Ray, s::Union{Ray, Line}) = Ray(r.a - s.a)
(+)(r::Ray, s::Ray) = Ray(r.a + s.a)
(+)(p::AbstractPoint, r::Ray) = p + coord(r)

for op in [:dot, :cross]
    @eval begin
        $op(x::Union{SymPoint,Ray,Line}, y) = $op(x.a, y)
        $op(x, y::Union{SymPoint,Ray,Line}) = $op(x, y.a)
        $op(x::Union{SymPoint,Ray,Line}, y::Union{SymPoint,Ray,Line}) = $op(x.a, y.a)
    end
end

(*)(α, r::T) where {T<:Union{SymPoint,Ray,Line}} = T(α * r.a)
(*)(r::T, α) where {T<:Union{SymPoint,Ray,Line}} = T(r.a * α)
/(r::T, α) where {T<:Union{SymPoint,Ray,Line}} = T(r.a / α)

Base.convert{T}(::Type{Vector{T}}, x::Union{SymPoint,Ray,Line}) = convert(Vector{T}, x.a)

Line{N,T}(a::MyVec) where {N,T} = Line{N,T}(myvec(T, a))

SymPoint(a::Union{Point,AbstractVector}) = SymPoint{fulldim(a), eltype(a)}(a)
Ray(a::Union{Vec,AbstractVector}) = Ray{fulldim(a), eltype(a)}(a)
Line(a::Union{Vec,AbstractVector}) = Line{fulldim(a), eltype(a)}(a)

for ElemT in (:SymPoint, :Ray, :Line)
    @eval begin
        $ElemT{N,Tout}(v::SymPoint{N,Tin}) where {N,Tin,Tout} = $ElemT{N,Tout}(vecconv(Tout, v.a))
        $ElemT{N,T}(v::SymPoint{N,T}) = v
    end
end

const AbstractRay{N, T} = Union{Ray{N, T}, Line{N, T}}
const FixedVRepElement{N,T} = Union{Point{N,T}, Ray{N,T}, Line{N,T}, SymPoint{N,T}}
const VRepElement{N,T} = Union{FixedVRepElement{N,T}, AbstractVector{T}}
const RepElement{N,T} = Union{HRepElement{N,T}, VRepElement{N,T}}
const FixedRepElement{N,T} = Union{HRepElement{N,T}, FixedVRepElement{N,T}}

FullDim(::Union{FixedRepElement{N}, Type{<:FixedRepElement{N}}}) where {N} = FullDim{N}()
MultivariatePolynomials.coefficienttype(::Union{FixedRepElement{N, T}, Type{<:FixedRepElement{N, T}}}) where {N, T} = T

islin(::Union{SymPoint, Line, Type{<:Union{SymPoint, Line}}}) = true
islin(::Union{AbstractPoint, Ray, Type{<:Union{AbstractPoint, Ray}}}) = false
ispoint(::Union{SymPoint, AbstractPoint, Type{<:Union{SymPoint, AbstractPoint}}}) = true
ispoint(::Union{Line, Ray, Type{<:Union{Line, Ray}}}) = false
isray(::Union{SymPoint, AbstractPoint, Type{<:Union{SymPoint, AbstractPoint}}}) = false
isray(::Union{Line, Ray, Type{<:Union{Line, Ray}}}) = true

coord(v::ElemT) where {ElemT<:Union{Point,AbstractVector}} = v
coord(v::ElemT) where {ElemT<:Union{HRepElement,SymPoint,Ray,Line}} = v.a

const VRepElementContainer{N,T} = Union{Ray{N,T}, Line{N,T}, SymPoint{N,T}}

function (*)(P::AbstractMatrix, v::ElemT) where {N, T, ElemT<:VRepElementContainer{N, T}}
      Tout = mypromote_type(T, eltype(P))
      ElemTout = similar_type(ElemT, FullDim{size(P, 1)}(), Tout)
      return ElemTout(P * v.a)
end
# FIXME there seem to be a Julia bug, with Vector, it does not recognize that the zeropad method
#       works
function _zeropad(v, n::Integer)
    if n == 0
        v
    else
        T = MultivariatePolynomials.coefficienttype(v)
        ElemTout = similar_type(typeof(v), FullDim(v) + FullDim{abs(n)}())
        if n < 0
            aout = [zeros(T, -n); coord(v)]
        else
            aout = [coord(v); zeros(T, n)]
        end
        ElemTout(aout)
    end
end
zeropad(v::ElemT, n::Integer) where {N, T, ElemT <: VRepElement{N, T}} = _zeropad(v, n)
zeropad(v::AbstractVector, n::Integer) = _zeropad(v, n)

for ElemT in [:HalfSpace, :HyperPlane, :SymPoint, :Ray, :Line]
    @eval begin
        similar_type(::Type{$ElemT{N,T,AT}}, dout::FullDim{Nout}, ::Type{Tout}) where {N,T,AT,Nout,Tout} = $ElemT{Nout,Tout,similar_type(AT, dout, Tout)}
    end
end

mydot(a::AbstractVector, r::Ray) = mydot(a, r.a)
mydot(a::AbstractVector, l::Line) = mydot(a, l.a)

ininterior{N}(r::Ray{N}, h::HalfSpace{N}) = myneg(mydot(h.a, r))
ininterior{N}(l::Line{N}, h::HalfSpace{N}) = myneg(mydot(h.a, l))
ininterior{N}(p::Point{N}, h::HalfSpace{N}) = mylt(mydot(h.a, p), h.β)
ininterior{N}(p::AbstractVector, h::HalfSpace{N}) = mylt(mydot(h.a, p), h.β)
ininterior{N}(p::SymPoint{N}, h::HalfSpace{N}) = mylt(mydot(h.a, p.p), h.β)

inrelativeinterior(p::VRepElement, h::HalfSpace) = ininterior(p, h)

Base.in(r::Ray{N}, h::HalfSpace{N}) where {N} = mynonpos(mydot(h.a, r))
Base.in(l::Line{N}, h::HalfSpace{N}) where {N} = mynonpos(mydot(h.a, l))
Base.in(p::Point{N}, h::HalfSpace{N}) where {N} = myleq(mydot(h.a, p), h.β)
Base.in(p::AbstractVector, h::HalfSpace{N}) where {N} = myleq(mydot(h.a, p), h.β)
Base.in(p::SymPoint{N}, h::HalfSpace{N}) where {N} = myleq(mydot(h.a, p.p), h.β)

ininterior(p::VRepElement, h::HyperPlane) = false
inrelativeinterior(p::VRepElement, h::HyperPlane) = p in h

Base.in(r::Ray{N}, h::HyperPlane{N}) where {N} = myeqzero(mydot(h.a, r))
Base.in(l::Line{N}, h::HyperPlane{N}) where {N} = myeqzero(mydot(h.a, l))
Base.in(p::Point{N}, h::HyperPlane{N}) where {N} = myeq(mydot(h.a, p), h.β)
Base.in(p::AbstractVector, h::HyperPlane{N}) where {N} = myeq(mydot(h.a, p), h.β)
Base.in(p::SymPoint{N}, h::HyperPlane{N}) where {N} = myeq(mydot(h.a, p.p), h.β)

import Base.vec
function Base.vec(x::FixedVector{N,T}) where {N,T}
    y = Vector{T}(N)
    for i in 1:N
        y[i] = x[i]
    end
    y
end
#Base.vec{ElemT::VRepElementContainer}(x::ElemT) = ElemT(vec(x.a))

function pushbefore(a::AbstractSparseVector{T}, β::T) where T
    b = spzeros(T, length(a)+1)
    b[1] = β
    b[2:end] = a
    b
end
pushbefore(a::AbstractVector, β) = [β; a]
function pushbefore(a::ElemT, β, ElemTout = similar_type(ElemT, FullDim{N+1}())) where {N, ElemT<:FixedVector{N}}
    ElemTout([β; vec(a)])
end

function lift(h::ElemT) where {N, T, ElemT <: HRepElement{N, T}}
    ElemT(pushbefore(h.a, -h.β), zero(T))
end
lift(h::Ray{N,T}) where {N,T} = Ray{N+1,T}(pushbefore(h.a, zero(T)))
lift(h::Line{N,T}) where {N,T} = Line{N+1,T}(pushbefore(h.a, zero(T)))
lift(h::Point{N,T}) where {N,T} = pushbefore(h, one(T))
lift(h::AbstractVector{T}) where {T} = Ray{length(h)+1,T}(pushbefore(h, one(T)))
lift(h::SymPoint{N,T}) where {N,T} = Line{N+1,T}(pushbefore(h.a, one(T)))

translate(p::Union{Point,AbstractVector,SymPoint}, v) = p + v
translate(r::Union{Ray,Line}, v) = r

translate(h::ElemT, p) where {ElemT<:HRepElement} = ElemT(h.a, h.β + mydot(h.a, p))

_simplify(a::AbstractVector) = a
_simplify(a::AbstractVector, β) = a, β

iszo(g) = iszero(g) || g == 1 # TODO isone in v0.7

function _simplify(a::AbstractVector{<:Integer})
    g = gcd(a)
    if !iszo(g)
        div.(a, g)
    else
        a
    end
end
function _simplify(a::AbstractVector{<:Integer}, β::Integer)
    g = gcd(gcd(a), β)
    if !iszo(g)
        div.(a, g), div(β, g)
    else
        a, β
    end
end
function _simplify(a::AbstractVector{<:Rational})
    g = gcd(numerator.(a))
    if !iszo(g)
        a = a ./ g
    end
    g = gcd(denominator.(a))
    if !iszo(g)
        a .* g
    else
        a
    end
end
function _simplify(a::AbstractVector{<:Rational}, β::Rational)
    g = gcd(gcd(numerator.(a)), numerator(β))
    if !iszo(g)
        a = a ./ g
        β = β / g
    end
    g = gcd(gcd(denominator.(a)), denominator(β))
    if !iszo(g)
        a .* g, β * g
    else
        a, β
    end
end

simplify(h::HalfSpace{N, T}) where {N, T} = HalfSpace{N, T}(_simplify(h.a, h.β)...)
simplify(h::HyperPlane{N, T}) where {N, T} = HyperPlane{N, T}(_simplify(h.a, h.β)...)
simplify(r::T) where {T<:Union{Ray,Line}} = T(_simplify(coord(r)))
# Cannot scale points
simplify(p::Union{Point, AbstractVector, SymPoint}) = p

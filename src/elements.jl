import GeometryTypes.Point
export Point
export HRepElement, HalfSpace, HyperPlane
export VRepElement, AbstractPoint, AbstractRay, Ray, Line
#export SymPoint
export islin, isray, ispoint, coord, lift, simplify

const MyVec{N,T} = Union{Vec{N,T},AbstractVector{T}}
_vec{T}(::Type{T}, a::AbstractVector) = AbstractArray{T}(a)
_vec{T}(::Type{T}, a::AbstractVector{T}) = a
function _vec{N,T}(::Type{T}, a::StaticArrays.SVector{N})
    StaticArrays.SVector{N,T}(a)
end
_vec{N,T}(::Type{T}, a::StaticArrays.SVector{N,T}) = a
_vec{N,T}(::Type{T}, a::Vec{N}) = Vec{N,T}(a)
_vec{N,T}(::Type{T}, a::Vec{N,T}) = a

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
HalfSpace{N, T}(a::MyVec, β) where {N, T} = HalfSpace{N, T}(_vec(T, a), T(β))

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
HyperPlane{N, T}(a::MyVec, β) where {N, T} = HyperPlane{N, T}(_vec(T, a), T(β))

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

Base.:-(h1::HRepElement, h2::HRepElement) = HalfSpace(h1.a - h2.a, h1.β - h2.β)
Base.:-(h1::HyperPlane, h2::HyperPlane) = HyperPlane(h1.a - h2.a, h1.β - h2.β)

Base.:(*)(h::HyperPlane, α::Real) = HyperPlane(h.a * α, h.β * α)
Base.:(*)(α::Real, h::HyperPlane) = HyperPlane(α * h.a, α * h.β)
Base.:(*)(h::HalfSpace, α::Real) = HalfSpace(h.a * α, h.β * α)
Base.:(*)(α::Real, h::HalfSpace) = HalfSpace(α * h.a, α * h.β)

function Base.:(/)(h::ElemT, P::Matrix) where {N, T, ElemT<:HRepElement{N, T}}
    Tout = _promote_type(T, eltype(P))
    ElemTout = similar_type(ElemT, FullDim{size(P, 2)}(), Tout)
    ElemTout(Matrix{Tout}(P) * _vec(Tout, h.a), Tout(h.β))
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

origin(::Type{<:SparseVector{T}}, ::FullDim{N}) where {N, T} = spzeros(T, N)
origin(::Type{Vector{T}}, ::FullDim{N}) where {N, T} = zeros(T, N)
origin(VT::Type{<:AbstractVector}, ::FullDim) = zeros(VT)
# Canonical basis vector
function basis(::Type{Vector{T}}, ::FullDim{N}, i::Int) where {N, T}
    v = zeros(T, N)
    v[i] = one(T)
    v
end
function basis(::Type{SVector{N, T}}, ::FullDim{N}, i::Int) where {N, T}
    SVector{N, T}(ntuple(j -> j == i ? one(T) : zero(T), Val{N}))
end

#"""
#    struct SymPoint{N, T, AT <: AbstractPoint{N, T}}
#        a::AT
#    end
#
#The convex hull of `a` and `-a`.
#"""
#struct SymPoint{N, T, AT <: AbstractPoint{N, T}}
#    a::AT
#    function SymPoint{N, T, AT}(a::AT) where {N, T, AT<:AbstractPoint{N, T}}
#        new{N, T, AT}(a)
#    end
#end
#SymPoint{N, T, AT}(sympoint::SymPoint) where {N, T, AT} = SymPoint{N, T, AT}(AT(sympoint.a))
#SymPoint{N, T}(a::AbstractPoint{N, T}) where {N, T} = SymPoint{N, T, typeof(a)}(a)
#SymPoint(a::AbstractPoint) = SymPoint{fulldim(a), eltype(a)}(a)

#const AnyPoint{N, T} = Union{SymPoint{N, T}, AbstractPoint{N, T}}
const AnyPoint{N, T} = AbstractPoint{N, T}

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
Base.convert(::Type{Ray{N, T, AT}}, ray::Ray) where {N, T, AT} = Ray{N, T, AT}(ray)
Ray{N, T}(a::AT) where {N, T, AT<:MyVec{N, T}} = Ray{N, T, AT}(a)
Ray(a::MyVec) = Ray{fulldim(a), eltype(a)}(a)

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
Base.convert(::Type{Line{N, T, AT}}, line::Line) where {N, T, AT} = Line{N, T, AT}(line)
Line{N, T}(a::AT) where {N, T, AT<:MyVec{N, T}} = Line{N, T, AT}(a)
Line(a::MyVec) = Line{fulldim(a), eltype(a)}(a)

#const VStruct{N, T} = Union{SymPoint{N, T}, Line{N, T}, Ray{N, T}}
const VStruct{N, T} = Union{Line{N, T}, Ray{N, T}}

Base.:(==)(a::T, b::T) where T<:VStruct = coord(a) == coord(b)
Base.getindex(x::VStruct, i) = x.a[i]
Base.vec(x::VStruct) = vec(x.a)

Base.:-(h::ElemT) where {ElemT<:Union{HyperPlane, HalfSpace}} = ElemT(-h.a, -h.β)
Base.:-(elem::ElemT) where ElemT<:VStruct = ElemT(-coord(elem))
# Used in remproj
Base.:-(p::AbstractVector, l::Line) = p - coord(l)
# Ray - Line is done in remproj
Base.:-(r::Ray, s::Union{Ray, Line}) = Ray(r.a - s.a)
Base.:+(r::Ray, s::Ray) = Ray(r.a + s.a)
Base.:+(p::AbstractPoint, r::Ray) = p + coord(r)

for op in [:dot, :cross]
    @eval begin
        Base.$op(x::VStruct, y) = $op(x.a, y)
        Base.$op(x, y::VStruct) = $op(x, y.a)
        Base.$op(x::VStruct, y::VStruct) = $op(x.a, y.a)
    end
end

for ElT in (:HyperPlane, :HalfSpace, :Line, :Ray)
    @eval begin
        Base.promote_rule(::Type{$ElT{N, T, AT}}, ::Type{$ElT{N, T, AT}}) where {N, T, AT} = $ElT{N, T, AT}
        Base.promote_rule(::Type{$ElT{N, T, AS}}, ::Type{$ElT{N, T, AT}}) where {N, T, AS, AT} = $ElT{N, T, promote_type(AS, AT)}
        function Base.promote_rule(::Type{$ElT{N, S, AS}}, ::Type{$ElT{N, T, AT}}) where {N, S, T, AS, AT}
            U = promote_type(S, T)
            $ElT{N, U, promote_type(similar_type(AS, U), similar_type(AT, U))}
        end
    end
end

Base.:(*)(α, r::T) where T<:VStruct = T(α * r.a)
Base.:(*)(r::T, α) where T<:VStruct = T(r.a * α)
Base.:(/)(r::T, α) where T<:VStruct = T(r.a / α)

const AbstractRay{N, T} = Union{Ray{N, T}, Line{N, T}}
const FixedVRepElement{N,T} = Union{Point{N, T}, VStruct{N, T}}
const VRepElement{N,T} = Union{FixedVRepElement{N,T}, AbstractVector{T}}
const RepElement{N,T} = Union{HRepElement{N,T}, VRepElement{N,T}}
const FixedRepElement{N,T} = Union{HRepElement{N,T}, FixedVRepElement{N,T}}
const Element{N, T} = Union{HRepElement{N, T}, VRepElement{N, T}}

FullDim(::Union{FixedRepElement{N}, Type{<:FixedRepElement{N}}}) where {N} = FullDim{N}()
MultivariatePolynomials.coefficienttype(::Union{FixedRepElement{N, T}, Type{<:FixedRepElement{N, T}}}) where {N, T} = T

#islin(::Union{SymPoint, Line, Type{<:Union{SymPoint, Line}}}) = true
islin(::Union{Line, Type{<:Line}}) = true
islin(::Union{AbstractPoint, Ray, Type{<:Union{AbstractPoint, Ray}}}) = false
ispoint(::Union{AnyPoint, Type{<:AnyPoint}}) = true
ispoint(::Union{Line, Ray, Type{<:Union{Line, Ray}}}) = false
isray(v) = !ispoint(v)

coord(v::ElemT) where {ElemT<:Union{Point,AbstractVector}} = v
coord(v::ElemT) where {ElemT<:Union{HRepElement,VStruct}} = v.a

function Base.:*(P::AbstractMatrix, v::ElemT) where {N, T, ElemT<:VStruct{N, T}}
      Tout = _promote_type(T, eltype(P))
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

for ElemT in [:HalfSpace, :HyperPlane, :Ray, :Line] # , :SymPoint
    @eval begin
        similar_type(::Type{$ElemT{N,T,AT}}, dout::FullDim{Nout}, ::Type{Tout}) where {N,T,AT,Nout,Tout} = $ElemT{Nout,Tout,similar_type(AT, dout, Tout)}
    end
end

ininterior{N}(r::Ray{N}, h::HalfSpace{N}) = _neg(h.a ⋅ r)
ininterior{N}(l::Line{N}, h::HalfSpace{N}) = _neg(h.a ⋅ l)
ininterior{N}(p::Point{N}, h::HalfSpace{N}) = _lt(h.a ⋅ p, h.β)
ininterior{N}(p::AbstractVector, h::HalfSpace{N}) = _lt(h.a ⋅ p, h.β)
#ininterior{N}(p::SymPoint{N}, h::HalfSpace{N}) = _lt(h.a ⋅ p.p, h.β)

inrelativeinterior(p::VRepElement, h::HalfSpace) = ininterior(p, h)

Base.in(r::Ray{N}, h::HalfSpace{N}) where N = _nonpos(h.a ⋅ r)
Base.in(l::Line{N}, h::HalfSpace{N}) where N = _nonpos(h.a ⋅ l)
Base.in(p::Point{N}, h::HalfSpace{N}) where N = _leq(h.a ⋅ p, h.β)
Base.in(p::AbstractVector, h::HalfSpace{N}) where N = _leq(h.a ⋅ p, h.β)
#function Base.in(p::SymPoint{N}, h::HalfSpace{N}) where N
#    ap = h.a ⋅ p.p
#    _leq(ap, h.β) && _leq(-ap, h.β)
#end

ininterior(p::VRepElement, h::HyperPlane) = false
inrelativeinterior(p::VRepElement, h::HyperPlane) = p in h

Base.in(r::Ray{N}, h::HyperPlane{N}) where {N} = isapproxzero(h.a ⋅ r)
Base.in(l::Line{N}, h::HyperPlane{N}) where {N} = isapproxzero(h.a ⋅ l)
Base.in(p::Point{N}, h::HyperPlane{N}) where {N} = _isapprox(h.a ⋅ p, h.β)
Base.in(p::AbstractVector, h::HyperPlane{N}) where {N} = _isapprox(h.a ⋅ p, h.β)
#Base.in(p::SymPoint{N}, h::HyperPlane{N}) where {N} = isapproxzero(h.β) && isapproxzero(h.a ⋅ p.a)

function Base.vec(x::FixedVector{N,T}) where {N,T}
    y = Vector{T}(N)
    for i in 1:N
        y[i] = x[i]
    end
    y
end
#Base.vec{ElemT::VStruct}(x::ElemT) = ElemT(vec(x.a))

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

function lift(h::HRepElement{N, T}) where {N, T}
    similar_type(typeof(h), FullDim{N+1}(), T)(pushbefore(h.a, -h.β), zero(T))
end
lift(h::Ray{N,T}) where {N,T} = Ray{N+1,T}(pushbefore(h.a, zero(T)))
lift(h::Line{N,T}) where {N,T} = Line{N+1,T}(pushbefore(h.a, zero(T)))
lift(h::Point{N,T}) where {N,T} = pushbefore(h, one(T))
lift(h::AbstractVector{T}) where {T} = pushbefore(h, one(T))
#lift(h::SymPoint{N,T}) where {N,T} = SymPoint{N+1,T}(pushbefore(h.a, one(T)))

#translate(p::Union{Point,AbstractVector,SymPoint}, v) = p + v
translate(p::Union{Point,AbstractVector}, v) = p + v
translate(r::Union{Ray,Line}, v) = r

translate(h::ElemT, p) where {ElemT<:HRepElement} = ElemT(h.a, h.β + h.a ⋅ p)

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
simplify(p::AnyPoint) = p

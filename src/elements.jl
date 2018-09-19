import GeometryTypes.Point
export Point
export HRepElement, HalfSpace, HyperPlane
export VRepElement, AbstractRay, Ray, Line
#export SymPoint
export islin, isray, ispoint, coord, lift, simplify

_vec(::Type{T}, a::AbstractVector) where {T} = AbstractArray{T}(a)
_vec(::Type{T}, a::AbstractVector{T}) where {T} = a
function _vec(::Type{T}, a::StaticArrays.SVector{N}) where {N, T}
    StaticArrays.SVector{N, T}(a)
end
_vec(::Type{T}, a::StaticArrays.SVector{N, T}) where {N, T} = a

abstract type HRepElement{T} end

"""
    struct HalfSpace{T, AT} <: HRepElement{T}
        a::AT
        β::T
    end

An halfspace defined by the set of points ``x`` such that ``\\langle a, x \\rangle \\le \\beta``.
"""
struct HalfSpace{T, AT <: AbstractVector{T}} <: HRepElement{T}
    a::AT
    β::T
    function HalfSpace{T, AT}(a::AT, β::T) where {T, AT<:AbstractVector{T}}
        new{T, AT}(a, β)
    end
end

HalfSpace{T}(a::AT, β::T) where {T, AT <: AbstractVector{T}} = HalfSpace{T, AT}(a, β)
HalfSpace{T}(a::AbstractVector, β) where {T} = HalfSpace{T}(_vec(T, a), T(β))

"""
    struct HyperPlane{T, AT} <: HRepElement{T}
        a::AT
        β::T
    end

An hyperplane defined by the set of points ``x`` such that ``\\langle a, x \\rangle = \\beta``.
"""
struct HyperPlane{T, AT<:AbstractVector{T}} <: HRepElement{T}
    a::AT
    β::T
    function HyperPlane{T, AT}(a::AT, β::T) where {T, AT<:AbstractVector{T}}
        new{T, AT}(a, β)
    end
end

HyperPlane{T}(a::AT, β::T) where {T, AT <: AbstractVector{T}} = HyperPlane{T, AT}(a, β)
HyperPlane{T}(a::AbstractVector, β) where {T} = HyperPlane{T}(_vec(T, a), T(β))

HalfSpace(a::AbstractVector{T}, β::T) where {T} = HalfSpace{T}(a, β)
HalfSpace(a::AbstractVector{S}, β::T) where {S, T} = HalfSpace{promote_type(S, T)}(a, β)
HyperPlane(a::AbstractVector{T}, β::T) where {T} = HyperPlane{T}(a, β)
HyperPlane(a::AbstractVector{S}, β::T) where {S, T} = HyperPlane{promote_type(S, T)}(a, β)

Base.convert(::Type{HalfSpace{T, AT}}, h::HalfSpace) where {T, AT} = HalfSpace{T, AT}(AT(h.a), T(h.β))
Base.convert(::Type{HyperPlane{T, AT}}, h::HyperPlane) where {T, AT} = HyperPlane{T, AT}(AT(h.a), T(h.β))
Base.convert(::Type{HalfSpace{T, AT}}, h::HalfSpace{T, AT}) where {T, AT} = h
Base.convert(::Type{HyperPlane{T, AT}}, h::HyperPlane{T, AT}) where {T, AT} = h

islin(::Union{HalfSpace, Type{<:HalfSpace}}) = false
islin(::Union{HyperPlane, Type{<:HyperPlane}}) = true

Base.:-(h1::HRepElement, h2::HRepElement) = HalfSpace(h1.a - h2.a, h1.β - h2.β)
Base.:-(h1::HyperPlane, h2::HyperPlane) = HyperPlane(h1.a - h2.a, h1.β - h2.β)

Base.:(*)(h::HyperPlane, α::Real) = HyperPlane(h.a * α, h.β * α)
Base.:(*)(α::Real, h::HyperPlane) = HyperPlane(α * h.a, α * h.β)
Base.:(*)(h::HalfSpace, α::Real) = HalfSpace(h.a * α, h.β * α)
Base.:(*)(α::Real, h::HalfSpace) = HalfSpace(α * h.a, α * h.β)

function Base.:(/)(h::ElemT, P::Matrix) where {T, ElemT<:HRepElement{T}}
    Tout = _promote_type(T, eltype(P))
    ElemTout = similar_type(ElemT, FullDim{size(P, 2)}(), Tout)
    ElemTout(Matrix{Tout}(P) * _vec(Tout, h.a), Tout(h.β))
end
function zeropad(a::Vector{T}, n::Integer)
    if n < 0
        return [zeros(T, -n); a]
    else
        return [a; zeros(T, n)]
    end
end
noparam_type(::HalfSpace) = HalfSpace
noparam_type(::HyperPlane) = HyperPlane
function zeropad(h::HRepElement, n::Integer)
    if n == 0
        return h
    else
        return noparam_type(zeropad(h.a, n), h.β)
    end
end

# Point: -> A same Rep should always return the same of the two types so that when points and sympoints will have different accessors it will be type stable
# AbstractVector{T}
# Linear Point:
# SymPoint{T}
# Ray:
# Ray{T}
# Linear Ray:
# Line{T}

origin(::Type{<:SparseVector{T}}, N::Int) where {T} = spzeros(T, N)
origin(::Type{Vector{T}}, N::Int) where {T} = zeros(T, N)
origin(VT::Type{<:AbstractVector}, ::FullDim) = zeros(VT)
# Canonical basis vector
function basis(::Type{Vector{T}}, N::Int, i::Int) where {T}
    v = zeros(T, N)
    v[i] = one(T)
    v
end
function basis(::Type{SVector{N, T}}, ::Int, i::Int) where {N, T}
    SVector{N, T}(ntuple(j -> j == i ? one(T) : zero(T), Val(N)))
end

"""
    struct Ray{T, AT <: AbstractVector{T}}
        a::AT
    end

The conic hull of `a`, i.e. the set of points `λa` where `λ` is any nonnegative real number.
"""
struct Ray{T, AT <: AbstractVector{T}}
    a::AT
    function Ray{T, AT}(a::AT) where {T, AT<:AbstractVector{T}}
        new{T, AT}(a)
    end
end
Ray{T, AT}(ray::Ray) where {T, AT} = Ray{T, AT}(AT(ray.a))
Base.convert(::Type{Ray{T, AT}}, ray::Ray) where {T, AT} = Ray{T, AT}(ray)
Ray{T}(a::AT) where {T, AT<:AbstractVector{T}} = Ray{T, AT}(a)
Ray(a::AbstractVector) = Ray{fulldim(a), eltype(a)}(a)

"""
    struct Line{T, AT <: AbstractVector{T}}
        a::AT
    end

The conic hull of `a` and `-a`, i.e. the set of points `λa` where `λ` is any real number.
"""
struct Line{T, AT<:AbstractVector{T}}
    a::AT
    function Line{T, AT}(a::AT) where {T, AT<:AbstractVector{T}}
        new{T, AT}(a)
    end
end
Line{T, AT}(line::Line) where {T, AT} = Line{T, AT}(AT(line.a))
Base.convert(::Type{Line{T, AT}}, line::Line) where {T, AT} = Line{T, AT}(line)
Line{T}(a::AT) where {T, AT<:AbstractVector{T}} = Line{T, AT}(a)
Line(a::AbstractVector) = Line{fulldim(a), eltype(a)}(a)

#const VStruct{T} = Union{SymPoint{T}, Line{T}, Ray{T}}
const VStruct{T} = Union{Line{T}, Ray{T}}

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
Base.:+(p::AbstractVector, r::Ray) = p + coord(r)

for op in [:dot, :cross]
    @eval begin
        Compat.LinearAlgebra.$op(x::VStruct, y) = $op(x.a, y)
        Compat.LinearAlgebra.$op(x, y::VStruct) = $op(x, y.a)
        Compat.LinearAlgebra.$op(x::VStruct, y::VStruct) = $op(x.a, y.a)
    end
end

for ElT in (:HyperPlane, :HalfSpace, :Line, :Ray)
    @eval begin
        Base.promote_rule(::Type{$ElT{T, AT}}, ::Type{$ElT{T, AT}}) where {T, AT} = $ElT{T, AT}
        # Allowing mixing e.g. sparse vector with vector would not be helpful as code meant
        # to use sparse polyhedra would lose sparsity silently. Same thing for SVector
        Base.promote_rule(::Type{$ElT{T, VT}}, ::Type{$ElT{T, WT}}) where {T, VT, WT} = error("Cannot mix Polyhedra elements of vector type $VT and $WT")
        function Base.promote_rule(::Type{$ElT{S, AS}}, ::Type{$ElT{T, AT}}) where {S, T, AS, AT}
            U = promote_type(S, T)
            promote_type($ElT{U, similar_type(AS, U)}, $ElT{U, similar_type(AT, U)})
        end
    end
end

Base.:(*)(α, r::T) where T<:VStruct = T(α * r.a)
Base.:(*)(r::T, α) where T<:VStruct = T(r.a * α)
Base.:(/)(r::T, α) where T<:VStruct = T(r.a / α)

const VRepElement{T} = Union{VStruct{T}, AbstractVector{T}}
const RepElement{T} = Union{HRepElement{T}, VRepElement{T}}

FullDim(el::Union{RepElement, Type{<:RepElement}}) = FullDim(el)
MultivariatePolynomials.coefficienttype(::Union{RepElement{T}, Type{<:RepElement{T}}}) where {T} = T

#islin(::Union{SymPoint, Line, Type{<:Union{SymPoint, Line}}}) = true
islin(::Union{Line, Type{<:Line}}) = true
islin(::Union{AbstractVector, Ray, Type{<:Union{AbstractVector, Ray}}}) = false
ispoint(::Union{AbstractVector, Type{<:AbstractVector}}) = true
ispoint(::Union{Line, Ray, Type{<:Union{Line, Ray}}}) = false
isray(v) = !ispoint(v)

coord(v::ElemT) where {ElemT<:Union{Point,AbstractVector}} = v
coord(v::ElemT) where {ElemT<:Union{HRepElement,VStruct}} = v.a

function Base.:*(P::AbstractMatrix, v::ElemT) where {T, ElemT<:VStruct{T}}
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
zeropad(v::ElemT, n::Integer) where {T, ElemT <: VRepElement{T}} = _zeropad(v, n)
zeropad(v::AbstractVector, n::Integer) = _zeropad(v, n)

for ElemT in [:HalfSpace, :HyperPlane, :Ray, :Line] # , :SymPoint
    @eval begin
        function similar_type(::Type{$ElemT{T, AT}}, dout::FullDim, ::Type{Tout}) where {T, AT, Tout}
            return $ElemT{Tout, similar_type(AT, dout, Tout)}
        end
    end
end

ininterior(r::Ray, h::HalfSpace) = _neg(h.a ⋅ r)
ininterior(l::Line, h::HalfSpace) = _neg(h.a ⋅ l)
ininterior(p::Point, h::HalfSpace) = _lt(h.a ⋅ p, h.β)
ininterior(p::AbstractVector, h::HalfSpace) = _lt(h.a ⋅ p, h.β)
#ininterior(p::SymPoint, h::HalfSpace) = _lt(h.a ⋅ p.p, h.β)

inrelativeinterior(p::VRepElement, h::HalfSpace) = ininterior(p, h)

Base.in(r::Ray, h::HalfSpace) = _nonpos(h.a ⋅ r)
Base.in(l::Line, h::HalfSpace) = _nonpos(h.a ⋅ l)
Base.in(p::Point, h::HalfSpace) = _leq(h.a ⋅ p, h.β)
Base.in(p::AbstractVector, h::HalfSpace) = _leq(h.a ⋅ p, h.β)
#function Base.in(p::SymPoint, h::HalfSpace)
#    ap = h.a ⋅ p.p
#    _leq(ap, h.β) && _leq(-ap, h.β)
#end

ininterior(p::VRepElement, h::HyperPlane) = false
inrelativeinterior(p::VRepElement, h::HyperPlane) = p in h

Base.in(r::Ray, h::HyperPlane) = isapproxzero(h.a ⋅ r)
Base.in(l::Line, h::HyperPlane) = isapproxzero(h.a ⋅ l)
Base.in(p::Point, h::HyperPlane) = _isapprox(h.a ⋅ p, h.β)
Base.in(p::AbstractVector, h::HyperPlane) = _isapprox(h.a ⋅ p, h.β)
#Base.in(p::SymPoint, h::HyperPlane) = isapproxzero(h.β) && isapproxzero(h.a ⋅ p.a)

#function Base.vec(x::FixedVector{T}) where {T}
#    y = Vector{T}(N)
#    for i in 1:N
#        y[i] = x[i]
#    end
#    y
#end
#Base.vec{ElemT::VStruct}(x::ElemT) = ElemT(vec(x.a))

function pushbefore(a::AbstractSparseVector{T}, β::T) where T
    b = spzeros(T, length(a)+1)
    b[1] = β
    b[2:end] = a
    b
end
pushbefore(a::AbstractVector, β) = [β; a]
function pushbefore(a::ElemT, β, ElemTout = similar_type(ElemT, increment(FullDim(a))))
    ElemTout([β; vec(a)])
end

function lift(h::HRepElement{T}) where {T}
    similar_type(typeof(h), increment(FullDim(h)), T)(pushbefore(h.a, -h.β), zero(T))
end
lift(h::Ray{T}) where {T} = Ray{T}(pushbefore(h.a, zero(T)))
lift(h::Line{T}) where {T} = Line{T}(pushbefore(h.a, zero(T)))
lift(h::Point{T}) where {T} = pushbefore(h, one(T))
lift(h::AbstractVector{T}) where {T} = pushbefore(h, one(T))
#lift(h::SymPoint{T}) where {T} = SymPoint{T}(pushbefore(h.a, one(T)))

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

simplify(h::HalfSpace{T}) where {T} = HalfSpace{T}(_simplify(h.a, h.β)...)
simplify(h::HyperPlane{T}) where {T} = HyperPlane{T}(_simplify(h.a, h.β)...)
simplify(r::T) where {T<:Union{Ray,Line}} = T(_simplify(coord(r)))
# Cannot scale points
simplify(p::AbstractVector) = p

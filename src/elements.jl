export HRepElement, HalfSpace, HyperPlane
export VRepElement, Ray, Line
export islin, isray, ispoint, coord, lift, simplify

_vec(::Type{T}, a::AbstractVector) where {T} = AbstractArray{T}(a)
_vec(::Type{T}, a::AbstractVector{T}) where {T} = a
function _vec(::Type{T}, a::StaticArrays.SVector{N}) where {N, T}
    StaticArrays.SVector{N, T}(a)
end
_vec(::Type{T}, a::StaticArrays.SVector{N, T}) where {N, T} = a

abstract type HRepElement{T, AT} end

"""
    struct HalfSpace{T, AT} <: HRepElement{T, AT}
        a::AT
        β::T
    end

An halfspace defined by the set of points ``x`` such that ``\\langle a, x \\rangle \\le \\beta``.
"""
struct HalfSpace{T, AT <: AbstractVector{T}} <: HRepElement{T, AT}
    a::AT
    β::T
    function HalfSpace{T, AT}(a::AT, β::T) where {T, AT<:AbstractVector{T}}
        new{T, AT}(a, β)
    end
end

HalfSpace{T}(a::AT, β::T) where {T, AT <: AbstractVector{T}} = HalfSpace{T, AT}(a, β)
HalfSpace{T}(a::AbstractVector, β) where {T} = HalfSpace{T}(_vec(T, a), T(β))

"""
    struct HyperPlane{T, AT} <: HRepElement{T, AT}
        a::AT
        β::T
    end

An hyperplane defined by the set of points ``x`` such that ``\\langle a, x \\rangle = \\beta``.
"""
struct HyperPlane{T, AT<:AbstractVector{T}} <: HRepElement{T, AT}
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

function Base.convert(::Type{HalfSpace{T, AT}}, h::HalfSpace) where {T, AT}
    return HalfSpace{T, AT}(convert(AT, h.a), convert(T, h.β))
end
function Base.convert(::Type{HyperPlane{T, AT}}, h::HyperPlane) where {T, AT}
    return HyperPlane{T, AT}(convert(AT, h.a), convert(T, h.β))
end
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

function Base.:(/)(h::ElemT, P::UniformScaling) where {T, ElemT<:HRepElement{T}}
    Tout = _promote_type(T, eltype(P))
    ElemTout = similar_type(ElemT, FullDim(h), Tout)
    ElemTout(P * _vec(Tout, h.a), Tout(h.β))
end
function Base.:(/)(h::ElemT, P::AbstractMatrix) where {T, ElemT<:HRepElement{T}}
    Tout = _promote_type(T, eltype(P))
    ElemTout = similar_type(ElemT, size(P, 2), Tout)
    ElemTout(AbstractMatrix{Tout}(P) * _vec(Tout, h.a), Tout(h.β))
end

# Point: -> A same Rep should always return the same of the two types so that when points and sympoints will have different accessors it will be type stable
# AbstractVector{T}
# Ray:
# Ray{T}
# Linear Ray:
# Line{T}

origin(::Type{<:SparseVector{T}}, N::Int) where {T} = spzeros(T, N)
origin(::Type{Vector{T}}, N::Int) where {T} = zeros(T, N)
origin(VT::Type{<:AbstractVector}, ::Int) = zeros(VT)
# Canonical basis vector
function basis(::Type{Vector{T}}, N::Int, i::Int) where {T}
    v = zeros(T, N)
    v[i] = one(T)
    v
end
function basis(::Type{StaticArrays.SVector{N, T}}, ::StaticArrays.Size, i::Int) where {N, T}
    StaticArrays.SVector{N, T}(ntuple(j -> j == i ? one(T) : zero(T), Val(N)))
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
Ray(a::AbstractVector) = Ray{eltype(a)}(a)

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
Line(a::AbstractVector) = Line{eltype(a)}(a)

Base.:(==)(a::Ray, b::Ray) = coord(a) == coord(b)
Base.:(==)(a::Line, b::Line) = coord(a) == coord(b)

const VStruct{T, AT} = Union{Line{T, AT}, Ray{T, AT}}
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
        LinearAlgebra.$op(x::VStruct, y) = $op(x.a, y)
        LinearAlgebra.$op(x, y::VStruct) = $op(x, y.a)
        LinearAlgebra.$op(x::VStruct, y::VStruct) = $op(x.a, y.a)
    end
end

for ElT in (:HyperPlane, :HalfSpace, :Line, :Ray)
    @eval begin
        Base.promote_rule(::Type{$ElT{T, AT}}, ::Type{$ElT{T, AT}}) where {T, AT} = $ElT{T, AT}
        # Allowing mixing e.g. sparse vector with vector would not be helpful as code meant
        # to use sparse polyhedra would lose sparsity silently. Same thing for StaticArrays.SVector
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
const StructElement{T, AT} = Union{VStruct{T, AT}, HRepElement{T, AT}}

vectortype(::Type{<:StructElement{T, AT}}) where {T, AT} = AT
vectortype(AT::Type{<:AbstractVector}) = AT

FullDim(::Type{<:StructElement{T, AT}}) where {T, AT} = FullDim(AT)
FullDim(el::StructElement) = FullDim(coord(el))
coefficient_type(::Union{RepElement{T}, Type{<:RepElement{T}}}) where {T} = T

islin(::Union{Line, Type{<:Line}}) = true
islin(::Union{AbstractVector, Ray, Type{<:Union{AbstractVector, Ray}}}) = false
ispoint(::Union{AbstractVector, Type{<:AbstractVector}}) = true
ispoint(::Union{Line, Ray, Type{<:Union{Line, Ray}}}) = false
isray(v) = !ispoint(v)

coord(v::AbstractVector) = v
coord(v::Union{HRepElement, VStruct}) = v.a

function Base.:*(P::AbstractMatrix, v::ElemT) where {T, ElemT<:VStruct{T}}
      Tout = _promote_type(T, eltype(P))
      ElemTout = similar_type(ElemT, size(P, 1), Tout)
      return ElemTout(P * v.a)
end

function zeropad(a::AbstractSparseVector{T}, n::Integer) where T
    if iszero(n)
        return a
    elseif n < 0
        # Add type assert to check that vcat does not make it dense
        return [spzeros(T, -n); a]::AbstractSparseVector{T}
    else
        return [a; spzeros(T, n)]::AbstractSparseVector{T}
    end
end
function (zeropad(a::StaticArrays.SVector{N1, T}, ::StaticArrays.Size{N2})::StaticArrays.SVector{N1+abs(N2[1]), T}) where {T, N1, N2}
    # TODO, should probably de a generated function to make it type stable
    if iszero(N2[1])
        return a
    else
        z = StaticArrays.@SVector zeros(T, abs(N2[1]))
        if N2[1] < 0
            return vcat(z, a)::StaticArrays.SVector{N1-N2[1], T}
        else
            return vcat(a, z)::StaticArrays.SVector{N1+N2[1], T}
        end
    end
end
function zeropad(a::Vector{T}, n::Integer) where T
    if iszero(n)
        return a
    elseif n < 0
        return [zeros(T, -n); a]
    else
        return [a; zeros(T, n)]
    end
end
zeropad(h::HRepElement, d::FullDim) = constructor(h)(zeropad(h.a, d), h.β)
zeropad(v::VStruct, d::FullDim)     = constructor(v)(zeropad(v.a, d))

for ElemT in [:HalfSpace, :HyperPlane, :Ray, :Line]
    @eval begin
        function similar_type(::Type{$ElemT{T, AT}}, dout::FullDim, ::Type{Tout}) where {T, AT, Tout}
            return $ElemT{Tout, similar_type(AT, dout, Tout)}
        end
    end
end

ininterior(r::Ray, h::HalfSpace) = _neg(h.a ⋅ r)
ininterior(l::Line, h::HalfSpace) = _neg(h.a ⋅ l)
ininterior(p::AbstractVector, h::HalfSpace) = _lt(h.a ⋅ p, h.β)

inrelativeinterior(p::VRepElement, h::HalfSpace) = ininterior(p, h)

Base.in(r::Ray, h::HalfSpace) = _nonpos(h.a ⋅ r)
Base.in(l::Line, h::HalfSpace) = _nonpos(h.a ⋅ l)
Base.in(p::AbstractVector, h::HalfSpace) = _leq(h.a ⋅ p, h.β)

ininterior(p::VRepElement, h::HyperPlane) = false
inrelativeinterior(p::VRepElement, h::HyperPlane) = p in h

Base.in(r::Ray, h::HyperPlane) = isapproxzero(h.a ⋅ r)
Base.in(l::Line, h::HyperPlane) = isapproxzero(h.a ⋅ l)
Base.in(p::AbstractVector, h::HyperPlane) = _isapprox(h.a ⋅ p, h.β)

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
function pushbefore(a::StaticArrays.SVector, β)
    StaticArrays.SVector(β, a...)
end

constructor(::HyperPlane) = HyperPlane
constructor(::HalfSpace) = HalfSpace
constructor(::Ray) = Ray
constructor(::Line) = Line

function lift(h::HRepElement{T}) where {T}
    constructor(h)(pushbefore(h.a, -h.β), zero(T))
end
lift(v::VStruct{T}) where {T} = constructor(v)(pushbefore(v.a, zero(T)))
lift(v::AbstractVector{T}) where {T} = pushbefore(v, one(T))

translate(p::AbstractVector, v) = p + v
translate(r::VStruct, v) = r

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
simplify(r::VStruct) = constructor(r)(_simplify(coord(r)))
# Cannot scale points
simplify(p::AbstractVector) = p

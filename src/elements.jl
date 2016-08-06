mydot(a, b) = sum(map(*, a, b))
mydot(a::AbstractVector, b::AbstractVector) = dot(f, v)
mydot(f::FixedVector, v::FixedVector) = dot(f, v)
#mydot{T<:Union{Ray,Line}}(a::T, b::T) = mydot(a.r, b.r)
#mydot{T<:Union{Ray,Line}}(a, b::T) = mydot(a, b.r)
#mydot{T<:Union{Ray,Line}}(a::T, b) = mydot(a.r, b)

abstract HRepElement{N,T}

# <a, x> <= β
type HalfSpace{N,T} <: HRepElement{N,T}
  a
  β
end

type HyperPlane{N,T} <: HRepElement{N,T}
  a
  β
end

islin(v::HalfSpace) = false
islin(v::HyperPlane) = true

function (*){ElemT<:HRepElement}(h::ElemT, P::Matrix)
  Tout = mypromote_type(eltype(ElemT), eltype(P))
  ElemTout = changeboth(ElemT, size(P, 2), Tout)
  ElemTout(h.a' * P, h.β)
end
function zeropad{ElemT<:HRepElement}(h::ElemT, n::Integer)
  if n == 0
    h
  else
    ElemTout = changefulldim(ElemT, fulldim(h) + abs(n))
    T = eltype(ElemT)
    if n < 0
      aout = [spzeros(T, n); h.a]
    else
      aout = [h.a; spzeros(T, n)]
    end
    ElemTout(aout, h.β)
  end
end
function Base.round{ElemT<:HRepElement}(h::ElemT, P::Matrix)
  ElemT(round(h.a), round(h.β))
end

# Point:
# Point{N, T} or AbstractVector{T}
# Linear Point:
# SymPoint{N,T}
# Ray:
# Vec{N, T} or Ray{N, T}
# Linear Ray:
# Line{N, T}

type SymPoint{N, T}
  a
end

type Ray{N, T}
  a
end

type Line{N,T}
  a
end

typealias FixedVRepElement{N,T} Union{Point{N,T},Vec{N,T},Ray{N,T},Line{N,T},SymPoint{N,T}}
typealias VRepElement{N,T} Union{FixedVRepElement{N,T},AbstractVector{T}}
typealias RepElement{N,T} Union{HRepElement{N,T}, VRepElement{N,T}}
typealias FixedRepElement{N,T} Union{HRepElement{N,T}, FixedVRepElement{N,T}}

fulldim{N}(e::FixedRepElement{N}) = N
eltype{N,T}(e::FixedRepElement{N, T}) = T
fulldim(v::AbstractVector) = length(v)
eltype{T}(v::AbstractVector{T}) = T

islin{T<:Union{Point,AbstractVector,Vec,Ray}}(v::T) = false
islin{T<:Union{SymPoint,Line}}(v::T) = true
isray{T<:Union{Point,AbstractVector,SymPoint}}(v::T) = false
isray{T<:Union{Vec,Ray,Line}}(v::T) = true

#coord{ElemT<:Union{Point,AbstractVector,Vec}}(v::ElemT) = v
#coord{ElemT<:Union{SymPoint,Ray,Line}}(v::ElemT) = v.a

typealias VRepElementContainer Union{SymPoint, Ray, Line}

@generated function (*){ElemT<:VRepElement}(P::Matrix, v::ElemT)
  if ElemT <: VRepElementContainer
    Tout = mypromote_type(eltype(ElemT), eltype(P))
    ElemTout = changeboth(ElemT, size(P, 1), Tout)
    return :(:ElemTout(P * v.a))
  else
    return :(P * v.a)
  end
end
function zeropad{ElemT<:VRepElement}(v::ElemT, n::Integer)
  if n == 0
    v
  else
    ElemTout = changefulldim(ElemT, fulldim(h) + abs(n))
    T = eltype(ElemT)
    if n < 0
      aout = [spzeros(T, n); coef(v)]
    else
      aout = [coef(v); spzeros(T, n)]
    end
    ElemTout(aout)
  end
end

for ElemT in [:HalfSpace, :HyperPlane, :Point, :SymPoint, :Vec, :Ray, :Line]
  @eval begin
    changeeltype{N,T,Tout}(::Type{$ElemT{N,T}}, ::Type{Tout}) = $ElemT{N,Tout}
    changefulldim{N,T}(::Type{$ElemT{N,T}}, Nout) = $ElemT{Nout,T}
    changeboth{N,T,Tout}(::Type{$ElemT{N,T}}, Nout, ::Type{Tout}) = $ElemT{Nout,Tout}
  end
end

changeeltype{VecT<:AbstractVector,Tout}(::Type{VecT}, ::Type{Tout}) = AbstractVector{Tout}
changefulldim{VecT<:AbstractVector}(::Type{VecT}, Nout) = VecT
changeboth{VecT<:AbstractVector,Tout}(::Type{VecT}, Nout, ::Type{Tout}) = AbstractVector{Tout}

Base.in{N}(h::HalfSpace{N}, r::Vec{N}) = mynonpos(mydot(h.a, r))
Base.in{N}(h::HalfSpace{N}, r::Ray{N}) = mynonpos(mydot(h.a, r.r))
Base.in{N}(h::HalfSpace{N}, l::Line{N}) = mynonpos(mydot(h.a, l.r))
Base.in{N}(h::HalfSpace{N}, p::Point{N}) = myleq(mydot(h.a, p), β)
Base.in{N}(h::HalfSpace{N}, p::AbstractVector) = myleq(mydot(h.a, p), β)
Base.in{N}(h::HalfSpace{N}, p::SymPoint{N}) = myleq(mydot(h.a, p.p), β)

Base.in{N}(h::HyperPlane{N}, r::Vec{N}) = myeqzero(mydot(h.a, r))
Base.in{N}(h::HyperPlane{N}, r::Ray{N}) = myeqzero(mydot(h.a, r.r))
Base.in{N}(h::HyperPlane{N}, l::Line{N}) = myeqzero(mydot(h.a, l.r))
Base.in{N}(h::HyperPlane{N}, p::Point{N}) = myeq(mydot(h.a, p), β)
Base.in{N}(h::HyperPlane{N}, p::AbstractVector) = myeq(mydot(h.a, p), β)
Base.in{N}(h::HyperPlane{N}, p::SymPoint{N}) = myeq(mydot(h.a, p.p), β)

import Base.vec
function Base.vec{N,T}(x::FixedVector{N,T})
  y = Vector{T}(N)
  for i in 1:N
    y[i] = x[i]
  end
  y
end
#Base.vec{ElemT::VRepElementContainer}(x::ElemT) = ElemT(vec(x.a))

pushbefore(a::AbstractVector, β) = [β, a]
function pushbefore{ElemT<:FixedVector}(a::ElemT, β, ElemTout = changefulldim(ElemT, fulldim(ElemT)+1))
  ElemTout([β, vec(a)])
end

function lift{ElemT <: HRepElement}(h::ElemT)
  ElemT(pushbefore(h.a, h.β), zero(eltype(ElemT)))
end
lift{N,T}(h::Vec{N,T}) = pushbefore(h, zero(T))
lift{N,T}(h::Ray{N,T}) = Ray{N+1,T}(pushbefore(h.a, zero(T)))
lift{N,T}(h::Line{N,T}) = Line{N+1,T}(pushbefore(h.a, zero(T)))
lift{N,T}(h::Point{N,T}) = pushbefore(h, one(T)), Vec{N+1,T}
lift{T}(h::AbstractVector{T}) = Ray{length(h)+1,T}(pushbefore(h.a, one(T)))
lift{N,T}(h::SymPoint{N,T}) = Line{N+1,T}(pushbefore(h.a, one(T)))

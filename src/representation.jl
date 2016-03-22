import Base.round, Base.eltype

export Representation, HRepresentation, VRepresentation, splitvertexrays!, fulldim

abstract Representation{T <: Real}

Base.eltype{T <: Real}(desc::Representation{T}) = T

# No copy since I do not modify anything and a copy is done when building a polyhedron

# myfree is for wrapper such as CDoubleRepresentation which use GMPRational

type HRepresentation{T <: Real} <: Representation{T}
  # Ax <= b
  A::Array{T, 2}
  b::Array{T, 1}
  linset::IntSet

  function HRepresentation(A::Array{T, 2}, b::Array{T, 1}, linset::IntSet)
    if size(A, 1) != length(b)
      error("The length of b must be equal to the number of rows of A")
    end
    if ~isempty(linset) && last(linset) > length(b)
      error("The elements of linset should be between 1 and the number of rows of A/length of b")
    end
    ine = new(A, b, linset)
    finalizer(ine, myfree)
    ine
  end
end

function myfree{T<:Real}(ine::HRepresentation{T})
  # Nothing to free
end

function HRepresentation{S <: Real, T <: Real}(A::Matrix{S}, b::Vector{T}, linset::IntSet=IntSet([]))
  U = promote_type(S, T)
  HRepresentation{U}(Matrix{U}(A), Vector{U}(b), linset)
end
HRepresentation{T <: Real}(A::Array{T, 2}, b::Array{T, 1}, linset::IntSet=IntSet([])) = HRepresentation{T}(A, b, linset)

fulldim(ine::HRepresentation) = size(ine.A, 2)

Base.round{T<:AbstractFloat}(ine::HRepresentation{T}) = HRepresentation{T}(Base.round(ine.A), Base.round(ine.b), ine.linset)

function Base.intersect{T}(ine1::HRepresentation{T}, ine2::HRepresentation{T})
  A = [ine1.A; ine2.A]
  b = [ine1.b; ine2.b]
  linset = copy(ine1.linset)
  for i in 1:size(ine2.A, 1)
    if i in ine2.linset
      push!(linset, size(ine1.A, 1) + i)
    end
  end
  HRepresentation{T}(A, b, linset)
end

type VRepresentation{T <: Real} <: Representation{T}
  V::Array{T, 2} # each row is a vertex/ray
  R::Array{T, 2} # rays
  vertex::IntSet # vertex or ray in V
  Vlinset::IntSet # linear or not
  Rlinset::IntSet # linear or not

  function VRepresentation(V::Array{T, 2}, R::Array{T, 2}, vertex::IntSet, Vlinset::IntSet, Rlinset::IntSet)
    if length(R) > 0 && length(V) > 0 && size(V, 2) != size(R, 2)
      error("The dimension of the vertices and rays should be the same")
    end
    if ~isempty(vertex) && last(vertex) > size(V, 1)
      error("The elements of vertex should be between 1 and the number of rows of V")
    end
    if ~isempty(Vlinset) && last(Vlinset) > size(V, 1)
      error("The elements of Vlinset should be between 1 and the number of rows of V")
    end
    if ~isempty(Rlinset) && last(Rlinset) > size(R, 1)
      error("The elements of Rlinset should be between 1 and the number of rows of R")
    end
    desc = new(V, R, vertex, Vlinset, Rlinset)
    finalizer(desc, myfree)
    desc
  end
end

function myfree{T<:Real}(desc::VRepresentation{T})
  # Nothing to free
end

VRepresentation{T}(V::Array{T, 2}, R::Array{T, 2}, vertex::IntSet=IntSet(1:size(V,1)), Vlinset::IntSet=IntSet([]), Rlinset::IntSet=IntSet([])) = VRepresentation{T}(V, R, vertex, Vlinset, Rlinset)

VRepresentation{T <: Real}(V::Array{T, 2}, vertex::IntSet=IntSet(1:size(V,1)), linset::IntSet=IntSet([])) = VRepresentation(V, Matrix{T}(0, size(V, 2)), vertex, linset, IntSet([]))

fulldim(ext::VRepresentation) = size(ext.V, 2)

Base.round{T<:AbstractFloat}(ext::VRepresentation{T}) = VRepresentation{T}(Base.round(ext.V), Base.round(ext.R), ext.vertex, ext.Vlinset, ext.Rlinset)

function (+){T<:Real}(ext1::VRepresentation{T}, ext2::VRepresentation{T})
  V = [ext1.V; ext2.V]
  R = [ext1.R; ext2.R]
  vertex = copy(ext1.vertex)
  Vlinset = copy(ext1.Vlinset)
  for i in 1:size(ext2.V, 1)
    if i in ext2.vertex
      push!(vertex, size(ine1.V, 1) + i)
    end
    if i in ext2.Vlinset
      push!(Vlinset, size(ine1.V, 1) + i)
    end
  end
  Rlinset = copy(ext1.Rlinset)
  for i in 1:size(ext2.R, 1)
    if i in ext2.Rlinset
      push!(Rlinset, size(ine1.R, 1) + i)
    end
  end
  VRepresentation{T}(V, R, vertex, Vlinset, Rlinset)
end

function splitvertexrays!{T<:Real}(ext::VRepresentation{T})
  nV = length(ext.vertex)
  if nV != size(ext.V, 1)
    nR = size(ext.R, 1) + size(ext.V, 1) - nV
    newV = Array(T, nV, size(ext.V, 2))
    newR = Array(T, nR, size(ext.V, 2))
    newR[1:size(ext.R, 1), :] = ext.R
    newVlinset = IntSet([])
    curV = 1
    curR = size(ext.R, 1) + 1
    for i = 1:size(ext.V, 1)
      if i in ext.vertex
        newV[curV, :] = ext.V[i, :]
        if i in ext.Vlinset
          push!(newVlinset, curV)
        end
        curV += 1
      else
        newR[curR, :] = ext.V[i, :]
        if i in ext.Vlinset
          push!(ext.Rlinset, curR)
        end
        curR += 1
      end
    end
    ext.V = newV
    ext.R = newR
    ext.vertex = IntSet(1:nV)
    ext.Vlinset = newVlinset
  end
end

# Representation -> Representation

Base.convert{T, S}(::Type{Representation{T}}, ine::HRepresentation{S}) = Base.convert(HRepresentation{T}, ine)
Base.convert{T, S}(::Type{Representation{T}}, ext::VRepresentation{S}) = Base.convert(VRepresentation{T}, ext)

Base.convert{T, S}(::Type{HRepresentation{T}}, ine::HRepresentation{S}) = HRepresentation{T}(Array{T}(ine.A), Array{T}(ine.b), ine.linset)

Base.convert{T, S}(::Type{VRepresentation{T}}, ext::VRepresentation{S}) = VRepresentation{T}(Array{T}(ext.V), Array{T}(ext.R), ext.vertex, ext.Vlinset, ext.Rlinset)

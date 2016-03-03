import Base.round, Base.eltype
abstract Description{T <: Real}

Base.eltype{T <: Real}(desc::Description{T}) = T

# No copy since I do not modify anything and a copy is done when building a polyhedron

# myfree is for wrapper such as CDoubleDescription which use GMPRational

type InequalityDescription{T <: Real} <: Description{T}
  # Ax <= b
  A::Array{T, 2}
  b::Array{T, 1}
  linset::IntSet

  function InequalityDescription(A::Array{T, 2}, b::Array{T, 1}, linset::IntSet)
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

function myfree{T<:Real}(ine::InequalityDescription{T})
  # Nothing to free
end

function InequalityDescription{S <: Real, T <: Real}(A::Matrix{S}, b::Vector{T}, linset::IntSet=IntSet([]))
  U = promote_type(S, T)
  InequalityDescription{U}(Matrix{U}(A), Vector{U}(b), linset)
end
InequalityDescription{T <: Real}(A::Array{T, 2}, b::Array{T, 1}, linset::IntSet=IntSet([])) = InequalityDescription{T}(A, b, linset)

InequalityDescription{T <: Real}(A::Array{T, 2}, b::Array{T, 1}) = InequalityDescription(A, b, IntSet([]))

Base.round{T<:AbstractFloat}(ine::InequalityDescription{T}) = InequalityDescription{T}(Base.round(ine.A), Base.round(ine.b), ine.linset)

type GeneratorDescription{T <: Real} <: Description{T}
  V::Array{T, 2} # each row is a vertex/ray
  R::Array{T, 2} # rays
  vertex::IntSet # vertex or ray in V
  Vlinset::IntSet # linear or not
  Rlinset::IntSet # linear or not

  function GeneratorDescription(V::Array{T, 2}, R::Array{T, 2}, vertex::IntSet, Vlinset::IntSet, Rlinset::IntSet)
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

function myfree{T<:Real}(desc::GeneratorDescription{T})
  # Nothing to free
end

GeneratorDescription{T}(V::Array{T, 2}, R::Array{T, 2}, vertex::IntSet, Vlinset::IntSet, Rlinset::IntSet) = GeneratorDescription{T}(V, R, vertex, Vlinset, Rlinset)

GeneratorDescription{T <: Real}(V::Array{T, 2}, vertex::IntSet, linset::IntSet) = GeneratorDescription(V, Array{T, 2}(0, size(V, 2)), vertex, linset, IntSet([]))
GeneratorDescription{T <: Real}(V::Array{T, 2}, vertex::IntSet) = GeneratorDescription(V, vertex, IntSet([]))
GeneratorDescription{T <: Real}(V::Array{T, 2}) = GeneratorDescription(V, IntSet(1:size(V,1)))

Base.round{T<:AbstractFloat}(ext::GeneratorDescription{T}) = GeneratorDescription{T}(Base.round(ext.V), Base.round(ext.R), ext.vertex, ext.Vlinset, ext.Rlinset)

function splitvertexrays!{T<:Real}(ext::GeneratorDescription{T})
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

# Description -> Description

Base.convert{T, S}(::Type{Description{T}}, ine::InequalityDescription{S}) = Base.convert(InequalityDescription{T}, ine)
Base.convert{T, S}(::Type{Description{T}}, ext::GeneratorDescription{S}) = Base.convert(GeneratorDescription{T}, ext)

Base.convert{T, S}(::Type{InequalityDescription{T}}, ine::InequalityDescription{S}) = InequalityDescription{T}(Array{T}(ine.A), Array{T}(ine.b), ine.linset)

Base.convert{T, S}(::Type{GeneratorDescription{T}}, ext::GeneratorDescription{S}) = GeneratorDescription{T}(Array{T}(ext.V), Array{T}(ext.R), ext.vertex, ext.Vlinset, ext.Rlinset)

export Description, InequalityDescription, GeneratorDescription, splitvertexrays!

import Base.round, Base.eltype

export Representation, HRepresentation, VRepresentation, SimpleHRepresentation, LiftedHRepresentation, SimpleVRepresentation, LiftedVRepresentation, fulldim

abstract Representation{N, T <: Real}
abstract HRepresentation{N,T} <: Representation{N,T}
abstract VRepresentation{N,T} <: Representation{N,T}

Base.eltype{N,T}(repr::Representation{N,T}) = T
fulldim{N}(repr::Representation{N}) = N

# myfree is for wrapper such as CDDLib which use GMPRational
function myfree{N,T<:Real}(repr::Representation{N,T})
  # Nothing to free
end

# No copy since I do not modify anything and a copy is done when building a polyhedron

type SimpleHRepresentation{N, T} <: HRepresentation{N, T}
  # Ax <= b
  A::Array{T, 2}
  b::Array{T, 1}
  linset::IntSet

  function SimpleHRepresentation(A::Array{T, 2}, b::Array{T, 1}, linset::IntSet=IntSet())
    if size(A, 1) != length(b)
      error("The length of b must be equal to the number of rows of A")
    end
    if ~isempty(linset) && last(linset) > length(b)
      error("The elements of linset should be between 1 and the number of rows of A/length of b")
    end
    if size(A, 2) != N
      error("dimension does not match")
    end
    ine = new(A, b, linset)
    finalizer(ine, myfree)
    ine
  end
end

function SimpleHRepresentation{S <: Real, T <: Real}(A::Matrix{S}, b::Vector{T}, linset::IntSet=IntSet())
  U = promote_type(S, T)
  SimpleHRepresentation{size(A,2),U}(Matrix{U}(A), Vector{U}(b), linset)
end
SimpleHRepresentation{T <: Real}(A::Array{T, 2}, b::Array{T, 1}, linset::IntSet=IntSet()) = SimpleHRepresentation{size(A,2),T}(A, b, linset)

Base.copy{N,T}(ine::SimpleHRepresentation{N,T}) = SimpleHRepresentation{N,T}(copy(ine.A), copy(ine.b), copy(ine.linset))

Base.round{N,T<:AbstractFloat}(ine::SimpleHRepresentation{N,T}) = SimpleHRepresentation{N,T}(Base.round(ine.A), Base.round(ine.b), copy(ine.linset))

function affinehull{N,T}(ine::SimpleHRepresentation{N,T})
  neqs = length(ine.linset)
  A = Matrix{T}(neqs, N)
  b = Vector{T}(neqs)
  cur = 1
  for i in 1:size(ine.A, 1)
    if i in ine.linset
      A[cur, :] = ine.A[i, :]
      b[cur] = ine.b[i]
      cur += 1
    end
  end
  SimpleHRepresentation{N,T}(A, b, IntSet(1:neqs))
end

function Base.intersect{N,T}(ine1::SimpleHRepresentation{N,T}, ine2::SimpleHRepresentation{N,T})
  A = [ine1.A; ine2.A]
  b = [ine1.b; ine2.b]
  linset = copy(ine1.linset)
  for i in 1:size(ine2.A, 1)
    if i in ine2.linset
      push!(linset, size(ine1.A, 1) + i)
    end
  end
  SimpleHRepresentation{N,T}(A, b, linset)
end

function (*){N1, N2, S, T}(ine1::SimpleHRepresentation{N1,S}, ine2::SimpleHRepresentation{N2,T})
  U = promote_type(S, T)
  A = [ine1.A zeros(U, size(ine1.A, 1), size(ine2.A, 2)); zeros(U, size(ine2.A, 1), size(ine1.A, 2)) ine2.A]
  linset = copy(ine1.linset)
  for lin in ine2.linset
    push!(linset, size(ine1.A, 1) + lin)
  end
  SimpleHRepresentation{N1+N2,U}(A, [ine1.b; ine2.b], linset)
end
function (*){N}(ext::SimpleHRepresentation{N}, P::Matrix)
  if size(P, 1) != N
    error("The number of rows of P must match the dimension of the H-representation")
  end
  SimpleHRepresentation(ine.A * P, ine.b, copy(ine.linset))
end

type LiftedHRepresentation{N, T} <: HRepresentation{N, T}
  # Ax >= 0, it is [b -A] * [z; x] where z = 1
  A::Array{T, 2}
  linset::IntSet

  function LiftedHRepresentation(A::Array{T, 2}, linset::IntSet=IntSet())
    if ~isempty(linset) && last(linset) > size(A, 1)
      error("The elements of linset should be between 1 and the number of rows of A")
    end
    if size(A, 2) != N+1
      error("dimension does not match")
    end
    ine = new(A, linset)
    finalizer(ine, myfree)
    ine
  end
end

LiftedHRepresentation{T <: Real}(A::Array{T, 2}, linset::IntSet=IntSet()) = LiftedHRepresentation{size(A,2)-1,T}(A, linset)
Base.copy{N,T}(ine::LiftedHRepresentation{N,T}) = LiftedHRepresentation{N,T}(copy(ine.A), copy(ine.linset))

Base.round{N,T<:AbstractFloat}(ine::LiftedHRepresentation{N,T}) = LiftedHRepresentation{N,T}(Base.round(ine.A), copy(ine.linset))

function affinehull{N,T}(ine::LiftedHRepresentation{N,T})
  neqs = length(ine.linset)
  A = Matrix{T}(neqs, N+1)
  cur = 1
  for i in 1:size(ine.A, 1)
    if i in ine.linset
      A[cur, :] = ine.A[i, :]
      cur += 1
    end
  end
  LiftedHRepresentation{N,T}(A, IntSet(1:neqs))
end

function Base.intersect{N,T}(ine1::LiftedHRepresentation{N,T}, ine2::LiftedHRepresentation{N,T})
  A = [ine1.A; ine2.A]
  linset = copy(ine1.linset)
  for i in 1:size(ine2.A, 1)
    if i in ine2.linset
      push!(linset, size(ine1.A, 1) + i)
    end
  end
  LiftedHRepresentation{N,T}(A, linset)
end

function (*){N1, N2, S, T}(ine1::LiftedHRepresentation{N1,S}, ine2::LiftedHRepresentation{N2,T})
  U = promote_type(S, T)
  A1 = ine1.A[:,2:end]
  A2 = ine2.A[:,2:end]
  b2 = ine2.A[:,1]
  A = [ine1.A zeros(U, size(A1, 1), size(A2, 2)); b2 zeros(U, size(A2, 1), size(A1, 2)) A2]
  linset = copy(ine1.linset)
  for lin in ine2.linset
    push!(linset, size(ine1.A, 1) + lin)
  end
  LiftedHRepresentation{N1+N2,U}(A, linset)
end
function (*){N}(ine::LiftedHRepresentation{N}, P::Matrix)
  if size(P, 1) != N
    error("The number of rows of P must match the dimension of the H-representation")
  end
  LiftedHRepresentation([ine.A[:,1] ine.A[:,2:end] * P], copy(ine.linset))
end

function Base.intersect{N,S,T}(ine1::HRepresentation{N,S}, ine2::HRepresentation{N,T})
  U = promote_type(S, T)
  intersect(Representation{N,U}(ine1), Representation{N,U}(ine2))
end

function Base.intersect{N}(ine1::LiftedHRepresentation{N}, ine2::SimpleHRepresentation{N})
  intersect(ine1, LiftedHRepresentation(ine2))
end
Base.intersect{N}(ine1::SimpleHRepresentation{N}, ine2::LiftedHRepresentation{N}) = intersect(ine2, ine1)


function Base.convert{N,T}(::Type{LiftedHRepresentation{N,T}}, ine::SimpleHRepresentation)
  LiftedHRepresentation{N,T}([ine.b -ine.A], copy(ine.linset))
end
LiftedHRepresentation{N,T}(ext::SimpleHRepresentation{N,T}) = LiftedHRepresentation{N,T}(ext)
function Base.convert{N,T}(::Type{SimpleHRepresentation{N,T}}, ine::LiftedHRepresentation{N,T})
  # The Lifted H-representation is the cone
  # [b -A] [z; x] >= 0
  # To get the Simple H-representation Ax <= 0, we add the equality z = 1
  # and then eliminate the variable z.
  # If we look at this elimination using the Block Elimination method,
  # we see that using z = 1, from the inequality
  # b_i z >= sum A_ij x_j
  # we get
  # b_i >= sum A_ij x_j
  # Note that if A_ij = 0 for all j then we should not keep this inequality:

  # if there is no rows, it returns a Vector{Union{}}
  filter = Vector{Bool}(map(i -> !myeqzero(maximum(abs(ine.A[i,2:end]))), 1:size(ine.A, 1)))
  SimpleHRepresentation{N,T}(-ine.A[filter,2:end], ine.A[filter,1], copy(ine.linset))
end
SimpleHRepresentation{N,T}(ine::LiftedHRepresentation{N,T}) = SimpleHRepresentation{N,T}(ine)

type SimpleVRepresentation{N,T} <: VRepresentation{N,T}
  V::Array{T, 2} # each row is a vertex
  R::Array{T, 2} # each row is a ray
  Vlinset::IntSet
  Rlinset::IntSet

  function SimpleVRepresentation(V::Array{T, 2}, R::Array{T, 2}, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet())
    if length(R) > 0 && size(R, 2) != N
      error("dimension does not match")
    end
    if length(V) > 0 && size(V, 2) != N
      error("dimension does not match")
    end
    if ~isempty(Vlinset) && last(Vlinset) > size(V, 1)
      error("The elements of Vlinset should be between 1 and the number of rows of V")
    end
    if ~isempty(Rlinset) && last(Rlinset) > size(R, 1)
      error("The elements of Rlinset should be between 1 and the number of rows of R")
    end
    repr = new(V, R, Vlinset, Rlinset)
    finalizer(repr, myfree)
    repr
  end
end

SimpleVRepresentation{T <: Real}(V::Array{T, 2}, R::Array{T, 2}, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet()) = SimpleVRepresentation{size(V,2),T}(V, R, Vlinset, Rlinset)

SimpleVRepresentation{T <: Real}(V::Array{T, 2}, linset::IntSet=IntSet()) = SimpleVRepresentation{size(V, 2),T}(V, Matrix{T}(0, size(V, 2)), linset, IntSet())

Base.copy{N,T}(ext::SimpleVRepresentation{N,T}) = SimpleVRepresentation{N,T}(copy(ext.V), copy(ext.R), copy(ext.Vlinset), copy(ext.Rlinset))

Base.round{N,T<:AbstractFloat}(ext::SimpleVRepresentation{N,T}) = SimpleVRepresentation{N,T}(Base.round(ext.V), Base.round(ext.R), copy(ext.Vlinset), copy(ext.Rlinset))

function (+){N,T<:Real}(ext1::SimpleVRepresentation{N,T}, ext2::SimpleVRepresentation{N,T})
  V = [ext1.V; ext2.V]
  R = [ext1.R; ext2.R]
  Vlinset = copy(ext1.Vlinset)
  for i in ext2.Vlinset
    push!(Vlinset, size(ext1.V, 1) + i)
  end
  Rlinset = copy(ext1.Rlinset)
  for i in ext2.Rlinset
    push!(Rlinset, size(ext1.R, 1) + i)
  end
  SimpleVRepresentation{N,T}(V, R, Vlinset, Rlinset)
end

function (*){N}(P::Matrix, ext::SimpleVRepresentation{N})
  if size(P, 2) != N
    error("The number of columns of P must match the dimension of the V-representation")
  end
  SimpleVRepresentation(ext.V * P', ext.R * P', copy(ext.Vlinset), copy(ext.Rlinset))
end

type LiftedVRepresentation{N,T} <: VRepresentation{N,T}
  R::Array{T,2} # each row is a vertex if the first element is 1 and a ray otherwise
  linset::IntSet

  function LiftedVRepresentation(R::Array{T, 2}, linset::IntSet=IntSet([]))
    if length(R) > 0 && size(R, 2) != N+1
      error("dimension does not match")
    end
    if ~isempty(linset) && last(linset) > size(R, 1)
      error("The elements of linset should be between 1 and the number of rows of R")
    end
    repr = new(R, linset)
    finalizer(repr, myfree)
    repr
  end
end

LiftedVRepresentation{T <: Real}(R::Array{T, 2}, linset::IntSet=IntSet()) = LiftedVRepresentation{size(R,2)-1,T}(R, linset)

Base.copy{N,T}(ext::LiftedVRepresentation{N,T}) = LiftedVRepresentation{N,T}(copy(ext.R), copy(ext.linset))

Base.round{N,T<:AbstractFloat}(ext::LiftedVRepresentation{N,T}) = LiftedVRepresentation{N,T}(Base.round(ext.R), copy(ext.linset))

function (+){N,T<:Real}(ext1::LiftedVRepresentation{N,T}, ext2::LiftedVRepresentation{N,T})
  R = [ext1.R; ext2.R]
  linset = copy(ext1.linset)
  for i in ext2.linset
    push!(linset, size(ext1.R, 1) + i)
  end
  LiftedVRepresentation{N,T}(R, linset)
end

function (*){N}(P::Matrix, ext::LiftedVRepresentation{N})
  if size(P, 2) != N
    error("The number of columns of P must match the dimension of the V-representation")
  end
  LiftedVRepresentation([ext.R[:,1] ext.R[:,2:end] * P'], copy(ext.linset))
end

function (+){N,S,T}(ext1::VRepresentation{N,S}, ext2::VRepresentation{N,T})
  U = promote_type(S, T)
  Representation{N,U}(ext1) + Representation{N,U}(ext2)
end

function (+){N}(ext1::LiftedVRepresentation{N}, ext2::SimpleVRepresentation{N})
  ext1 + LiftedVRepresentation(ext2)
end
(+){N}(ext1::SimpleVRepresentation{N}, ext2::LiftedVRepresentation{N}) = ext2 + ext1

function Base.convert{N,T}(::Type{LiftedVRepresentation{N,T}}, ext::SimpleVRepresentation{N,T})
  R = [ones(T, size(ext.V, 1)) ext.V; zeros(T, size(ext.R, 1)) ext.R]
  linset = copy(ext.Vlinset)
  for i in ext.Rlinset
    push!(linset, size(ext.V, 1) + i)
  end
  LiftedVRepresentation{N,T}(R, linset)
end
LiftedVRepresentation{N,T}(ext::SimpleVRepresentation{N,T}) = LiftedVRepresentation{N,T}(ext)

function Base.convert{N,T}(::Type{SimpleVRepresentation{N,T}}, ext::LiftedVRepresentation{N,T})
  rays = IntSet()
  m = size(ext.R, 1)
  for i in 1:m
    if myeqzero(ext.R[i,1])
      push!(rays, i)
    end
  end
  rays = collect(rays)
  R = ext.R[rays, 2:end]
  verts = collect(setdiff(IntSet(1:m), rays))
  V = ext.R[verts, 2:end]
  z = ext.R[verts, 1]
  # z is the lifted coordinate. It should be one.
  # Since we have a cone, it is homogeneous so we can
  # scale the row so that z[i] equals one.
  for i in 1:length(verts)
    if !myeq(z[i], one(T))
      V[i,:] = V[i,:] / z[i]
    end
  end
  Rlinset = IntSet()
  for i in 1:length(rays)
    if rays[i] in ext.linset
      push!(Rlinset, i)
    end
  end
  Vlinset = IntSet()
  for i in 1:length(verts)
    if verts[i] in ext.linset
      push!(Vlinset, i)
    end
  end
  SimpleVRepresentation{N,T}(V,R,Vlinset,Rlinset)
end
SimpleVRepresentation{N,T}(ext::LiftedVRepresentation{N,T}) = SimpleVRepresentation{N,T}(ext)

# Representation -> Representation

Base.convert{N, T, S}(::Type{Representation{N, T}}, ine::SimpleHRepresentation{N, S}) = Base.convert(SimpleHRepresentation{N, T}, ine)
Base.convert{N, T, S}(::Type{Representation{N, T}}, ext::LiftedHRepresentation{N, S}) = Base.convert(LiftedHRepresentation{N, T}, ext)
Base.convert{N, T, S}(::Type{Representation{N, T}}, ext::SimpleVRepresentation{N, S}) = Base.convert(SimpleVRepresentation{N, T}, ext)
Base.convert{N, T, S}(::Type{Representation{N, T}}, ext::LiftedVRepresentation{N, S}) = Base.convert(LiftedVRepresentation{N, T}, ext)

Base.convert{N, T, S}(::Type{SimpleHRepresentation{N, T}}, ine::SimpleHRepresentation{N, S}) = SimpleHRepresentation{N, T}(Array{T}(ine.A), Array{T}(ine.b), ine.linset)
Base.convert{N, T, S}(::Type{LiftedHRepresentation{N, T}}, ine::LiftedHRepresentation{N, S}) = LiftedHRepresentation{N, T}(Array{T}(ine.A), ine.linset)

Base.convert{N, T, S}(::Type{LiftedVRepresentation{N, T}}, ext::LiftedVRepresentation{N, S}) = LiftedVRepresentation{N, T}(Array{T}(ext.R), ext.linset)
Base.convert{N, T, S}(::Type{SimpleVRepresentation{N, T}}, ext::SimpleVRepresentation{N, S}) = SimpleVRepresentation{N, T}(Array{T}(ext.V), Array{T}(ext.R), ext.Vlinset, ext.Rlinset)

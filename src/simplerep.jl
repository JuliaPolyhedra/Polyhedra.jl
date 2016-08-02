export SimpleHRepresentation, SimpleVRepresentation

# H-Representation

# No copy since I do not modify anything and a copy is done when building a polyhedron
type SimpleHRepresentation{N, T} <: HRepresentation{N, T}
  # Ax <= b
  A::AbstractMatrix{T}
  b::AbstractVector{T}
  linset::IntSet

  function SimpleHRepresentation(A::AbstractMatrix{T}, b::AbstractVector{T}, linset::IntSet=IntSet())
    if size(A, 1) != length(b)
      error("The length of b must be equal to the number of rows of A")
    end
    if ~isempty(linset) && last(linset) > length(b)
      error("The elements of linset should be between 1 and the number of rows of A/length of b")
    end
    if size(A, 2) != N
      error("dimension does not match")
    end
    new(A, b, linset)
  end
end

function SimpleHRepresentation{S <: Real, T <: Real}(A::AbstractMatrix{S}, b::AbstractVector{T}, linset::IntSet=IntSet())
  U = promote_type(S, T)
  SimpleHRepresentation{size(A,2),U}(AbstractMatrix{U}(A), AbstractVector{U}(b), linset)
end
SimpleHRepresentation{T <: Real}(A::AbstractMatrix{T}, b::AbstractVector{T}, linset::IntSet=IntSet()) = SimpleHRepresentation{size(A,2),T}(A, b, linset)

function SimpleHRepresentation{N, T}(it::HRepIterator{N, T})
  A = Matrix{T}(length(it), N)
  b = Vector{T}(length(it))
  linset = IntSet()
  for (i, hrep) in enumerate(it)
    A[i,:] = hrep[1]
    b[i] = hrep[2]
    if hrep[3]
      push!(linset, i)
    end
  end
  new(A, b, linset)
end
function SimpleHRepresentation{N, T}(;eqs::Nullable{EqIterator{N, T}}=nothing, ineq::Nullable{IneqIterator{N, T}}=nothing)
  neq = isnull(eqs) ? 0 : length(eqs)
  nineq = isnull(ineqs) ? 0 : length(ineqs)
  nhrep = neq + nineq
  A = Matrix{T}(nhrep, N)
  b = Vector{T}(nhrep)
  linset = IntSet(1:neq)
  if !(eqs === nothing)
    for (i, eq) in enumerate(get(eqs))
      A[i,:] = eq[1]
      b[i] = eq[2]
    end
  end
  if !(ineqs === nothing)
    for (i, ineq) in enumerate(get(ineqs))
      A[neq+i,:] = ineq[1]
      b[neq+i] = ineq[2]
    end
  end
  new(A, b, linset)
end

Base.length(ine::SimpleHRepresentation) = size(ine.A, 1)

Base.copy{N,T}(ine::SimpleHRepresentation{N,T}) = SimpleHRepresentation{N,T}(copy(ine.A), copy(ine.b), copy(ine.linset))

starthrep(ine::SimpleHRepresentation) = 1
donehrep(ine::SimpleHRepresentation, state) = state > length(ine)
nexthrep(ine::SimpleHRepresentation, state) = ((ine.A[state,:], ine.b[state], state in ine.linset), state+1)

neqs(ine::SimpleHRepresentation) = length(ine.linset)
starteq(ine::SimpleHRepresentation) = start(ine.linset)
doneeq(ine::SimpleHRepresentation, state) = done(ine.linset)
function nexteq(ine::SimpleHRepresentation, state)
  (i, nextstate) = next(ine.linset)
  ((ine.A[i,:], ine.b[i]), nextstate)
end

function nextz(is::IntSet, i)
  while i in is
    i += 1
  end
  i
end
nineqs(ine::SimpleHRepresentation) = length(ine) - neqs(ine)
startineq(ine::SimpleHRepresentation) = nextz(ine.linset, 1)
doneineq(ine::SimpleHRepresentation, state) = state > length(ine)
nextineq(ine::SimpleHRepresentation, state) = ((ine.A[state,:], ine.b[state]), nextz(state+1))

# V-Representation

type SimpleVRepresentation{N,T} <: VRepresentation{N,T}
  V::AbstractMatrix{T} # each row is a vertex
  R::AbstractMatrix{T} # each row is a ray
  Vlinset::IntSet
  Rlinset::IntSet

  function SimpleVRepresentation(V::AbstractMatrix{T}, R::AbstractMatrix{T}, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet())
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
    new(V, R, Vlinset, Rlinset)
  end
end

SimpleVRepresentation{T <: Real}(V::AbstractMatrix{T}, R::AbstractMatrix{T}, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet()) = SimpleVRepresentation{size(V,2),T}(V, R, Vlinset, Rlinset)

SimpleVRepresentation{T <: Real}(V::AbstractMatrix{T}, linset::IntSet=IntSet()) = SimpleVRepresentation{size(V, 2),T}(V, similar(V, 0, size(V, 2)), linset, IntSet())

Base.copy{N,T}(ext::SimpleVRepresentation{N,T}) = SimpleVRepresentation{N,T}(copy(ext.V), copy(ext.R), copy(ext.Vlinset), copy(ext.Rlinset))

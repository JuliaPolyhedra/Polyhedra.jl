# No copy since I do not modify anything and a copy is done when building a polyhedron
type LPHRepresentation{N, T} <: HRepresentation{N, T}
  # Ax <= b
  A::AbstractMatrix{T}
  l::AbstractVector{T}
  u::AbstractVector{T}
  colleqs::IntSet
  colgeqs::IntSet
  coleqs::IntSet
  lb::AbstractVector{T}
  ub::AbstractVector{T}
  rowleqs::IntSet
  rowgeqs::IntSet
  roweqs::IntSet

  function LPHRepresentation(A::AbstractMatrix{T}, l::AbstractVector{T}, u::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T})
    if length(l) != length(u) || size(A, 1) != length(l)
      error("The length of l and u must be equal to the number of rows of A")
    end
    if length(lb) != length(ub) || size(A, 2) != length(lb)
      error("The length of lb and ub must be equal to the number of columns of A")
    end
    if size(A, 2) != N
      error("Type dimension does not match the number of rows of A")
    end
    colleqs = IntSet()
    colgeqs = IntSet()
    coleqs = IntSet()
    for i in 1:N
      if l[i] > typemin(T)
        push!(colgeqs, i)
      end
      if u[i] < typemax(T)
        push!(colleqs, i)
      end
      if l[i] > typemin(T) && u[i] < typemax(T) && myeq(l[i], u[i])
        push!(coleqs, i)
      end
    end
    rowleqs = IntSet()
    rowgeqs = IntSet()
    roweqs = IntSet()
    for i in 1:size(A, 1)
      if lb[i] > typemin(T)
        push!(rowgeqs, i)
      end
      if ub[i] < typemax(T)
        push!(rowleqs, i)
      end
      if lb[i] > typemin(T) && ub[i] < typemax(T) && myeq(lb[i], ub[i])
        push!(roweqs, i)
      end
    end
    new(A, l, u, colleqs, colgeqs, coleqs, lb, ub, rowleqs, rowgeqs, roweqs)
  end
end

LPHRepresentation{T <: Real}(A::AbstractMatrix{T}, l::AbstractVector{T}, u::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T}) = LPHRepresentation{size(A,2),T}(A, l, u, lb, ub)
function LPHRepresentation(A::AbstractMatrix, l::AbstractVector, u::AbstractVector, lb::AbstractVector, ub::AbstractVector)
  T = promote_type(eltype(A), eltype(l), eltype(u), eltype(lb), eltype(ub))
  LPHRepresentation{size(A,2),T}(AbstractMatrix{T}(A), AbstractVector{T}(l), AbstractVector{T}(u), AbstractVector{T}(lb), AbstractVector{T}(ub))
end

function LPHRepresentation{N, T}(it::HRepIterator{N, T})
  A = Matrix{T}(length(it), N)
  lb = Vector{T}(length(it))
  ub = Vector{T}(length(it))
  MathProgBase.warn_no_inf(T)
  l = fill(typemin(T), N)
  u = fill(typemax(T), N)
  for (i, hrep) in enumerate(it)
    A[i,:] = hrep[1]
    ub[i] = hrep[2]
    if hrep[3]
      lb[i] = ub[i]
    else
      lb[i] = typemin(T)
    end
  end
  LPHRepresentation{N, T}(A, l, u, lb, ub)
end
function LPHRepresentation{N, T}(;eqs::Nullable{EqIterator{N, T}}=nothing, ineq::Nullable{IneqIterator{N, T}}=nothing)
  neq = isnull(eqs) ? 0 : length(eqs)
  nineq = isnull(ineqs) ? 0 : length(ineqs)
  nhrep = neq + nineq
  A = Matrix{T}(nhrep, N)
  lb = Vector{T}(nhrep)
  ub = Vector{T}(nhrep)
  MathProgBase.warn_no_inf(T)
  l = fill(typemin(T), N)
  u = fill(typemax(T), N)
  if !(eqs === nothing)
    for (i, eq) in enumerate(get(eqs))
      A[i,:] = eq[1]
      lb[i] = eq[2]
      ub[i] = eq[2]
    end
  end
  if !(ineqs === nothing)
    for (i, ineq) in enumerate(get(ineqs))
      A[neq+i,:] = ineq[1]
      lb[neq+i] = typemin(T)
      ub[neq+i] = ineq[2]
    end
  end
  new(A, l, u, lb, ub)
end

Base.copy{N,T}(lp::LPHRepresentation{N,T}) = LPHRepresentation{N,T}(copy(A), copy(l), copy(u), copy(colleqs), copy(colgeqs), copy(coleqs), copy(lb), copy(ub), copy(rowleqs), copy(rowgeqs), copy(roweqs))

function checknext(lp::LPHRepresentation, colrow, i, lgeq, allowed)
  lgeq += 1
  ok = false
  while colrow <= 2 && !ok
    if colrow == 1
      lgeqs = (lp.colleqs, lp.colgeqs, lp.coleqs)
    else
      lgeqs = (lp.rowleqs, lp.rowgeqs, lp.roweqs)
    end
    while i <= (colrow == 1 ? size(lp.A, 2) : size(lp.A, 1)) && !ok
      while lgeq <= 3 && !ok
        if allowed(i) && i in lgeqs[lgeq]
          ok = true
        else
          lgeq += 1
        end
      end
      if !ok
        i += 1
        lgeq = 1
      end
    end
    if !ok
      colrow += 1
      i = 1
      lgeq = 1
    end
  end
  (colrow, i, lgeq)
end

neqs(lp::LPHRepresentation) = length(lp.coleqs) + length(lp.roweqs)
nineqs(lp::LPHRepresentation) = length(lp.colleqs) + length(lp.colgeqs) + length(lp.rowleqs) + length(lp.rowgeqs)

starthrep(lp::LPHRepresentation) = checknext(1, 0, 3, (i) -> true)
donehrep(lp::LPHRepresentation, state) = state[1] > 2
function nexthrep{N,T}(lp::LPHRepresentation{N,T}, state)
  colrow, i, lgeq = state[1], state[2], state[3]
  if state[1] == 1
    a = spzeros(T, N)
    a[i] = lgeq == 1 ? -one(T) : one(T)
    β = lgeq == 1 ? -lp.lb[i,:] : lp.ub[i,:]
  elseif state[1] == 1
    a = lgeq == 1 ? -lp.A[i,:] : lp.A[i,:]
    β = lgeq == 1 ? -lp.lb[i,:] : lp.ub[i,:]
  else
    error("The iterator is done")
  end
  ((a, β, lgeq == 3), checknext(colrow, i, lgeq, (i) -> true))
end

starteq(lp::LPHRepresentation) = checknext(1, 0, 3, (i) -> i == 3)
doneeq(lp::LPHRepresentation, state) = state[1] > 2
function nexteq{N,T}(lp::LPHRepresentation{N,T}, state)
  colrow, i, lgeq = state[1], state[2], state[3]
  if state[1] == 1
    a = spzeros(T, N)
    a[i] = lgeq == 1 ? -one(T) : one(T)
    β = lgeq == 1 ? -lp.lb[i,:] : lp.ub[i,:]
  elseif state[1] == 1
    a = lgeq == 1 ? -lp.A[i,:] : lp.A[i,:]
    β = lgeq == 1 ? -lp.lb[i,:] : lp.ub[i,:]
  else
    error("The iterator is done")
  end
  ((a, β, lgeq == 3), checknext(colrow, i, lgeq, (i) -> i == 3))
end

startineq(lp::LPHRepresentation) = checknext(1, 0, 3, (i) -> i <= 2)
doneineq(lp::LPHRepresentation, state) = state[1] > 2
function nextineq{N,T}(lp::LPHRepresentation{N,T}, state)
  colrow, i, lgeq = state[1], state[2], state[3]
  if state[1] == 1
    a = spzeros(T, N)
    a[i] = lgeq == 1 ? -one(T) : one(T)
    β = lgeq == 1 ? -lp.lb[i,:] : lp.ub[i,:]
  elseif state[1] == 1
    a = lgeq == 1 ? -lp.A[i,:] : lp.A[i,:]
    β = lgeq == 1 ? -lp.lb[i,:] : lp.ub[i,:]
  else
    error("The iterator is done")
  end
  ((a, β, lgeq == 3), checknext(colrow, i, lgeq, (i) -> i <= 2))
end

# Mandatory
export polyhedron, getinequalities, getgenerators, eliminate, detecthlinearities!, detectvlinearities!, removeredundantinequalities!, removeredundantgenerators!, isredundantinequality, isredundantgenerator, isstronglyredundantinequality, isstronglyredundantgenerator, inequalitiesarecomputed, generatorsarecomputed, loadpolyhedron!

if VERSION < v"0.5-"
  export normalize
  normalize(v,p=2) = v / norm(v,p)
end

polyhedron(repr::Representation) = polyhedron(repr, getlibraryfor(eltype(repr)))
Base.copy(p::Polyhedron)                                                             = error("not implemented")
Base.push!{N}(p::Polyhedron{N}, ine::HRepresentation{N})                             = error("not implemented")
Base.push!{N}(p::Polyhedron{N}, ext::VRepresentation{N})                             = error("not implemented")
inequalitiesarecomputed(p::Polyhedron)                                               = error("not implemented")
getinequalities(p::Polyhedron)                                                       = error("not implemented")
generatorsarecomputed(p::Polyhedron)                                                 = error("not implemented")
getgenerators(p::Polyhedron)                                                         = error("not implemented")
implementseliminationmethod(p::Polyhedron, ::Type{Val{:FourierMotzkin}})             = false
eliminate(p::Polyhedron, delset::IntSet, ::Type{Val{:FourierMotzkin}})               = error("not implemented")
implementseliminationmethod(p::Polyhedron, ::Type{Val{:BlockElimination}})           = false
eliminate(p::Polyhedron, delset::IntSet, ::Type{Val{:BlockElimination}})             = error("not implemented")
detecthlinearities!(p::Polyhedron)                                                   = error("not implemented")
detectvlinearities!(p::Polyhedron)                                                   = error("not implemented")
removeredundantinequalities!(p::Polyhedron)                                          = error("not implemented")
removeredundantgenerators!(p::Polyhedron)                                            = error("not implemented")
isredundantinequality(p::Polyhedron, i::Integer)                                     = error("not implemented")
isredundantgenerator(p::Polyhedron, i::Integer)                                      = error("not implemented")
isstronglyredundantinequality(p::Polyhedron, i::Integer)                             = error("not implemented")
isstronglyredundantgenerator(p::Polyhedron, i::Integer)                              = error("not implemented")
loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Type{Val{:ine}}) = error("not implemented")
loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Type{Val{:ext}}) = error("not implemented")

# These can optionally be reimplemented for speed by a library
export numberofinequalities, numberofgenerators, dim, affinehull, getredundantinequalities, getstronglyredundantinequalities, getredundantgenerators, getstronglyredundantgenerators, transforminequalities, transformgenerators, project, radialprojectoncut

loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Symbol) = loadpolyhedron!(p, filename, Val{extension})

function loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::AbstractString)
  s = findfirst(["ext", "ine"], filename)
  if s == 0
    error("Invalid extension $extension, please give 'ext' for V-representation or 'ine' for H-representation")
  end
  loadpolyhedron!(p, filename, [:ext, :ine][s])
end

eliminate(p::Polyhedron, method::Symbol) = eliminate(p, Val{method})
eliminate(p::Polyhedron, delset::IntSet, method::Symbol) = eliminate(p, delset::IntSet, Val{method})

eliminate{N}(p::Polyhedron{N}, method::Type{Val{:ProjectGenerators}}) = eliminate(p, IntSet(N), method)

function eliminate{N}(p::Polyhedron{N}, delset::IntSet=IntSet(N))
  fm = implementseliminationmethod(p, Val{:FourierMotzkin})
  be = implementseliminationmethod(p, Val{:BlockElimination})
  if (!fm && !be) || generatorsarecomputed(p)
    method = :ProjectGenerators
  elseif fm && (!be || delset == IntSet(N))
    method = :FourierMotzkin
  else
    method = :BlockElimination
  end
  eliminate(p, delset, Val{method})
end

function eliminate{N}(p::Polyhedron{N}, delset::IntSet, ::Type{Val{:ProjectGenerators}})
  ext = getgenerators(p)
  I = eye(Int, N)
  polyhedron(I[setdiff(IntSet(1:N), collect(delset)),:] * ext, getlibrary(p))
end

function call{N, S, T}(::Type{Polyhedron{N, S}}, p::Polyhedron{N, T})
  if !inequalitiesarecomputed(p) && generatorsarecomputed(p)
    repr = VRepresentation{N,S}(getgenerators(p))
  else
    repr = HRepresentation{N,S}(getinequalities(p))
  end
  polyhedron(repr, getlibraryfor(p, S))
end

function numberofinequalities(p::Polyhedron)
  size(getinequalities(p).A, 1)
end

function numberofgenerators(p::Polyhedron)
  ext = getgenerators(p)
  size(ext.V, 1) + size(ext.R, 1)
end

# eliminate the last dimension by default
eliminate{N,T}(p::Polyhedron{N,T})  = eliminate(p::Polyhedron, IntSet([N]))

function transformgenerators{N}(p::Polyhedron{N}, P::AbstractMatrix)
  # Each generator x is transformed to P * x
  # If P is orthogonal, the new axis are the rows of P.
  if size(P, 2) != N
    error("The number of columns of P must match the dimension of the polyhedron")
  end
  ext = P * getgenerators(p)
  polyhedron(ext, getlibraryfor(p, eltype(ext)))
end

function transforminequalities(p::Polyhedron, P::AbstractMatrix)
  # The new axis are the column of P.
  # Let y be the coordinates of a point x in these new axis.
  # We have x = P * y so y = P \ x.
  # We have
  # b = Ax = A * P * (P \ x) = (A * P) * y
  ine = getinequalities(p) * P
  polyhedron(ine, getlibraryfor(p, eltype(ine)))
end

function project{N,T}(p::Polyhedron{N,T}, P::Array)
  # Function to make x orthogonal to an orthonormal basis in Q
  # We first make the columns of P orthonormal
  if size(P, 1) != N
    error("The columns of P should have the same dimension than the polyhedron")
  end
  m = size(P, 2)
  if m > N
    error("P should have more columns than rows")
  end
  Q = Array{Float64}(P) # normalize will make it nonrational
  Proj = zeros(eltype(P), N, N)
  for i = 1:m
    Q[:,i] = normalize(Q[:,i] - Proj * Q[:,i])
    Proj += Q[:,i] * Q[:,i]'
  end
  if m == N
    basis = Q
  else
    # For the rest, we take the canonical basis and we look at
    # I - Proj * I
    I = eye(Float64, N)
    R = I - Proj
    # We take the n-m that have highest norm
    order = sortperm([dot(R[:,i], R[:,i]) for i in 1:N])
    R = I[:,order[m+1:N]]
    for i in 1:N-m
      R[:,i] = normalize(R[:,i] - Proj * R[:,i])
      Proj += R[:,i] * R[:,i]'
    end
    basis = [Q R]
  end
  eliminate(transforminequalities(p, basis), IntSet(m+1:N))
end

function radialprojectoncut{N}(p::Polyhedron{N}, cut::Vector, at)
  if myeqzero(at)
    error("at is zero")
  end
  if length(cut) != N
    error("The dimensions of the cut and of the polyhedron do not match")
  end
  ext = SimpleVRepresentation(getgenerators(p))
  V = copy(ext.V)
  R = copy(ext.R)
  for i in 1:size(V, 1)
    v = vec(ext.V[i,:])
    if !myeq(dot(cut, v), at)
      error("The nonhomogeneous part should be in the cut")
    end
  end
  for i in 1:size(R, 1)
    v = vec(ext.R[i,:])
    if myeqzero(v)
      # It can happen since I do not necessarily have removed redundancy
      v = zeros(eltype(v), length(v))
    elseif !myeq(dot(cut, v), at)
      if myeqzero(dot(cut, v))
        error("A ray is parallel to the cut") # FIXME is ok if some vertices are on the cut ? (i.e. at == 0, cut is not needed)
      end
      v = v * at / dot(cut, v)
    end
    R[i,:] = v
  end
  # no more rays nor linearity since at != 0
  ext2 = SimpleVRepresentation([V; R])
  polyhedron(ext2, getlibraryfor(p, eltype(ext2)))
end

function fulldim{N,T}(p::Polyhedron{N,T})
  N
end

function dim(p::Polyhedron)
  detecthlinearities!(p)
  ine = getinequalities(p)
  d = size(ine.A, 2) - length(ine.linset)
end

function affinehull(p::Polyhedron)
  detecthlinearities!(p)
  typeof(p)(affinehull(getinequalities(p)))
end

function getredundantinequalities(p::Polyhedron)
  red = IntSet([])
  for i in 1:numberofinequalities(p)
    if isredundantinequality(p, i)[1]
      push!(red, i)
    end
  end
  red
end
function getstronglyredundantinequalities(p::Polyhedron)
  red = IntSet([])
  for i in 1:numberofinequalities(p)
    if isstronglyredundantinequality(p, i)[1]
      push!(red, i)
    end
  end
  red
end

function getredundantgenerators(p::Polyhedron)
  red = IntSet([])
  for i in 1:numberofgenerators(p)
    if isredundantgenerator(p, i)[1]
      push!(red, i)
    end
  end
  red
end
function getstronglyredundantgenerators(p::Polyhedron)
  red = IntSet([])
  for i in 1:numberofgenerators(p)
    if isstronglyredundantgenerator(p, i)[1]
      push!(red, i)
    end
  end
  red
end
function isredundantgenerator(p::Polyhedron, x::Vector, vertex::Bool)
  ine = getinequalities(p)
  for i in 1:size(ine.A, 1)
    if i in ine.linset
      if !myeq(dot(ine.A[i, :], x), vertex ? ine.b[i] : 0)
        return (false, Nullable{Vector{eltype(ine)}}(ine.A[i,:]))
      end
    else
      if mygt(dot(ine.A[i, :], x), vertex ? ine.b[i] : 0)
        return (false, Nullable{Vector{eltype(ine)}}(ine.A[i,:]))
      end
    end
  end
  return (true, Nullable{Vector{eltype(ine)}}(nothing))
end

function Base.intersect{N}(p::Polyhedron{N}, ine::HRepresentation{N})
  inter = Base.intersect(getinequalities(p), ine)
  polyhedron(inter, getlibraryfor(p, eltype(inter)))
end
Base.intersect{N}(ine::HRepresentation{N}, p::Polyhedron{N}) = Base.intersect(p, ine)
Base.intersect{N}(p1::Polyhedron{N}, p2::Polyhedron{N}) = Base.intersect(p1, getinequalities(p2))

function (+){N}(p::Polyhedron{N}, ext::VRepresentation{N})
  sum = getgenerators(p) + ine
  polyhedron(sum, getlibraryfor(p, eltype(sum)))
end
(+){N}(ext::VRepresentation{N}, p::Polyhedron{N}) = Base.intersect(p, ext)
(+){N}(p1::Polyhedron{N}, p2::Polyhedron{N}) = Base.intersect(p1, getgenerators(p2))

function (*)(p1::Polyhedron, p2::Polyhedron)
  # iac1 = inequalitiesarecomputed(p1) ? 1 : 0
  # iac2 = inequalitiesarecomputed(p2) ? 1 : 0
  # gac1 = genratorsarecomputed(p1) ? 1 : 0
  # gac2 = genratorsarecomputed(p2) ? 1 : 0
  # if iac1 + iac2 >= gac1 + gac2
  repr = getinequalities(p1) * getinequalities(p2)
  # else
  #   repr = getgenerators(p1) * getgenerators(p2)
  # end
  polyhedron(repr, getlibraryfor(p1, eltype(repr)))
end

using MathProgBase

import MathProgBase.LinearQuadraticModel, MathProgBase.loadproblem!, MathProgBase.optimize!, MathProgBase.status, MathProgBase.getobjval, MathProgBase.getsolution, MathProgBase.getunboundedray, MathProgBase.linprog

export LPPolyhedron, SimpleLPPolyhedron

abstract LPPolyhedron{N, T}

type SimpleLPPolyhedron{N, T} <: LPPolyhedron{N, T}
  ext::SimpleVRepresentation{N, T}
  obj::Vector{T}
  sense::Symbol

  objval::Nullable{T}
  solution::Nullable{Vector{T}}
  status::Symbol
end

function LinearQuadraticModel{N, T}(p::Polyhedron{N, T})
  SimpleLPPolyhedron{N, T}(getgenerators(p), zeros(T, N), :None, nothing, nothing, :Undecided)
end
function loadproblem!(lpm::SimpleLPPolyhedron, obj, sense)
  if !(sense in [:Max, :Min])
    error("sense should be :Max or :Min")
  end
  if sum(abs(obj)) != 0
    lpm.obj = copy(obj)
    lpm.sense = sense
  end
end
function optimize!{N, T}(lpm::SimpleLPPolyhedron{N, T})
  if size(lpm.ext.V, 1) + size(lpm.ext.R, 1) == 0
    lpm.status = :Infeasible
  else
    if lpm.sense in [:Max, :Min]
      if lpm.sense == :Max
        better(a, b) = a > b
        mybetter(a, b) = mygt(a, b)
      else
        better(a, b) = a < b
        mybetter(a, b) = mylt(a, b)
      end
      lpm.status = :Infeasible
      for i in 1:size(lpm.ext.R, 1)
        v = vec(lpm.ext.R[i,:])
        objval = dot(lpm.obj, v)
        if lpm.status != :Unbounded && mybetter(objval, zero(T))
          lpm.status = :Unbounded
          lpm.objval = lpm.sense == :Max ? typemax(T) : typemin(T)
          lpm.solution = v
        end
      end
      if status != :Unbounded
        for i in 1:size(lpm.ext.V, 1)
          v = vec(lpm.ext.V[i,:])
          objval = dot(lpm.obj, v)
          if lpm.status == :Undecided || (lpm.status == :Optimal && better(objval, get(lpm.objval)))
            lpm.status = :Optimal
            lpm.objval = objval
            lpm.solution = v
          end
        end
      end
    else
      lpm.status = :Optimal
      lpm.objval = zero(T)
      lpm.solution = zeros(T,N)
    end
  end
end

function MathProgBase.status(lpm::SimpleLPPolyhedron)
  lpm.status
end
function MathProgBase.getobjval(lpm::SimpleLPPolyhedron)
  get(lpm.objval)
end
function MathProgBase.getsolution(lpm::SimpleLPPolyhedron)
  copy(get(lpm.solution))
end
function MathProgBase.getunboundedray(lpm::SimpleLPPolyhedron)
  copy(get(lpm.solution))
end

type LinprogSolution
    status
    objval
    sol
    attrs
end

function MathProgBase.linprog{N}(c::Vector, p::Polyhedron{N})
  if N != length(c)
    println("length of objective does not match dimension of polyhedron")
  end
  m = LinearQuadraticModel(p)
  loadproblem!(m, c, :Min)
  optimize!(m)
  stat = status(m)
  if stat == :Optimal
    return LinprogSolution(stat, getobjval(m), getsolution(m), Dict())
  elseif stat == :Unbounded
    attrs = Dict()
    attrs[:unboundedray] = getunboundedray(m)
    return LinprogSolution(stat, nothing, [], attrs)
  else
    return LinprogSolution(stat, nothing, [], Dict())
  end
end

function Base.isempty{N,T}(p::Polyhedron{N,T})
  linprog(zeros(T, N), p).status == :Infeasible
end

function isredundantinequalityaux(p::Polyhedron, a::Vector, b)
  sol = linprog(-a, p)
  if sol.status == :Unbounded
    (false, sol.attrs[:unboundedray], false)
  elseif sol.status == :Optimal && mygt(sol.objval, b)
    (false, sol.sol, true)
  else
    (true, nothing, false)
  end
end
function isredundantinequality(p::Polyhedron, a::Vector, b, eq::Bool)
  if eq
    sol = isredundantinequalityaux(p, a, b)
    if !sol[1]
      sol
    else
      isredundantinequalityaux(p, -a, -b)
    end
  else
    isredundantinequalityaux(p, a, b)
  end
end

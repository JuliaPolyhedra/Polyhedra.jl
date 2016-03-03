# Mandatory
import Base.copy, Base.push!
export polyhedron, getlibrary, getinequalities, getgenerators, eliminate!, detectlinearities!, removeredundantinequalities!, removeredundantgenerators!, isredundantinequality, isredundantgenerator, isstronglyredundantinequality, isstronglyredundantgenerator

polyhedron(desc::Description)                            = error("please specify a library")
getlibrary(p::Polyhedron)                                = error("not implemented")
Base.copy(p::Polyhedron)                                 = error("not implemented")
Base.push!(p::Polyhedron, ine::InequalityDescription)    = error("not implemented")
Base.push!(p::Polyhedron, ext::GeneratorDescription)     = error("not implemented")
inequalitiesarecomputed(p::Polyhedron)                   = error("not implemented")
getinequalities(p::Polyhedron)                           = error("not implemented")
generatorsarecomputed(p::Polyhedron)                     = error("not implemented")
getgenerators(p::Polyhedron)                             = error("not implemented")
eliminate!(p::Polyhedron, delset::IntSet)                = error("not implemented")
detectlinearities!(p::Polyhedron)                        = error("not implemented")
removeredundantinequalities!(p::Polyhedron)              = error("not implemented")
removeredundantgenerators!(p::Polyhedron)                = error("not implemented")
isredundantinequality(p::Polyhedron, i::Integer)         = error("not implemented")
isredundantgenerator(p::Polyhedron, i::Integer)          = error("not implemented")
isstronglyredundantinequality(p::Polyhedron, i::Integer) = error("not implemented")
isstronglyredundantgenerator(p::Polyhedron, i::Integer)  = error("not implemented")

# These can optionally be reimplemented for speed by a library
import Base.isempty
export numberofinequalities, numberofgenerators, eliminate, fulldim, dim, affinehull, getredundantinequalities, getstronglyredundantinequalities, getredundantgenerators, getstronglyredundantgenerators, transform, project

function numberofinequalities(p::Polyhedron)
  size(getinequalities(p).A, 1)
end

function numberofgenerators(p::Polyhedron)
  size(getinequalities(p).A, 1)
end

Base.isempty(p::Polyhedron) = numberofgenerators(p) == 0

function eliminate(p::Polyhedron, delset::IntSet)
  pcopy = copy(p)
  eliminate!(pcopy, delset)
  pcopy
end

# eliminate the last dimension by default
eliminate!(p::Polyhedron) = eliminate!(p::Polyhedron, IntSet([fulldim(p)]))
eliminate(p::Polyhedron)  = eliminate(p::Polyhedron, IntSet([fulldim(p)]))

function transform(p::Polyhedron, P::Array)
  # The new axis are the column of P.
  # Let y be the coordinates of a point x in these new axis.
  # We have x = P * y so y = P \ x.
  # We have
  # b = Ax = A * P * (P \ x) = (A * P) * y
  ine = getinequalities(p)
  newine = InequalityDescription(ine.A * P, ine.b, ine.linset)
  polyhedron(newine, getlibrary(p))
end

function project(p::Polyhedron, P::Array)
  # Function to make x orthogonal to an orthonormal basis in Q
  # We first make the columns of P orthonormal
  n = size(P, 1)
  m = size(P, 2)
  if m > n
    error("P should have more columns than rows")
  end
  Q = Array{Float64}(P) # normalize will make it nonrational
  Proj = zeros(eltype(P), n, n)
  for i = 1:m
    Q[:,i] = normalize(Q[:,i] - Proj * Q[:,i])
    Proj += Q[:,i] * Q[:,i]'
  end
  if m == n
    basis = Q
  else
    # For the rest, we take the canonical basis and we look at
    # I - Proj * I
    I = eye(Float64, n)
    R = I - Proj
    # We take the n-m that have highest norm
    order = sortperm([dot(R[:,i], R[:,i]) for i in 1:n])
    R = I[:,order[m+1:n]]
    for i in 1:n-m
      R[:,i] = normalize(R[:,i] - Proj * R[:,i])
      Proj += R[:,i] * R[:,i]'
    end
    basis = [Q R]
  end
  p2 = transform(p, basis)
  eliminate!(p2, IntSet(m+1:n))
  p2
end

function fulldim(p::Polyhedron)
  size(getinequalities(p).A, 2)
end

function dim(p::Polyhedron)
  detectlinearities!(p)
  ine = getinequalities(p)
  d = size(ine.A, 2) - length(ine.linset)
end

function affinehull(p::Polyhedron)
  detectlinearities!(p)
  ine = getinequalities(p)
  neqs = length(ine.linset)
  A = Matrix{eltype(ine)}(neqs, size(ine.A, 2))
  b = Vector{eltype(ine)}(neqs)
  cur = 1
  for i in 1:size(ine.A, 1)
    if i in ine.linset
      A[cur, :] = ine.A[i, :]
      b[cur] = ine.b[i]
      cur += 1
    end
  end
  typeof(p)(InequalityDescription(A, b, IntSet(1:neqs)))
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

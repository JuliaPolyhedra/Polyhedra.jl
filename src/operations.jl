# Mandatory
import Base.copy, Base.push!
export polyhedron, getinequalities, getgenerators, eliminate, detectlinearities!, removeredundantinequalities!, removeredundantgenerators!, isredundantinequality, isredundantgenerator, isstronglyredundantinequality, isstronglyredundantgenerator

if VERSION < v"0.5-"
  export normalize
  normalize(v,p=2) = v / norm(v,p)
end

polyhedron(desc::Representation) = polyhedron(desc, getlibraryfor(eltype(desc)))
Base.copy(p::Polyhedron)                                 = error("not implemented")
Base.push!(p::Polyhedron, ine::HRepresentation)          = error("not implemented")
Base.push!(p::Polyhedron, ext::VRepresentation)          = error("not implemented")
inequalitiesarecomputed(p::Polyhedron)                   = error("not implemented")
getinequalities(p::Polyhedron)                           = error("not implemented")
generatorsarecomputed(p::Polyhedron)                     = error("not implemented")
getgenerators(p::Polyhedron)                             = error("not implemented")
eliminate(p::Polyhedron, delset::IntSet)                 = error("not implemented")
detectlinearities!(p::Polyhedron)                        = error("not implemented")
removeredundantinequalities!(p::Polyhedron)              = error("not implemented")
removeredundantgenerators!(p::Polyhedron)                = error("not implemented")
isredundantinequality(p::Polyhedron, i::Integer)         = error("not implemented")
isredundantgenerator(p::Polyhedron, i::Integer)          = error("not implemented")
isstronglyredundantinequality(p::Polyhedron, i::Integer) = error("not implemented")
isstronglyredundantgenerator(p::Polyhedron, i::Integer)  = error("not implemented")

# These can optionally be reimplemented for speed by a library
import Base.isempty
export numberofinequalities, numberofgenerators, dim, affinehull, getredundantinequalities, getstronglyredundantinequalities, getredundantgenerators, getstronglyredundantgenerators, transforminequalities, transformgenerators, project, radialprojectoncut

function numberofinequalities(p::Polyhedron)
  size(getinequalities(p).A, 1)
end

function numberofgenerators(p::Polyhedron)
  size(getinequalities(p).A, 1)
end

Base.isempty(p::Polyhedron) = numberofgenerators(p) == 0

# eliminate the last dimension by default
eliminate{N,T}(p::Polyhedron{N,T})  = eliminate(p::Polyhedron, IntSet([N]))

function transformgenerators(p::Polyhedron, P::Matrix)
  # Each generator x is transformed to P * x
  # If P is orthogonal, the new axis are the rows of P.
  ext = getgenerators(p)
  newext = VRepresentation(ext.V * P', ext.R * P', ext.vertex, ext.Vlinset, ext.Rlinset)
  polyhedron(newext, getlibraryfor(p, eltype(P)))
end

function transforminequalities(p::Polyhedron, P::Matrix)
  # The new axis are the column of P.
  # Let y be the coordinates of a point x in these new axis.
  # We have x = P * y so y = P \ x.
  # We have
  # b = Ax = A * P * (P \ x) = (A * P) * y
  ine = getinequalities(p)
  newine = HRepresentation(ine.A * P, ine.b, ine.linset)
  polyhedron(newine, getlibraryfor(p, eltype(P)))
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
  ext = getgenerators(p)
  V = copy(ext.V)
  R = copy(ext.R)
  for i in 1:(size(V, 1)+size(R, 1))
    v = vec((i <= size(V, 1)) ? ext.V[i,:] : ext.R[i-size(V, 1),:])
    if i in ext.vertex
      if !myeq(dot(cut, v), at)
        error("The nonhomogeneous part should be in the cut")
      end
    else
      if myeqzero(v)
        # It can happen since I do not have remove redundancy
        v = zeros(eltype(v), length(v))
      elseif !myeq(dot(cut, v), at)
        if myeqzero(dot(cut, v))
          error("A ray is parallel to the cut")
        end
        v = v * at / dot(cut, v)
      end
    end
    if i <= size(V, 1)
      V[i,:] = v
    else
      R[i-size(V, 1),:] = v
    end
  end
  # no more rays nor linearity since at != 0
  ext2 = VRepresentation(V, R)
  polyhedron(ext2, getlibraryfor(p, eltype(ext2)))
end

function fulldim{N,T}(p::Polyhedron{N,T})
  N
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
  typeof(p)(HRepresentation(A, b, IntSet(1:neqs)))
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

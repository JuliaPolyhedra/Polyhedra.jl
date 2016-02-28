# Mandatory
import Base.copy, Base.push!
export getinequalities, getgenerators, eliminate!, detectlinearities!, removeredundantinequalities!, removeredundantgenerators!, isredundantinequality, isredundantgenerator, isstronglyredundantinequality, isstronglyredundantgenerator

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
export numberofinequalities, numberofgenerators, eliminate, fulldim, dim, affinehull, getredundantinequalities, getstronglyredundantinequalities, getredundantgenerators, getstronglyredundantgenerators

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

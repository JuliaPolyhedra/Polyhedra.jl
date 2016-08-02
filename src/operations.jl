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
hrepiscomputed(p::Polyhedron)                                                        = error("not implemented")
gethrep(p::Polyhedron)                                                               = error("not implemented")
vrepiscomputed(p::Polyhedron)                                                        = error("not implemented")
getvrep(p::Polyhedron)                                                               = error("not implemented")
implementseliminationmethod(p::Polyhedron, ::Type{Val{:FourierMotzkin}})             = false
eliminate(p::Polyhedron, delset::IntSet, ::Type{Val{:FourierMotzkin}})               = error("not implemented")
implementseliminationmethod(p::Polyhedron, ::Type{Val{:BlockElimination}})           = false
eliminate(p::Polyhedron, delset::IntSet, ::Type{Val{:BlockElimination}})             = error("not implemented")
#loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Type{Val{:ine}}) = error("not implemented")
#loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Type{Val{:ext}}) = error("not implemented") # FIXME ExtFileVRepresentation or just ExtFile

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
    if decomposedvfast(p)
      f2 = (i,x) -> (AbstractVector{eltype(Tout)}(x[1]), x[2])
      polyhedron(PointIterator(p, f2), RayIterator(p, f2), getlibraryfor(p, S))
    else
      f3 = (i,x) -> (AbstractVector{eltype(Tout)}(x[1]), x[2], x[3])
      polyhedron(VRepIterator(p, f3), getlibraryfor(p, S))
    end
  else
    if decomposedvfast(p)
      f2 = (i,x) -> (AbstractVector{eltype(Tout)}(x[1]), S(x[2]))
      polyhedron(IneqIterator(p, f2), EqIterator(p, f2), getlibraryfor(p, S))
    else
      f3 = (i,x) -> (AbstractVector{eltype(Tout)}(x[1]), S(x[2]), x[3])
      polyhedron(HRepIterator(p, f3), getlibraryfor(p, S))
    end
  end
end

# eliminate the last dimension by default
eliminate{N,T}(p::Polyhedron{N,T})  = eliminate(p::Polyhedron, IntSet([N]))

# function transformgenerators{N}(p::Polyhedron{N}, P::AbstractMatrix)
#   # Each generator x is transformed to P * x
#   # If P is orthogonal, the new axis are the rows of P.
#   if size(P, 2) != N
#     error("The number of columns of P must match the dimension of the polyhedron")
#   end
#   ext = P * getgenerators(p)
#   polyhedron(ext, getlibraryfor(p, eltype(ext)))
# end
#
# function transforminequalities(p::Polyhedron, P::AbstractMatrix)
#   # The new axis are the column of P.
#   # Let y be the coordinates of a point x in these new axis.
#   # We have x = P * y so y = P \ x.
#   # We have
#   # b = Ax = A * P * (P \ x) = (A * P) * y
#   ine = getinequalities(p) * P
#   polyhedron(ine, getlibraryfor(p, eltype(ine)))
# end

# function (*){N,S}(A::AbstractMatrix{S}, p::Polyhedron{N})
#   if size(A, 2) != N
#     error("Incompatible dimension")
#   end
#   if generatorsarecomputed(p)
#     transformgenerators(p, A)
#   else # FIXME not wokring
#     ine = SimpleHRepresentation(getinequalities(p))
#     nnew = size(A, 1)
#     if false
#       # CDD works with delset not at the end ?
#       newA = [ine.A spzeros(S, size(ine.A, 1), nnew);
#                   A  -speye(S, nnew, nnew)]
#       delset = IntSet(nnew+(1:N))
#     else
#       newA = [spzeros(S, size(ine.A, 1), nnew) ine.A;
#                -speye(S, nnew, nnew) A]
#       delset = IntSet(1:N)
#     end
#     newb = [ine.b; spzeros(S, nnew)]
#     newlinset = ine.linset ∪ IntSet(N+(1:nnew))
#     newine = SimpleHRepresentation(newA, newb, newlinset)
#     newpoly = polyhedron(newine, getlibraryfor(p, eltype(newine)))
#     eliminate(newpoly, IntSet(nnew+(1:N)))
#   end
# end

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
  eliminate(p * basis, IntSet(m+1:N))
end

# TODO rewrite, it is just cutting a cone with a half-space, nothing more
# function radialprojectoncut{N}(p::Polyhedron{N}, cut::Vector, at)
#   if myeqzero(at)
#     error("at is zero")
#   end
#   if length(cut) != N
#     error("The dimensions of the cut and of the polyhedron do not match")
#   end
#   ext = SimpleVRepresentation(getgenerators(p))
#   V = copy(ext.V)
#   R = copy(ext.R)
#   for i in 1:size(V, 1)
#     v = vec(ext.V[i,:])
#     if !myeq(dot(cut, v), at)
#       error("The nonhomogeneous part should be in the cut")
#     end
#   end
#   for i in 1:size(R, 1)
#     v = vec(ext.R[i,:])
#     if myeqzero(v)
#       # It can happen since I do not necessarily have removed redundancy
#       v = zeros(eltype(v), length(v))
#     elseif !myeq(dot(cut, v), at)
#       if myeqzero(dot(cut, v))
#         error("A ray is parallel to the cut") # FIXME is ok if some vertices are on the cut ? (i.e. at == 0, cut is not needed)
#       end
#       v = v * at / dot(cut, v)
#     end
#     R[i,:] = v
#   end
#   # no more rays nor linearity since at != 0
#   ext2 = SimpleVRepresentation([V; R])
#   polyhedron(ext2, getlibraryfor(p, eltype(ext2)))
# end

#function fulldim{N,T}(p::Polyhedron{N,T})
#  N
#end

function dim(p::Polyhedron)
  detecthlinearities!(p)
  fulldim(p) - neqs(p)
end

# function affinehull(p::Polyhedron)
#   detecthlinearities!(p)
#   typeof(p)(affinehull(getinequalities(p)))
# end

function isredundantgenerator(p::Polyhedron, x::Vector, vertex::Bool, cert=false)
  for (a, β) in eqs(p)
    if !myeq(dot(a, x), vertex ? β : 0)
      return cert ? (false, Nullable{Vector{eltype(p)}}(a)) : false
    end
  end
  for (a, β) in ineqs(p)
    if mygt(dot(a, x), vertex ? β : 0)
      return cert ? (false, Nullable{Vector{eltype(p)}}(a)) : false
    end
  end
  cert ? (true, Nullable{Vector{eltype(ine)}}(nothing)) : true
end

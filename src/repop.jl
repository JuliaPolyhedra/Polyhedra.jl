function affinehull{T<:HRep}(p::T)
  T(eqs(p))
end

# Always type of first arg
function Base.intersect{T1<:HRep, T2<:HRep}(p1::T1, p2::T2)
  if eltype(T1) != eltype(T2)
    error("Cannot take the intersection of polyhedra of different element type")
  end
  if fulldim(T1) != fulldim(T2)
    error("Cannot take the intersection of polyhedra of different dimension")
  end
  T1(HRepIterator([p1, p2]))
end

# Always type of first arg
function (+){T1<:VRep, T2<:VRep}(p1::T1, p2::T2)
  if eltype(T1) != eltype(T2)
    error("Cannot take the intersection of polyhedra of different element type")
  end
  if fulldim(T1) != fulldim(T2)
    error("Cannot take the intersection of polyhedra of different dimension")
  end
  T1(VRepIterator([p1, p2]))
end

# p1 has priority
function usehrep(p1::Polyhedron, p2::Polyhedron)
  hrepiscomputed(p1) && (!vrepiscomputed(p1) || hrepiscomputed(p2))
end

# Always type of first arg
@generated function (*){T1<:Rep, T2<:Rep}(p1::T1, p2::T2)
  if eltype(T1) != eltype(T2)
    error("Cannot take the cartesian product between polyhedra of different element type")
  end
  T = eltype(T1)
  N1 = fulldim(T1)
  N2 = fulldim(T2)
  hashrep = T1 <: HRepresentation || T2 <: HRepresentation
  hasvrep = T1 <: VRepresentation || T2 <: VRepresentation
  Tout = changefulldim(T1, fulldim(T1)+fulldim(T2))
  if hashrep && hasvrep
    error("Cannot take the cartesian product between a H-Representation and a V-Representation")
  elseif hashrep || (!hasvrep && usehrep(p1, p2))
    # TODO fastdecompose
    :(Tout(HRepIterator([p1, p2], (i, (a, β, lin)) -> (i == 1 ? [a; spzeros(T, N2)] : [spzeros(T, N1); a], β, lin))))
  else
    # TODO fastdecompose
    :(Tout(VRepIterator([p1, p2], (i, (x, point, lin)) -> (i == 1 ? [x; spzeros(T, N2)] : [spzeros(T, N1); x], point, lin))))
  end
end

function (*){RepT<:HRep}(hrep::RepT, P::AbstractMatrix)
  if size(P, 1) != fulldim(T)
    error("The number of rows of P must match the dimension of the H-representation")
  end
  Tout = mypromote_type(eltype(RepT), eltype(P))
  RepTout = changeboth(RepT, size(P, 2), Tout)
  if decomposedhfast(hrep)
    f = (i,x) -> (Vector{Tout}(x[1] * P), Tout(x[2]))
    eqs = EqIterator(hrep, f)
    ineqs = IneqIterator(hrep, f)
    if RepT <: HRepresentation
      Tout(ineqs=ineqs, eqs=eqs)
    else
      polyhedron(ineqs, eqs, getlibraryfor(p, Tout))
    end
  else
    f = (i,x) -> (Vector{Tout}(x[1] * P), Tout(x[2]), x[3])
    hreps = HRepIterator(hrep, f)
    if RepT <: HRepresentation
      Tout(hreps)
    else
      polyhedron(hreps, getlibraryfor(p, Tout))
    end
  end
end
function (*){RepT<:VRep}(P::AbstractMatrix, vrep::RepT)
  if size(P, 2) != fulldim(T)
    error("The number of rows of P must match the dimension of the H-representation")
  end
  Tout = mypromote_type(eltype(RepT), eltype(P))
  RepTout = changeboth(RepT, size(P, 1), Tout)
  if decomposedvfast(vrep)
    f = (i,x) -> (Vector{Tout}(P * x[1]), x[2])
    points = PointIterator(vrep, f)
    rays = RayIterator(vrep, f)
    if RepT <: VRepresentation
      Tout(points=points, rays=rays)
    else
      polyhedron(points, rays, getlibraryfor(p, Tout))
    end
  else
    f = (i,x) -> (Vector{Tout}(P * x[1]), x[2], x[3])
    vreps = VRepIterator(vrep, f)
    if RepT <: VRepresentation
      Tout(vreps)
    else
      polyhedron(vreps, getlibraryfor(p, Tout))
    end
  end
end

function Base.round{N,T<:AbstractFloat}(rep::HRepresentation{N,T})
  f2round = (i,x) -> (round(x[1]), round(x[2]))
  f3round = (i,x) -> (round(x[1]), round(x[2]), x[3])
  if decomposedfast(rep)
    typeof(rep)(eqs=EqIterator(rep, f2round), ineqs=IneqIterator(rep, f2round))
  else
    typeof(rep)(HRepIterator(rep, f3round))
  end
end
function Base.round{N,T<:AbstractFloat}(rep::VRepresentation{N,T})
  f2round = (i,x) -> (round(x[1]), x[2])
  f3round = (i,x) -> (round(x[1]), x[2], x[3])
  if decomposedfast(rep)
    typeof(rep)(eqs=EqIterator(rep, f2round), ineqs=IneqIterator(rep, f2round))
  else
    typeof(rep)(VRepIterator(rep, f3round))
  end
end

function gethredundantindices(hrep::HRep)
  red = IntSet([])
  for i in 1:nhreps(hrep)
    if ishredundant(hrep, i)
      push!(red, i)
    end
  end
  red
end
function gethstronglyredundantindices(hrep::HRep)
  red = IntSet([])
  for i in 1:nhreps(hrep)
    if ishstronglyredundant(p, i)[1]
      push!(red, i)
    end
  end
  red
end
function getvredundantindices(vrep::VRep)
  red = IntSet([])
  for i in 1:nvreps(vrep)
    if isvredundant(p, i)
      push!(red, i)
    end
  end
  red
end
function getvstronglyredundantindices(hrep::HRep)
  red = IntSet([])
  for i in 1:nvreps(p)
    if isvstronglyredundant(p, i)
      push!(red, i)
    end
  end
  red
end

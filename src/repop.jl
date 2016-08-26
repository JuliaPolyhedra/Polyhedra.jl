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
  f = (i, x) -> zeropad(x, i == 1 ? N2 : -N1)
  if hashrep && hasvrep
    error("Cannot take the cartesian product between a H-Representation and a V-Representation")
  elseif hashrep || (!hasvrep && usehrep(p1, p2))
    # TODO fastdecompose
    :(Tout(HRepIterator([p1, p2], f)))
  else
    # TODO fastdecompose
    :(Tout(VRepIterator([p1, p2], f)))
  end
end

function (*){RepT<:HRep}(hrep::RepT, P::AbstractMatrix)
  if size(P, 1) != fulldim(T)
    error("The number of rows of P must match the dimension of the H-representation")
  end
  Tout = mypromote_type(eltype(RepT), eltype(P))
  RepTout = changeboth(RepT, size(P, 2), Tout)
  f = (i, h) -> h * P
  if decomposedhfast(hrep)
    eqs = eqs(hrep, f)
    ineqs = ineqs(hrep, f)
    if RepT <: HRepresentation
      Tout(ineqs=ineqs, eqs=eqs)
    else
      polyhedron(ineqs, eqs, getlibraryfor(p, Tout))
    end
  else
    hreps = hrep(hrep, f)
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
  f = (i, v) -> P * v
  if decomposedvfast(vrep)
    points = points(vrep, f)
    rays = rays(vrep, f)
    if RepT <: VRepresentation
      Tout(points=points, rays=rays)
    else
      polyhedron(points, rays, getlibraryfor(p, Tout))
    end
  else
    vreps = vrep(vrep, f)
    if RepT <: VRepresentation
      Tout(vreps)
    else
      polyhedron(vreps, getlibraryfor(p, Tout))
    end
  end
end

function Base.round{N,T<:AbstractFloat}(rep::HRepresentation{N,T})
  f = (i, h) -> round(h)
  if decomposedfast(rep)
    typeof(rep)(eqs=eqs(rep, f), ineqs=ineqs(rep, f))
  else
    typeof(rep)(hrep(rep, f))
  end
end
function Base.round{N,T<:AbstractFloat}(rep::VRepresentation{N,T})
  f = (i, v) -> round(v)
  if decomposedfast(rep)
    typeof(rep)(eqs=eqs(rep, f), ineqs=ineqs(rep, f))
  else
    typeof(rep)(vrep(rep, f))
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

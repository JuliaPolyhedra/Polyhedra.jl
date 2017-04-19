export affinehull

function affinehull{T<:HRep}(p::T)
    T(eqs=eqs(p))
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
        error("Cannot take the Minkowski sum of polyhedra of different element type")
    end
    if fulldim(T1) != fulldim(T2)
        error("Cannot take the Minkowski sum of polyhedra of different dimension")
    end
    T1(VRepIterator([p1, p2]))
end

# p1 has priority
function usehrep(p1::Polyhedron, p2::Polyhedron)
    hrepiscomputed(p1) && (!vrepiscomputed(p1) || hrepiscomputed(p2))
end

# Always type of first arg
#@generated
function (*){RepT1<:Rep, RepT2<:Rep}(p1::RepT1, p2::RepT2)
    if eltype(RepT1) != eltype(RepT2)
        error("Cannot take the cartesian product between polyhedra of different element type")
    end
    T = eltype(RepT1)
    N1 = fulldim(RepT1)
    N2 = fulldim(RepT2)
    Nout = N1 + N2
    hashrep = RepT1 <: HRepresentation || RepT2 <: HRepresentation
    hasvrep = RepT1 <: VRepresentation || RepT2 <: VRepresentation
    RepTout = changefulldim(RepT1, Nout)
    f = (i, x) -> zeropad(x, i == 1 ? N2 : -N1)
    if hashrep && hasvrep
        error("Cannot take the cartesian product between a H-Representation and a V-Representation")
    elseif hashrep || (!hasvrep && usehrep(p1, p2))
        # TODO fastdecompose
        # FIXME Nin, Tin are only the N and T of p1. This does not make sense.
        #       Do we really need these 2 last parameters ? I guess we should remove them
        RepTout(HRepIterator{Nout, T, N1, T}([p1, p2], f))
    else
        # TODO fastdecompose
        #:($RepTout(VRepIterator{$Nout, $T, $N1, $T}([p1, p2], $f)))
        RepTout(VRepIterator{Nout, T, N1, T}([p1, p2], f))
    end
end

function (*){RepT<:HRep}(rep::RepT, P::AbstractMatrix)
    Nin = fulldim(rep)
    Tin = eltype(rep)
    if size(P, 1) != Nin
        error("The number of rows of P must match the dimension of the H-representation")
    end
    Nout = size(P, 2)
    Tout = mypromote_type(eltype(RepT), eltype(P))
    RepTout = lazychangeboth(RepT, Nout, Tout)
    f = (i, h) -> h * P
    if decomposedhfast(rep)
        eqs = EqIterator{Nout,Tout,Nin,Tin}([rep], f)
        ineqs = IneqIterator{Nout,Tout,Nin,Tin}([rep], f)
        if RepT <: HRepresentation
            RepTout(ineqs=ineqs, eqs=eqs)
        else
            polyhedron(ineqs, eqs, getlibraryfor(rep, Nout, Tout))
        end
    else
        hreps = HRepIterator{Nout,Tout,Nin,Tin}([rep], f)
        if RepT <: HRepresentation
            RepTout(hreps)
        else
            polyhedron(hreps, getlibraryfor(rep, Nout, Tout))
        end
    end
end
function (*){RepT<:VRep}(P::AbstractMatrix, rep::RepT)
    Nin = fulldim(rep)
    Tin = eltype(rep)
    if size(P, 2) != Nin
        error("The number of rows of P must match the dimension of the H-representation")
    end
    Nout = size(P, 1)
    Tout = mypromote_type(eltype(RepT), eltype(P))
    RepTout = changeboth(RepT, Nout, Tout)
    f = (i, v) -> P * v
    if decomposedvfast(rep)
        points = PointIterator{Nout,Tout,Nin,Tin}([rep], f)
        rays = RayIterator{Nout,Tout,Nin,Tin}([rep], f)
        if RepT <: VRepresentation
            RepTout(points=points, rays=rays)
        else
            polyhedron(points, rays, getlibraryfor(rep, Nout, Tout))
        end
    else
        vreps = VRepIterator{Nout,Tout,Nin,Tin}([rep], f)
        if RepT <: VRepresentation
            RepTout(vreps)
        else
            polyhedron(vreps, getlibraryfor(rep, Nout, Tout))
        end
    end
end

function Base.round{N,T<:AbstractFloat}(rep::HRepresentation{N,T})
    f = (i, h) -> round(h)
    if decomposedfast(rep)
        typeof(rep)(eqs=eqs(rep, f), ineqs=ineqs(rep, f))
    else
        typeof(rep)(hreps(rep, f))
    end
end
function Base.round{N,T<:AbstractFloat}(rep::VRepresentation{N,T})
    f = (i, v) -> round(v)
    if decomposedfast(rep)
        typeof(rep)(eqs=eqs(rep, f), ineqs=ineqs(rep, f))
    else
        typeof(rep)(vreps(rep, f))
    end
end

export gethredundantindices, getvredundantindices

function gethredundantindices(hrep::HRep; strongly=false, solver = defaultLPsolverfor(hrep))
    red = IntSet()
    for i in 1:nhreps(hrep)
        if ishredundant(hrep, i; strongly, solver)
            push!(red, i)
        end
    end
    red
end
function getvredundantindices(vrep::VRep; strongly=false, solver = defaultLPsolverfor(vrep))
    red = IntSet()
    for i in 1:nvreps(vrep)
        if isvredundant(p, i; strongly, solver)
            push!(red, i)
        end
    end
    red
end

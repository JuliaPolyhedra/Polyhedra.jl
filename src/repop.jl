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
            polyhedron(ineqs, eqs, getlibraryfor(rep, Tout))
        end
    else
        hreps = HRepIterator{Nout,Tout,Nin,Tin}([rep], f)
        if RepT <: HRepresentation
            RepTout(hreps)
        else
            polyhedron(hreps, getlibraryfor(rep, Tout))
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
            polyhedron(points, rays, getlibraryfor(rep, Tout))
        end
    else
        vreps = VRepIterator{Nout,Tout,Nin,Tin}([rep], f)
        if RepT <: VRepresentation
            RepTout(vreps)
        else
            polyhedron(vreps, getlibraryfor(rep, Tout))
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

export convexhull
# Always type of first arg
function Base.intersect{N, T1, T2, RepTin<:HRep{N, T1}}(p1::RepTin, p2::HRep{N, T2})
    Tout = promote_type(T1, T2)
    RepTout = lazychangeeltype(RepTin, Tout)
    RepTout(HRepIterator([HRep{N, Tout}(p1), HRep{N, Tout}(p2)]))
end

"""
    convexhull(P1::VRep, P2::VRep)

Takes the convex hull of `P1` and `P2` ``\\{\\, \\lambda x + (1-\\lambda) y : x \\in P1, y \\in P2 \\,\\}``.
It is very efficient between two V-representations or between two polyhedron for which the V-representation has already been computed.
However, if `P1` (resp. `P2`) is a polyhedron for which the V-representation has not been computed yet, it will trigger a representation conversion which is costly.

The type of the result will be chosen closer to the type of `P1`. For instance, if `P1` is a polyhedron (resp. V-representation) and `P2` is a V-representation (resp. polyhedron), `convexhull(P1, P2)` will be a polyhedron (resp. V-representation).
If `P1` and `P2` are both polyhedra (resp. V-representation), the resulting polyhedron type (V-representation type) will be computing according to the type of `P1`.
The coefficient type however, will be promoted as required taking both the coefficient type of `P1` and `P2` into account.
```
"""
function convexhull{N, T1, T2, RepTin<:VRep{N, T1}}(p1::RepTin, p2::VRep{N, T2})
    Tout = promote_type(T1, T2)
    # Always type of first arg
    RepTout = lazychangeeltype(RepTin, Tout)
    RepTout(VRepIterator([VRep{N, Tout}(p1), VRep{N, Tout}(p2)]))
end

function (+){N, T1, T2, RepTin<:VRep{N, T1}}(p1::RepTin, p2::VRep{N, T2})
    Tout = promote_type(T1, T2)
    # Always type of first arg
    RepTout = lazychangeeltype(RepTin, Tout)
    ps = PointsHull{N, Tout}([po1 + po2 for po1 in points(p1) for po2 in points(p2)])
    rs = RaysHull(AbstractRay{N, Tout}[collect(rays(p1)); collect(rays(p2))])
    RepTout(points(ps), rays(rs))
end

# p1 has priority
function usehrep(p1::Polyhedron, p2::Polyhedron)
    hrepiscomputed(p1) && (!vrepiscomputed(p1) || hrepiscomputed(p2))
end

function hcartesianproduct{N1, N2, T, RepT1<:HRep{N1, T}, RepT2<:HRep{N2, T}}(p1::RepT1, p2::RepT2)
    Nout = N1 + N2
    # Always type of first arg
    RepTout = changefulldim(RepT1, Nout)
    f = (i, x) -> zeropad(x, i == 1 ? N2 : -N1)
    # TODO fastdecompose
    # FIXME Nin, Tin are only the N and T of p1. This does not make sense.
    #       Do we really need these 2 last parameters ? I guess we should remove them
    RepTout(HRepIterator{Nout, T, N1, T}([p1, p2], f))
end
function vcartesianproduct{N1, N2, T, RepT1<:VRep{N1, T}, RepT2<:VRep{N2, T}}(p1::RepT1, p2::RepT2)
    Nout = N1 + N2
    # Always type of first arg
    RepTout = changefulldim(RepT1, Nout)
    f1 = (i, x) -> zeropad(x, N2)
    f2 = (i, x) -> zeropad(x, -N1)
    # TODO fastdecompose
    # FIXME Nin, Tin are only the N and T of p1. This does not make sense.
    #       Do we really need these 2 last parameters ? I guess we should remove them
    q1 = changefulldim(RepT1, Nout)(VRepIterator{Nout, T, N1, T}([p1], f1))
    q2 = changefulldim(RepT2, Nout)(VRepIterator{Nout, T, N2, T}([p2], f2))
    q1 + q2
end
(*)(p1::HRep, p2::HRep) = hcartesianproduct(p1, p2)
(*)(p1::VRep, p2::VRep) = vcartesianproduct(p1, p2)
function (*)(p1::Polyhedron, p2::Polyhedron)
    if usehrep(p1, p2)
        hcartesianproduct(p1, p2)
    else
        vcartesianproduct(p1, p2)
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
    if RepT <: HRepresentation
        RepTout = lazychangeboth(RepT, Nout, Tout)
    end
    f = (i, h) -> h * P
    if decomposedhfast(rep)
        eqs = EqIterator{Nout,Tout,Nin,Tin}([rep], f)
        ineqs = IneqIterator{Nout,Tout,Nin,Tin}([rep], f)
        if RepT <: HRepresentation
            RepTout(eqs, ineqs)
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
            RepTout(points, rays)
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

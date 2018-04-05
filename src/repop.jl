export convexhull, convexhull!, conichull

const HAny{N, T} = Union{HRep{N, T}, HRepElement{N, T}}
const VAny{N, T} = Union{VRep{N, T}, VRepElement{N, T}}

_promote_reptype(P1::Type{<:HAffineSpace}, ::Type{<:HAffineSpace}) = P1
_promote_reptype(P1::Type{<:HAffineSpace}, ::Type{<:HRep}) = hreptype(P1)
_promote_reptype(P1::Type{<:HRep}, ::Type{<:HRep}) = P1

_promote_reptype(P1::Type{<:VPolytope}, ::Type{<:VPolytope}) = P1
_promote_reptype(P1::Type{<:VPolytope}, ::Type{<:VRep}) = vreptype(P1)
_promote_reptype(P1::Type{<:VLinearSpace}, ::Type{<:VLinearSpace}) = P1
_promote_reptype(P1::Type{<:VLinearSpace}, ::Type{<:VCone}) = conetype(P1)
_promote_reptype(P1::Type{<:VLinearSpace}, ::Type{<:VRep}) = vreptype(P1)
_promote_reptype(P1::Type{<:VCone}, ::Type{<:VCone}) = P1
_promote_reptype(P1::Type{<:VCone}, ::Type{<:VRep}) = vreptype(P1)
_promote_reptype(P1::Type{<:VRep}, ::Type{<:VRep}) = P1

function promote_reptype(P1::Type{<:Rep{N, S}}, P2::Type{<:Rep{N, T}}) where {N, S, T}
    U = promote_type(S, T)
    # P2 might not support coefficient type U so we avoid calling similar_type on P2
    _promote_reptype(similar_type(P1, U), P2)
end
promote_reptype(P1::Type{<:Rep{N}}, P2::Type{<:Rep{N}}, P::Type{<:Rep{N}}...) where N = promote_reptype(promote_reptype(P1, P2), P...)

"""
    intersect(P1::HRep, P2::HRep)

Takes the intersection of `P1` and `P2` ``\\{\\, x : x \\in P_1, x \\in P_2 \\,\\}``.
It is very efficient between two H-representations or between two polyhedron for which the H-representation has already been computed.
However, if `P1` (resp. `P2`) is a polyhedron for which the H-representation has not been computed yet, it will trigger a representation conversion which is costly.
See the [Polyhedral Computation FAQ](http://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node25.html) for a discussion on this operation.

The type of the result will be chosen closer to the type of `P1`. For instance, if `P1` is a polyhedron (resp. H-representation) and `P2` is a H-representation (resp. polyhedron), `intersect(P1, P2)` will be a polyhedron (resp. H-representation).
If `P1` and `P2` are both polyhedra (resp. H-representation), the resulting polyhedron type (resp. H-representation type) will be computed according to the type of `P1`.
The coefficient type however, will be promoted as required taking both the coefficient type of `P1` and `P2` into account.
"""
function Base.intersect(p::HRep{N}...) where N
    RepTout = promote_reptype(typeof.(p)...)
    T = MultivariatePolynomials.coefficienttype(RepTout)
    RepTout(hmap((i, x) -> similar_type(typeof(x), T)(x), FullDim{N}(), T, p...)...)
end
Base.intersect(p::Rep, el::HRepElement) = p ∩ intersect(el)

Base.intersect(hps::HyperPlane...) = hrep([hps...])
Base.intersect(hss::HalfSpace...) = hrep([hss...])
Base.intersect(h1::HyperPlane, h2::HalfSpace) = hrep([h1], [h2])
Base.intersect(h1::HalfSpace, h2::HyperPlane) = h2 ∩ h1
Base.intersect(h1::Union{HRep{N}, HRepElement{N}}, h2::Union{HRep{N}, HRepElement{N}}, hs::Union{HRep{N}, HRepElement{N}}...) where N = intersect(h1 ∩ h2, hs...)

"""
    intersect!(p1::VRep, p2::VRep)

Same as [`intersect`](@ref) except that `p1` is modified to be equal to the intersection.
"""
Base.intersect!(p::HRep{N}, ine::HRepresentation{N}) where {N} = error("intersect! not implemented for $(typeof(p)). It probably does not support in-place modification, try `intersect` (without the `!`) instead.")

"""
    convexhull(P1::VRep, P2::VRep)

Takes the convex hull of `P1` and `P2` ``\\{\\, \\lambda x + (1-\\lambda) y : x \\in P_1, y \\in P_2 \\,\\}``.
It is very efficient between two V-representations or between two polyhedron for which the V-representation has already been computed.
However, if `P1` (resp. `P2`) is a polyhedron for which the V-representation has not been computed yet, it will trigger a representation conversion which is costly.

The type of the result will be chosen closer to the type of `P1`. For instance, if `P1` is a polyhedron (resp. V-representation) and `P2` is a V-representation (resp. polyhedron), `convexhull(P1, P2)` will be a polyhedron (resp. V-representation).
If `P1` and `P2` are both polyhedra (resp. V-representation), the resulting polyhedron type (resp. V-representation type) will be computed according to the type of `P1`.
The coefficient type however, will be promoted as required taking both the coefficient type of `P1` and `P2` into account.
"""
function convexhull(p::VRep{N}...) where N
    RepTout = promote_reptype(typeof.(p)...)
    T = MultivariatePolynomials.coefficienttype(RepTout)
    RepTout(vmap((i, x) -> similar_type(typeof(x), T)(x), FullDim{N}(), T, p...)...)::RepTout # FIXME without this type annotation even convexhull(::PointsHull{2,Int64,Array{Int64,1}}, ::PointsHull{2,Int64,Array{Int64,1}}) is not type stable, why ?
end
convexhull(p::Rep, el::VRepElement) = convexhull(p, convexhull(el))

convexhull(ps::AbstractPoint...) = vrep([ps...])
convexhull(ls::Line...) = vrep([ls...])
convexhull(rs::Ray...) = vrep([rs...])
convexhull(p::AbstractPoint{N, T}, r::Union{Line{N, T}, Ray{N, T}}) where {N, T} = vrep([p], [r])
convexhull(r::Union{Line{N, T}, Ray{N, T}}, p::AbstractPoint{N, T}) where {N, T} = conichull(p, r)
convexhull(l::Line{N, T}, r::Ray{N, T}) where {N, T} = vrep([l], [r])
convexhull(r::Ray{N, T}, l::Line{N, T}) where {N, T} = conichull(l, r)
convexhull(p1::VAny{N, T}, p2::VAny{N, T}, ps::VAny{N, T}...) where {N, T} = convexhull(convexhull(p1, p2), ps...)
function convexhull(p::VAny{N}...) where N
    T = promote_type(MultivariatePolynomials.coefficienttype.(p)...)
    f(p) = similar_type(typeof(p), T)(p)
    convexhull(f.(p)...)
end

"""
    convexhull!(p1::VRep, p2::VRep)

Same as [`convexhull`](@ref) except that `p1` is modified to be equal to the convex hull.
"""
convexhull!(p::VRep{N}, ine::VRepresentation{N}) where {N} = error("convexhull! not implemented for $(typeof(p)). It probably does not support in-place modification, try `convexhull` (without the `!`) instead.")

# conify: same than conichull except that conify(::VRepElement) returns a VRepElement and not a V-representation
conify(v::VRep) = vrep(lines(v), [collect(rays(v)); collect(Ray.(points(v)))])
conify(v::VCone) = v
conify(p::AbstractPoint) = Ray(p)
conify(r::Union{Line, Ray}) = r

conichull(p...) = convexhull(conify.(p)...)

function sumpoints(::FullDim{N}, ::Type{T}, p1, p2) where {N, T}
    _tout(p) = similar_type(typeof(p), T)(p)
    ps = [_tout(po1 + po2) for po1 in points(p1) for po2 in points(p2)]
    tuple(ps)
end
sumpoints(::FullDim{N}, ::Type{T}, p1::Rep, p2::VCone) where {N, T} = RepIterator{N, T}.(preps(p1))
sumpoints(::FullDim{N}, ::Type{T}, p1::VCone, p2::Rep) where {N, T} = RepIterator{N, T}.(preps(p2))

function Base.:+(p1::VRep{N, T1}, p2::VRep{N, T2}) where {N, T1, T2}
    Tout = promote_type(T1, T2)
    RepTout = promote_reptype(typeof(p1), typeof(p2))
    RepTout(sumpoints(FullDim{N}(), Tout, p1, p2)..., RepIterator{N, Tout}.(rreps(p1, p2))...)
end
Base.:+(p::Rep, el::Union{Line, Ray}) = p + vrep([el])
Base.:+(el::Union{Line, Ray}, p::Rep) = p + el

# p1 has priority
function usehrep(p1::Polyhedron, p2::Polyhedron)
    hrepiscomputed(p1) && (!vrepiscomputed(p1) || hrepiscomputed(p2))
end

function hcartesianproduct(p1::RepT1, p2::RepT2) where {N1, N2, T1, T2, RepT1<:HRep{N1, T1}, RepT2<:HRep{N2, T2}}
    dout = FullDim{N1+N2}()
    Tout = promote_type(T1, T2)
    # Always type of first arg
    RepTout = similar_type(RepT1, dout, Tout)
    f = (i, x) -> zeropad(x, i == 1 ? N2 : -N1)
    RepTout(hmap(f, dout, Tout, p1, p2)...)
end
function vcartesianproduct(p1::RepT1, p2::RepT2) where {N1, N2, T1, T2, RepT1<:VRep{N1, T1}, RepT2<:VRep{N2, T2}}
    dout = FullDim{N1+N2}()
    Tout = promote_type(T1, T2)
    # Always type of first arg
    RepTout = similar_type(RepT1, dout)
    f1 = (i, x) -> zeropad(x, N2)
    f2 = (i, x) -> zeropad(x, -N1)
    q1 = similar_type(RepT1, dout, Tout)(vmap(f1, dout, Tout, p1)...)
    q2 = similar_type(RepT2, dout, Tout)(vmap(f2, dout, Tout, p2)...)
    q1 + q2
end
cartesianproduct(p1::HRep, p2::HRep) = hcartesianproduct(p1, p2)
cartesianproduct(p1::VRep, p2::VRep) = vcartesianproduct(p1, p2)

function cartesianproduct(p1::Polyhedron, p2::Polyhedron)
    if usehrep(p1, p2)
        hcartesianproduct(p1, p2)
    else
        vcartesianproduct(p1, p2)
    end
end

"""
    *(p1::Rep, p2::Rep)

Cartesian product between the polyhedra `p1` and `p2`.
"""
Base.:(*)(p1::Rep, p2::Rep) = cartesianproduct(p1, p2)

"""
    \\(P::AbstractMatrix, p::HRep)

Transform the polyhedron represented by ``p`` into ``P^{-1} p`` by transforming each halfspace ``\\langle a, x \\rangle \\le \\beta`` into ``\\langle P^\\top a, x \\rangle \\le \\beta`` and each hyperplane ``\\langle a, x \\rangle = \\beta`` into ``\\langle P^\\top a, x \\rangle = \\beta``.
"""
Base.:(\)(P::AbstractMatrix, rep::HRep) = rep / P'

"""
    /(p::HRep, P::AbstractMatrix)

Transform the polyhedron represented by ``p`` into ``P^{-T} p`` by transforming each halfspace ``\\langle a, x \\rangle \\le \\beta`` into ``\\langle P a, x \\rangle \\le \\beta`` and each hyperplane ``\\langle a, x \\rangle = \\beta`` into ``\\langle P a, x \\rangle = \\beta``.
"""
function Base.:(/)(p::RepT, P::AbstractMatrix) where {Nin, Tin, RepT<:HRep{Nin, Tin}}
    if size(P, 2) != Nin
        throw(DimensionMismatch("The number of rows of P must match the dimension of the H-representation"))
    end
    f = (i, h) -> h / P
    # For a matrix P of StaticArrays, `dout` should be type stable
    dout = FullDim{size(P, 1)}()
    Tout = _promote_type(Tin, eltype(P))
    RepTout = similar_type(RepT, dout, Tout)
    RepTout(hmap(f, dout, Tout, p)...)
end

function Base.:(*)(rep::HRep, P::AbstractMatrix)
    warn("`*(p::HRep, P::AbstractMatrix)` is deprecated. Use `P \\ p` or `p / P'` instead.")
    P \ rep
end

"""
    *(P::AbstractMatrix, p::VRep)

Transform the polyhedron represented by ``p`` into ``P p`` by transforming each element of the V-representation (points, symmetric points, rays and lines) `x` into ``P x``.
"""
function Base.:(*)(P::AbstractMatrix, p::RepT) where {Nin, Tin, RepT<:VRep{Nin, Tin}}
    if size(P, 2) != Nin
        throw(DimensionMismatch("The number of rows of P must match the dimension of the V-representation"))
    end
    f = (i, v) -> P * v
    # For a matrix P of StaticArrays, `dout` should be type stable
    dout = FullDim{size(P, 1)}()
    Tout = _promote_type(Tin, eltype(P))
    RepTout = similar_type(RepT, dout, Tout)
    RepTout(vmap(f, dout, Tout, p)...)
end

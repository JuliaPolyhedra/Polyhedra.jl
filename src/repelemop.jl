export ininterior, inrelativeinterior

hyperplane(h::HyperPlane) = h
hyperplane(h::HalfSpace) = HyperPlane(h.a, h.β)

halfspace(h::HyperPlane) = HalfSpace(h.a, h.β)
halfspace(h::HalfSpace) = h

line(r::Ray) = Line(coord(r))
linetype(::Type{Ray{N, T, AT}}) where {N, T, AT} = Line{N, T, AT}

###############
# TRANSLATION #
###############

export translate

function htranslate(p::HRep{N, T}, v::Union{AbstractPoint{N, S}}) where {N, S, T}
    f = (i, h) -> translate(h, v)
    Tout = Base.promote_op(+, T, S)
    similar_type(typeof(p), Tout)(hmap(f, FullDim{N}(), Tout, p)...)
end

function vtranslate{N, S, T}(p::VRep{N, T}, v::Union{AbstractPoint{N, S}})
    f = (i, u) -> translate(u, v)
    Tout = Base.promote_op(+, T, S)
    similar_type(typeof(p), Tout)(vmap(f, FullDim{N}(), Tout, p)...)
end

translate(p::HRep, v) = htranslate(p, v)
translate(p::VRep, v) = vtranslate(p, v)
function translate(p::Polyhedron, v)
    if hrepiscomputed(p)
        htranslate(p, v)
    else
        vtranslate(p, v)
    end
end

######
# IN #
######
function _vinh(v::VRepElement{N}, it::ElemIt{<:HRepElement{N}}, infun) where N
    for h in it
        if !infun(v, h)
            return false
        end
    end
    return true
end
function _vinh(v::VRepElement{N}, hr::HRep{N}, infun) where N
    all(map(it -> _vinh(v, it, infun), hreps(hr)))
end

"""
    in(p::VRepElement, h::HRepElement)

Returns whether `p` is in `h`.
If `h` is an hyperplane, it returns whether ``\\langle a, x \\rangle \\approx \\beta``.
If `h` is an halfspace, it returns whether ``\\langle a, x \\rangle \\le \\beta``.

    in(p::VRepElement, h::HRep)

Returns whether `p` is in `h`, e.g. in all the hyperplanes and halfspaces supporting `h`.
"""
Base.in(v::VRepElement{N}, hr::HRep{N}) where {N} = _vinh(v, hr, Base.in)

"""
    ininterior(p::VRepElement, h::HRepElement)

Returns whether `p` is in the interior of `h`.
If `h` is an hyperplane, it always returns `false`.
If `h` is an halfspace ``\\langle a, x \\rangle \\leq \\beta``, it returns whether `p` is in the open halfspace ``\\langle a, x \\rangle < \\beta``

    ininterior(p::VRepElement, h::HRep)

Returns whether `p` is in the interior of `h`, e.g. in the interior of all the hyperplanes and halfspaces supporting `h`.
"""
ininterior(v::VRepElement{N}, hr::HRep{N}) where {N} = _vinh(v, hr, ininterior)

"""
    inrelativeinterior(p::VRepElement, h::HRepElement)

Returns whether `p` is in the relative interior of `h`.
If `h` is an hyperplane, it is equivalent to `p in h` since the relative interior of an hyperplane is itself.
If `h` is an halfspace, it is equivalent to `ininterior(p, h)`.

    inrelativeinterior(p::VRepElement, h::HRep)

Returns whether `p` is in the relative interior of `h`, e.g. in the relative interior of all the hyperplanes and halfspaces supporting `h`.
"""
inrelativeinterior(v::VRepElement{N}, hr::HRep{N}) where {N} = _vinh(v, hr, inrelativeinterior)

function _hinv(h::HRepElement{N}, vr::ElemIt{<:VRepElement{N}}) where N
    all(in.(vr, h))
end
function _hinv(h::HRepElement{N}, vr::VRep{N}) where N
    all(_hinv.(h, vreps(vr)))
end
function _hinh(h::HalfSpace{N}, hr::HRep{N}, solver) where N
    # ⟨a, x⟩ ≦ β -> if β < max ⟨a, x⟩ then h is outside
    sol = linprog(-h.a, hr, solver)
    if sol.status == :Unbounded
        false
    elseif sol.status == :Infeasible
        true
    elseif sol.status == :Optimal
        _leq(-sol.objval, h.β)
    else
        error("Solver returned with status $(sol.status)")
    end
end
function _hinh(h::HyperPlane{N}, hr::HRep{N}, solver) where N
    _hinh(halfspace(h), hr, solver) && _hinh(halfspace(-h), hr, solver)
end

Base.issubset(vr::VRepresentation{N}, h::HRepElement{N}) where {N} = _hinv(h, vr)
Base.issubset(hr::HRepresentation{N}, h::HRepElement{N}) where {N} = _hinh(h, hr)

"""
    issubset(p::Rep, h::HRepElement)

Returns whether `p` is a subset of `h`, i.e. whether `h` supports the polyhedron `p`.
"""
function Base.issubset(p::Polyhedron{N}, h::HRepElement{N}, solver=Polyhedra.solver(p)) where N
    if vrepiscomputed(p)
        _hinv(h, p)
    else
        _hinh(h, p, solver)
    end
end

function Base.in(h::HRepElement, p::Rep)
    warn("in(h::HRepElement, p::Rep) is deprecated. Use issubset(p, h) instead.")
    issubset(p, h)
end

################
# INTERSECTION #
################
_intres(h::HyperPlane, ins, inp) = inp
_intres(h::HalfSpace, ins, inp) = [ins; inp]

# Unused as allpoints and allrays are used
#function _pushinout!(ins, out, l::Line, h::HalfSpace)
#    r1 = Ray(l.a)
#    r2 = Ray(-l.a)
#    if r1 in h
#        push!(ins, r1)
#        push!(out, r2)
#    else
#        push!(ins, r2)
#        push!(out, r1)
#    end
#end
#function _pushinout!(ins, out, p::SymPoint, h::HalfSpace)
#    p1 = Point(p.a)
#    p2 = Point(-p.a)
#    if p1 in h
#        push!(ins, p1)
#        push!(out, p2)
#    else
#        push!(ins, p2)
#        push!(out, p1)
#    end
#end

function _pushinout!(ins, out, pr::Union{AbstractPoint, Ray}, h::HalfSpace)
    if pr in h
        push!(ins, pr)
    else
        push!(out, pr)
    end
end

"""
    intersect(v::VRep{N, T}, h::HRepElement)

Compute the intersection of `v` with an halfspace or hyperplane `h`.
The method used by default is to keep the V-representation element of `v`
that are in `h` and add new ones generated as the intersection between
the hyperplane defining `h` and the segment between two adjacent
V-representation elements of `v` that are in either sides of the hyperplane.
See Lemma 3 of [1] for more detail on the method.

[1] Fukuda, K. and Prodon, A.
Double description method revisited
*Combinatorics and computer science*, *Springer*, **1996**, 91-111
"""
function Base.intersect(v::VRep{N, T}, h::HRepElement) where {N, T}
    PointT = pointtype(v)
    pins = PointT[] # Inside
    pinp = PointT[] # In plane
    pout = PointT[] # Outside
    hp = hyperplane(h)
    hs = halfspace(h)
    for p in points(v) # we want sympoints to be treated as 2 points
        if p in hp
            push!(pinp, p)
        else
            _pushinout!(pins, pout, p, hs)
        end
    end
    RayT = raytype(v)
    rins = RayT[]
    rinp = RayT[]
    rout = RayT[]
    for r in allrays(v) # we want lines to be treated as 2 rays
        if r in hp
            push!(rinp, r)
        else
            _pushinout!(rins, rout, r, hs)
        end
    end
    # \    For every pair of point inside/outside,
    # _.__ add the intersection of the segment with the hyperplane.
    #   \  In CDD and ConvexHull.jl, they only consider segments that
    #    \ belongs to the boundary which is way faster than what we do here
    for p in pins
        ap = h.a ⋅ p
        for q in pout
            aq = h.a ⋅ q
            λ = (aq - h.β) / (aq - ap)
            @assert 0 <= λ <= 1
            push!(pinp, simplify(λ * p + (1 - λ) * q))
        end
    end
    # Similar but with rays
    for r in rins
        ar = h.a ⋅ r
        for s in rout
            as = h.a ⋅ s
            # should take
            # λ = as / (as - ar)
            @assert 0 <= as / (as - ar) <= 1
            # By homogeneity we can avoid the division and do
            #newr = as * r - ar * s
            # but this can generate very large numbers (see JuliaPolyhedra/Polyhedra.jl#48)
            # so we still divide
            newr = (as * r - ar * s) / (as - ar)
            # In CDD, it does as * r - ar * s but then it normalize the ray
            # by dividing it by its smallest nonzero entry (see dd_CreateNewRay)
            if !isapproxzero(coord(newr))
                push!(rinp, simplify(newr))
            end
        end
    end
    # Similar but with one point and one ray
    for (ps, rs) in ((pins, rout), (pout, rins))
        for p in ps
            ap = h.a ⋅ p
            for r in rs
                ar = h.a ⋅ r
                λ = (h.β - ap) / ar
                push!(pinp, p + λ * r)
            end
        end
    end
    psout = _intres(h, pins, pinp)
    rsout = _intres(h, rins, rinp)
    if isempty(psout)
        # Empty polyhedron, there may be rays left,
        # Example 1: for 0x_1 + x_2 = -1 ∩ 0x_1 + x_2 = 1, the line (0, 1) is detected as correct
        # Example 2: for 0x_1 + 0x_2 = 1, the lines (1, 0) and (0, 1) are detected as correct
        # but since there is no point, the polyhedron is empty and we should drop all rays/lines
        empty!(rsout)
    end
    typeof(v)(psout, linetype(eltype(rsout))[], rsout)
end

export ininterior, inrelativeinterior

hyperplane(h::HyperPlane) = h
hyperplane(h::HalfSpace) = HyperPlane(h.a, h.β)

halfspace(h::HyperPlane) = HalfSpace(h.a, h.β)
halfspace(h::HalfSpace) = h

line(r::Ray) = Line(coord(r))

###############
# TRANSLATION #
###############

export translate

function translate(p::HRep{N, T}, v::Union{AbstractArray{S}, Point{N, S}}) where {N, S, T}
    f = (i, h) -> translate(h, v)
    Tout = Base.promote_op(+, T, S)
    lazychangeeltype(typeof(p), Tout)(HRepIterator{N, Tout, N, T}([p], f))
end

function translate{N, S, T}(p::VRep{N, T}, v::Union{AbstractArray{S}, Point{N, S}})
    f = (i, u) -> translate(u, v)
    Tout = Base.promote_op(+, T, S)
    lazychangeeltype(typeof(p), Tout)(VRepIterator{N, Tout, N, T}([p], f))
end

######
# IN #
######
function _vinh(v::VRepElement{N}, hr::HRep{N}, infun) where N
    for h in hreps(hr)
        if !infun(v, h)
            return false
        end
    end
    return true
end
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

function _hinv(h::HRepElement{N}, vr::VRep{N}) where N
    for v in vreps(vr)
        if !(v in h)
            return false
        end
    end
    return true
end
function _hinh(h::HalfSpace{N}, hr::HRep{N}, solver) where N
    # ⟨a, x⟩ ≦ β -> if β < max ⟨a, x⟩ then h is outside
    sol = linprog(-h.a, hr, solver)
    if sol.status == :Unbounded
        false
    elseif sol.status == :Infeasible
        true
    elseif sol.status == :Optimal
        myleq(-sol.objval, h.β)
    else
        error("Solver returned with status $(sol.status)")
    end
end
function _hinh(h::HyperPlane{N}, hr::HRep{N}, solver) where N
    _hinh(halfspace(h), hr, solver) && _hinh(halfspace(-h), hr, solver)
end

Base.in(h::HRepElement{N}, vr::VRepresentation{N}) where {N} = _hinv(h, vr)
Base.in(h::HRepElement{N}, hr::HRepresentation{N}) where {N} = _hinh(h, hr)
function Base.in(h::HRepElement{N}, p::Polyhedron{N}, solver = defaultLPsolverfor(p)) where N
    if vrepiscomputed(p)
        _hinv(h, p)
    else
        _hinh(h, p, solver)
    end
end


################
# INTERSECTION #
################
function _intres(v::T, h::HyperPlane, pins, pinp, rins, rinp) where T
    T(pinp, rinp)
end

function _intres(v::T, h::HalfSpace, pins, pinp, rins, rinp) where T
    T([pins; pinp], [rins; rinp])
end

function _pushinout!(ins, out, l::Line, h::HalfSpace)
    r1 = Ray(l.a)
    r2 = Ray(-l.a)
    if r1 in h
        push!(ins, r1)
        push!(out, r2)
    else
        push!(ins, r2)
        push!(out, r1)
    end
end
function _pushinout!(ins, out, p::SymPoint, h::HalfSpace)
    p1 = Point(p.a)
    p2 = Point(-p.a)
    if p1 in h
        push!(ins, p1)
        push!(out, p2)
    else
        push!(ins, p2)
        push!(out, p1)
    end
end
function _pushinout!(ins, out, pr::Union{AbstractPoint, AbstractRay}, h::HalfSpace)
    if pr in h
        push!(ins, pr)
    else
        push!(out, pr)
    end
end

function Base.intersect(v::VRep{N, T}, h::HRepElement) where {N, T}
    pins = AbstractPoint{N, T}[] # Inside
    pinp = AbstractPoint{N, T}[] # In plane
    pout = AbstractPoint{N, T}[] # Outside
    hp = hyperplane(h)
    hs = halfspace(h)
    for p in points(v)
        if p in hp
            push!(pinp, p)
        else
            _pushinout!(pins, pout, p, hs)
        end
    end
    if !haspoints(v)
        # The origin
        if !myeqzero(h.β)
            _pushinout!(pins, pout, zeros(T, N), hs)
        end
    end
    rins = AbstractRay{N, T}[]
    rinp = AbstractRay{N, T}[]
    rout = AbstractRay{N, T}[]
    for r in rays(v)
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
        ap = mydot(h.a, p)
        for q in pout
            aq = mydot(h.a, q)
            λ = (aq - h.β) / (aq - ap)
            @assert 0 <= λ <= 1
            push!(pinp, λ * p + (1 - λ) * q)
        end
    end
    # Similar but with rays
    for r in rins
        ar = mydot(h.a, r)
        for s in rout
            as = mydot(h.a, s)
            # should take
            # λ = as / (as - ar)
            @assert 0 <= as / (as - ar) <= 1
            # but by homogeneity we can avoid the division
            newr = as * r - ar * s
            if !myeqzero(coord(newr))
                push!(rinp, newr)
            end
        end
    end
    # Similar but with one point and one ray
    for (ps, rs) in ((pins, rout), (pout, rins))
        for p in ps
            ap = mydot(h.a, p)
            for r in rs
                ar = mydot(h.a, r)
                λ = (h.β - ap) / ar
                push!(pinp, p + λ * r)
            end
        end
    end
    _intres(v, h, pins, pinp, rins, rinp)
end

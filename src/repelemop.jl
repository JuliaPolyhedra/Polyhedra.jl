hyperplane(h::HyperPlane) = h
hyperplane(h::HalfSpace) = HyperPlane(h.a, h.β)

halfspace(h::HyperPlane) = HalfSpace(h.a, h.β)
halfspace(h::HalfSpace) = h

line(r::Ray) = Line(coord(r))

###############
# TRANSLATION #
###############

export translate

function translate{N, S, T}(p::HRep{N, T}, v::Union{AbstractArray{S}, Point{N, S}})
    f = (i, h) -> translate(h, v)
    Tout = Base.promote_op(+, T, S)
    lazychangeeltype(typeof(p), Tout)(HRepIterator{N, Tout, N, T}([p], f))
end


######
# IN #
######
function Base.in{N}(v::VRepElement{N}, hr::HRep{N})
    for h in hreps(hr)
        if !(v in h)
            return false
        end
    end
    return true
end

function _hinv{N}(h::HRepElement{N}, vr::VRep{N})
    for v in vreps(vr)
        if !(v in h)
            return false
        end
    end
    return true
end
function _hinh{N}(h::HalfSpace{N}, hr::HRep{N}, solver)
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
function _hinh{N}(h::HyperPlane{N}, hr::HRep{N}, solver)
    _hinh(halfspace(h), hr, solver) && _hinh(halfspace(-h), hr, solver)
end

Base.in{N}(h::HRepElement{N}, vr::VRepresentation{N}) = _hinv(h, vr)
Base.in{N}(h::HRepElement{N}, hr::HRepresentation{N}) = _hinh(h, hr)
function Base.in{N}(h::HRepElement{N}, p::Polyhedron{N}, solver = defaultLPsolverfor(p))
    if vrepiscomputed(p)
        _hinv(h, p)
    else
        _hinh(h, p, solver)
    end
end


################
# INTERSECTION #
################
function _intres{T}(v::T, h::HyperPlane, pins, pinp, rins, rinp)
    T(pinp, rinp)
end

function _intres{T}(v::T, h::HalfSpace, pins, pinp, rins, rinp)
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

function Base.intersect{N, T}(v::VRep{N, T}, h::HRepElement)
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

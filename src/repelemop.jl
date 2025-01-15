export ininterior, inrelativeinterior, support_function

hyperplane(h::HyperPlane) = h
hyperplane(h::HalfSpace) = HyperPlane(h.a, h.β)

halfspace(h::HyperPlane) = HalfSpace(h.a, h.β)
halfspace(h::HalfSpace) = h

line(r::Ray) = Line(coord(r))
linetype(::Type{Ray{T, AT}}) where {T, AT} = Line{T, AT}

###############
# TRANSLATION #
###############

export translate

function htranslate(p::HRep{T}, v::Union{AbstractVector{S}}) where {S, T}
    f = (i, h) -> translate(h, v)
    d = FullDim(p)
    Tout = Base.promote_op(+, T, S)
    similar(p, d, Tout, hmap(f, d, Tout, p)...)
end

function vtranslate(p::VRep{T}, v::Union{AbstractVector{S}}) where {S, T}
    f = (i, u) -> translate(u, v)
    d = FullDim(p)
    Tout = Base.promote_op(+, T, S)
    similar(p, d, Tout, vmap(f, d, Tout, p)...)
end

"""
    translate(p::Polyhedra.Rep, v::AbstractVector)

Computes translation of the polyhedron `p` with the vector `v`.
That is, computes
```math
\\{\\, x + v \\mid x \\in p \\,\\}.
```
By default, if the H-representation, it simply translates every hyperplanes and halfspace,
otherwise, it translates every points of the V-representation.
That is, this operation can be achieved both in the H-representation and V-representation
hence does not trigger any representation conversion.
"""
function translate end

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
function _vinh(v::VRepElement, it::ElemIt{<:HRepElement}, infun; tol)
    for h in it
        if !infun(v, h; tol)
            return false
        end
    end
    return true
end

function _vinh(v::VRepElement, hr::HRep, infun; kws...)
    # The line below fails on Julia v0.7. On the simplextest,
    # map(...) returns (false,) and then all((false,)) returns true
    # because (false,)[1] returns true...
    # all(map(it -> _vinh(v, it, infun), hreps(hr)))
    for it in hreps(hr)
        if !_vinh(v, it, infun; kws...)
            return false
        end
    end
    return true
end

"""
    in(p::VRepElement, h::HRepElement)

Returns whether `p` is in `h`.
If `h` is an hyperplane, it returns whether ``\\langle a, x \\rangle \\approx \\beta``.
If `h` is an halfspace, it returns whether ``\\langle a, x \\rangle \\le \\beta``.

    in(p::VRepElement, h::HRep)

Returns whether `p` is in `h`, e.g. in all the hyperplanes and halfspaces supporting `h`.
"""
function Base.in(v::VRepElement{T}, hr::HRep{U}; tol = _default_tol(T, U)) where {T,U}
    _vinh(v, hr, Base.in; tol)
end

"""
    ininterior(p::VRepElement, h::HRepElement)

Returns whether `p` is in the interior of `h`.
If `h` is an hyperplane, it always returns `false`.
If `h` is an halfspace ``\\langle a, x \\rangle \\leq \\beta``, it returns whether `p` is in the open halfspace ``\\langle a, x \\rangle < \\beta``

    ininterior(p::VRepElement, h::HRep)

Returns whether `p` is in the interior of `h`, e.g. in the interior of all the hyperplanes and halfspaces supporting `h`.
"""
function ininterior(v::VRepElement{T}, hr::HRep{U}; tol = _default_tol(T, U)) where {T,U}
    return _vinh(v, hr, ininterior; tol)
end

"""
    inrelativeinterior(p::VRepElement, h::HRepElement; tol)

Returns whether `p` is in the relative interior of `h`.
If `h` is an hyperplane, it is equivalent to `p in h` since the relative interior of an hyperplane is itself.
If `h` is an halfspace, it is equivalent to `ininterior(p, h)`.

    inrelativeinterior(p::VRepElement, h::HRep; tol)

Returns whether `p` is in the relative interior of `h`, e.g. in the relative interior of all the hyperplanes and halfspaces supporting `h`.
"""
function inrelativeinterior(v::VRepElement{T}, hr::HRep{U}; tol = _default_tol(T, U)) where {T,U}
    return _vinh(v, hr, inrelativeinterior; tol)
end

structural_nonzero_indices(a::SparseArrays.SparseVector) = SparseArrays.nonzeroinds(a)
structural_nonzero_indices(a::AbstractVector) = eachindex(a)
function _dot_terms(a, x::Vector{MOI.VariableIndex}, T::Type)
    return MOI.ScalarAffineTerm{T}[
        MOI.ScalarAffineTerm{T}(a[i], x[i])
        for i in structural_nonzero_indices(a) if !iszero(a[i])
    ]
end
function _dot(a, x::Vector{MOI.VariableIndex}, T::Type)
    return MOI.ScalarAffineFunction(_dot_terms(a, x, T), zero(T))
end

function support_function_model(h::AbstractVector, rep::Rep, solver)
    length(h) != fulldim(rep) && throw(DimensionMismatch())
    model, T = layered_optimizer(solver)
    x, cx = MOI.add_constrained_variables(model, PolyhedraOptSet(rep))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), _dot(h, x, T))
    return model
end

"""
   support_function(h::AbstractVector, rep::Rep, solver=Polyhedra.linear_objective_solver(p))

Return the value of the support function of `rep` at `h`. See [Section 13, R15] for more details.

[R15] Rockafellar, R.T. *Convex analysis*. Princeton university press, **2015**.
"""
function support_function(h::AbstractVector, rep::Rep, solver=Polyhedra.linear_objective_solver(rep))
    model = support_function_model(h, rep, solver)
    MOI.optimize!(model)
    term = MOI.get(model, MOI.TerminationStatus())
    if term == MOI.OPTIMAL
        return MOI.get(model, MOI.ObjectiveValue())
    else
        error("Cannot determine the value of the support function at $h",
              " because the linear program terminated with status $term.")
    end
end

function _hinv(h::HRepElement, vr::ElemIt{<:VRepElement}; kws...)
    return all(v -> in(v, h; kws...), vr)
end

function _hinv(h::HRepElement, vr::VRep; kws...)
    return all(v -> _hinv(h, v; kws...), vreps(vr))
end

function _hinh(h::HalfSpace, hr::HRep, solver; tol)
    # ⟨a, x⟩ ≦ β -> if β < max ⟨a, x⟩ then h is outside
    model = support_function_model(h.a, hr, solver)
    MOI.optimize!(model)
    term = MOI.get(model, MOI.TerminationStatus())
    if term == MOI.DUAL_INFEASIBLE
        return false
    elseif term == MOI.INFEASIBLE
        return true
    elseif term == MOI.OPTIMAL
        obj = MOI.get(model, MOI.ObjectiveValue())
        # If `h` and `hr` have rational representations,
        # the default tolerance is `MA.Zero`.
        # However, if the solver only supports `Float64` then
        # `obj` is `Float64` so we'll have to use some tolerance
        if obj isa AbstractFloat && tol isa MA.Zero
            tol = _default_tol(typeof(obj))
        end
        return _leq(obj, h.β; tol)
    else
        error("Cannot determine whether the polyhedron is contained in the",
              " halfspace or not because the linear program terminated with",
              " status $term.")
    end
end
function _hinh(h::HyperPlane, hr::HRep, solver; kws...)
    _hinh(halfspace(h), hr, solver; kws...) && _hinh(halfspace(-h), hr, solver; kws...)
end

function Base.issubset(vr::VRepresentation{T}, h::HRepElement{U}; tol = _default_tol(T, U)) where {T,U}
    return _hinv(h, vr; tol)
end

function Base.issubset(hr::HRepresentation{T}, h::HRepElement{U}; tol = _default_tol(T, U)) where {T,U}
    return _hinh(h, hr; tol)
end

"""
    issubset(p::Rep, h::HRepElement)

Returns whether `p` is a subset of `h`, i.e. whether `h` supports the polyhedron `p`.
"""
function Base.issubset(p::Polyhedron{T}, h::HRepElement{U}, solver=Polyhedra.linear_objective_solver(p); tol = _default_tol(T, U)) where {T,U}
    if vrepiscomputed(p)
        _hinv(h, p; tol)
    else
        _hinh(h, p, solver; tol)
    end
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

function _pushinout!(ins, out, pr::Union{AbstractVector, Ray}, h::HalfSpace)
    if pr in h
        push!(ins, pr)
    else
        push!(out, pr)
    end
end

"""
    intersect(v::VRepresentation{T}, h::HRepElement)

Compute the intersection of `v` with an halfspace or hyperplane `h`.
The method used by default is to keep the V-representation element of `v`
that are in `h` and add new ones generated as the intersection between
the hyperplane defining `h` and the segment between two adjacent
V-representation elements of `v` that are in either sides of the hyperplane.
See Lemma 3 of [FP96] for more detail on the method.

[FP96] Fukuda, K. and Prodon, A.
**Double description method revisited**
*Combinatorics and computer science*, *Springer*, **1996**, 91-111
"""
function Base.intersect(v::VRepresentation{T}, h::HRepElement) where {T}
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
        ap = _dot(h.a, p)
        for q in pout
            aq = _dot(h.a, q)
            λ = (aq - h.β) / (aq - ap)
            @assert 0 <= λ <= 1
            push!(pinp, simplify(λ * p + (1 - λ) * q))
        end
    end
    # Similar but with rays
    for r in rins
        ar = _dot(h.a, r)
        for s in rout
            as = _dot(h.a, s)
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
            ap = _dot(h.a, p)
            for r in rs
                ar = _dot(h.a, r)
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
    similar(v, psout, linetype(eltype(rsout))[], rsout)
end

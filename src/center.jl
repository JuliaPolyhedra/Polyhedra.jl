export maximum_radius_with_center, hchebyshevcenter, vchebyshevcenter, chebyshevcenter

"""
    maximum_radius_with_center(h::HRep, center)

Return the maximum radius `r` such that the Euclidean ball of center `center`
and radius `r` is included in the polyhedron `h`.
"""
function maximum_radius_with_center(h::HRep{T}, center) where T
    if any(hp -> !isapproxzero(coord(hp)), hyperplanes(h))
        return zero(T)
    else
        radius = nothing
        for hs in halfspaces(h)
            n = norm(hs.a, 2)
            if iszero(n)
                if hs.β < 0
                    # The polyhedron is empty
                    return zero(T)
                end
            else
                new_radius = (hs.β - _dot(center, hs.a)) / n
                if radius === nothing
                    radius = new_radius
                else
                    radius = min(radius, new_radius)
                end
            end
        end
        if radius === nothing
            error("The polyhedron is the full space, its maximum radius is infinite.")
        elseif radius < 0
            error("The polyhedron does not contain the provided center $center.")
        end
        return radius
    end
end

_shrink(h::HyperPlane, radius, T::Type) = convert(similar_type(typeof(h), T), h)
_shrink(h::HalfSpace, radius, T::Type) = HalfSpace{T}(h.a, h.β - norm(h.a, 2) * radius)
function _shrink(p::HRep, radius)
    # Chebyshev center only works with `Float64` since it uses a `Float64`
    # solver (we may change that when there is interest in using a `BigFloat`
    # solver. `radius` is of type `Float64` so `T = typeof(radius)` is adequate.
    T = typeof(radius)
    f = (i, h) -> _shrink(h, radius, T)
    d = FullDim(p)
    return similar(p, d, T, hmap(f, d, T, p)...)
end

# TODO reference thesis where proper chebyshev center is defined
"""
    hchebyshevcenter(p::HRep[, solver]; linearity_detected=false, proper=true)

Return a tuple with the center and radius of the largest euclidean ball contained in the polyhedron `p`.
Throws an error if the polyhedron is empty or if the radius is infinite.
Linearity is detected first except if `linearity_detected`.

Note that a polytope may have several Chebyshev center.
In general, the set of Chebyshev center of a polytope `p` is a polytope which has a lower dimension than `p` if `p` has a positive dimension.
For instance, if `p` is the rectangle `[-2, 2] x [-1, 1]`, the Chebyshev radius of `p` is 1
and the set of Chebyshev centers is `[-1, 1] x {0}`.
The *proper* Chebyshev center is `(0, 0)`, the Chebyshev center of `[-1, 1] x {0}`.
If `!proper` then any Chebyshev center is returned (the one returned depends on the solver).
Otherwise the proper Chebyshev center is computed.
The proper Chebyshev center is defined by induction on the dimension of `p`.
If `p` has dimension 0 then it is a singleton and its proper Chebyshev center
  is the only element of `p`.
Otherwise, the dimension of the set `q` of Chebyshev centers of `p` is smaller than
the dimension of `p` and the proper Chebyshev center of `p` is the proper Chebyshev center of `q`.
"""
function hchebyshevcenter(p::HRepresentation, solver=default_solver(p; T=Float64); # Need Float64 for `norm(a, 2)`
                          linearity_detected=false, # workaround for https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/25
                          proper=true, verbose=1)
    if !linearity_detected
        p = detecthlinearity(p, solver)
    end
    model, T = layered_optimizer(solver)
    c = MOI.add_variables(model, fulldim(p))
    _constrain_in(model, hyperplanes(p), c, T)
    rs, _ = MOI.add_constrained_variables(model, MOI.Nonnegatives(1))
    r = rs[1]
    for (i, hs) in enumerate(halfspaces(p))
        func, set = _constrain_in_func_set(hs, c, T)
        push!(func.terms, MOI.ScalarAffineTerm{T}(norm(hs.a, 2), r))
        MOI.add_constraint(model, func, set)
    end
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.VariableIndex}(), r)
    MOI.optimize!(model)
    term = MOI.get(model, MOI.TerminationStatus())
    if term ∉ [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        if term == MOI.INFEASIBLE
            error("An empty polyhedron has no H-Chebyshev center.")
        elseif term == MOI.DUAL_INFEASIBLE
            error("The polyhedron contains a Euclidean ball of arbitrarily large radius.")
        else
            _unknown_status(model, term, "computing the H-Chebyshev center.")
        end
    end
    radius = MOI.get(model, MOI.VariablePrimal(), r)
    if proper
        q = detecthlinearity(_shrink(p, radius), solver)
        p_dim = dim(p, true)
        q_dim = dim(q, true)
        if !iszero(q_dim)
            if q_dim >= p_dim
                error("The dimension of the set of Chebyshev centers `$q` is `$q_dim` while we expect it to have a smaller dimension than the original polyhedron which has dimension `$p_dim`.")
            end
            if verbose >= 1
                println("The set `q` of Chebyshev centers has dimension `$q_dim` which is nonzero but lower than the previous dimension `$p_dim`: we now compute the set of Chebyshev center of `q`.")
            end
            center, _ = hchebyshevcenter(q, solver, linearity_detected=true)
            return center, radius
        end
    end
    center = MOI.get(model, MOI.VariablePrimal(), c)
    return center, radius
end
function hchebyshevcenter(p::Polyhedron, solver=default_solver(p; T=Float64); linearity_detected=false, kws...)
    if !linearity_detected
        # We are going to detect it so we might as well do it at the level
        # of the `Polyhedron` so that it is saved and not recomputed if it is needed later.
        detecthlinearity!(p, solver) # FIXME `solver` was forced to be for `T=Float64` but it might not be the best choice here
    end
    return hchebyshevcenter(hrep(p), solver; linearity_detected=true, kws...)
end


# TODO solver here should not be VRepOptimizer
"""
    vchebyshevcenter(p::VRep[, solver])

Return a tuple with the center and radius of the smallest euclidean ball containing the polyhedron `p`.
Throws an error if the polyhedron is empty or if the radius is infinite (i.e. `p` is not a polytope, it contains rays).
"""
function vchebyshevcenter(p::VRep, solver=Polyhedra.default_solver(p))
    error("TODO")
end

"""
    chebyshevcenter(p::Rep[, solver])

If `p` is a H-representation or is a polyhedron for which the H-representation has already been computed, calls `hchebyshevcenter`, otherwise, call `vchebyshevcenter`.
"""
function chebyshevcenter(p::Polyhedron, solver=Polyhedra.default_solver(p; T=Float64); kws...)
    if hrepiscomputed(p)
        hchebyshevcenter(p, solver; kws...)
    else
        vchebyshevcenter(p, solver; kws...)
    end
end
chebyshevcenter(p::HRepresentation, solver=Polyhedra.default_solver(p; T=Float64)) = hchebyshevcenter(p, solver)
chebyshevcenter(p::VRepresentation, solver=Polyhedra.default_solver(p; T=Float64)) = vchebyshevcenter(p, solver)

#######################
# Linearity detection #
#######################

export dim, detecthlinearity!, detecthlinearity, detectvlinearity!, detectvlinearity, detect_new_linearities

"""
    dim(h::HRep, current=false)

Returns the dimension of the affine hull of the polyhedron.
That is the number of non-redundant hyperplanes that define it.
If `current` is `true` then it simply returns the dimension according the current number of hyperplanes, assuming that the H-linearity has already been detected.
Otherwise, it first calls [`detecthlinearity!`](@ref).
"""
function dim(h::HRep, current=false; kws...)
    if !current
        detecthlinearity!(h; kws...)
    end
    fulldim(h) - nhyperplanes(h)
end

"""
    detecthlinearity!(p::HRep)

Detects all the hyperplanes contained in the H-representation and remove all redundant hyperplanes.

## Examples
The representation
```julia
h = HalfSpace([1, 1], 1]) ∩ HalfSpace([-1, -1], -1)
```
contains the hyperplane `HyperPlane([1, 1], 1)`.
"""
detecthlinearity!(p::HRep) = error("detecthlinearity! not implemented for $(typeof(p))")
function detecthlinearity!(p::Polyhedron, solver=default_solver(p); kws...)
    if hredundancy(p) == UNKNOWN_REDUNDANCY
        sethrep!(p, detecthlinearity(hrep(p), solver; kws...), LINEARITY_DETECTED)
    end
end

"""
    detectvlinearity!(p::VRep)

Detects all the lines contained in the V-representation and remove all redundant lines.

## Examples
The representation
```julia
v = conichull([1, 1], [-1, -1])
```
contains the line `Line([1, 1])`.
"""
detectvlinearity!(p::VRep) = error("detectvlinearity! not implemented for $(typeof(p))")
function detectvlinearity!(p::Polyhedron, solver=default_solver(p); kws...)
    if vredundancy(p) == UNKNOWN_REDUNDANCY
        setvrep!(p, detectvlinearity(vrep(p), solver; kws...), LINEARITY_DETECTED)
    end
end

# No ray so no line
function detectvlinearity!(p::VPolytope) end
detectvlinearity(p::VPolytope) = p

linearize(h::HalfSpace) = hyperplane(h)
linearize(r::Ray) = line(r)
_aff_push!(aff::HAffineSpace, h) = intersect!(aff, h)
_aff_push!(aff::VLinearSpace, l) = convexhull!(aff, l)

function _detect_opposite_element(aff, non_opposite, element; kws...)
    element = remproj(element, aff)
    if !isapproxzero(element; kws...) && !any(el -> el ≈ element, non_opposite)
        lin = linearize(element)
        i = findfirst(el -> linearize(el) ≈ lin, non_opposite)
        if i === nothing
            push!(non_opposite, element)
        else
            deleteat!(non_opposite, i)
            _aff_push!(aff, lin)
        end
        return i !== nothing
    else
        return false
    end
end
function _detect_opposite_elements(aff, non_opposite, elements; kws...)
    newlin = true
    for i in 1:fulldim(aff) # could use `while newlin` but `for`-loop is safer.
        newlin || break
        newlin = false
        empty!(non_opposite)
        # Project each ray/halfspace orthogonal to the line hull/affine hull.
        # Remove rays/halfspaces that become zero and detect new lines/hyperplanes
        # with rays/halfspaces that becomes opposite to each other.
        for element in elements
            newlin |= _detect_opposite_element(aff, non_opposite, element; kws...)
        end
    end
end

"""
    detect_new_linearities(rep::HRepresentation, solver; verbose=0)

Given a polyhedron with H-representation `rep`, detect whether a new hyperplane can be generated from the halfspaces in `halfspaces` using an linear program solved by `solver`.
The method is similar to the method used for lines described as follows.
This function is automatically called by `removehredundancy` if a solver is provided.

    detect_new_linearities(rep::VRepresentation, solver; verbose=0)

Given a cone defined by the V-representation `rep` (ignoring the points in the representation if any), detect whether a new line can be generated from the rays in `rays` using an linear program solved by `solver`.
The method is as follows (suppose `lines` is empty for simplicity).
This function is automatically called by `removevredundancy` if a solver is provided.

If there was a line `l` in the cone, it would mean that there exist `μ >= 0` and `ν >= 0` such that
`Σ μ_i r_i = l` and `Σ ν_i r_i = -l`. We deduce from this that `Σ λ_i r_i = 0` where `λ = μ + ν`.

Conversely, if there are `λ >= 0` such that `Σ λ_i r_i = 0` then let `j` be the index of `λ` with largest magnitude
(to make sure it is nonzero).
We have `Σ_{i != j} λ_i/λ_j r_i = -r_j`. As both `r_j` and `-r_j` are in the cone,
`r_j` generates a line in the cone.
However, this means that we now have `Σ_{i != j} λ_i/λ_j r_i ≡ 0 (mod r_j)` so if there is another `λ_i`
with nonzero value, we can transform it to a line as well.
In summary, we have a line `r_i` for each `i` such that `λ_i != 0`.

The dual program is:
max z
    r_i'x >= z
When the primal is feasible, the dual program may still be feasible.
We know that `z = 0` by strong duality as the objective value needs to be equal to the objective of the primal which is zero.
So the constraints are `r_i'x >= 0`. If we have `r_i'x > 0` for some `i`, it means that `-r_i` does not belong to the cone
hence `r_i` can be dropped for the purpose of searching for lines.

## Note

In CDDLib, the dual program is solved, if the objective value is zero then linearity are found
by, for each `i` such that `r_i'x = 0`, solve an LP to find whether `-r_i` belongs to the cone.
CDDLib ignores the primal results provided in `λ` which directly gives linearity without the need
to solve an LP for each ray.
This method is therefore significantly more efficient as it's complexity is `O(dimension of linespace)` which is upper
bounded by `O(fulldim)` while the method of CDDLib is `O(number of rays)`.
"""
function detect_new_linearities(rep::Representation, solver; verbose=0, kws...)
    lins = _linearity(rep)
    nonlins = _nonlinearity(rep)
    isempty(nonlins) && return Int[]
    model, T = layered_optimizer(solver)
    verbose >= 2 && (_model = JuMP.direct_model(model))
    is_lin = falses(length(nonlins))
    active = trues(length(nonlins))
    # We pass `true` as we break homogeneity of the problem with `sum(λ) == 1`.
    hull = _zero_hull(rep, T)
    λ, cλ, sum_con = _hull(model, T, hull, rep, eachindex(nonlins), true)
    if verbose >= 2
        for (i, λ) in enumerate(λ)
            MOI.set(model, MOI.VariableName(), λ, "λ[$i]")
        end
    end
    if !isempty(lins)
        η, _, _ = _hull(model, T, hull, rep, eachindex(lins))
    end
    if rep isa HRepresentation
        MOI.add_constraint.(model, hull[1:end-1], MOI.EqualTo(zero(T)))
        β_con = MOI.add_constraint(model, hull[end], MOI.LessThan(zero(T)))
    else
        MOI.add_constraint.(model, hull, MOI.EqualTo(zero(T)))
    end
    for i in 1:fulldim(rep) # safer than `while ...`
        (count(is_lin) == length(is_lin) || iszero(count(active))) && break
        verbose >= 2 && println(_model)
        # FIXME stopping when we have enough hyperplanes to prove that it's empty
        #       should also resolve the issue with presolve.
        is_feasible(model, "detecting new linearity (you may need to activate presolve for some solvers, e.g. by replacing `GLPK.Optimizer` with `optimizer_with_attributes(GLPK.Optimizer, \"presolve\" => GLPK.GLP_ON)`).") || break
        if rep isa HRepresentation && !isapproxzero((β_primal = MOI.get(model, MOI.ConstraintPrimal(), β_con);); kws...)
            verbose >= 1 && @info("The polyhedron is empty as $β_primal is negative.")
            return nothing
        end
        ray_to_line = Int[]
        ray_to_drop = Int[]
        # If two halfspaces make a hyperplane with a small error, then their `λ` can be chosen very large
        # so as to cancel each others and create another halfspace with the error.
        # However, the `λ` for this other halfspace will be quite small in comparison so we need to
        # scale to detect that.
        λ_l∞ = maximum(abs(MOI.get(model, MOI.VariablePrimal(), λ[i])) for i in eachindex(λ) if active[i])
        for i in eachindex(λ)
            if active[i] && !is_lin[i]
                if !isapproxzero((primal = MOI.get(model, MOI.VariablePrimal(), λ[i]);) / λ_l∞; kws...)
                    # `λ_i > 0`, we know that `-r_i` belongs to the cone so we transform the ray into a line.
                    verbose >= 1 && @info("$(i)th element is linear as $primal is positive.")
                    is_lin[i] = true
                    push!(ray_to_line, i)
                elseif !isapproxzero((dual = MOI.get(model, MOI.ConstraintDual(), cλ[i]);); kws...)
                    # `r_i'x > 0`, we know that `r_i` does not belong to the cone so we just drop the ray from the search for lines.
                    verbose >= 1 && @info("$(i)th element is nonlinear as $dual is positive.")
                    active[i] = false
                    push!(ray_to_drop, i)
                end # otherwise, we are still uncertain about this ray and we'll do a new solve.
            end
        end
        # Query `primal` and `dual` before doing any `MOI.delete` and `MOI.modify` as they
        # won't be available afterwards.
        for i in ray_to_line
            MOI.delete(model, cλ[i])
            MOI.modify(model, sum_con, MOI.ScalarCoefficientChange(λ[i], zero(T)))
        end
        for i in ray_to_drop
            MOI.delete(model, λ[i])
        end
    end
    return findall(is_lin)
end

_no_nonlin_word(::HRepresentation) = "affine"
_no_nonlin_word(::VRepresentation) = "bounded"
_lin_word(::HRepresentation) = "halfspaces"
_lin_word(::VRepresentation) = "lines"
_nonlin_word(::HRepresentation) = "hyperplanes"
_nonlin_word(::VRepresentation) = "rays"
_nonlin_type(h::HRepresentation) = halfspacetype(h)
_nonlin_type(v::VRepresentation) = raytype(v)
_hasnonlinearity(h::HRepresentation) = hashalfspaces(h)
_hasnonlinearity(v::VRepresentation) = hasrays(v)
_nonlinearity(h::HRepresentation) = halfspaces(h)
_nonlinearity(v::VRepresentation) = rays(v)
_linearity_space(h::HRepresentation, current) = affinehull(h, current)
_linearity_space(v::VRepresentation, current) = linespace(v, current)

struct OppositeMockOptimizer end
function _detect_linearity(rep::Representation, solver; kws...)
    aff = _linearity_space(rep, true)
    if _hasnonlinearity(rep)
        if solver === nothing
            @warn("""
Cannot detect exact linearity as no solver was provided and the polyhedron is not $(_no_nonlin_word(rep)).
As fallback, we will only detect $(_lin_word(rep)) from opposite $(_nonlin_word(rep)) but that may not detect all $(_lin_word(rep)).
Set a solver if you believe that the polyhedron may have more linearity.
""" * NO_SOLVER_HELP)
            solver = OppositeMockOptimizer
        end
        if solver == OppositeMockOptimizer
            els = _nonlin_type(rep)[]
            _detect_opposite_elements(aff, els, _nonlinearity(rep); kws...)
        else
            new_lins = detect_new_linearities(rep, solver; kws...)
            if new_lins === nothing
                empty!(aff)
                els = eltype(_nonlinearity(rep))[]
            else
                els = _nonlinearity(rep)
                if !isempty(new_lins)
                    for idx in collect(eachindex(els))[new_lins]
                        _aff_push!(aff, linearize(get(rep, idx)))
                    end
                end
                els = elements_without(els, new_lins)
            end
        end
    else
        els = _nonlinearity(rep)
    end
    return removeduplicates(aff), els
end
function detecthlinearity(hr::HRepresentation, solver; kws...)
    aff, hs = _detect_linearity(hr, solver; kws...)
    typeof(hr)(FullDim(hr), aff.hyperplanes, hs)
end
function detectvlinearity(vr::VRepresentation, solver; kws...)
    aff, rays = _detect_linearity(vr, solver; kws...)
    typeof(vr)(FullDim(vr), preps(vr)..., aff.lines, rays)
end

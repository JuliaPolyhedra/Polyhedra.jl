######################
# Redundancy removal #
######################

# Redundancy
export detecthlinearity!, detectvlinearity!, dim
export isredundant, removevredundancy!, removevredundancy, removehredundancy!, gethredundantindices, getvredundantindices

"""
    detecthlinearity!(p::VRep)

Detects all the hyperplanes contained in the H-representation and remove all redundant hyperplanes.

## Examples
The representation
```julia
h = HalfSpace([1, 1], 1]) ∩ HalfSpace([-1, -1], -1)
```
contains the hyperplane `HyperPlane([1, 1], 1)`.
"""
detecthlinearity!(p::HRep) = error("detecthlinearity! not implemented for $(typeof(p))")
detecthlinearity!(p::Polyhedron) = sethrep!(p, removeduplicates(hrep(p)))

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
detectvlinearity!(p::Polyhedron) = setvrep!(p, removeduplicates(vrep(p)))
function detectvlinearity!(p::VPolytope) end # No ray so no line

"""
    dim(h::HRep, current=false)

Returns the dimension of the affine hull of the polyhedron.
That is the number of non-redundant hyperplanes that define it.
If `current` is `true` then it simply returns the dimension according the current number of hyperplanes, assuming that the H-linearity has already been detected.
Otherwise, it first calls [`detecthlinearity!`](@ref).
"""
function dim(h::HRep, current=false)
    if !current
        detecthlinearity!(h)
    end
    fulldim(h) - nhyperplanes(h)
end

"""
    isredundant(p::Rep, idx::Index; strongly=false)

Return a `Bool` indicating whether the element with index `idx` can be removed without changing the polyhedron represented by `p`.
If `strongly` is `true`,
* if `idx` is an H-representation element `h`, it returns `true` only if no V-representation element of `p` is in the hyperplane of `h`.
* if `idx` is a V-representation element `v`, it returns `true` only if `v` is in the relative interior of `p`.
"""
function isredundant end

"""
    removehredundancy!(p::HRep)

Removes the elements of the H-representation of `p` that can be removed without changing the polyhedron represented by `p`. That is, it only keeps the halfspaces corresponding to facets of the polyhedron.
"""
function removehredundancy! end
function removehredundancy!(p::Polyhedron)
    detectvlinearity!(p)
    detecthlinearity!(p)
    sethrep!(p, removehredundancy(hrep(p), vrep(p)))
end

function is_feasible(model, message)
    @assert JuMP.objective_sense(model) == MOI.FEASIBILITY_SENSE
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    # The objective sense is `MOI.FEASIBILITY_SENSE` so
    # `INFEASIBLE_OR_UNBOUNDED` means infeasible because it cannot be
    # unbounded
    if status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
        return false
    elseif status == MOI.OPTIMAL
        return true
    else
        error("Solver returned $term when ", mesage)
    end
end

"""
    removevredundancy!(p::VRep; strongly=false)

Removes the elements of the V-representation of `p` that can be removed without
changing the polyhedron represented by `p`. That is, it only keeps the extreme
points and rays. This operation is often called "convex hull" as the remaining
points are the extreme points of the convex hull of the initial set of points.
If `strongly=true`, weakly redundant points, i.e., points that are not extreme
but are not in the relative interior either, may be kept.
"""
function removevredundancy! end
function removevredundancy!(p::Polyhedron; strongly=false)
    solver = nothing
    if !strongly && !hrepiscomputed(p) && supportssolver(typeof(p))
        solver = default_solver(p)
        if solver === nothing
            @warn("`removevredundancy` will trigger the computation of the" *
                  " H-representation, which is computationally demanding" *
                  " because no solver was provided to the library. If this is" *
                  " expected, call `computehrep!` explicitely before calling" *
                  " this function to remove this warning.")
        end
    end
    if solver === nothing
        detecthlinearity!(p)
        detectvlinearity!(p)
        setvrep!(p, removevredundancy(vrep(p), hrep(p), strongly=strongly))
    else
        setvrep!(p, removevredundancy(vrep(p), solver))
    end
end

function nonredundant_elements(vr::VRepresentation, model, hull, z, elements, indices, λ)
    _fill_hull(hull, λ, vr, indices)
    crefs = @constraint(model, z .== hull)
    for idx in indices
        if JuMP.has_lower_bound(λ[idx])
            lb = JuMP.lower_bound(λ[idx])
            JuMP.delete_lower_bound(λ[idx])
        else
            lb = nothing
        end
        JuMP.fix(λ[idx], 0.0)
        element = get(vr, idx)
        for (zi, xi) in zip(z, coord(element))
            JuMP.fix(zi, xi)
        end
        if !isone(length(indices)) &&
            is_feasible(model, "attempting to determine whether $element of index $index is redundant.")
            lb = nothing
            JuMP.delete(model, λ[idx])
            for aff in hull
                # TODO Add a function in JuMP.AffExpr API for that instead of using this hack
                delete!(aff.terms, λ[idx])
            end
        else
            push!(elements, element)
            JuMP.unfix(λ[idx])
        end
        if lb !== nothing
            JuMP.set_lower_bound(λ[idx], 0.0)
        end
    end
    for cref in crefs
        JuMP.delete(model, cref)
    end
end
function nonredundant_lines(vr::VPolytope, model, hull, z)
    return tuple()
end
function nonredundant_lines(vr::VRepresentation, model, hull, z)
    nr_lines = linetype(vr)[]
    if haslines(vr)
        lines_idx = collect(eachindex(lines(vr)))
        λ = @variable(model, [lines_idx])
        nonredundant_elements(vr, model, hull, z, nr_lines, lines_idx, λ)
    end
    return (nr_lines,)
end
function nonredundant_rays(vr::VLinearSpace, model, hull, z)
    return tuple()
end
function nonredundant_rays(vr::VPolytope, model, hull, z)
    return tuple()
end
function nonredundant_rays(vr::VRepresentation, model, hull, z)
    nr_rays = raytype(vr)[]
    if hasrays(vr)
        rays_idx = collect(eachindex(rays(vr)))
        λ = @variable(model, [rays_idx], lower_bound = 0.0)
        nonredundant_elements(vr, model, hull, z, nr_rays, rays_idx, λ)
    end
    return (nr_rays,)
end
function nonredundant_points(vr::VCone, model, hull, z)
    return tuple()
end
function nonredundant_points(vr::VRepresentation, model, hull, z)
    nr_points = pointtype(vr)[]
    if haspoints(vr)
        points_idx = collect(eachindex(points(vr)))
        λ = @variable(model, [points_idx], lower_bound = 0.0)
        @constraint(model, sum(λ) == 1)
        nonredundant_elements(vr, model, hull, z, nr_points, points_idx, λ)
    end
    return (nr_points,)
end


removevredundancy(vr::VEmptySpace, solver::Solver) = vr

"""
    removevredundancy(vr::VRepresentation)

Return a V-representation of the polyhedron represented by `vr` all the
elements of `vr` except the redundant ones, i.e. the elements that can
be expressed as convex combination of other ones.
"""
function removevredundancy(vr::VRepresentation, solver::Solver)
    model = ParameterJuMP.ModelWithParams(solver)
    z = ParameterJuMP.add_parameters(model, zeros(fulldim(vr)))
    hull = [zero(JuMP.AffExpr) for i in 1:fulldim(vr)]
    nr_lines = nonredundant_lines(vr, model, hull, z)
    nr_rays = nonredundant_rays(vr, model, hull, z)
    nr_points = nonredundant_points(vr, model, hull, z)
    return vrep(nr_points..., nr_lines..., nr_rays...; d=FullDim(vr))
end
function _fill_hull(hull::Vector{JuMP.AffExpr}, λ, p, idxs)
    for idx in idxs
        x = get(p, idx)
        c = coord(x)
        for i in eachindex(hull)
            JuMP.add_to_expression!(hull[i], c[i], λ[idx])
        end
    end
end

function _filter(f, it)
    # FIXME returns a Vector{Any}
    #collect(Iterators.filter(f, it)) # filter does not implement length so we need to collect
    ret = eltype(it)[]
    for el in it
        if f(el)
            push!(ret, el)
        end
    end
    ret
end
function removevredundancy(vrepit::VIt, hrep::HRep; nl=nlines(hrep), kws...)
    _filter(v -> !isredundant(hrep, v; nl=nl, kws...), vrepit)
end

# Remove redundancy in the V-representation using the H-representation
# There shouldn't be any duplicates in hrep for this to work
function removevredundancy(vrep::VRep, hrep::HRep; kws...)
    nl = nlines(vrep)
    typeof(vrep)(
        FullDim(vrep),
        removevredundancy.(vreps(vrep), hrep; nl=nl, kws...)...
    )::typeof(vrep)
end

function removehredundancy(hrepit::HIt, vrep::VRep; strongly=false, d=dim(vrep))
    _filter(h -> !isredundant(vrep, h, strongly=strongly, d=d), hrepit)
end

# Remove redundancy in the H-representation using the V-representation
# There shouldn't be any duplicates in vrep for this to work
function removehredundancy(hrep::HRep, vrep::VRep; strongly=false)
    R = BitSet()
    d = dim(hrep, true) # TODO dim(hrep)
    typeof(hrep)(FullDim(hrep),
                 removehredundancy.(hreps(hrep), vrep,
                                    strongly=strongly, d=d)...)
end


#function gethredundantindices(hrep::HRep; strongly=false, solver=Polyhedra.solver(hrep))
#    red = BitSet()
#    for (i, h) in enumerate(hreps(hrep))
#        if ishredundant(hrep, h; strongly=strongly, solver=solver)
#            push!(red, i)
#        end
#    end
#    red
#end
#function getvredundantindices(vrep::VRep; strongly = true, solver=Polyhedra.solver(vrep))
#    red = BitSet()
#    for (i, v) in enumerate(vreps(vrep))
#        if isvredundant(vrep, v; strongly=strongly, solver=solver)
#            push!(red, i)
#        end
#    end
#    red
#end

# V-redundancy
# If p is an H-representation, nl needs to be given otherwise if p is a Polyhedron, it can be asked to p.
# TODO nlines should be the number of non-redundant lines so something similar to dim
function isredundant(p::HRep{T}, v::Union{AbstractVector, Line, Ray}; strongly = false, nl::Int=nlines(p), solver=nothing) where {T}
    # v is in every hyperplane otherwise it would not be valid
    hcount = nhyperplanes(p) + count(h -> v in hyperplane(h), halfspaces(p))
    strong = (isray(v) ? fulldim(p)-1 : fulldim(p)) - nl
    return hcount < (strongly ? min(strong, 1) : strong)
end
# A line is never redundant but it can be a duplicate
isredundant(p::HRep{T}, v::Line; strongly = false, nl::Int=nlines(p), solver=nothing) where {T} = false

# H-redundancy
# If p is a V-representation, nl needs to be given otherwise if p is a Polyhedron, it can be asked to p.
function isredundant(p::VRep{T}, h::HRepElement; strongly = false, d::Int=dim(p), solver=nothing) where {T}
    checkvconsistency(p)
    hp = hyperplane(h)
    pcount = count(p -> p in hp, points(p))
    # every line is in h, otherwise it would not be valid
    rcount = nlines(p) + count(r -> r in hp, rays(p))
    pcount < min(d, 1) || (!strongly && pcount + rcount < d)
end
# An hyperplane is never redundant but it can be a duplicate
isredundant(p::VRep{T}, h::HyperPlane; strongly = false, d::Int=dim(p), solver=nothing) where {T} = false

# H-redundancy
#function ishredundantaux(p::HRep, a, β, strongly, solver)
#    sol = MPB.linprog(-a, p, solver)
#    if sol.status == :Unbounded
#        false
#    elseif sol.status == :Optimal
#        if _gt(sol.objval, β)
#            false
#        elseif _geq(sol.objval, β)
#            if strongly
#                false
#            else
#                true
#            end
#        else
#            true
#        end
#    end
#end
#function ishredundant(p::Rep, h::HRepElement; strongly = false, solver=Polyhedra.solver(p))
#    if islin(h)
#        sol = ishredundantaux(p, h.a, h.β, strongly, solver)
#        if !sol[1]
#            sol
#        else
#            ishredundantaux(p, -h.a, -h.β, strongly, solver)
#        end
#    else
#        ishredundantaux(p, h.a, h.β, strongly, solver)
#    end
#end

######################
# Duplicates removal #
######################
export removeduplicates

"""
    removeduplicates(rep::Representation)

Removes the duplicates in the Representation.

* In an H-representation, it removes the redundant hyperplanes and it remove an halfspace when it is equal to another halfspace in the affine hull.
  For instance, `HalfSpace([1, 1], 1)` is equal to `HalfSpace([1, 0], 0)` in the affine hull generated by `HyperPlane([0, 1], 1])`.
* In a V-representation, it removes the redundant lines and it remove a point (resp. ray) when it is equal to another point (resp. ray) in the line hull.
  For instance, in the line hull generated by `Line([0, 1])`, `[1, 1]` is equal to `[1, 0]` and `Ray([2, 2])` is equal to `Ray([1, 0])`.
"""
function removeduplicates end

# H-duplicates
# Separate function so that it is compiled with a concrete type for p
#function vpupdatedup!(aff, points, sympoints, p::SymPoint)
#    found = false
#    # sympoints are treated before points so there shouldn't be any
#    @assert isempty(points)
##    for (i, q) in enumerate(points)
##        if (coord(p) - q) in aff || (coord(p) + q) in aff
##            found = true
##            deleteat!(points, i)
##            push!(sympoints, p)
##            break
##        end
##    end
#    if !found && !any(sp -> (coord(sp) - coord(p)) in aff || (coord(sp) + coord(p)) in aff, sympoints)
#        push!(sympoints, p)
#    end
#end
function vpupdatedup!(aff, points, p::AbstractVector)
    if !any(point -> (point - p) in aff, points)
        push!(points, p)
    end
end
#function vrupdatedup!(rays, lines, l::Line)
#    if !any(isapprox.(lines, [l]))
#        for (i, r) in enumerate(rays)
#            if line(r) ≈ l
#                found = true
#                deleteat!(rays, i)
#                break
#            end
#        end
#        push!(lines, l)
#    end
#end
function vrupdatedup!(aff::VLinearSpace, rays::Vector{<:Ray}, r)
    r = remproj(r, aff)
    if !isapproxzero(r) && !any(ray -> remproj(ray, aff) ≈ r, rays)
        l = line(r)
        found = false
        for (i, s) in enumerate(rays)
            if line(s) ≈ l
                deleteat!(rays, i)
                found = true
                break
            end
        end
        if found
            convexhull!(aff, l)
        else
            push!(rays, r)
        end
        found
    else
        false
    end
end
function premovedups(vrep::VRepresentation, aff::VLinearSpace)
    ps = pointtype(vrep)[]
    for p in points(vrep)
        vpupdatedup!(aff, ps, p)
    end
    tuple(ps)
end
function removeduplicates(vrep::VPolytope)
    typeof(vrep)(FullDim(vrep), premovedups(vrep, emptyspace(vrep))...)
end
function removeduplicates(vrep::VRepresentation)
    aff = linespace(vrep, true)
    newlin = true
    rs = raytype(vrep)[]
    while newlin
        newlin = false
        empty!(rs)
        for r in rays(vrep)
            newlin |= vrupdatedup!(aff, rs, r)
        end
    end
    typeof(vrep)(FullDim(vrep), premovedups(vrep, aff)..., aff.lines, rs)
end

# H-duplicates
#function hupdatedup!(hp, hs, h::HyperPlane)
#    if !any(isapprox.(hp, [h]))
#        for (i, s) in enumerate(hs)
#            if hyperplane(s) ≈ h
#                deleteat!(hs, i)
#                break # There should be no duplicate in hp so no need to continue
#            end
#        end
#        push!(hp, h)
#    end
#end
function hupdatedup!(aff::HAffineSpace, hss, h::HalfSpace)
    h = remproj(h, aff)
    if !isapproxzero(h) && !any(hs -> _isapprox(remproj(hs, aff), h), hss)
        hp = hyperplane(h)
        found = false
        for (i, hs) in enumerate(hss)
            # TODO Not enough, e.g.
            # x <= 1
            # y <= 1
            # x + y >= 2
            if hyperplane(hs) ≈ hp
                deleteat!(hss, i)
                found = true
                break # There should be no duplicate in hp so no need to continue
            end
        end
        if found
            intersect!(aff, hp)
        else
            push!(hss, h)
        end
        found
    else
        false
    end
end
function removeduplicates(hrep::HRepresentation{T}) where {T}
    aff = affinehull(hrep, true)
    newlin = true
    hs = halfspacetype(hrep)[]
    while newlin
        newlin = false
        aff = removeduplicates(aff)
        empty!(hs)
        for h in halfspaces(hrep)
            newlin |= hupdatedup!(aff, hs, h)
        end
    end
    typeof(hrep)(FullDim(hrep), aff.hyperplanes, hs)
end

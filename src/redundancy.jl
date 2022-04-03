######################
# Redundancy removal #
######################

# Redundancy
export isredundant, removevredundancy!, removevredundancy, removehredundancy!, removehredundancy

"""
    @enum Redundancy UNKNOWN_REDUNDANCY LINEARITY_DETECTED NO_REDUNDANCY

Redundancy of a H-representation or V-representation.

* `UNKNOWN_REDUNDANCY`: It is unknown whether there are any undetected linearity or redundancy.
* `LINEARITY_DETECTED`: There are no undetected linearity.
* `NO_REDUNDANCY`: There are no undetected linearity not redundancy.

An *undetected linearity* for a V-representation is a [`Line`](@ref) that is implicit in a conic hull of [`Ray`](@ref)s.
For instance, `conichull([1, 1], [-1, -1])` has undetected linearity because it contains
the [`Line`](@ref) `Line([1, 1])`.
Some undetected linearity is less obvious, e.g., `conichull([1, 0, -1], [0, 1, 1], [-1, -1, 0])`
contains the [`Line`](@ref) `Line([1, 1, 0])` as the sum of the first two [`Ray`](@ref)s is `Ray([1, 1, 0])`.

An *undetected linearity* for a H-representation is a [`HyperPlane`](@ref) that is implicit in an intersection of [`HalfSpace`](@ref)s.
For instance, `HalfSpace([1, 1], 1) ∩ HalfSpace([-1, -1], 1)` has undetected linearity because it contains
the [`HyperPlane`](@ref) `HyperPlane([1, 1], 1)`.
Some undetected linearity is less obvious, e.g., `HalfSpace([1, 0], -1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([-1, -1], 0)`
contains the [`HyperPlane`](@ref) `HyperPlane([1, 1], 0)` as the sum of the first two [`HalfSpace`](@ref)s is `HalfSpace([1, 1], 0)`.
"""
@enum Redundancy UNKNOWN_REDUNDANCY LINEARITY_DETECTED NO_REDUNDANCY

"""
    hredundancy(p::Polyhedron)

Return the [`Redundancy`](@ref) of [`hrep(p)`](@ref `hrep`).
"""
hredundancy(::Polyhedron) = UNKNOWN_REDUNDANCY

"""
    vredundancy(p::Polyhedron)

Return the [`Redundancy`](@ref) of [`vrep(p)`](@ref `hrep`).
"""
vredundancy(::Polyhedron) = UNKNOWN_REDUNDANCY

"""
    isredundant(p::Rep, idx::Index; strongly=false)

Return a `Bool` indicating whether the element with index `idx` can be removed without changing the polyhedron represented by `p`.
If `strongly` is `true`,
* if `idx` is an H-representation element `h`, it returns `true` only if no V-representation element of `p` is in the hyperplane of `h`.
* if `idx` is a V-representation element `v`, it returns `true` only if `v` is in the relative interior of `p`.
"""
function isredundant end

# `isredundant(::Polyhedron, ::HIndex)` would cause ambiguity with LRSLib
# which defines `isredundant(::LRSLib.Polyhedron, ::Index)` so we add a
# `_fallback` layer to circumvent this.
function isredundant(p::Polyhedron, idx::Index; kws...)
    return isredundant_fallback(p, idx; kws...)
end

function isredundant_fallback(p::Polyhedron, idx::HIndex; kws...)
    hredundancy(p) == NO_REDUNDANCY && return false
    solver = _solver_warn(p, false, false)
    if solver !== nothing
        @warn("""
`isredundant` with a solver is not supported yet so it triggers the computation of the V-representation, which
is computationally demanding because no solver was provided to the library.
If this is expected, call `computevrep!` explicitely before calling this
function to remove this warning.
""")
    end
    # `detecthlinearity(p)` might remove redundancy but we get `h` first here.
    h = get(p, idx)
    detectvlinearity!(p)
    detecthlinearity!(p)
    return isredundant(vrep(p), h; d = dim(p), kws...)
end

function isredundant_fallback(p::Polyhedron, idx::VIndex; strongly=false, kws...)
    vredundancy(p) == NO_REDUNDANCY && return false
    solver = _solver_warn(p, true, strongly)
    if solver !== nothing
        @warn("""
`isredundant` with a solver is not supported yet so it triggers the computation of the H-representation, which
is computationally demanding because no solver was provided to the library.
If this is expected, call `computehrep!` explicitely before calling this
function to remove this warning.
""")
    end
    # `detectvlinearity(p)` might remove redundancy but we get `v` first here.
    v = get(p, idx)
    detecthlinearity!(p)
    detectvlinearity!(p)
    return isredundant(hrep(p), h; strongly=strongly, kws...)
end

function _solver_warn(p::Polyhedron, v_or_h::Bool, strongly::Bool)
    computed = v_or_h ? hrepiscomputed(p) : vrepiscomputed(p)
    v = v_or_h ? 'v' : 'h'
    V = uppercase(v)
    h = v_or_h ? 'h' : 'v'
    solver = nothing
    if !strongly && !computed && supportssolver(typeof(p))
        solver = default_solver(p)
        if solver === nothing
            @warn("""
`remove$(v)redundancy!` will trigger the computation of the $V-representation, which
is computationally demanding because no solver was provided to the library.
If this is expected, call `compute$(h)rep!` explicitely before calling this
function to remove this warning.
""" * NO_SOLVER_HELP)
        end
    end
    return solver

end

"""
    removehredundancy!(p::HRep)

Removes the elements of the H-representation of `p` that can be removed without changing the polyhedron represented by `p`. That is, it only keeps the halfspaces corresponding to facets of the polyhedron.
"""
function removehredundancy! end
function removehredundancy!(p::Polyhedron)
    hredundancy(p) == NO_REDUNDANCY && return
    solver = _solver_warn(p, false, false)
    if solver === nothing
        detectvlinearity!(p)
        detecthlinearity!(p)
        nonred = removehredundancy(hrep(p), vrep(p))
    else
        detecthlinearity!(p)
        nonred = removehredundancy(hrep(p), solver)
    end
    sethrep!(p, nonred, NO_REDUNDANCY)
end

"""
    removevredundancy!(p::VRep; strongly=false, planar=true)

Removes the elements of the V-representation of `p` that can be removed without
changing the polyhedron represented by `p`. That is, it only keeps the extreme
points and rays. This operation is often called "convex hull" as the remaining
points are the extreme points of the convex hull of the initial set of points.
If `strongly=true`, weakly redundant points, i.e., points that are not extreme
but are not in the relative interior either, may be kept.
If `fulldim(p)` is 2, `strongly` is `false` and `planar` is `true`, a planar convex hull algorithm
is used.
"""
function removevredundancy! end
function removevredundancy!(p::Polyhedron; strongly=false, planar=true, kws...)
    vredundancy(p) == NO_REDUNDANCY && return
    if fulldim(p) == 2 && !strongly && planar
        setvrep!(p, planar_hull(vrep(p)), NO_REDUNDANCY)
    else
        solver = _solver_warn(p, true, strongly)
        if solver === nothing
            detecthlinearity!(p; kws...)
            detectvlinearity!(p; kws...)
            nonred = removevredundancy(vrep(p), hrep(p), strongly=strongly)
        else
            detectvlinearity!(p; kws...)
            nonred = removevredundancy(vrep(p), solver)
        end
        # If `strongly` then we only remove strongly redundant elements
        # henwe we cannot say that the redundancy is `NO_REDUNDANCY`.
        setvrep!(p, nonred, strongly ? LINEARITY_DETECTED : NO_REDUNDANCY)
    end
end

# Interpolating the string can take more time than the LP solve:
# https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/279
# This `struct` allows to interpolate it lazily only when an error is actually thrown.
struct _RedundantMessage{ET, IT}
    element::ET
    index::IT
end
function Base.show(io::IO, message::_RedundantMessage)
    print(io, "attempting to determine whether `")
    print(io, message.element)
    print(io, "` of index `")
    print(io, message.index)
    print(io, "` is redundant.")
end

function _redundant_indices(rep::Representation, model::MOI.ModelLike, T::Type,
                            hull, indices, clean)
    red_indices = Int[]
    isempty(indices) && return red_indices
    λ, cλ, sum_one = _hull(model, T, hull, rep, indices)
    hull_con = MOI.add_constraint.(model, hull, MOI.EqualTo(zero(T)))
    for (i, idx) in enumerate(indices)
        if cλ === nothing
            fix_con = MOI.add_constraint(model, λ[i], MOI.EqualTo(zero(T)))
        else
            fix_con = MOI.transform(model, cλ[i], MOI.EqualTo(zero(T)))
        end
        element = get(rep, idx)
        for (ci, ai) in zip(hull_con, _coord_for_hull(element))
            MOI.set(model, MOI.ConstraintSet(), ci, MOI.EqualTo{T}(ai))
        end
        #if !isone(length(indices)) &&
        if  is_feasible(model, _RedundantMessage(element, idx))
            MOI.delete(model, λ[i])
            push!(red_indices, i)
        else
            if cλ === nothing
                MOI.delete(model, fix_con)
            else
                cλ[i] = MOI.transform(model, fix_con, MOI.GreaterThan(zero(T)))
            end
        end
    end
    if clean
        MOI.delete(model, hull_con)
        if !isempty(red_indices)
            rm = Set(λ[red_indices])
            hull = map(h -> MOI.Utilities.filter_variables(vi -> !(vi in rm), h), hull)
        end
    end
    return red_indices
end

function elements_without(elements, rm)
    if isempty(rm)
        return elements
    else
        return deleteat!(collect(elements), rm)
    end
end
function _nonredundant(rep, model, T, hull, elements, clean)
    red = _redundant_indices(rep, model, T, hull, eachindex(elements), clean)
    return (elements_without(elements, red),)
end

function nonredundant_hyperplanes(vr::HRepresentation, model, T, hull, clean)
    return _nonredundant(vr, model, T, hull, hyperplanes(vr), clean)
end
nonredundant_halfspaces(vr::HAffineSpace, model, T, hull, clean) = tuple()
function nonredundant_halfspaces(vr::HRepresentation, model, T, hull, clean)
    return _nonredundant(vr, model, T, hull, halfspaces(vr), clean)
end

"""
    removehredundancy(hr::HRepresentation, solver)

Return a H-representation of the polyhedron represented by `hr` with all the
elements of `hr` except the redundant ones, i.e. the elements that can
be expressed as convex combination of other ones.
"""
function removehredundancy(hr::HRep, solver)
    model, T = layered_optimizer(solver)
    hull = _zero_hull(hr, T)
    # TODO It's much more efficient to remove redundant hyperplanes with `removeduplicates`.
    #      Moreover, linearity may have already been detected.
    #      The only advantage over `removeduplicates` is that it does not alter
    #      the hyperplanes, it just remove the redundant ones.
    nr_hyperplanes = nonredundant_hyperplanes(hr, model, T, hull, true)
    nr_halfspaces = nonredundant_halfspaces(hr, model, T, hull, false)
    return hrep(nr_hyperplanes..., nr_halfspaces...; d=FullDim(hr))
end

nonredundant_lines(vr::VPolytope, model, T, hull, clean) = tuple()
function nonredundant_lines(vr::VRepresentation, model, T, hull, clean)
    return _nonredundant(vr, model, T, hull, lines(vr), clean)
end
nonredundant_rays(vr::VLinearSpace, model, T, hull, clean) = tuple()
nonredundant_rays(vr::VPolytope, model, T, hull, clean) = return tuple()
function nonredundant_rays(vr::VRepresentation, model, T, hull, clean)
    return _nonredundant(vr, model, T, hull, rays(vr), clean)
end
nonredundant_points(vr::VCone, model, T, hull, clean) = tuple()
function nonredundant_points(vr::VRepresentation, model, T, hull, clean)
    return _nonredundant(vr, model, T, hull, points(vr), clean)
end

removevredundancy(vr::VEmptySpace, solver) = vr

"""
    removevredundancy(vr::VRepresentation, solver)

Return a V-representation of the polyhedron represented by `vr` with all the
elements of `vr` except the redundant ones, i.e. the elements that can
be expressed as convex combination of other ones.
"""
function removevredundancy(vr::VRepresentation, solver)
    model, T = layered_optimizer(solver)
    hull = _zero_hull(vr, T)
    # TODO It's much more efficient to remove redundant lines with `removeduplicates`.
    #      Moreover, linearity may have already been detected.
    #      The only advantage over `removeduplicates` is that it does not alter
    #      the lines, it just remove the redundant ones.
    nr_lines = nonredundant_lines(vr, model, T, hull, true)
    nr_rays = nonredundant_rays(vr, model, T, hull, true)
    nr_points = nonredundant_points(vr, model, T, hull, false)
    return vrep(nr_points..., nr_lines..., nr_rays...; d=FullDim(vr))
end

_dim_for_hull(h::HRepresentation) = fulldim(h) + 1
_dim_for_hull(v::VRepresentation) = fulldim(v)
function _zero_hull(rep::Representation, T::Type)
    return [zero(MOI.ScalarAffineFunction{T}) for i in 1:_dim_for_hull(rep)]
end
_coord_for_hull(h::HRepElement) = [h.a; h.β]
_coord_for_hull(r::VRepElement) = coord(r)
function _fill_hull(hull::Vector{MOI.ScalarAffineFunction{T}}, λ, rep, idxs) where T
    for (λ, idx) in zip(λ, idxs)
        a = _coord_for_hull(get(rep, idx))
        for i in eachindex(a)
            push!(hull[i].terms, MOI.ScalarAffineTerm{T}(a[i], λ))
        end
    end
end
function _hull(model::MOI.ModelLike, ::Type{T}, hull::Vector{MOI.ScalarAffineFunction{T}}, rep, idxs, sum_one = idxs isa PointIndices) where T
    λ = MOI.add_variables(model, length(idxs))
    if !(idxs isa Union{HyperPlaneIndices, LineIndices})
        cλ = [MOI.add_constraint(model, λ, MOI.GreaterThan(zero(T))) for λ in λ]
    else
        cλ = nothing
    end
    if sum_one
        sum_func = MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.(one(T), λ),
            zero(T)
        )
        sum_con = MOI.add_constraint(model, sum_func, MOI.EqualTo(one(T)))
    else
        sum_con = nothing
    end
    _fill_hull(hull, λ, rep, idxs)
    return λ, cλ, sum_con
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
function _removevred_withhred(vrep::VRep, hrep::HRep; kws...)
    nl = nlines(vrep)
    typeof(vrep)(
        FullDim(vrep),
        removevredundancy.(vreps(vrep), hrep; nl=nl, kws...)...
    )::typeof(vrep)
end
# Split in two methods to avoid ambiguity with `(::VRepresentation, solver::Any)`
function removevredundancy(vrep::Polyhedron, hrep::HRep; kws...)
    _removevred_withhred(vrep, hrep; kws...)
end
function removevredundancy(vrep::VRepresentation, hrep::HRep; kws...)
    _removevred_withhred(vrep, hrep; kws...)
end

function removehredundancy(hrepit::HIt, vrep::VRep; strongly=false, d=dim(vrep))
    _filter(h -> !isredundant(vrep, h, strongly=strongly, d=d), hrepit)
end

# Remove redundancy in the H-representation using the V-representation
# There shouldn't be any duplicates in vrep for this to work
function _removehred_withvred(hrep::HRep, vrep::VRep; strongly=false)
    R = BitSet()
    d = dim(hrep, true) # TODO dim(hrep)
    typeof(hrep)(FullDim(hrep),
                 removehredundancy.(hreps(hrep), vrep,
                                    strongly=strongly, d=d)...)
end
# Split in two methods to avoid ambiguity with `(::HRepresentation, solver::Any)`
function removehredundancy(hrep::Polyhedron, vrep::VRep; kws...)
    _removehred_withvred(hrep, vrep; kws...)
end
function removehredundancy(hrep::HRepresentation, vrep::VRep; kws...)
    _removehred_withvred(hrep, vrep; kws...)
end

function coord_matrix!(A, elements, condition, offset)
    for el in elements
        if condition(el)
            offset += 1
            A[offset, :] = coord(el)
        end
    end
    return offset
end

# V-redundancy
# If p is an H-representation, nl needs to be given otherwise if p is a Polyhedron, it can be asked to p.
# TODO nlines should be the number of non-redundant lines so something similar to dim
function isredundant(p::HRep{T}, v::Union{AbstractVector, Line, Ray}; strongly = false, nl::Int=nlines(p), solver=nothing) where {T}
    # v is in every hyperplane otherwise it would not be valid
    hcount = nhyperplanes(p) + count(h -> v in hyperplane(h), halfspaces(p))
    strong = (isray(v) ? fulldim(p) - 1 : fulldim(p)) - nl
    bound = strongly ? min(strong, 1) : strong
    if hcount < bound
        return true
    else
        A = Matrix{T}(undef, hcount, fulldim(p))
        offset = coord_matrix!(A, hyperplanes(p), _ -> true, 0)
        offset = coord_matrix!(A, halfspaces(p), h -> v in hyperplane(h), offset)
        @assert offset == size(A, 1)
        return rank(A) < bound
    end
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
    if pcount < min(d, 1) || (!strongly && pcount + rcount < d)
        return true
    else
        strongly && return false
        A = Matrix{T}(undef, max(0, pcount - 1) + rcount, fulldim(p))
        offset = 0
        if pcount > 1
            orig = nothing
            for x in points(p)
                if x in hp
                    if orig === nothing
                        orig = x
                    else
                        offset += 1
                        @. A[offset, :] = x - orig
                    end
                end
            end
        end
        if !iszero(rcount)
            offset = coord_matrix!(A, rays(p), r -> r in hp, offset)
        end
        @assert offset == size(A, 1)
        return rank(A) + min(1, pcount) < d
    end
end
# An hyperplane is never redundant but it can be a duplicate
isredundant(p::VRep{T}, h::HyperPlane; strongly = false, d::Int=dim(p), solver=nothing) where {T} = false

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

# V-duplicates
# Function barrier approach:
# Separate function so that it is compiled with a concrete type for p
function vpupdatedup!(aff, points, p::AbstractVector)
    if !any(point -> (point - p) in aff, points)
        push!(points, p)
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
removeduplicates(vrep::VLinearSpace, solver = nothing) = detectvlinearity(vrep, solver)
function removeduplicates(vrep::VRepresentation, solver = nothing)
    aff, rs = _detect_linearity(vrep, solver)
    typeof(vrep)(FullDim(vrep), premovedups(vrep, aff)..., aff.lines, rs)
end

# H-duplicates
removeduplicates(hr::HRepresentation, solver = nothing) = detecthlinearity(hr, solver)

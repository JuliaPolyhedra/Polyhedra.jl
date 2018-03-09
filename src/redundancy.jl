######################
# Redundancy removal #
######################

# Redundancy
export isredundant, removevredundancy!, removehredundancy!, gethredundantindices, getvredundantindices

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

"""
    removevredundancy!(p::VRep)

Removes the elements of the V-representation of `p` that can be removed without changing the polyhedron represented by `p`. That is, it only keeps the extreme points and rays. This operation is often called "convex hull" as the remaining points are the extreme points of the convex hull of the initial set of points.
"""
function removevredundancy! end

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
function removevredundancy(vrepit::VIt, hrep::HRep; strongly=true, nl=nlines(hrep))
    _filter(v -> !isredundant(hrep, v, strongly=strongly, nl=nl), vrepit)
end

# Remove redundancy in the V-representation using the H-representation
# There shouldn't be any duplicates in hrep for this to work
function removevredundancy(vrep::VRep, hrep::HRep; strongly=true)
    nl = nlines(vrep)
    typeof(vrep)(removevredundancy.(vreps(vrep), hrep, strongly=strongly, nl=nl)...)::typeof(vrep) # FIXME return type annotation needed in Julia v0.6.2
end

function removehredundancy(hrepit::HIt, vrep::VRep; strongly=true, d=dim(vrep))
    _filter(h -> !isredundant(vrep, h, strongly=strongly, d=d), hrepit)
end

# Remove redundancy in the H-representation using the V-representation
# There shouldn't be any duplicates in vrep for this to work
function removehredundancy(hrep::HRep, vrep::VRep; strongly=true)
    R = IntSet()
    d = dim(hrep, true) # TODO dim(hrep)
    typeof(hrep)(removehredundancy.(hreps(hrep), vrep, strongly=strongly, d=d)...)
end


#function gethredundantindices(hrep::HRep; strongly=false, solver = defaultLPsolverfor(hrep))
#    red = IntSet()
#    for (i, h) in enumerate(hreps(hrep))
#        if ishredundant(hrep, h; strongly=strongly, solver=solver)
#            push!(red, i)
#        end
#    end
#    red
#end
#function getvredundantindices(vrep::VRep; strongly = true, solver = defaultLPsolverfor(vrep))
#    red = IntSet()
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
function isredundant(p::HRep{N,T}, v::VRepElement; strongly = true, nl::Int=nlines(p), solver = JuMP.UnsetSolver) where {N,T}
    count = nhyperplanes(p) # v is in every hyperplane otherwise it would not be valid
    for h in halfspaces(p)
        if v in hyperplane(h)
            count += 1
        end
    end
    strong = (isray(v) ? N-1 : N) - nl
    count < (strongly ? strong : min(strong, 1))
end
# A line is never redundant but it can be a duplicate
isredundant(p::HRep{N,T}, v::Line; strongly = true, nl::Int=nlines(p), solver = JuMP.UnsetSolver) where {N,T} = false

# H-redundancy
# If p is a V-representation, nl needs to be given otherwise if p is a Polyhedron, it can be asked to p.
function isredundant(p::VRep{N,T}, h::HRepElement; strongly = true, d::Int=dim(p), solver = JuMP.UnsetSolver) where {N,T}
    checkvconsistency(p)
    pcount = 0
    hp = hyperplane(h)
    for v in allpoints(p)
        if v in hp
            pcount += 1
        end
    end
    rcount = nlines(p) # every line is in h, otherwise it would not be valid
    for v in rays(p)
        if v in hp
            rcount += 1
        end
    end
    pcount < min(d, 1) || (strongly && pcount + rcount < d)
end
# An hyperplane is never redundant but it can be a duplicate
isredundant(p::VRep{N,T}, h::HyperPlane; strongly = true, d::Int=dim(p), solver = JuMP.UnsetSolver) where {N,T} = false

function ishredundant(args...; kws...)
    warn("ishredundant is deprecated, use isredundant intead")
    isredundant(args...; kws...)
end
function isvredundant(args...; kws...)
    warn("isvredundant is deprecated, use isredundant intead")
    isredundant(args...; kws...)
end

# H-redundancy
#function ishredundantaux(p::HRep, a, β, strongly, solver)
#    sol = linprog(-a, p, solver)
#    if sol.status == :Unbounded
#        false
#    elseif sol.status == :Optimal
#        if mygt(sol.objval, β)
#            false
#        elseif mygeq(sol.objval, β)
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
#function ishredundant(p::Rep, h::HRepElement; strongly = false, solver = defaultLPsolverfor(p))
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

# H-duplicates
# Separate function so that it is compiled with a concrete type for p
function vpupdatedup!(aff, points, sympoints, p::SymPoint)
    found = false
    for (i, q) in enumerate(points)
        if (coord(p) - q) in aff || (coord(p) + q) in aff
            found = true
            deleteat!(points, i)
            push!(sympoints, p)
            break
        end
    end
    if !found && !any(sp -> (sp - p) in aff || (sp + p) in aff, sympoints)
        push!(sympoints, p)
    end
end
function vpupdatedup!(aff, points, sympoints, p::AbstractPoint)
    if !any(point -> (point - p) in aff, points) && !any(sp -> (coord(sp) - p) in aff || (coord(sp) + p) in aff, sympoints)
        found = false
        for (i, q) in enumerate(points)
            if p + q in aff
                found = true
                deleteat!(points, i)
                push!(sympoints, SymPoint(p))
                break
            end
        end
        if !found
            push!(points, p)
        end
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
function vrupdatedup!(aff::VAffineSpace, rays::Vector{<:Ray}, r)
    r = remproj(r, aff)
    if !myeqzero(r) && !any(ray -> remproj(ray, aff) ≈ r, rays)
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
function premovedups(vrep::VRepresentation, aff::VAffineSpace)
    ps = pointtype(vrep)[]
    sps = sympointtype(vrep)[]
    for p in sympoints(vrep)
        vpupdatedup!(aff, ps, sps, p)
    end
    for p in points(vrep)
        vpupdatedup!(aff, ps, sps, p)
    end
    sps, ps
end
function removeduplicates(vrep::VPolytope)
    typeof(vrep)(premovedups(vrep, emptyspace(vrep))...)
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
    typeof(vrep)(premovedups(vrep, aff)..., aff.lines, rs)
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
    if !myeqzero(h) && !any(hs -> myeq(remproj(hs, aff), h), hss)
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
function removeduplicates(hrep::HRepresentation{N, T}) where {N, T}
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
    typeof(hrep)(aff.hps, hs)
end

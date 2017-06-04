######################
# Redundancy removal #
######################

# Redundancy
export removevredundancy!, removehredundancy!, isvredundant, ishredundant, gethredundantindices, getvredundantindices
removevredundancy!(p::VRep) = error("removevredundancy! not implemented for $(typeof(p))")
removehredundancy!(p::HRep) = error("removehredundancy! not implemented for $(typeof(p))")

# Remove redundancy in the V-representation using the H-representation
# There shouldn't be any duplicates in hrep for this to work
function removevredundancy(vrep::VRep, hrep::HRep; strongly=true)
    R = IntSet()
    nl = nlines(vrep)
    for (i, v) in enumerate(vreps(vrep))
        if isvredundant(hrep, v, strongly=strongly, nl=nl)
            push!(R, i)
        end
    end
    vrep[setdiff(1:nvreps(vrep), R)]
end

# Remove redundancy in the H-representation using the V-representation
# There shouldn't be any duplicates in vrep for this to work
function removehredundancy(hrep::HRep, vrep::VRep; strongly=true)
    R = IntSet()
    d = dim(hrep, true) # TODO dim(hrep)
    for (i, h) in enumerate(hreps(hrep))
        if ishredundant(vrep, h, strongly=strongly, d=d)
            push!(R, i)
        end
    end
    hrep[setdiff(1:nhreps(hrep), R)]
end


function gethredundantindices(hrep::HRep; strongly=false, solver = defaultLPsolverfor(hrep))
    red = IntSet()
    for (i, h) in enumerate(hreps(hrep))
        if ishredundant(hrep, h; strongly=strongly, solver=solver)
            push!(red, i)
        end
    end
    red
end
function getvredundantindices(vrep::VRep; strongly = true, solver = defaultLPsolverfor(vrep))
    red = IntSet()
    for (i, v) in enumerate(vreps(vrep))
        if isvredundant(vrep, v; strongly=strongly, solver=solver)
            push!(red, i)
        end
    end
    red
end

# V-redundancy
# If p is an H-representation, nl needs to be given otherwise if p is a Polyhedron, it can be asked to p.
# TODO nlines should be the number of non-redundant lines so something similar to dim
function isvredundant{N,T}(p::HRep{N,T}, v::VRepElement; strongly = true, nl::Int=nlines(p), solver = JuMP.UnsetSolver)
    count = 0
    for h in hreps(p)
        if v in hyperplane(h)
            count += 1
        end
    end
    strong = (isray(v) ? N-1 : N) - nl
    count < (strongly ? strong : min(strong, 1))
end
# A line is never redundant but it can be a duplicate
isvredundant{N,T}(p::HRep{N,T}, v::Line; strongly = true, nl::Int=nlines(p), solver = JuMP.UnsetSolver) = false

# H-redundancy
# If p is a V-representation, nl needs to be given otherwise if p is a Polyhedron, it can be asked to p.
function ishredundant{N,T}(p::VRep{N,T}, h::HRepElement; strongly = true, d::Int=dim(p), solver = JuMP.UnsetSolver)
    rcount = 0
    pcount = 0
    hp = hyperplane(h)
    for v in vreps(p)
        if v in hp
            if isray(v)
                rcount += 1
            else
                pcount += 1
            end
        end
    end
    if myeqzero(h.β) && !haspoints(p)
        # The origin, see #28
        pcount += 1
    end
    pcount < min(d, 1) || (strongly && pcount + rcount < d)
end
# An hyperplane is never redundant but it can be a duplicate
ishredundant{N,T}(p::VRep{N,T}, h::HyperPlane; strongly = true, d::Int=dim(p), solver = JuMP.UnsetSolver) = false

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
function vpupdatedup!(aff, points, sympoints, p)
    if !any(point -> (point - p) in aff, points) && !any(sp -> (coord(sp) - p) in aff || (coord(sp) + p) in aff, sympoints)
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
function vrupdatedup!(aff, rays, r)
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
            push!(aff, l)
        else
            push!(rays, r)
        end
        found
    else
        false
    end
end
function removeduplicates{N, T}(vrep::VRepresentation{N, T})
    aff = linespace(vrep, true)
    newlin = true
    while newlin
        newlin = false
        rs = Union{Vec{N, T}, Ray{N, T}}[]
        for r in rays(vrep)
            if !islin(r)
                newlin |= vrupdatedup!(aff, rs, r)
            end
        end
    end
    ps = Union{Point{N, T}, AbstractVector{T}}[]
    sympoints = SymPoint{N, T}[]
    for p in points(vrep)
        vpupdatedup!(aff, ps, sympoints, p)
    end
    typeof(vrep)([sympoints; ps], [aff.lines; rs])
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
            push!(aff, hp)
        else
            push!(hss, h)
        end
        found
    else
        false
    end
end
function removeduplicates{N, T}(hrep::HRepresentation{N, T})
    aff = affinehull(hrep, true)
    newlin = true
    while newlin
        newlin = false
        aff = removeduplicates(aff)
        hs = HalfSpace{N, T}[]
        for h in ineqs(hrep)
            newlin |= hupdatedup!(aff, hs, h)
        end
    end
    typeof(hrep)(aff.hps, hs)
end

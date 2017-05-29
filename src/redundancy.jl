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
        if isvredundant(hrep, v, strongly=strongly, nl=nlines(vrep))
            push!(R, i)
        end
    end
    vrep[setdiff(1:nvreps(vrep), R)]
end

function gethredundantindices(hrep::HRep; strongly=false, solver = defaultLPsolverfor(hrep))
    red = IntSet()
    for i in 1:nhreps(hrep)
        if ishredundant(hrep, i; strongly, solver)
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

# H-redundancy
function ishredundantaux(p::Polyhedron, a, β, strongly, cert, solver)
    sol = linprog(-a, p, solver)
    if sol.status == :Unbounded
        cert ?  (false, sol.attrs[:unboundedray], :UnboundedRay) : false
    elseif sol.status == :Optimal
        if mygt(sol.objval, β)
            cert ? (false, sol.sol, :ExteriorPoint) : false
        elseif mygeq(sol.objval, β)
            if strongly
                cert ? (false, sol.sol, :BoundaryPoint) : false
            else
                cert ? (true, sol.sol, :BoundaryPoint) : true
            end
        else
            cert ? (true, nothing, :NotApplicable) : true
        end
    end
end
function ishredundant(p::Rep, h::HRepElement; strongly=false, cert=false, solver = defaultLPsolverfor(p))
    if islin(h)
        sol = ishredundantaux(p, h.a, h.β, strongly, cert, solver)
        if !sol[1]
            sol
        else
            ishredundantaux(p, -h.a, -h.β, strongly, cert, solver)
        end
    else
        ishredundantaux(p, h.a, h.β, strongly, cert, solver)
    end
end

######################
# Duplicates removal #
######################

# Linearity
export detecthlinearities!, detectvlinearities!
detectvlinearities!(p::VRep) = error("detectvlinearities! not implemented for $(typeof(p))")
detecthlinearities!(p::HRep) = error("detecthlinearities! not implemented for $(typeof(p))")

# H-duplicates
# Separate function so that it is compiled with a concrete type for p
function vpupdatedup!(points, sympoints, p::SymPoint)
    found = false
    mp = -p
    for (i, q) in enumerate(points)
        if coord(p) ≈ q || coord(mp) ≈ q
            found = true
            deleteat!(points, i)
            push!(sympoints, p)
            break
        end
    end
    if !found && !any(isapprox.(sympoints, [p]))
        push!(sympoints, p)
    end
end
function vpupdatedup!(points, sympoints, p)
    mp = -p
    if !any(isapprox.(points, [p])) && !any(isapprox.(coord.(sympoints), [p])) && !any(isapprox.(coord.(sympoints), [mp]))
        push!(points, p)
    end
end
function vrupdatedup!(rays, lines, l::Line)
    if !any(isapprox.(lines, [l]))
        for (i, r) in enumerate(rays)
            if line(r) ≈ l
                found = true
                deleteat!(rays, i)
                break
            end
        end
        push!(lines, l)
    end
end
function vrupdatedup!(rays, lines, r)
    l = line(r)
    if !any(isapprox.(rays, [r])) && !any(isapprox.(lines, [l]))
        mr = -r
        found = false
        for (i, s) in enumerate(rays)
            if line(s) ≈ l
                deleteat!(rays, i)
                push!(lines, l)
                found = true
                break
            end
        end
        if !found
            push!(rays, r)
        end
    end
end
function removeduplicates{N, T}(vrep::VRepresentation{N, T})
    ps = Union{Point{N, T}, AbstractVector{T}}[]
    sympoints = SymPoint{N, T}[]
    for p in points(vrep)
        vpupdatedup!(ps, sympoints, p)
    end
    rs = Union{Vec{N, T}, Ray{N, T}}[]
    lines = Line{N, T}[]
    for r in rays(vrep)
        vrupdatedup!(rs, lines, r)
    end
    typeof(vrep)([sympoints; ps], [lines; rs])
end

# H-duplicates
function hupdatedup!(hp, hs, h::HyperPlane)
    if !any(isapprox.(hp, [h]))
        for (i, s) in enumerate(hs)
            if hyperplane(s) ≈ h
                deleteat!(hs, i)
                break # There should be no duplicate in hp so no need to continue
            end
        end
        push!(hp, h)
    end
end
function hupdatedup!(hp, hs, h::HalfSpace)
    if !any(isapprox.(hp, [hyperplane(h)])) && !any(isapprox.(hs, [h]))
        push!(hs, h)
    end
end
function removeduplicates{N, T}(vrep::HRepresentation{N, T})
    hp = HyperPlane{N, T}[]
    hs = HalfSpace{N, T}[]
    for p in hreps(vrep)
        hupdatedup!(hp, hs, p)
    end
    typeof(vrep)(hp, hs)
end

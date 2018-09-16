export hchebyshevcenter, vchebyshevcenter, chebyshevcenter
using JuMP
"""
    hchebyshevcenter(p::HRep[, solver])

Return a tuple with the center and radius of the largest euclidean ball contained in the polyhedron `p`.
Throws an error if the polyhedron is empty or if the radius is infinite.
"""
function hchebyshevcenter(p::HRep, solver=Polyhedra.solver(p))
    m = Model(solver=solver)
    c = @variable m [1:fulldim(p)]
    for hp in hyperplanes(p)
        a = Vector{Float64}(hp.a)
        β = Float64(hp.β)
        @constraint m dot(a, c) == β
    end
    @variable m r[1:nhalfspaces(p)] >= 0
    for (i, hs) in enumerate(halfspaces(p))
        a = Vector{Float64}(hs.a)
        β = Float64(hs.β)
        @constraint m dot(a, c) + r[i] * norm(a, 2) <= β
    end
    @variable m minr >= 0
    @constraint m minr .<= r
    @objective m Max minr
    status = solve(m)
    if status != :Optimal
        if status == :Infeasible
            error("An empty polyhedron has no H-Chebyshev center")
        elseif status == :Unbounded
            error("The polyhedron contains euclidean ball of arbitrary large radius")
        else
            error("Solver returned $status when computing the H-Chebyshev center")
        end
    end
    @constraint m minr == getvalue(minr)
    @variable m maxr >= 0
    @constraint m maxr .>= r
    @objective m Min maxr
    status = solve(m)
    @assert status == :Optimal
    (getvalue(c), getvalue(minr))
end

# TODO solver here should not be VRepSolver
"""
    vchebyshevcenter(p::VRep[, solver])

Return a tuple with the center and radius of the smallest euclidean ball containing the polyhedron `p`.
Throws an error if the polyhedron is empty or if the radius is infinite (i.e. `p` is not a polytope, it contains rays).
"""
function vchebyshevcenter(p::VRep, solver=Polyhedra.solver(p))
    error("TODO")
end

"""
    chebyshevcenter(p::Rep[, solver])

If `p` is a H-representation or is a polyhedron for which the H-representation has already been computed, calls `hchebyshevcenter`, otherwise, call `vchebyshevcenter`.
"""
function chebyshevcenter(p::Polyhedron, solver=Polyhedra.solver(p))
    if hrepiscomputed(p)
        hchebyshevcenter(p, solver)
    else
        vchebyshevcenter(p, solver)
    end
end
chebyshevcenter(p::HRepresentation, solver=Polyhedra.solver(p)) = hchebyshevcenter(p, solver)
chebyshevcenter(p::VRepresentation, solver=Polyhedra.solver(p)) = vchebyshevcenter(p, solver)

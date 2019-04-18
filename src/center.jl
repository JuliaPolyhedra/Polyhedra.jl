export hchebyshevcenter, vchebyshevcenter, chebyshevcenter
using JuMP

"""
    hchebyshevcenter(p::HRep[, solver])

Return a tuple with the center and radius of the largest euclidean ball contained in the polyhedron `p`.
Throws an error if the polyhedron is empty or if the radius is infinite.
"""
function hchebyshevcenter(p::HRep, solver=Polyhedra.default_solver(p; T=Float64)) # Need Float64 for `norm(a, 2)`
    model = JuMP.Model(solver)
    c = JuMP.@variable(model, [1:fulldim(p)])
    for hp in hyperplanes(p)
        a = convert(Vector{Float64}, hp.a)
        β = convert(Float64, hp.β)
        JuMP.@constraint(model, dot(a, c) == β)
    end
    JuMP.@variable(model, r[1:nhalfspaces(p)] >= 0)
    for (i, hs) in enumerate(halfspaces(p))
        a = convert(Vector{Float64}, hs.a)
        β = convert(Float64, hs.β)
        JuMP.@constraint(model, dot(a, c) + r[i] * norm(a, 2) <= β)
    end
    JuMP.@variable(model, minr >= 0)
    JuMP.@constraint(model, minr .<= r)
    JuMP.@objective(model, Max, minr)
    JuMP.optimize!(model)
    term = JuMP.termination_status(model)
    if term != MOI.OPTIMAL
        if term == MOI.INFEASIBLE
            error("An empty polyhedron has no H-Chebyshev center.")
        elseif term == MOI.DUAL_INFEASIBLE
            error("The polyhedron contains euclidean ball of arbitrary large radius.")
        else
            error("Solver returned $term when computing the H-Chebyshev center.")
        end
    end
    JuMP.@constraint(model, minr == JuMP.value(minr))
    JuMP.@variable(model, maxr >= 0)
    JuMP.@constraint(model, maxr .>= r)
    JuMP.@objective(model, Min, maxr)
    JuMP.optimize!(model)
    term = JuMP.termination_status(model)
    if term == MOI.OPTIMAL
        return (JuMP.value.(c), JuMP.value(minr))
    else
        error("Solver returned $term when computing the H-Chebyshev center.")
    end
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
function chebyshevcenter(p::Polyhedron, solver=Polyhedra.default_solver(p; T=Float64))
    if hrepiscomputed(p)
        hchebyshevcenter(p, solver)
    else
        vchebyshevcenter(p, solver)
    end
end
chebyshevcenter(p::HRepresentation, solver=Polyhedra.default_solver(p; T=Float64)) = hchebyshevcenter(p, solver)
chebyshevcenter(p::VRepresentation, solver=Polyhedra.default_solver(p; T=Float64)) = vchebyshevcenter(p, solver)

module PolyhedraJuMPExt

import JuMP
import Polyhedra: hrep, LPHRep, polyhedron, _optimize!
using Polyhedra: Rep, Projection, _moi_set, fulldim, dimension_names, PolyhedraToLPBridge, ProjectionBridge

"""
    hrep(model::JuMP.Model)

Builds an H-representation from the feasibility set of the JuMP model `model`.
Note that if non-linear constraint are present in the model, they are ignored.
"""
hrep(model::JuMP.Model) = LPHRep(model)
LPHRep(model::JuMP.Model) = LPHRep(JuMP.backend(model))
polyhedron(model::JuMP.Model, args...) = polyhedron(hrep(model), args...)
_optimize!(model::JuMP.Model) = JuMP.optimize!(model)

struct VariableInSet{V <: JuMP.ScalarVariable, S <: Union{Rep, Projection}} <: JuMP.AbstractVariable
    variables::Vector{V}
    set::S
end
function JuMP.build_variable(error_fun::Function, variables::Vector{<:JuMP.ScalarVariable}, set::Union{Rep, Projection})
    if length(variables) != fulldim(set)
        error("Number of variables ($(length(variables))) does not match the full dimension of the polyhedron ($(fulldim(set))).")
    end
    return VariableInSet(variables, set)
end
function JuMP.add_variable(model::JuMP.AbstractModel, v::VariableInSet, names)
    dim_names = dimension_names(v.set)
    if dim_names !== nothing
        names = copy(names)
        for i in eachindex(names)
            if isempty(names[i]) && !isempty(dim_names[i])
                names[i] = dim_names[i]
            end
        end
    end
    JuMP.add_bridge(model, PolyhedraToLPBridge)
    JuMP.add_bridge(model, ProjectionBridge)
    return JuMP.add_variable(model, JuMP.VariablesConstrainedOnCreation(v.variables, _moi_set(v.set)), names)
end
function JuMP.build_constraint(error_fun::Function, func, set::Rep)
    return JuMP.BridgeableConstraint(
        JuMP.build_constraint(error_fun, func, _moi_set(set)),
        PolyhedraToLPBridge)
end
function JuMP.build_constraint(error_fun::Function, func, set::Projection)
    return JuMP.BridgeableConstraint(JuMP.BridgeableConstraint(
        JuMP.build_constraint(error_fun, func, _moi_set(set)),
        ProjectionBridge), PolyhedraToLPBridge)
end

end

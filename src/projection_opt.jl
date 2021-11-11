struct ProjectionOptSet{T, RepT<:Rep{T}, I} <: MOI.AbstractVectorSet
    p::Projection{RepT, I}
end
MOI.dimension(set::ProjectionOptSet) = fulldim(set.p)
Base.copy(set::ProjectionOptSet) = set

struct ProjectionBridge{T, F, RepT, I} <: MOI.Bridges.Constraint.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    constraint::MOI.ConstraintIndex{F, PolyhedraOptSet{T, RepT}}
    dimensions::I
end
function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ProjectionBridge{T, F, RepT, I}}, model::MOI.ModelLike,
    f::MOI.AbstractVectorFunction, p::ProjectionOptSet) where {T, F, RepT, I}
    vf = MOIU.eachscalar(f)
    func = Vector{eltype(vf)}(undef, fulldim(p.p.set))
    for (i, j) in enumerate(p.p.dimensions)
        func[j] = vf[i]
    end
    N = fulldim(p.p.set)
    variables = MOI.add_variables(model, N - fulldim(p.p))
    for (i, j) in enumerate(setdiff(1:N, p.p.dimensions))
        func[j] = variables[i]
    end
    constraint = MOI.add_constraint(model, MOI.Utilities.vectorize(func), PolyhedraOptSet(p.p.set))
    return ProjectionBridge{T, F, RepT, I}(variables, constraint, p.p.dimensions)
end

MOI.supports_constraint(::Type{ProjectionBridge{T}},
                        ::Type{<:MOI.AbstractVectorFunction},
                        ::Type{<:ProjectionOptSet{T}}) where {T} = true
function MOI.Bridges.added_constrained_variable_types(::Type{<:ProjectionBridge})
    return Tuple{DataType}[]
end
function MOI.Bridges.added_constraint_types(::Type{ProjectionBridge{T, F, RepT, I}}) where {T, F, RepT, I}
    return [(F, PolyhedraOptSet{T, RepT})]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ProjectionBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{ProjectionOptSet{T, RepT, I}}) where {T, RepT, I}
    return ProjectionBridge{T, F, RepT, I}
end

# Attributes, Bridge acting as an model
function MOI.get(b::ProjectionBridge{T, F, RepT, I}, ::MOI.NumberOfConstraints{F, ProjectionOptSet{RepT, I}}) where {T, F, RepT, I}
    return 1
end
function MOI.get(b::ProjectionBridge{T, F, RepT, I}, ::MOI.ListOfConstraintIndices{F, ProjectionOptSet{RepT, I}}) where {T, F, RepT, I}
    return [b.constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, b::ProjectionBridge)
    MOI.delete(model, b.constraint)
    for vi in b.variables
        MOI.delete(model, vi)
    end
end
_moi_set(set::Projection) = ProjectionOptSet(set)
function JuMP.build_constraint(error_fun::Function, func, set::Projection)
    return JuMP.BridgeableConstraint(JuMP.BridgeableConstraint(
        JuMP.build_constraint(error_fun, func, _moi_set(set)),
        ProjectionBridge), PolyhedraToLPBridge)
end

using JuMP
struct DummyOptimizer <: MOI.AbstractOptimizer end
MOI.is_empty(::DummyOptimizer) = true
MOI.supports_constraint(::DummyOptimizer, ::Type{<:MOI.AbstractFunction}, ::Type{<:MOI.AbstractSet}) = true
function MOI.optimize!(::DummyOptimizer, model::MOI.ModelLike)
    return MOI.Utilities.identity_index_map(model), false
end
MOI.get(::DummyOptimizer, ::MOI.TerminationStatus) = MOI.OTHER_ERROR
MOI.get(::DummyOptimizer, ::MOI.RawStatusString) = "Dummy status"

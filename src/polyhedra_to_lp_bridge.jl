"""
    PolyhedraToLPBridge{T}

The `PolyhedraToLPBridge` converts a constraint `VF`-in-`PolyhedraOptSet`
into the constraints `F`-in-`EqualTo` for the hyperplanes and `F`-to-`LessThan`
for halfspaces.
"""
struct PolyhedraToLPBridge{T, F} <: MOI.Bridges.AbstractBridge
    hyperplanes::Vector{MOI.ConstraintIndex{F, MOI.EqualTo{T}}}
    halfspaces::Vector{MOI.ConstraintIndex{F, MOI.LessThan{T}}}
end
function PolyhedraToLPBridge{T, F}(model, f::MOI.AbstractVectorFunction, p::PolyhedraOptSet) where {T, F}
    vf = MOIU.eachscalar(f)
    hps = [
        MOIU.add_scalar_constraint(model, func(T, h.a, vf, F), MOI.EqualTo(convert(T, h.β)))
        for h in hyperplanes(p.rep)
    ]
    hss = [
        MOIU.add_scalar_constraint(model, func(T, h.a, vf, F), MOI.LessThan(convert(T, h.β)))
        for h in halfspaces(p.rep)
    ]
    return PolyhedraToLPBridge{T, F}(hps, hss)
end

function func(T::Type, a::AbstractVector, vf, F::Type)
    func = zero(F)
    for (α, f) in zip(a, vf)
        MOIU.operate!(+, T, func, MOIU.operate!(*, T, f, convert(T, α)))
    end
    return func
end

# start allowing everything (scalar)
MOI.supports_constraint(::Type{PolyhedraToLPBridge{T}},
                        ::Type{<:MOI.AbstractVectorFunction},
                        ::Type{<:PolyhedraOptSet}) where {T} = true
function MOI.Bridges.added_constraint_types(::Type{PolyhedraToLPBridge{T, F}}) where {T, F}
    return [(F, MOI.EqualTo{T}), (F, MOI.LessThan{T})]
end
function MOI.Bridges.concrete_bridge_type(::Type{<:PolyhedraToLPBridge{T}},
                              VF::Type{<:MOI.AbstractVectorFunction},
                              ::Type{<:PolyhedraOptSet}) where T
    SF = MOIU.scalar_type(VF)
    TermType = MOIU.promote_operation(*, T, T, SF)
    F = MOIU.promote_operation(+, T, TermType, TermType)
    return PolyhedraToLPBridge{T, F}
end

# Attributes, Bridge acting as an model
function MOI.get(b::PolyhedraToLPBridge{T, F}, ::MOI.NumberOfConstraints{F, MOI.EqualTo{T}}) where {T, F}
    return length(b.hyperplanes)
end
function MOI.get(b::PolyhedraToLPBridge{T, F}, ::MOI.ListOfConstraintIndices{F, MOI.EqualTo{T}}) where {T, F}
    return b.hyperplanes
end
function MOI.get(b::PolyhedraToLPBridge{T, F}, ::MOI.NumberOfConstraints{F, MOI.LessThan{T}}) where {T, F}
    return length(b.halfspaces)
end
function MOI.get(b::PolyhedraToLPBridge{T, F}, ::MOI.ListOfConstraintIndices{F, MOI.LessThan{T}}) where {T, F}
    return b.halfspaces
end

# Indices
function MOI.delete(model::MOI.ModelLike, b::PolyhedraToLPBridge)
    for h in b.hyperplanes
        MOI.delete(model, h)
    end
    for h in b.halfspaces
        MOI.delete(model, h)
    end
end

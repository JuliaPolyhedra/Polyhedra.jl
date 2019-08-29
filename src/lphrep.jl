export LPHRep

MOI.Utilities.@model(_MOIModel,
                     (), (MOI.EqualTo, MOI.LessThan,), (), (),
                     (), (MOI.ScalarAffineFunction,), (), ())
# We need the `SingleVariable` constraints to be bridged so we should say that
# they are not supported. We notably exclude `Integer` as we just ignore
# integrality constraints. Binary constraint should be bridged to integrality
# once https://github.com/JuliaOpt/MathOptInterface.jl/issues/704 is done.
function MOI.supports_constraint(
    ::_MOIModel{T}, ::Type{MOI.SingleVariable},
    ::Type{<:Union{MOI.EqualTo{T}, MOI.GreaterThan{T}, MOI.LessThan{T},
                   MOI.Interval{T}, MOI.ZeroOne}}) where T
    return false
end


mutable struct LPHRep{T} <: HRepresentation{T}
    model::_MOIModel{T}
    hyperplane_indices::Union{Nothing,
                              Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},
                                                         MOI.EqualTo{T}}}}
    halfspace_indices::Union{Nothing,
                             Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},
                                                        MOI.LessThan{T}}}}
end
function LPHRep(model::_MOIModel{T}) where T
    return LPHRep(model, nothing, nothing)
end
function LPHRep(model::MOI.ModelLike, T::Type = Float64)
    _model = _MOIModel{T}()
    bridged = MOI.Bridges.LazyBridgeOptimizer(_model)
    # Only enable constraint bridges that don't create variables and don't add
    # any variable bridge so that there is an identity mapping betwenen
    # variables of `model` and polyhedra dimensions.
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.GreaterToLessBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.LessToGreaterBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.NonnegToNonposBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.NonposToNonnegBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.ScalarizeBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.VectorizeBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.ScalarFunctionizeBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.VectorFunctionizeBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.SplitIntervalBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.NormInfinityBridge{T})
    # This one creates variables so the user need to consider the polyhedra as the
    # feasible set of the extended formulation.
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.NormOneBridge{T})
    MOI.copy_to(bridged, model)
    return LPHRep(_model)
end
# returns `Int64` so need to convert for 32-bit system
FullDim(rep::LPHRep) = convert(Int, MOI.get(rep.model, MOI.NumberOfVariables()))

hvectortype(::Type{LPHRep{T}}) where {T} = SparseVector{T, Int}
similar_type(::Type{LPHRep{S}}, ::FullDim, ::Type{T}) where {S, T} = LPHRep{T}
fulltype(::Type{LPHRep{T}}) where {T} = LPHRep{T}

LPHRep(h::HRep{T}) where {T} = LPHRep{T}(h)
LPHRep{T}(h::HRep) where {T} = convert(LPHRep{T}, h)

function _saf(a::AbstractVector{T}, vars) where T
    func = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{T}[], zero(T))
    for i in eachindex(a)
        coef = a[i]
        if !iszero(coef)
            push!(func.terms, MOI.ScalarAffineTerm(coef, vars[i]))
        end
    end
    return func
end

function LPHRep{T}(d::FullDim,
                   hyperplanes::ElemIt{<:HyperPlane{T}},
                   halfspaces::ElemIt{<:HalfSpace{T}}) where {T}
    model = _MOIModel{T}()
    vars = MOI.add_variables(model, fulldim(d))
    for hyperplane in hyperplanes
        func = _saf(hyperplane.a, vars)
        set = MOI.EqualTo(hyperplane.β)
        MOI.add_constraint(model, func, set)
    end
    for halfspace in halfspaces
        func = _saf(halfspace.a, vars)
        set = MOI.LessThan(halfspace.β)
        MOI.add_constraint(model, func, set)
    end
    return LPHRep(model)
end

function Base.copy(lp::LPHRep{T}) where {T}
    model = _MOIModel{T}()
    MOI.copy_to(model, lp.model)
    return LPHRep(model)
end

function constraint_indices(rep::LPHRep{T},
                            ::Union{HyperPlaneIndex, HyperPlaneIndices}) where T
    if rep.hyperplane_indices === nothing
        rep.hyperplane_indices = MOI.get(
            rep.model,
            MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{T},
                                        MOI.EqualTo{T}}())
    end
    return rep.hyperplane_indices
end
function MOI.add_constraint(rep::LPHRep{T},
                            func::MOI.ScalarAffineFunction{T},
                            set::MOI.EqualTo{T}) where T
    rep.hyperplane_indices = nothing
    return MOI.add_constraint(rep.model, func, set)
end

function constraint_indices(rep::LPHRep{T},
                            ::Union{HalfSpaceIndex, HalfSpaceIndices}) where T
    if rep.halfspace_indices === nothing
        rep.halfspace_indices = MOI.get(
            rep.model,
            MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{T},
                                        MOI.LessThan{T}}())
    end
    return rep.halfspace_indices
end
function MOI.add_constraint(rep::LPHRep{T},
                            func::MOI.ScalarAffineFunction{T},
                            set::MOI.LessThan{T}) where T
    rep.halfspace_indices = nothing
    return MOI.add_constraint(rep.model, func, set)
end

function Base.length(idxs::HIndices{T, LPHRep{T}}) where {T}
    return length(constraint_indices(idxs.rep, idxs))
end
function Base.isempty(idxs::HIndices{T, LPHRep{T}}) where {T}
    return isempty(constraint_indices(idxs.rep, idxs))
end
function startindex(idxs::HIndices{T, LPHRep{T}}) where {T}
    if isempty(idxs)
        return nothing
    else
        return eltype(idxs)(1)
    end
end
function Base.get(rep::LPHRep{T}, idx::HIndex{T}) where {T}
    ci = constraint_indices(rep, idx)[idx.value]
    func = MOI.get(rep.model, MOI.ConstraintFunction(), ci)::MOI.ScalarAffineFunction
    indices = [t.variable_index.value for t in func.terms]
    values = [t.coefficient for t in func.terms]
    a = sparsevec(indices, values, FullDim(rep))
    set = MOI.get(rep.model, MOI.ConstraintSet(), ci)
    β = MOI.constant(set) - func.constant
    if idx isa HyperPlaneIndex
        return HyperPlane(a, β)
    else
        return HalfSpace(a, β)
    end
end
function nextindex(rep::LPHRep{T}, idx::HIndex{T}) where {T}
    if idx.value >= length(constraint_indices(rep, idx))
        return nothing
    else
        return typeof(idx)(idx.value + 1)
    end
end

function Base.isvalid(lp::LPHRep{T}, idx::HIndex{T}) where {T}
    return 1 ≤ idx.value ≤ length(constraint_indices(rep, idx))
end

dualtype(::Type{LPHRep{T}}, ::Type{AT}) where {T, AT} = dualtype(Intersection{T, AT, Int}, AT)
dualfullspace(h::LPHRep, d::FullDim, ::Type{T}, ::Type{AT}) where {T, AT} = dualfullspace(Intersection{T, AT, Int}, d, T, AT)

export LPHRep

MOI.Utilities.@model(_MOIModel,
                     (), (MOI.EqualTo, MOI.LessThan,), (), (),
                     (), (MOI.ScalarAffineFunction,), (), ())
# We need the `SingleVariable` constraints to be bridged so we should say that
# they are not supported. We notably exclude `Integer` as we just ignore
# integrality constraints. Binary constraint should be bridged to integrality
# once https://github.com/jump-dev/MathOptInterface.jl/issues/704 is done.
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
    MOI.Bridges.add_bridge(bridged, PolyhedraToLPBridge{T})
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

supports_names(::Type{<:LPHRep}) = true
function dimension_names(h::LPHRep)
    if MOI.VariableName() in MOI.get(h.model, MOI.ListOfVariableAttributesSet())
        names = [
            MOI.get(h.model, MOI.VariableName(), MOI.VariableIndex(i))
            for i in 1:fulldim(h)
        ]
    else
        return nothing
    end
end

function _constrain_in_func_set(h::HyperPlane, vars, T)
    return _dot(h.a, vars, T), MOI.EqualTo{T}(h.β)
end
function _constrain_in_func_set(h::HalfSpace, vars, T)
    return _dot(h.a, vars, T), MOI.LessThan{T}(h.β)
end
function _constrain_in(model, hs::HIt, vars, T)
    for h in hs
        MOI.add_constraint(model, _constrain_in_func_set(h, vars, T)...)
    end
end

function LPHRep{T}(d::FullDim,
                   hyperplanes::ElemIt{<:HyperPlane{T}},
                   halfspaces::ElemIt{<:HalfSpace{T}};
                   dimension_names = nothing) where {T}
    model = _MOIModel{T}()
    vars = MOI.add_variables(model, fulldim(d))
    if dimension_names !== nothing
        if length(dimension_names) !== fulldim(d)
            throw(DimensionMismatch("Length of dimension_names ($(length(dimension_names))) does not match the full dimension of the polyhedron ($(fulldim(d)))."))
        end
        for (i, name) in enumerate(dimension_names)
            if !isempty(name)
                MOI.set(model, MOI.VariableName(), MOI.VariableIndex(i), name)
            end
        end
    end
    _constrain_in(model, hyperplanes, vars, T)
    _constrain_in(model, halfspaces, vars, T)
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
function MOI.delete(rep::LPHRep{T},
                    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}}) where T
    rep.hyperplane_indices = nothing
    return MOI.delete(rep.model, ci)
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
function MOI.delete(rep::LPHRep{T},
                    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.LessThan{T}}) where T
    rep.halfspace_indices = nothing
    return MOI.delete(rep.model, ci)
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
    func = MOI.get(rep.model, MOI.ConstraintFunction(), ci)::MOI.ScalarAffineFunction{T}
    # MOI uses `Int64` but `SparseArrays` uses `Int32` by default so `Int64` will create
    # issues with, e.g. preimages with `spzeros(d, n)`, etc...
    indices = Int[t.variable_index.value for t in func.terms]
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

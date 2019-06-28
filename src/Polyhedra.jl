module Polyhedra

using SparseArrays
using LinearAlgebra

export Polyhedron

abstract type Library end
abstract type Polyhedron{T} end

import JuMP, ParameterJuMP
const Solver = JuMP.OptimizerFactory
const SolverOrNot = Union{Nothing, Solver}

coefficient_type(::Union{AbstractVector{T}, Type{<:AbstractVector{T}}}) where T = T
similar_type(::Type{<:Vector}, ::Int, ::Type{T}) where T = Vector{T}
similar_type(::Type{SparseVector{S, IT}}, ::Int, ::Type{T}) where {S, IT, T} = SparseVector{T, IT}

emptymatrix(::Type{MT}, m, n) where {MT<:AbstractMatrix} = MT(undef, m, n)
emptymatrix(::Type{SparseMatrixCSC{T, Int}}, m, n) where T = spzeros(T, m, n)
similar_type(::Type{<:Matrix}, ::Type{T}) where T = Matrix{T}
similar_type(::Type{SparseMatrixCSC{S, I}}, ::Type{T}) where {S, I, T} = SparseMatrixCSC{T, I}

# Interface/Definitions
include("dimension.jl")
include("elements.jl")
include("comp.jl")
include("representation.jl")
include("indices.jl")
include("incidence.jl")
include("iterators.jl")
include("polyhedron.jl")

vvectortype(rep::Rep) = vvectortype(typeof(rep))
hvectortype(rep::Rep) = hvectortype(typeof(rep))
# TODO Only define vectortype for the type for Polyhedron v0.4
function vectortype(RepT::Type{<:Polyhedron})
    @assert hvectortype(RepT) == vvectortype(RepT)
    hvectortype(RepT)
end
vectortype(HRepT::Type{<:HRepresentation}) = hvectortype(HRepT)
vectortype(VRepT::Type{<:VRepresentation}) = vvectortype(VRepT)
vectortype(p::Rep) = vectortype(typeof(p))

vectortype(::Type{<:AbstractSparseArray{T}}) where T = SparseVector{T, Int}
vectortype(::Type{<:AbstractMatrix{T}}) where T = Vector{T}

hmatrixtype(RepT::Type{<:HRep}, T::Type) = matrixtype(similar_type(hvectortype(RepT), T))
vmatrixtype(RepT::Type{<:VRep}, T::Type) = matrixtype(similar_type(vvectortype(RepT), T))
matrixtype(::Type{<:AbstractVector{T}}) where T = Matrix{T}
matrixtype(::Type{<:AbstractSparseVector{T}}) where T = SparseMatrixCSC{T, Int}

function similar_type(::Type{ET}, ::Type{Tout}) where {Tout, ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, FullDim(ET), Tout)
end
function similar_type(::Type{ET}, d::FullDim) where {ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, d, coefficient_type(ET))
end

# Operations
include("repop.jl")
include("center.jl")
include("repelemop.jl")
include("aff.jl")
include("redundancy.jl")
include("projection.jl")

# Implementations of representations
include("vecrep.jl")
include("mixedrep.jl")
include("lphrep.jl")
include("jump.jl")
include("matrep.jl")
include("liftedrep.jl")
include("doubledescription.jl") # FIXME move it after projection.jl once it stops depending on LiftedRep
include("interval.jl") # 1D polyhedron
include("defaultlibrary.jl")

include("fulldim.jl")

# Optimization
include("opt.jl")
include("polyhedra_to_lp_bridge.jl")
include("vrep_optimizer.jl")
include("default.jl")

# Visualization
include("show.jl")
include("recipe.jl")
include("decompose.jl")

end # module

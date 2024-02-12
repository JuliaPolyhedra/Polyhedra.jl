import StaticArrays
function similar_type(SAT::Type{<:StaticArrays.SVector},
                      size::StaticArrays.Size, T::Type)
    StaticArrays.similar_type(SAT, T, size)
end
# Called when mixing polyhedra with one using StaticArrays and the other not.
similar_type(::Type{<:StaticArrays.SVector}, ::Int, T::Type) = Vector{T}
similar_type(::Type{<:Vector}, ::StaticArrays.Size, T::Type) = Vector{T}

# The dimension is Int for Vector, SparseVector and
# StaticArrays.Size for StaticArrays.SVector
# This allows similar_type to be compiled as only one method for Vector
# and SparseVector and be type stable for StaticArrays.SVector

const FullDim = Union{Int, StaticArrays.Size}

"""
    typed_fulldim(p)::FullDim

Similar to [`fulldim`](@ref) but used for type stability with the vector type.
If the vector type is `StaticArrays.SVector` then it returns a
`StaticArrays.Size`.
"""
typed_fulldim(::Type{<:AbstractVector}) = -1 # Shouldn't hurt as it will not be used
typed_fulldim(v::AbstractVector) = length(v)
typed_fulldim(v::Union{StaticArrays.SVector{N}, Type{<:StaticArrays.SVector{N}}}) where N = StaticArrays.Size(v)

function FullDim_convert(::Type{StaticArrays.Size{N}}, d::Int) where N
    @assert N[1] == d
    return StaticArrays.Size{N}()
end
function FullDim_convert(::Type{StaticArrays.Size{N}},
                         d::StaticArrays.Size{N}) where N
    return d
end
FullDim_convert(::Type{Int}, d::Int) = d
function FullDim_convert(::Type{Int}, d::StaticArrays.Size{N}) where N
    return N[1]
end

"""
    fulldim(rep::Rep)::Int

Returns the dimension of the space in which polyhedron, representation,
element or vector is defined. That is, a straight line in a 3D space has
`fulldim` 3 even if its dimension is 1.
"""
fulldim(p) = fulldim(typed_fulldim(p))

function fulldim end

fulldim(N::Integer) = convert(Int, N)
fulldim(::StaticArrays.Size{N}) where N = N[1]

neg_fulldim(N::Int) = -N
# FIXME May need generated function to make it type stable
neg_fulldim(::StaticArrays.Size{N}) where N = StaticArrays.Size{(-N[1],)}()

map_fulldim(f::Function, N::Integer) = f(N)
map_fulldim(f::Function, ::StaticArrays.Size{N}) where N = StaticArrays.Size{(f(N[1]),)}()

sum_fulldim(N1::Int, N2::Int) = N1 + N2
# FIXME May need generated function to make it type stable
sum_fulldim(::StaticArrays.Size{N1}, ::StaticArrays.Size{N2}) where {N1, N2} = StaticArrays.Size{(N1[1] + N2[1],)}()
# Called in cartesian product of two polyhedra, only the second using static arrays.
sum_fulldim(N1::Int, ::StaticArrays.Size{N2}) where N2 = N1 + N2[1]
sum_fulldim(::StaticArrays.Size{N1}, N2::Int) where N1 = N1[1] + N2

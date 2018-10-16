FullDim_rec() = error("Cannot infer dimension of polyhedron constructed from no element.")
function FullDim_rec(it::It, its::Union{Rep, It}...)
    if vectortype(eltype(it)) <: StaticArrays.SVector
        return FullDim(eltype(it))
    elseif isempty(it)
        return FullDim_rec(its...)
    else
        return FullDim(first(it))
    end
end

FullDim_rep(rep::Rep, other_reps::Union{Nothing, Rep}...) = FullDim(rep)
FullDim_rep(rep::Nothing, other_reps::Union{Nothing, Rep}...) = FullDim_rep(other_reps...)

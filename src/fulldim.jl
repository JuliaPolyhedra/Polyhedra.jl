function FullDim_rec()
    return error(
        "Cannot infer dimension of polyhedron constructed from no element.",
    )
end
function FullDim_rec(it::It, its::Union{Rep,It}...)
    if vectortype(eltype(it)) <: StaticArrays.SVector
        return typed_fulldim(eltype(it))
    elseif isempty(it)
        return FullDim_rec(its...)
    else
        return typed_fulldim(first(it))
    end
end

FullDim_rep(rep::Rep, other_reps::Union{Nothing,Rep}...) = typed_fulldim(rep)
function FullDim_rep(rep::Nothing, other_reps::Union{Nothing,Rep}...)
    return FullDim_rep(other_reps...)
end

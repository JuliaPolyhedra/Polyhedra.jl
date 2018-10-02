# -1 is the dimension of an empty polyhedron, here it is used as the
# *full* dimension of a polyhedron with no element
FullDim_rec() = -1
function FullDim_rec(it::ElemIt{<:Union{AT, VStruct{T, AT},
                                        HRepElement{T, AT}}},
                     its::Union{Rep, It}...) where {T,
                                                    AT <: StaticArrays.SVector}
    return FullDim(AT)
end
function FullDim_rec(it::It, its::Union{Rep, It}...)
    if isempty(it)
        return FullDim_rec(its...)
    else
        return FullDim(first(it))
    end
end

FullDim_rep(rep::Rep, other_reps::Union{Nothing, Rep}...) = FullDim(rep)
FullDim_rep(rep::Nothing, other_reps::Union{Nothing, Rep}...) = FullDim(other_reps...)

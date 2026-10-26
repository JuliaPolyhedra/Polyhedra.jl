struct InconsistentVRep{T,AT,D<:Polyhedra.FullDim} <: VRepresentation{T}
    points::Polyhedra.PointsHull{T,AT,D}
    rays::Polyhedra.RaysHull{T,AT,D}
    function InconsistentVRep{T,AT,D}(
        d::Polyhedra.FullDim,
        points,
        lines,
        rays,
    ) where {T,AT,D}
        return new{T,AT,D}(
            Polyhedra.PointsHull(d, points),
            Polyhedra.RaysHull(d, lines, rays),
        )
    end
end
Polyhedra.FullDim(rep::InconsistentVRep) = Polyhedra.FullDim(rep.points)
function Polyhedra.dualtype(
    ::Type{InconsistentVRep{T,AT,D}},
    ::Type{AT},
) where {T,AT,D}
    return Polyhedra.Intersection{T,AT,D}
end
Polyhedra.hvectortype(::Type{<:InconsistentVRep{T,AT}}) where {T,AT} = AT
Polyhedra.vvectortype(::Type{<:InconsistentVRep{T,AT}}) where {T,AT} = AT
function Polyhedra.similar_type(
    PT::Type{<:InconsistentVRep},
    d::Polyhedra.FullDim,
    ::Type{T},
) where {T}
    return InconsistentVRep{
        T,
        Polyhedra.similar_type(Polyhedra.hvectortype(PT), d, T),
        typeof(d),
    }
end
function Polyhedra.fulltype(::Type{InconsistentVRep{T,AT,D}}) where {T,AT,D}
    return InconsistentVRep{T,AT,D}
end
#Polyhedra.@subrepelem InconsistentVRep SymPoint points
Polyhedra.@subrepelem InconsistentVRep Point points
Polyhedra.@subrepelem InconsistentVRep Line rays
Polyhedra.@subrepelem InconsistentVRep Ray rays

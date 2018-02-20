export hrep, vrep

"""
    hrep(p::Polyhedron)

Returns an H-representation for the polyhedron `p`.
"""
hrep(p::Polyhedron) = error("`hrep` not implemented for `$(eltype(p))`")

"""
    vrep(p::Polyhedron)

Returns a V-representation for the polyhedron `p`.
"""
vrep(p::Polyhedron) = error("`vrep` not implemented for `$(eltype(p))`")

# Mandatory
export polyhedron, loadpolyhedron!
export hrepiscomputed, hrep, vrepiscomputed, vrep
export volume, surface

"""
    polyhedron(rep::Representation{T})

Creates a polyhedron from the representation `rep` using the default library included in the Polyhedra package.
"""
polyhedron(rep::Representation{T}) where {T} = polyhedron(rep, default_library(FullDim(rep), T))

"""
    hrepiscomputed(p::Polyhedron)

Returns whether the H-representation of this polyhedron has been computed.
"""
function hrepiscomputed end

"""
    hrep(p::Polyhedron)

Returns an H-representation for the polyhedron `p`.
"""
hrep(p::Polyhedron) = error("`hrep` not implemented for `$(eltype(p))`")

"""
    vrepiscomputed(p::Polyhedron)

Returns whether the V-representation of this polyhedron has been computed.
"""
function vrepiscomputed end

"""
    vrep(p::Polyhedron)

Returns a V-representation for the polyhedron `p`.
"""
vrep(p::Polyhedron) = error("`vrep` not implemented for `$(eltype(p))`")

"""
    volume(p::Polyhedron{T}) where {T}

Returns the `fulldim(p)`-dimensional hyper-volume of the polyhedron `p`.
Returns `Inf` or `-one(T)` if it is infinite depending on whether the type `T` has an infinite value.
"""
function volume end

"""
    surface(p::Polyhedron{T}) where {T}

Returns the `fulldim(p)-1`-dimensional hyper-volume of the surface of the polyhedron `p`.
Returns `Inf` or `-one(T)` if it is infinite depending on whether the type `T` has an infinite value.
"""
function surface end

loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Type{Val{:ine}}) = error("loadpolyhedron! not implemented for .ine")
loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Type{Val{:ext}}) = error("loadpolyhedron! not implemented for .ext") # FIXME ExtFileVRepresentation or just ExtFile

loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Symbol) = loadpolyhedron!(p, filename, Val{extension})

function loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::AbstractString)
    s = something(findfirst(isequal(extension), ["ext", "ine"]), 0)
    if s == 0
        error("Invalid extension $extension, please give 'ext' for V-representation or 'ine' for H-representation")
    end
    loadpolyhedron!(p, filename, [:ext, :ine][s])
end

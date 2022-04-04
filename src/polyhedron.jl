# Mandatory
export polyhedron, loadpolyhedron!
export hrepiscomputed, hrep, vrepiscomputed, vrep
export volume, volume_simplex, unscaled_volume_simplex, surface, center_of_mass

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
    Polyhedra.resethrep!(p::Polyhedron, h::HRepresentation, redundancy::Redundancy = UNKNOWN_REDUNDANCY)

Reset the H-representation of `p` to `h`.
The redundancy of `p` is assumed to be `redundancy`; see [`Polyhedra.Redundancy`](@ref).

!!! info
    The representation is not assumed to be a valid representation for `p`
    so it invalidates the V-representation of `p`.
    Use [`Polyhedra.sethrep!`](@ref) if `h` is known to be a valid representation for `p`.
"""
function resethrep! end

"""
    Polyhedra.sethrep!(p::Polyhedron, h::HRepresentation, redundancy::Redundancy = UNKNOWN_REDUNDANCY)

Reset the H-representation of `p` to `h`.
The redundancy of `p` is assumed to be `redundancy`; see [`Polyhedra.Redundancy`](@ref).

!!! warning
    The representation is assumed to be a valid representation for `p`
    so it does not invalidate the V-representation of `p` if it was already computed previously.
    Use [`Polyhedra.resethrep!`](@ref) if `h` is not known to be a valid representation for `p`.
"""
function sethrep! end

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
    Polyhedra.resetvrep!(p::Polyhedron, v::VRepresentation, redundancy::Redundancy = UNKNOWN_REDUNDANCY)

Reset the V-representation of `p` to `v`.
The redundancy of `p` is assumed to be `redundancy`; see [`Polyhedra.Redundancy`](@ref).

!!! info
    The representation is not assumed to be a valid representation for `p`
    so it invalidates the H-representation of `p`.
    Use [`Polyhedra.setvrep!`](@ref) if `v` is known to be a valid representation for `p`.
"""
function resetvrep! end

"""
    Polyhedra.setvrep!(p::Polyhedron, v::HRepresentation, redundancy::Redundancy = UNKNOWN_REDUNDANCY)

Reset the H-representation of `p` to `v`.
The redundancy of `p` is assumed to be `redundancy`; see [`Polyhedra.Redundancy`](@ref).

!!! warning
    The representation is assumed to be a valid representation for `p`
    so it does not invalidate the H-representation of `p` if it was already computed previously.
    Use [`Polyhedra.resetvrep!`](@ref) if `v` is not known to be a valid representation for `p`.
"""
function setvrep! end

"""
    volume(p::Polyhedron{T}) where {T}

Returns the `fulldim(p)`-dimensional hyper-volume of the polyhedron `p`.
Returns `Inf` or `-one(T)` if it is infinite depending on whether the type `T` has an infinite value.
"""
function volume(p::Polyhedron{T}) where T
    return reduce(+, unscaled_volume_simplex(Δ) for Δ in triangulation(p);
                  init=zero(T)) / factorial(fulldim(p))
end

function unscaled_volume_simplex(s)
    A = Matrix{coefficient_type(s)}(undef, fulldim(s), npoints(s) - 1)
    v = first(points(s))
    for (i, p) in enumerate(Iterators.drop(points(s), 1))
        A[:, i] = p - v
    end
    return abs(det(A))
end
function volume_simplex(s)
    return unscaled_volume_simplex(s) / factorial(fulldim(s))
end

"""
    center_of_mass(p::Polyhedron{T}) where {T}

Returns the center of mass of `p`, represented as a `Vector{T}` of length `fulldim(p)`.
Throws an error if `p` is degenerate.
"""
function center_of_mass(p::Polyhedron{T}) where T
    # Implementation strategy: For a simplex, the center of mass coincides with
    # the centroid which is easy to compute.  So we triangulate `p` into
    # simplices and compute a volume-weighted average of these simplices'
    # centers of mass.
    simplices = triangulation(p)
    isempty(simplices) && error("Tried to compute center of mass of a degenerate polyhedron")
    unscaled_center = zeros(T, fulldim(p))
    unscaled_vol = zero(T)
    for Δ in simplices
        v = unscaled_volume_simplex(Δ)
        unscaled_center += v * sum(points(Δ))
        unscaled_vol += v
    end
    return unscaled_center / (unscaled_vol * (fulldim(p) + 1))
end

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

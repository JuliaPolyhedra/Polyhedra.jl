# Mandatory
export polyhedron, hrep, vrep, hrepiscomputed, vrepiscomputed, loadpolyhedron!
export volume, surface

polyhedron(rep::Representation{N, T}) where {N, T} = polyhedron(rep, default_library(FullDim{N}(), T))

"""
    volume(p::Polyhedron{N, T}) where {N, T}

Returns the `N`-dimensional hyper-volume of the polyhedron `p`.
Returns `Inf` or `-one(T)` if it is infinite depending on whether the type `T` has an infinite value.
"""
function volume end

"""
    surface(p::Polyhedron{N, T}) where {N, T}

Returns the `N-1`-dimensional hyper-volume of the surface of the polyhedron `p`.
Returns `Inf` or `-one(T)` if it is infinite depending on whether the type `T` has an infinite value.
"""
function surface end

"""
    hrepiscomputed(p::Polyhedron)

Returns whether the H-representation of this polyhedron has been computed.
"""
function hrepiscomputed end

function hrep end

"""
    vrepiscomputed(p::Polyhedron)

Returns whether the V-representation of this polyhedron has been computed.
"""
function vrepiscomputed(p::Polyhedron) end

function vrep(p::Polyhedron) end

loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Type{Val{:ine}}) = error("loadpolyhedron! not implemented for .ine")
loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Type{Val{:ext}}) = error("loadpolyhedron! not implemented for .ext") # FIXME ExtFileVRepresentation or just ExtFile

# These can optionally be reimplemented for speed by a library
export dim, radialprojectoncut

loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::Symbol) = loadpolyhedron!(p, filename, Val{extension})

function loadpolyhedron!(p::Polyhedron, filename::AbstractString, extension::AbstractString)
    s = findfirst(["ext", "ine"], filename)
    if s == 0
        error("Invalid extension $extension, please give 'ext' for V-representation or 'ine' for H-representation")
    end
    loadpolyhedron!(p, filename, [:ext, :ine][s])
end

# function transformgenerators{N}(p::Polyhedron{N}, P::AbstractMatrix)
#   # Each generator x is transformed to P * x
#   # If P is orthogonal, the new axis are the rows of P.
#   if size(P, 2) != N
#     error("The number of columns of P must match the dimension of the polyhedron")
#   end
#   ext = P * getgenerators(p)
#   polyhedron(ext, getlibraryfor(p, eltype(ext)))
# end
#
# function transforminequalities(p::Polyhedron, P::AbstractMatrix)
#   # The new axis are the column of P.
#   # Let y be the coordinates of a point x in these new axis.
#   # We have x = P * y so y = P \ x.
#   # We have
#   # b = Ax = A * P * (P \ x) = (A * P) * y
#   ine = getinequalities(p) * P
#   polyhedron(ine, getlibraryfor(p, eltype(ine)))
# end

# function (*){N,S}(A::AbstractMatrix{S}, p::Polyhedron{N})
#   if size(A, 2) != N
#     error("Incompatible dimension")
#   end
#   if generatorsarecomputed(p)
#     transformgenerators(p, A)
#   else # FIXME not wokring
#     ine = SimpleHRepresentation(getinequalities(p))
#     nnew = size(A, 1)
#     if false
#       # CDD works with delset not at the end ?
#       newA = [ine.A spzeros(S, size(ine.A, 1), nnew);
#                   A  -speye(S, nnew, nnew)]
#       delset = IntSet(nnew+(1:N))
#     else
#       newA = [spzeros(S, size(ine.A, 1), nnew) ine.A;
#                -speye(S, nnew, nnew) A]
#       delset = IntSet(1:N)
#     end
#     newb = [ine.b; spzeros(S, nnew)]
#     newlinset = ine.linset ∪ IntSet(N+(1:nnew))
#     newine = SimpleHRepresentation(newA, newb, newlinset)
#     newpoly = polyhedron(newine, getlibraryfor(p, eltype(newine)))
#     eliminate(newpoly, IntSet(nnew+(1:N)))
#   end
# end

#function fulldim{N,T}(p::Polyhedron{N,T})
#  N
#end

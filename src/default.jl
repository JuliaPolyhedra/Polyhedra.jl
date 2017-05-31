export defaultLPsolverfor

getlibraryfor{T}(n::Int, ::Type{T}) = getlibraryfor(Val{n}, T)
getlibraryfor(::Type{Val{N}}, ::Type{T}) where {N, T} = SimplePolyhedraLibrary()

# Libraries can implement this function to be the default for a certain type T
# getlibraryfor{T<:Real}(p::LibraryPolyhedronType, n, ::Type{T})
# This function should be implemented by libraries
# or it would be like not implementing `similar` for AbstractArray
getlibraryfor{T}(p::Polyhedron, n::Int, ::Type{T}) = getlibraryfor(n, T)
getlibraryfor{N, T}(p::Polyhedron{N}, ::Type{T}) = getlibraryfor(p, N, T)
getlibraryfor{N, T}(p::Polyhedron{N, T}, n::Int) = getlibraryfor(p, n, T)
getlibrary{N, T}(p::Polyhedron{N, T}) = getlibraryfor(p, N, T)

function defaultLPsolverfor{N,T}(p::Rep{N,T}, solver=JuMP.UnsetSolver())
    if vrepiscomputed(p)
        SimpleVRepSolver()
    else
        solver
    end
end

defaultLPsolverfor(::VRepresentation, solver=JuMP.UnsetSolver()) = SimpleVRepSolver()
defaultLPsolverfor(::HRepresentation, solver=JuMP.UnsetSolver()) = solver

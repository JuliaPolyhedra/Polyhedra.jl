export defaultLPsolverfor

getlibraryfor(n::Int, ::Type{T}) where {T} = getlibraryfor(Val{n}, T)
getlibraryfor(::Type{Val{N}}, ::Type{T}) where {N, T} = SimplePolyhedraLibrary{T}()
getlibraryfor(::Type{Val{1}}, ::Type{T}) where T = IntervalLibrary{T}()

# Libraries can implement this function to be the default for a certain type T
# getlibraryfor{T<:Real}(p::LibraryPolyhedronType, n, ::Type{T})
# This function should be implemented by libraries
# or it would be like not implementing `similar` for AbstractArray
getlibraryfor(p::Polyhedron, n::Int, ::Type{T}) where {T} = getlibraryfor(n, T)
getlibraryfor(p::Polyhedron{N}, ::Type{T}) where {N, T} = getlibraryfor(p, N, T)
getlibraryfor(p::Polyhedron{N, T}, n::Int) where {N, T} = getlibraryfor(p, n, T)
getlibrary(p::Polyhedron{N, T}) where {N, T} = getlibraryfor(p, N, T)

function defaultLPsolverfor(p::Rep{N,T}, solver=JuMP.UnsetSolver()) where {N,T}
    if vrepiscomputed(p)
        SimpleVRepSolver()
    else
        solver
    end
end

defaultLPsolverfor(::VRepresentation, solver=JuMP.UnsetSolver()) = SimpleVRepSolver()
defaultLPsolverfor(::HRepresentation, solver=JuMP.UnsetSolver()) = solver

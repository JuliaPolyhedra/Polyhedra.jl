export defaultLPsolverfor
# these are the default library packages in order of decreasing precedence
# for each supported polyhedra type

const ExactLibraries = [(:CDDLib, :CDDLibrary, [:exact]),
                        (:LRSLib, :LRSLibrary, Symbol[]),
                        (:ConvexHull, :ConvexHullLib, [:exact])]
const FloatLibraries = [(:CDDLib, :CDDLibrary, [:float]),
                        (:QHull, :QHullLib, Symbol[]),
                        (:ConvexHull, :ConvexHullLib, [:float])]

# This function is inspired from MathProgBase/src/defaultsolvers.jl
for (librarytype, AbstractType) in [("Exact", Real), ("Float", AbstractFloat)]
    libraries = Symbol(librarytype * "Libraries")
    @eval function getlibraryfor{T<:$AbstractType}(::Type, ::Type{T})
        for (pkgname, libraryname, args) in $libraries
            try
                eval(Expr(:import, pkgname))
            catch
                if isdir(Pkg.dir((string(pkgname))))
                    warn("Package ",string(pkgname),
                         " is installed but couldn't be loaded. ",
                         "You may need to run `Pkg.build(\"$pkgname\")`")
                end
                continue
            end
            #ex = Expr(:(=), $(quot(defaultname)),
            #          Expr(:call, Expr(:., pkgname, quot(solvername))))
            #return eval(Expr(:call, Expr(:., pkgname, quot(libraryname))))
            return eval(Expr(:call, Expr(:., pkgname, Base.Meta.quot(libraryname)), Base.Meta.quot.(args)...))
        end
        return SimplePolyhedraLibrary()
    end
end

getlibraryfor{T}(n::Int, ::Type{T}) = getlibraryfor(Val{n}, T)

# Libraries can implement this function to be the default for a certain type T
# getlibraryfor{T<:Real}(p::LibraryPolyhedronType, n, ::Type{T})
# This function should be implemented by libraries
# or it would be like not implementing `similar` for AbstractArray
getlibraryfor{T}(p::Polyhedron, n::Int, ::Type{T}) = getlibraryfor(n, T)
getlibraryfor{N, T}(p::Polyhedron{N}, ::Type{T}) = getlibraryfor(p, N, T)
getlibraryfor{N, T}(p::Polyhedron{N, T}, n::Int) = getlibraryfor(p, n, T)
getlibrary{N, T}(p::Polyhedron{N, T}) = getlibraryfor(p, N, T)

function defaultLPsolverfor{N,T}(p::Rep{N,T}, solver=nothing)
    if vrepiscomputed(p)
        SimpleVRepSolver()
    else
        if solver === nothing
            JuMP.UnsetSolver
        else
            solver
        end
    end
end

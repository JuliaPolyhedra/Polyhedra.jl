# these are the default library packages in order of decreasing precedence
# for each supported polyhedra type

const ExactLibraries = [(:CDDLib, :CDDLib),
                        (:LRSLib, :LRSLib),
                        (:ConvexHull, :ConvexHullLib)]
const FloatLibraries = [(:CDDLib, :CDDLib),
                        (:QHull, :QHullLib),
                        (:ConvexHull, :ConvexHullLib)]

# This function is inspired from MathProgBase/src/defaultsolvers.jl
function getlibraryfor{T<:AbstractFloat}(n, ::Type{T})
    for (pkgname, libraryname) in $libraries
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
        ex = Expr(:(=), $(quot(defaultname)),
                  Expr(:call, Expr(:., pkgname, quot(solvername))))
        eval(ex)
        ex = Expr(:call, $(quot(t)), $(quot(defaultname)))
        return eval(ex)
    end
end

# Libraries can implement this function to be the default for a certain type T
# getlibraryfor{T<:Real}(p::LibraryPolyhedronType, n, ::Type{T})
# This function should be implemented by libraries
# or it would be like not implementing `similar` for AbstractArray
getlibraryfor{T}(p::Polyhedron, ::Type{T}) = getlibraryfor(T)
getlibrary(p::Polyhedron) = getlibraryfor(p, eltype(p))

function defaultLPsolverfor{N,T}(p::Rep{N,T})
    if vrepiscomputed(p)
        SimpleVRepSolver()
    else
        MathProgBase.defaultLPsolver
    end
end

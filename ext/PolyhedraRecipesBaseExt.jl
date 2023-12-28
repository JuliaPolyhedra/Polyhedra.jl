module PolyhedraRecipesBaseExt

using LinearAlgebra
using RecipesBase
using Polyhedra
using Polyhedra: basis, _semi_hull

function planar_contour(p::Polyhedron)
    if fulldim(p) != 2
        if fulldim(p) == 3
            error("Plotting 3-dimensional polyhedron with Plots is not supported, use Makie or MeshCat.")
        else
            error("Plotting $(fulldim(p))-dimensional polyhedron with Plots is not supported.")
        end
    end
    removevredundancy!(p)
    if hasallrays(p)
        error("Rays not supported yet in the 2D plotting recipe.")
    end
    ps = collect(points(p))
    if isempty(ps)
        error("Plotting empty polyhedron is not supported.")
    end
    sort!(ps, by = x -> x[1])
    counterclockwise = (p1, p2) -> dot(cross([p1; 0], [p2; 0]), [0, 0, 1])
    sweep_norm = basis(eltype(ps), fulldim(p), 1)
    top = _semi_hull(ps,  1, counterclockwise, sweep_norm)
    bot = _semi_hull(ps, -1, counterclockwise, sweep_norm)
    if !isempty(top) && !isempty(bot)
        @assert top[end] == bot[1]
        pop!(top)
    end
    if !isempty(bot) && !isempty(top)
        @assert bot[end] == top[1]
        pop!(bot)
    end
    hull = [top; bot]
    push!(hull, hull[1])  # ensure shape is closed
    return [p[1] for p in hull], [p[2] for p in hull]
end

@recipe function f(p::Polyhedron)
    seriestype --> :shape
    legend --> false
    planar_contour(p)
end

end

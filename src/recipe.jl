using RecipesBase

function getsemihull(ps::Vector{PT}, sign_sense, counterclockwise, yray = nothing) where PT
    hull = PT[]
    if length(ps) == 0
        return hull
    end
    prev = sign_sense == 1 ? first(ps) : last(ps)
    cur = prev
    for j in (sign_sense == 1 ? (2:length(ps)) : ((length(ps)-1):-1:1))
        while prev != cur && counterclockwise(cur - prev, ps[j] - prev) >= 0
            cur = prev
            pop!(hull)
            if !isempty(hull)
                prev = last(hull)
            end
        end
        if yray !== nothing && counterclockwise(ps[j] - cur, yray) >= 0
            break
        else
            push!(hull, cur)
            prev = cur
            cur = ps[j]
        end
    end
    push!(hull, cur)
    hull
end


@recipe function f(p::Polyhedron)
    if fulldim(p) != 2
        error("Plotting 3-dimensional polyhedron with Plots is not supported, use Makie, MeshCat or DrakeVisualizer")
    end
    removevredundancy!(p)
    if hasrays(p)
        warn("Rays not supported yet in the 2D plotting recipe")
    end
    ps = collect(points(p))
    sort!(ps, by = x -> x[1])
    counterclockwise(p1, p2) = dot(cross([p1; 0], [p2; 0]), [0, 0, 1])
    hull = [getsemihull(ps, 1, counterclockwise); getsemihull(ps, -1, counterclockwise)]
    seriestype --> :shape
    legend --> false
    [p[1] for p in hull], [p[2] for p in hull]
end

function _semi_hull(ps::Vector{PT}, sign_sense, counterclockwise, sweep_norm, yray::Nothing = nothing) where PT
    hull = PT[]
    if isempty(ps)
        return hull
    end
    prev = sign_sense == 1 ? first(ps) : last(ps)
    cur = prev
    # Invariant:
    # We either have:
    # * `hull` is empty and `cur == prev` or
    # * `hull` is nonempty and `prev = last(hull)`.
    # In any case, the semihull of `ps[1:(j-1)]` is given by `[hull; cur]`.
    for j in (sign_sense == 1 ? (2:length(ps)) : ((length(ps)-1):-1:1))
        skip = false
        while prev != cur
            cur_vec = cur - prev
            psj_vec = ps[j] - prev
            cc = counterclockwise(cur_vec, psj_vec)
            if isapproxzero(cc)
                # `prev`, `cur` and `ps[j]` are on the same line
                # Two cases here:
                # 1) `prev`, `cur`, `ps[j]` are on a line perpendicular to `sweep_norm`
                #    The one that is not clockwise is redundant.
                # 2) `cur`, `ps[j]` and on the same ray starting from `prev` of direction
                #    `sweep_norm * sign_sense`
                #    The one that is closer to `prev` is redundant.
                dcur = dot(cur_vec, sweep_norm)
                dpsj = dot(psj_vec, sweep_norm)
                if isapproxzero(dcur) && isapproxzero(dpsj)
                    # Case 1
                    if sign_sense * counterclockwise(cur_vec, sweep_norm) < sign_sense * counterclockwise(psj_vec, sweep_norm)
                        skip = true
                        break
                    end
                elseif sign_sense * dcur > sign_sense * dpsj
                    skip = true
                    break
                end
            elseif cc < 0
                break
            end
            # `cur` is redundant so `cur` is dropped and `prev` becomes `cur`
            cur = prev
            pop!(hull)
            if !isempty(hull)
                prev = last(hull)
            end
        end
        if !skip
            push!(hull, cur)
            prev = cur
            cur = ps[j]
        end
    end
    push!(hull, cur)
    return hull
end

function _semi_hull(ps::Vector{PT}, sign_sense, counterclockwise, sweep_norm, yray) where PT
    hull = _semi_hull(ps, sign_sense, counterclockwise, sweep_norm)
    while length(hull) >= 2 && counterclockwise(hull[end] - hull[end - 1], yray) >= 0
        pop!(hull)
    end
    return hull
end

_colinear(counterclockwise, a, b, c) = isapproxzero(counterclockwise(b - a, c - a))

function _planar_hull(d::FullDim, points, lines, rays, counterclockwise, rotate)
    line = nothing
    lineleft = false
    lineright = false
    function checkleftright(r::Union{Ray, Line})
        cc = counterclockwise(r, line)
        if !isapproxzero(cc)
            if cc < 0 || islin(r)
                lineleft = true
            end
            if cc > 0 || islin(r)
                lineright = true
            end
        end
    end
    for l in lines
        if !isapproxzero(l)
            if line === nothing
                line = l
            else
                checkleftright(l)
            end
        end
    end
    xray = yray = nothing
    for r in rays
        isapproxzero(r) && continue
        if line === nothing
            if xray === nothing
                xray = yray = coord(r)
            else
                if Line(xray) ≈ linearize(r) && !(Ray(xray) ≈ r)
                    line = Line(xray)
                    checkleftright(Ray(yray))
                elseif Line(yray) ≈ linearize(r) && !(Ray(yray) ≈ r)
                    line = Line(yray)
                    checkleftright(Ray(xray))
                else
                    x_right = counterclockwise(r, xray) > 0
                    y_left = counterclockwise(r, yray) < 0
                    if x_right
                        if y_left
                            line = Line(xray)
                            lineleft = lineright = true
                        else
                            xray = coord(r)
                        end
                    else
                        if y_left
                            yray = coord(r)
                        end
                    end
                end
            end
        else
            checkleftright(r)
        end
    end
    if line === nothing
        if xray === nothing
            sweep_norm = rotate(basis(eltype(points), d, 1))
            if iszero(sum(abs, sweep_norm))
                sweep_norm = rotate(basis(eltype(points), d, 2))
            end
        else
            sweep_norm = rotate(xray)
        end
    else
        sweep_norm = rotate(coord(line))
    end
    sort!(points, by = Base.Fix2(dot, sweep_norm))

    _points = eltype(points)[]
    _lines = eltype(lines)[]
    _rays = eltype(rays)[]
    if line === nothing
        half_points = _semi_hull(points, 1, counterclockwise, sweep_norm, yray)
        if yray === nothing
            other_half = _semi_hull(points, -1, counterclockwise, sweep_norm, yray)[2:end-1]
            if !isempty(other_half) && length(half_points) > 1 &&
                _colinear(counterclockwise, half_points[1], half_points[2], other_half[end])
                start = 2
            else
                start = 1
            end
            if !isempty(other_half) && length(half_points) > 1 &&
                _colinear(counterclockwise, half_points[end], half_points[end - 1], other_half[1])
                stop = length(half_points) - 1
            else
                stop = length(half_points)
            end
            append!(_points, half_points[start:stop])
            append!(_points, other_half)
        else
            append!(_points, half_points)
            push!(_rays, Ray(xray))
            if !(Ray(xray) ≈ Ray(yray))
                push!(_rays, Ray(yray))
            end
        end
    else
        push!(_lines, line)
        push!(_points, first(points))
        if lineleft
            if lineright
                push!(_lines, Line(sweep_norm))
            else
                push!(_rays, Ray(sweep_norm))
            end
        elseif lineright
            push!(_rays, Ray(-sweep_norm))
        else
            if !(dot(first(points), sweep_norm) ≈ dot(last(points), sweep_norm))
                push!(_points, last(points))
            end
        end
    end
    return _points, _lines, _rays
end

counterclockwise(x, y) = x[1] * y[2] - x[2] * y[1]
rotate(x) = convert(typeof(x), StaticArrays.SVector(-x[2], x[1]))

function planar_hull(vr::VRepresentation)
    d = typed_fulldim(vr)
    vrep(_planar_hull(typed_fulldim(vr), collect(points(vr)), lines(vr), rays(vr), counterclockwise, rotate)...; d = d)
end

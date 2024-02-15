function _semi_hull(ps::Vector{PT}, sign_sense, counterclockwise, sweep_norm, yray::Nothing = nothing) where PT
    hull = PT[]
    if isempty(ps)
        return hull
    end
    prev = sign_sense == 1 ? first(ps) : last(ps)
    cur = prev
    flat_start = true
    # Invariant:
    # We either have:
    # * `hull` is empty and `cur == prev` or
    # * `hull` is nonempty and `prev = last(hull)`.
    # In any case, the semihull of `ps[1:(j-1)]` is given by `[hull; cur]`.
    for j in (sign_sense == 1 ? (2:length(ps)) : ((length(ps)-1):-1:1))
        skip = false
        flat = false
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
                flat = isapproxzero(dcur) && isapproxzero(dpsj)
                if flat
                    # Case 1
                    # There can be two flat plateaus. The first one at the start and then one at the end
                    # `flat_start` indicates in which case we are. We need to a different direction in both cases
                    sense = flat_start ? sign_sense : -sign_sense
                    ccjsweep = sense * counterclockwise(psj_vec, sweep_norm)
                    cccursweep = sense * counterclockwise(cur_vec, sweep_norm)
                    if cccursweep < 0 < ccjsweep
                        # In this case, `prev` -> `cur` was in the right direction but now that
                        # we discover `ps[j]`, we realize that `prev` is in the interval
                        # `[ps[j], cur]`. Even if `ps[j]` -> `cur` seems to be in the right direction, we cannot
                        # keep it because it is in the wrong order in `ps`, it will be picked up by the other
                        # call to `_semi_hull` in which case `ps[j]` should be skipped (this is typically the case for the first flat plateau)
                        # For the second (and last flat plateau), it will rather be the case that `cur` should be removed.
                        # The only thing that is sure here is that `prev` should be removed so we remove `prev` and we `continue`
                        # as the code should then automatically decide which of `ps[j]` or `cur` should be removed
                        pop!(hull)
                        prev = cur
                        if !isempty(hull)
                            prev = last(hull)
                        end
                        continue
                    elseif cccursweep < ccjsweep
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
        if prev != cur && !flat
            flat_start = false
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
    d = FullDim(vr)
    vrep(_planar_hull(FullDim(vr), collect(points(vr)), lines(vr), rays(vr), counterclockwise, rotate)...; d = d)
end

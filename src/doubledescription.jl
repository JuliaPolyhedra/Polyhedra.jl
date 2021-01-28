export doubledescription
# Naive implementation of the double description Algorithm
# See JuliaPolyhedra/ConvexHull.jl for an efficient implementation

polytypefor(::Type{T}) where {T <: Integer} = Rational{BigInt}
polytypefor(::Type{Float32}) = Float64
polytypefor(::Type{T}) where {T} = T

polyvectortype(a) = a
# TODO sparse halfspaces does not mean sparse points
polyvectortype(::Type{<:SparseVector{T}}) where T = Vector{T}

dualtype(RepT::Type{<:Representation}) = dualtype(RepT, polyvectortype(vectortype(RepT)))
function dualfullspace(rep::Representation, d::FullDim, ::Type{T}) where T
    dualfullspace(rep, d, T, polyvectortype(similar_type(vectortype(rep), d, T)))
end
function dualfullspace(rep::Representation{T}) where T
    dualfullspace(rep, FullDim(rep), polytypefor(T))
end

"""
    doubledescription(h::HRepresentation)

Computes the V-representation of the polyhedron represented by `h` using the Double-Description algorithm [1, 2].

    doubledescription(V::VRepresentation)

Computes the H-representation of the polyhedron represented by `v` using the Double-Description algorithm [1, 2].

[1] Motzkin, T. S., Raiffa, H., Thompson, G. L. and Thrall, R. M.
The double description method
*Contribution to the Theory of Games*, *Princeton University Press*, **1953**

[2] Fukuda, K. and Prodon, A.
Double description method revisited
*Combinatorics and computer science*, *Springer*, **1996**, 91-111
"""
function doubledescription end

function print_v_summary(v::VRep)
    print("$(npoints(v)) points, $(nrays(v)) rays and $(nlines(v)) lines")
end

function intersect_and_remove_redundancy(v, hs, h; verbose=0)
    if eltype(hs) <: HalfSpace
        str = "halfspace"
    else
        str = "hyperplane"
    end
    for (i, hel) in enumerate(hs)
        if verbose >= 1
            println("Intersecting $str $i/$(length(hs))")
        end
        v_int = v ∩ hel
        if verbose >= 3
            print("Removing duplicates: ")
            print_v_summary(v_int)
            println(".")
        end
        # removeduplicates is cheaper than removevredundancy since the latter
        # needs to go through all the hrep element
        # FIXME not sure what to do here but it will be revamped by
        #       https://github.com/JuliaPolyhedra/Polyhedra.jl/pull/195 anyway
        v_dup = removeduplicates(v_int, OppositeMockOptimizer)
        if verbose >= 3
            print("Removing redundancy: ")
            print_v_summary(v_dup)
            println(".")
        end
        v = removevredundancy(v_dup, h)
        if verbose >= 2
            print("After intersection:  ")
            print_v_summary(v)
            println(".")
        end
    end
    return v
end

function slow_doubledescription(h::HRepresentation, solver = nothing; kws...)
    # The redundancy of V-elements are measured using
    # the number of hyperplane they are in. If there are
    # redundant hyperplanes, it does not matter since no point
    # will be inside them but if there are duplicates it is a problem

    # FIXME Note that removevredundancy, uses `h` which contains all hyperplanes and halfspaces
    # already taken into account but also all the other ones. We should check that this
    # is the right idea.

    # FIXME not sure what to do here but it will be revamped by
    #       https://github.com/JuliaPolyhedra/Polyhedra.jl/pull/195 anyway
    h = removeduplicates(h, solver === nothing ? OppositeMockOptimizer : solver)
    v = dualfullspace(h)
    checkvconsistency(v)
    v = intersect_and_remove_redundancy(v, hyperplanes(h), h; kws...)
    v = intersect_and_remove_redundancy(v, halfspaces(h), h; kws...)
    return v
end

struct CutoffPointIndex
    cutoff::Int
    index::Int
end
Base.show(io::IO, p::CutoffPointIndex) = print(io, "p[$(p.cutoff), $(p.index)]")
struct CutoffRayIndex
    cutoff::Int
    index::Int
end
Base.show(io::IO, r::CutoffRayIndex) = print(io, "r[$(r.cutoff), $(r.index)]")
struct DoubleDescriptionData{PointT, RayT, LineT, HST}
    fulldim::Int
    halfspaces::Vector{HST}
    # Elements ordered by first halfspace cutting it off
    points::Vector{PointT}
    pz::Vector{BitSet}
    cutpoints::Vector{Vector{PointT}}
    cutpz::Vector{Vector{BitSet}}
    pin::Vector{Vector{CutoffPointIndex}}
    rays::Vector{RayT}
    rz::Vector{BitSet}
    cutrays::Vector{Vector{RayT}}
    cutrz::Vector{Vector{BitSet}}
    rin::Vector{Vector{CutoffRayIndex}}
    lines::Vector{LineT}
    cutline::Vector{Union{Nothing, LineT}}
    lineray::Vector{Union{Nothing, CutoffRayIndex}}
    nlines::Vector{Int}
end
function Base.show(io::IO, data::DoubleDescriptionData)
    println(io, "DoubleDescriptionData in $(data.fulldim) dimension:")
    println(io, data.points)
    println(io, data.rays)
    println(io, data.lines)
    for i in reverse(eachindex(data.cutpoints))
        println(io, " Halfspace $i: $(data.halfspaces[i]):")
        if !isempty(data.cutpoints[i])
            println(io, "  Cut points:")
            for j in eachindex(data.cutpoints[i])
                println(io, "  $j: ", data.cutpoints[i][j], " zero at: ", data.cutpz[i][j])
            end
        end
        if !isempty(data.pin[i])
            println(io, "  In: ", data.pin[i])
        end
        if !isempty(data.cutrays[i])
            println(io, "  Cut rays:")
            for j in eachindex(data.cutrays[i])
                println(io, "  $j: ", data.cutrays[i][j], " zero at: ", data.cutrz[i][j])
            end
        end
        if !isempty(data.rin[i])
            println(io, "  In: ", data.rin[i])
        end
        if data.cutline[i] !== nothing
            println(io, "  Cut line: ", data.cutline[i])
            if data.lineray[i] !== nothing
                println(io, "  Line ray: ", data.lineray[i])
            end
        end
        if !iszero(data.nlines[i])
            println(io, "  $(data.nlines[i]) uncut lines left")
        end
    end
end

function DoubleDescriptionData{PointT, RayT, LineT}(fulldim::Integer, hyperplanes, halfspaces) where {PointT, RayT, LineT}
    n = length(halfspaces)
    m = length(hyperplanes)
    return DoubleDescriptionData{PointT, RayT, LineT, eltype(halfspaces)}(
        fulldim,
        halfspaces,
        PointT[],
        BitSet[],
        [PointT[] for i in 1:n],
        [BitSet[] for i in 1:n],
        [CutoffPointIndex[] for i in 1:n],
        RayT[],
        BitSet[],
        [RayT[] for i in 1:n],
        [BitSet[] for i in 1:n],
        [CutoffRayIndex[] for i in 1:n],
        LineT[],
        Union{Nothing, LineT}[nothing for i in 1:(m + n)],
        Union{Nothing, CutoffRayIndex}[nothing for i in 1:n],
        zeros(Int, m + n)
    )
end
function tight_halfspace_indices(data::DoubleDescriptionData, p::CutoffPointIndex)
    if iszero(p.cutoff)
        return data.pz[p.index]
    else
        return data.cutpz[p.cutoff][p.index]
    end
end
function tight_halfspace_indices(data::DoubleDescriptionData, r::CutoffRayIndex)
    if iszero(r.cutoff)
        return data.rz[r.index]
    else
        return data.cutrz[r.cutoff][r.index]
    end
end
function Base.getindex(data::DoubleDescriptionData, p::CutoffPointIndex)
    if iszero(p.cutoff)
        return data.points[p.index]
    else
        return data.cutpoints[p.cutoff][p.index]
    end
end
function Base.getindex(data::DoubleDescriptionData, r::CutoffRayIndex)
    if iszero(r.cutoff)
        return data.rays[r.index]
    else
        return data.cutrays[r.cutoff][r.index]
    end
end

function _bitdot_range(b1::BitSet, b2::BitSet, i, n)
    count = 1 # They share the hyperplance `i`
    for j in (i + 1):n
        if j in b1 && j in b2
            count += 1
        end
    end
    return count
end
_ray_pair(p1, p2) = false
_ray_pair(::CutoffRayIndex, ::CutoffRayIndex) = true
# Necessary condition for adjacency.
# See Proposition 9 (NC1) of [FP96].
function is_adjacent_nc1(data, i, p1, p2)
    rhs = 1 + _ray_pair(p1, p2) + data.nlines[i] + _bitdot_range(
        tight_halfspace_indices(data, p1),
        tight_halfspace_indices(data, p2),
        i, length(data.halfspaces)
    )
    return data.fulldim >= rhs
end
function is_adjacency_breaker(data, i, p, p1, p2)
    return p != p1 && p != p2 && all((i + 1):length(data.halfspaces)) do j
        isin(data, j, p) || !(isin(data, j, p1) && isin(data, j, p2))
    end
end
function isadjacent(data, i::Integer, p1, p2)
    return is_adjacent_nc1(data, i, p1, p2) &&
        # According to Proposition 7 (c) of [FP96], we need to check
        # that there is not other point or ray that is in the same hyperplanes as `p1` and `p2`.
        # If it's the case, it is in hyperplane `i` in particular,
        # we can use it to restrict our attention to points and
        # rays in `pin[i]` and `rin[i]`.
        (_ray_pair(p1, p2) || !any(p -> is_adjacency_breaker(data, i, p, p1, p2), data.pin[i])) &&
        !any(r -> is_adjacency_breaker(data, i, r, p1, p2), data.rin[i])
end
isin(data, i, p) = i in tight_halfspace_indices(data, p)

resized_bitset(data) = sizehint!(BitSet(), length(data.halfspaces))
function add_index!(data, cutoff::Nothing, p::AbstractVector, tight::BitSet)
    push!(data.points, p)
    push!(data.pz, tight)
    return CutoffPointIndex(0, length(data.points))
end
function add_index!(data, cutoff::Integer, p::AbstractVector, tight::BitSet)
    push!(data.cutpoints[cutoff], p)
    push!(data.cutpz[cutoff], tight)
    return CutoffPointIndex(cutoff, length(data.cutpoints[cutoff]))
end
function add_index!(data, cutoff::Nothing, r::Polyhedra.Ray, tight::BitSet)
    push!(data.rays, r)
    push!(data.rz, tight)
    return CutoffRayIndex(0, length(data.rays))
end
function add_index!(data, cutoff::Integer, r::Polyhedra.Ray, tight::BitSet)
    push!(data.cutrays[cutoff], r)
    push!(data.cutrz[cutoff], tight)
    return CutoffRayIndex(cutoff, length(data.cutrays[cutoff]))
end

function add_in!(data, i, index::CutoffPointIndex)
    push!(data.pin[i], index)
end
function add_in!(data, i, index::CutoffRayIndex)
    push!(data.rin[i], index)
end
function set_in!(data, I, index)
    for i in I
        if isin(data, i, index)
            add_in!(data, i, index)
        end
    end
end
function add_element!(data, k, el, tight)
    cutoff = nothing
    for i in reverse(1:k)
        if data.cutline[i] !== nothing
            el = line_project(el, data.cutline[i], data.halfspaces[i])
            # The line is in all halfspaces from `i+1` up so projecting with it does not change it.
            push!(tight, i)
            index = add_adjacent_element!(data, i - 1, el, data.lineray[i], tight)
            set_in!(data, i:k, index)
            return index
        end
        # Could avoid computing `dot` twice between `el` and the halfspace here.
        if !(el in data.halfspaces[i])
            cutoff = i
            break
        end
        if el in hyperplane(data.halfspaces[i])
            push!(tight, i)
        end
    end
    index = add_index!(data, cutoff, el, tight)
    set_in!(data, (index.cutoff + 1):k, index)
    return index
end
function project_onto_affspace(data, offset, el, hyperplanes)
    for i in reverse(eachindex(hyperplanes))
        line = data.cutline[offset + i]
        h = hyperplanes[i]
        if line !== nothing
            el = line_project(el, line, h)
        elseif !(el in h)
            # e.g. 0x1 + 0x2 = 1 or -1 = x1 = 1, ...
            return nothing
        end
    end
    return el
end
function add_adjacent_element!(data, k, el, parent, tight)
    index = add_element!(data, k, el, tight)
    addintersection!(data, index, parent, nothing)
    return index
end

function combine(β, p1::AbstractVector, value1, p2::AbstractVector, value2)
    λ = (value2 - β) / (value2 - value1)
    return λ * p1 + (1 - λ) * p2
end
function combine(β, p::AbstractVector, pvalue, r::Polyhedra.Ray, rvalue)
    λ = (β - pvalue) / rvalue
    return p + λ * r
end
combine(β, r::Polyhedra.Ray, rvalue, p::AbstractVector, pvalue) = combine(β, p, pvalue, r, rvalue)
function combine(β, r1::Polyhedra.Ray, value1, r2::Polyhedra.Ray, value2)
    # should take
    # λ = value2 / (value2 - value1)
    @assert 0 <= value2 / (value2 - value1) <= 1
    # By homogeneity we can avoid the division and do
    #newr = value2 * r1 - value1 * r2
    # but this can generate very large numbers (see JuliaPolyhedra/Polyhedra.jl#48)
    # so we still divide
    newr = (value2 * r1 - value1 * r2) / (value2 - value1)
    # In CDD, it does value2 * r1 - value1 * r2 but then it normalize the ray
    # by dividing it by its smallest nonzero entry (see dd_CreateNewRay)
    return Polyhedra.simplify(newr)
end

combine(h, el1, el2) = combine(h.β, el1, h.a ⋅ el1, el2, h.a ⋅ el2)

"""
    addintersection!(data, idx1, idx2, hp_idx)

`hp_idx === nothing` means inherited adjacency, otherwise, it is
an index such that `idx1` and `idx2` are both in the
hyperplane `hp_idx`.
"""
function addintersection!(data, idx1, idx2, hp_idx)
    if idx1.cutoff > idx2.cutoff
        return addintersection!(data, idx2, idx1, hp_idx)
    end
    i = idx2.cutoff
       # Condition (c_k) in [FP96]
    if idx1.cutoff == idx2.cutoff ||
        isin(data, i, idx1) ||
        # If it's in both at some lower `j` then we'll
        # add it then otherwise, it will be added twice.
        any(j -> isin(data, j, idx1) && isin(data, j, idx2), (idx2.cutoff + 1):(i - 1)) ||
        (hp_idx !== nothing && !isadjacent(data, hp_idx, idx1, idx2))
        return
    end
    newel = combine(data.halfspaces[i], data[idx1], data[idx2])
    # `newel` and `idx1` have inherited adjacency, see 3.2 (i) of [FP96].
    # TODO are we sure that they are adjacent ?
    #      or should we check rank or combinatorial check ?
    #      In `DDMethodVariation2` of [FP96], `Adj` is not called
    #      in case `rj, rj'` have inherited adjacency
    tight = tight_halfspace_indices(data, idx1) ∩ tight_halfspace_indices(data, idx2)
    push!(tight, i)
    index = add_adjacent_element!(data, i - 1, newel, idx1, tight)
    set_in!(data, i:i, index)
    return index
end

_shift(el::AbstractVector, line::Line) = el + Polyhedra.coord(line)
_shift(el::Line, line::Line) = el + line
_shift(el::Ray, line::Line) = el + Polyhedra.Ray(Polyhedra.coord(line))
function _λ_proj(r::VStruct, line::Line, h::HRepElement)
    # Line or ray `r`:
    # (r + λ * line) ⋅ h.a == 0
    # λ = -r ⋅ h.a / (line ⋅ h.a)
    return -r ⋅ h.a / (line ⋅ h.a)
end
function _λ_proj(x::AbstractVector, line::Line, h::HRepElement)
    # Point `x`:
    # (x + λ * line) ⋅ h.a == h.β
    # λ = (h.β - x ⋅ h.a) / (line ⋅ h.a)
    return (h.β - x ⋅ h.a) / (line ⋅ h.a)
end
function line_project(el, line, h)
    λ = _λ_proj(el, line, h)
    return Polyhedra.simplify(_shift(el, λ * line))
end
function hline(data, line::Line, i, h)
    value = h.a ⋅ line
    if !Polyhedra.isapproxzero(value)
        if data.cutline[i] === nothing
            # Make `lineray` point inward
            data.cutline[i] = value > 0 ? -line : line
            return true, line
        else
            line = line_project(line, data.cutline[i], h)
        end
    end
    data.nlines[i] += 1
    return false, line
end
# TODO remove solver arg `_`, it is kept to avoid breaking code
function doubledescription(hr::HRepresentation, _ = nothing)
    v = Polyhedra.dualfullspace(hr)
    hps = Polyhedra.lazy_collect(hyperplanes(hr))
    hss = Polyhedra.lazy_collect(halfspaces(hr))
    data = DoubleDescriptionData{pointtype(v), raytype(v), linetype(v)}(fulldim(hr), hps, hss)
    for line in lines(v)
        cut = false
        for i in reverse(eachindex(hps))
            cut, line = hline(data, line, nhalfspaces(hr) + i, hps[i])
            cut && break
        end
        if !cut
            for i in reverse(eachindex(hss))
                cut, line = hline(data, line, i, hss[i])
                cut && break
            end
        end
        if !cut
            push!(data.lines, line)
        end
    end
    # Add line rays after all lines are added so that the rays can be `line_project`ed.
    # We only do that for halfspaces, hyperplanes do not create rays from cutoff lines.
    # We use increasing index order since higher index may need the `lineray` of lower index.
    for i in eachindex(hss)
        line = data.cutline[i]
        if line !== nothing
            ray = Polyhedra.Ray(Polyhedra.coord(line))
            data.lineray[i] = add_element!(data, i - 1, ray, BitSet((i + 1):length(hss)))
        end
    end
    @assert isone(npoints(v))
    # Add the origin
    orig = project_onto_affspace(data, nhalfspaces(hr), first(points(v)), hps)
    if orig !== nothing
        add_element!(data, nhalfspaces(hr), orig, resized_bitset(data))
    end
    for i in reverse(eachindex(hss))
        data.cutline[i] !== nothing && continue
        if isempty(data.cutpoints[i]) && isempty(data.cutrays[i])
            # Redundant, remove its contribution to avoid incorrect `isadjacent`
            for p in data.pin[i]
                delete!(tight_halfspace_indices(data, p), i)
            end
            for r in data.rin[i]
                delete!(tight_halfspace_indices(data, r), i)
            end
            continue
        end
        if i > 1
            # Catches new adjacent rays, see 3.2 (ii) of [FP96]
            for p1 in data.pin[i], p2 in data.pin[i]
                if p1.cutoff < p2.cutoff
                    addintersection!(data, p1, p2, i)
                end
            end
            for p in data.pin[i], r in data.rin[i]
                addintersection!(data, p, r, i)
            end
        end
        deleteat!(data.cutpoints, i)
        deleteat!(data.cutpz, i)
        if i > 1
            for r1 in data.rin[i], r2 in data.rin[i]
                # We encounter both `r1, r2` and `r2, r1`.
                # Break this symmetry with:
                if r1.cutoff < r2.cutoff
                    addintersection!(data, r1, r2, i)
                end
            end
        end
        deleteat!(data.cutrays, i)
        deleteat!(data.cutrz, i)
        deleteat!(data.pin, i)
        deleteat!(data.rin, i)
    end
    if isempty(data.points)
        # Empty polyhedron, there may be rays left,
        # Example 1: for 0x_1 + x_2 = -1 ∩ 0x_1 + x_2 = 1, the line (0, 1) is detected as correct
        # Example 2: for 0x_1 + 0x_2 = 1, the lines (1, 0) and (0, 1) are detected as correct
        # but since there is no point, the polyhedron is empty and we should drop all rays/lines
        empty!(data.lines)
        empty!(data.rays)
    end
    return similar(v, data.points, data.lines, data.rays)
end

function doubledescription(v::VRepresentation{T}, solver = nothing; kws...) where {T}
    checkvconsistency(v)
    lv = convert(LiftedVRepresentation{T, Matrix{T}}, v)
    R = -lv.R
    vl = doubledescription(MixedMatHRep{T}(R, zeros(T, size(R, 1)), lv.linset), solver; kws...)
    LiftedHRepresentation{T}(vl.R, vl.Rlinset)
end

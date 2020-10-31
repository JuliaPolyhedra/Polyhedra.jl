export triangulation, triangulation_indices

# Implementation of Cohen & Hickey algorithm for triangulation; see [CH79] and [BEK00].
#
# [CH79] Cohen, J., & Hickey, T. (1979). Two algorithms for determining volumes of convex polyhedra. Journal of the ACM (JACM), 26(3), 401-414.
# [BEK00] Büeler, B., Enge, A., & Fukuda, K. (2000). Exact volume computation for polytopes: a practical study. In Polytopes—combinatorics and computation (pp. 131-154). Birkhäuser, Basel.

# TODO parallelize: this is inherently parallelizable
function _triangulation(Δs, Δ, v_idx, h_idx, incident_idx, is_weak_adjacent, codim)
    @assert codim >= 0
    isempty(v_idx) && return
    v = first(v_idx)
    Δ = push!(copy(Δ), v)
    if isone(length(v_idx))
        # We should have `codim == 0` whenever we reach this point.  Due to
        # numerical issues, in cases where there is very near but not exact
        # vertex redundancy, we sometimes have `codim > 0`.  In these cases,
        # mathematically, `Δ` is a near-degenerate simplex (I think?), and our
        # numerically-imperfect version of `Δ` has too few vertices.  Rather
        # than cause dimension errors downstream, we simply omit such a `Δ`.
        if codim == 0
            push!(Δs, Δ)
        end
        return
    end
    tail = true
    for h in h_idx
        if !(v in incident_idx[h])
            tail = false
            # The adjacency may be outside the current face under scrutiny but that's ok,
            # it will simply end up calling `_triangulation` with an empty `v_idx`.
            weak_adjacent = [hj for hj in h_idx if hj != h && is_weak_adjacent[(h, hj)]]
            active = [point for point in v_idx if point in incident_idx[h]]
            _triangulation(Δs, Δ, active, weak_adjacent, incident_idx, is_weak_adjacent, codim - 1)
        end
    end
end
function triangulation_indices(p::Polyhedron)
    hasrays(p) && error("Triangulation only supported for polytope.")
    v_idx = eachindex(points(p))
    h_idx = eachindex(halfspaces(p))
    Δ = eltype(v_idx)[]
    Δs = typeof(Δ)[]
    incident_idx = Dict(h => Set(incidentpointindices(p, h)) for h in h_idx)
    is_weak_adjacent = Dict((hi, hj) => !isempty(incident_idx[hi] ∩ incident_idx[hj]) for hi in h_idx for hj in h_idx)
    _triangulation(Δs, Δ, v_idx, h_idx, incident_idx, is_weak_adjacent, fulldim(p))
    return Δs
end
function triangulation(p::Polyhedron)
    return map(Δ -> vrep(get.(p, Δ)), triangulation_indices(p))
end

export supportselimination, eliminate, project, fixandeliminate
export DefaultElimination, FourierMotzkin, BlockElimination, ProjectGenerators

abstract type EliminationAlgorithm end

struct DefaultElimination <: EliminationAlgorithm end

"""
    FourierMotzkin

Computation of the projection by computing the H-representation and applying the
Fourier-Motzkin elimination algorithm to it.

The resulting polyhedron may contain duplicates and/or redundancy. These can be
removed with [`removehredundancy!`](@ref).
"""
struct FourierMotzkin <: EliminationAlgorithm end

"""
    BlockElimination

Computation of the projection by computing the H-representation and applying
the block elimination algorithm to it.

The idea is as follows.
Consider a H-representation given by ``Ax + By \\le c``.
For any ``d``, the projection on the `x` variables
belongs to the halfspace ``d^\\top x \\le \\delta`` where ``\\delta``
is
```math
\\begin{aligned}
\\max \\;     & d^\\top x \\\\
\\text{s.t.} & A x + By \\le c
\\end{aligned}
```
By duality, this is equal to
```math
\\begin{aligned}
\\min \\;     & c^\\top \\nu \\\\
\\text{s.t.} & A^\\top \\nu = d \\\\
& B^\\top \\nu = 0 \\\\
& \\nu \\ge 0
\\end{aligned}
```
Let `R` be the matrix for which the columns are the extreme rays of the
polyhedral cone defined by
``B^\\top \\nu = 0, \\nu \\ge 0``.
The program can be rewritten as
```math
\\begin{aligned}
\\min \\;     & c^\\top \\nu \\\\
\\text{s.t.} & A^\\top \\nu = d \\\\
& \\nu = R\\lambda \\\\
& \\lambda \\ge 0
\\end{aligned}
```
which is equivalent to
```math
\\begin{aligned}
\\min \\;     & c^\\top R \\lambda \\\\
\\text{s.t.} & A^\\top R \\lambda = d \\\\
& \\lambda \\ge 0
\\end{aligned}
```
Dualizing, we obtain
```math
\\begin{aligned}
\\max \\;     & d^\\top x \\\\
\\text{s.t.} & R^\\top A x \\le R^\\top c
\\end{aligned}
```
where we see that ``R^\\top A x \\le R^\\top c``
is the H-representation of the projection.
The hyperplanes ``A x + B y = c`` are handled similarly, their multipliers
``\\nu`` are free so that the cone ``B^\\top \\nu = 0`` may have lines.
Each line ``\\ell`` of the cone gives the hyperplane
``\\ell^\\top A x = \\ell^\\top c`` of the projection.

The extreme rays and lines of the cone are computed by the library of the
polyhedron. The resulting polyhedron may contain redundancy. It can be removed
with [`removehredundancy!`](@ref).
"""
struct BlockElimination <: EliminationAlgorithm end

"""
    ProjectGenerators

Computation of the projection by computing the V-representation and projecting them.
"""
struct ProjectGenerators <: EliminationAlgorithm end

"""
    project(p::Polyhedron, pset, algo)

Equivalent to `eliminate(p, setdiff(1:fulldim(p), pset), algo)`.
"""
function project end

"""
    eliminate(p::Polyhedron, delset, algo::EliminationAlgorithm)

Eliminate the dimensions in `delset` by projecting the polyhedron onto the remaining dimension.
"""
function eliminate end

project(p::Polyhedron, pset) = eliminate(p, setdiff(1:fulldim(p), pset))
function project(p::Polyhedron, pset, algo)
    return eliminate(p, setdiff(1:fulldim(p), pset), algo)
end

supportselimination(p::Polyhedron, ::FourierMotzkin) = false
function eliminate(p::Polyhedron, delset, ::FourierMotzkin)
    return error("Fourier-Motzkin elimination not implemented for $(typeof(p))")
end
# The generic implementation of `eliminate` below can be used with any library
# but it is not used by `DefaultElimination` since the resulting H-representation
# may contain much more redundant halfspaces than the one obtained by
# `ProjectGenerators`. A library can opt in by returning `true`.
supportselimination(p::Polyhedron, ::BlockElimination) = false

function _combine_rows(::Type{U}, ν, elements, keep) where {U}
    a = zeros(U, length(keep))
    β = zero(U)
    for (νi, el) in zip(ν, elements)
        a .+= νi .* el.a[keep]
        β += νi * el.β
    end
    return a, β
end
_isapproxzero(a::AbstractVector; tol) = isempty(a) || isapproxzero(a; tol)

function eliminate(p::Polyhedron{T}, delset, ::BlockElimination) where {T}
    h = hrep(p)
    N = fulldim(h)
    keep = setdiff(1:N, delset)
    del = setdiff(1:N, keep)
    hps = collect(hyperplanes(h))
    hss = collect(halfspaces(h))
    elements = [hps; hss]
    m = length(elements)
    d = length(keep)
    if iszero(m)
        return similar(
            p,
            d,
            T,
            HyperPlane{T,Vector{T}}[],
            HalfSpace{T,Vector{T}}[];
            dimension_map = i -> (1, keep[i]),
        )
    end
    # Cone of the multipliers `ν` such that `B' ν = 0` and `ν ≥ 0` for the
    # multipliers corresponding to halfspaces
    cone = hrep(
        HyperPlane{T,Vector{T}}[
            HyperPlane(T[el.a[j] for el in elements], zero(T)) for j in del
        ],
        HalfSpace{T,Vector{T}}[
            HalfSpace(T[-(i == length(hps) + k) for i in 1:m], zero(T)) for
            k in eachindex(hss)
        ];
        d = m,
    )
    v = vrep(polyhedron(cone, similar_library(p, m, T)))
    U = coefficient_type(v)
    AT = Vector{U}
    tol = _default_tol(U)
    new_hps = HyperPlane{U,AT}[]
    for l in lines(v)
        a, β = _combine_rows(U, coord(l), elements, keep)
        # `0 = 0` is trivial, `0 = β` with `β ≠ 0` means that `p` is empty
        if !(_isapproxzero(a; tol) && isapproxzero(β; tol))
            push!(new_hps, simplify(HyperPlane(a, β)))
        end
    end
    new_hss = HalfSpace{U,AT}[]
    for r in rays(v)
        a, β = _combine_rows(U, coord(r), elements, keep)
        # `0 ≤ β` with `β ≥ 0` is trivial, otherwise it means that `p` is empty
        if !(_isapproxzero(a; tol) && _nonneg(β; tol))
            push!(new_hss, simplify(HalfSpace(a, β)))
        end
    end
    return similar(p, d, U, new_hps, new_hss; dimension_map = i -> (1, keep[i]))
end

function eliminate(p::Polyhedron, algo::EliminationAlgorithm)
    return eliminate(p, BitSet(fulldim(p)), algo)
end

# eliminate the last dimension by default
function eliminate(p::Polyhedron, delset = BitSet(fulldim(p)))
    return eliminate(p, delset, DefaultElimination())
end

function eliminate(p::Polyhedron, delset, ::DefaultElimination)
    fm = supportselimination(p, FourierMotzkin())
    be = supportselimination(p, BlockElimination())
    if (!fm && !be) || vrepiscomputed(p)
        algo = ProjectGenerators()
    elseif fm && (!be || (length(delset) == 1 && fulldim(p) in delset))
        algo = FourierMotzkin()
    else
        algo = BlockElimination()
    end
    return eliminate(p, delset, algo)
end

function eliminate(p::Polyhedron, delset, ::ProjectGenerators)
    N = fulldim(p)
    Id = Matrix(1I, N, N)
    Iproj = Id[collect(setdiff(BitSet(1:N), delset)), :]
    return Iproj * p
end

"""
    project(p::Polyhedron, P::AbstractMatrix)

Projects the polyhedron `p` into the `size(P, 2)`-dimensional linear subspace spanned by the columns of `P`.
The new orthonormal basis for this subspace is computed by applying the Gram–Schmidt process to the columns of `P`, starting with the first column.
"""
function project(p::Polyhedron, P::AbstractMatrix)
    # Function to make x orthogonal to an orthonormal basis in Q
    # We first make the columns of P orthonormal
    N = fulldim(p)
    if size(P, 1) != N
        throw(
            DimensionMismatch(
                "The columns of P should have the same dimension than the polyhedron",
            ),
        )
    end
    m = size(P, 2)
    if m > N
        throw(DimensionMismatch("P should have more rows than columns"))
    end
    Q = Array{Float64}(P) # normalize will make it nonrational
    Proj = zeros(eltype(P), N, N)
    for i in 1:m
        Q[:, i] = normalize(Q[:, i] - Proj * Q[:, i])
        Proj += Q[:, i] * Q[:, i]'
    end
    if m == N
        basis = Q
    else
        # For the rest, we take the canonical basis and we look at
        # I - Proj * I
        Id = Matrix(1.0I, N, N)
        R = Id - Proj
        # We take the n-m that have highest norm
        order = sortperm([dot(R[:, i], R[:, i]) for i in 1:N])
        R = Id[:, order[(m+1):N]]
        for i in 1:(N-m)
            R[:, i] = normalize(R[:, i] - Proj * R[:, i])
            Proj += R[:, i] * R[:, i]'
        end
        basis = [Q R]
    end
    return eliminate(basis \ p, BitSet((m+1):N))
end

function _fixelim(h::ElemT, I, J, v, dout::FullDim) where {ElemT<:HRepElement}
    return similar_type(ElemT, dout)(h.a[J], h.β - dot(h.a[I], v))
end

"""
    fixandeliminate(p::HRep{T}, I, v)

Fix the variables with indices in `I` to the corresponding value in `v`. This is equivalent to doing the following:
```julia
function ei(i)
    a = zeros(T, fulldim(p))
    a[i] = one(T)
    a
end
eliminate(p ∩ HyperPlane(ei(I[1]), v[1]) ∩ ... ∩ HyperPlane(ei(I[n]), v[n]))
```
where `n` is the length of `I` (and `v`), but it is much more efficient. The code above does a polyhedral projection while this function simply replaces
each halfspace `⟨a, x⟩ ≤ β` (resp. each hyperplane `⟨a, x⟩ = β`) by the halfspace `⟨a_J, x⟩ ≤ β - ⟨a_I, v⟩`
(resp. the hyperplane `⟨a_J, x⟩ = β - ⟨a_I, v⟩`) where `J = setdiff(1:fulldim(p), I)`.
"""
function fixandeliminate(p::HRep{S}, I, v) where {S}
    J = setdiff(1:fulldim(p), I)
    d = length(J)
    f = (i, h) -> _fixelim(h, I, J, v, d)
    T = promote_type(S, eltype(v))
    return similar(p, d, T, hmap(f, d, T, p)...)
end

# TODO rewrite, it is just cutting a cone with a half-space, nothing more
# export radialprojectoncut
# function radialprojectoncut(p::Polyhedron, cut::Vector, at)
#   if isapproxzero(at)
#     error("at is zero")
#   end
#   if length(cut) != fulldim(p)
#     error("The dimensions of the cut and of the polyhedron do not match")
#   end
#   ext = MixedMatVRep(getgenerators(p))
#   V = copy(ext.V)
#   R = copy(ext.R)
#   for i in 1:size(V, 1)
#     v = vec(ext.V[i,:])
#     if !myeq(dot(cut, v), at)
#       error("The nonhomogeneous part should be in the cut")
#     end
#   end
#   for i in 1:size(R, 1)
#     v = vec(ext.R[i,:])
#     if isapproxzero(v)
#       # It can happen since I do not necessarily have removed redundancy
#       v = zeros(eltype(v), length(v))
#     elseif !myeq(dot(cut, v), at)
#       if isapproxzero(dot(cut, v))
#         error("A ray is parallel to the cut") # FIXME is ok if some vertices are on the cut ? (i.e. at == 0, cut is not needed)
#       end
#       v = v * at / dot(cut, v)
#     end
#     R[i,:] = v
#   end
#   # no more rays nor linearity since at != 0
#   ext2 = MixedMatVRep([V; R])
#   polyhedron(ext2, getlibraryfor(p, eltype(ext2)))
# end

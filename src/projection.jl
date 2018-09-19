export supportselimination, eliminate, project, fixandeliminate
export DefaultElimination, FourierMotzkin, BlockElimination, ProjectGenerators

abstract type EliminationAlgorithm end

struct DefaultElimination <: EliminationAlgorithm end

"""
    FourierMotzkin

Computation of the projection by computing the H-representation and applying the Fourier-Motzkin elimination algorithm to it.
"""
struct FourierMotzkin <: EliminationAlgorithm end

"""
    BlockElimination

Computation of the projection by computing the H-representation and applying the block elimination algorithm to it.
"""
struct BlockElimination <: EliminationAlgorithm end

"""
    ProjectGenerators

Computation of the projection by computing the V-representation and projecting them.
"""
struct ProjectGenerators <: EliminationAlgorithm end

"""
    project(p::Polyhedron, pset, algo)

    Equivalent to `eliminate(p, setdiff(1:fulldim(p), pset), algo).
"""
function project end

"""
    eliminate(p::Polyhedron, delset, algo::EliminationAlgorithm)

Eliminate the dimensions in `delset` by projecting the polyhedron onto the remaining dimension.
"""
function eliminate end

project(p::Polyhedron, pset) = eliminate(p, setdiff(1:fulldim(p), pset))
project(p::Polyhedron, pset, algo) = eliminate(p, setdiff(1:fulldim(p), pset), algo)

supportselimination(p::Polyhedron, ::FourierMotzkin)   = false
eliminate(p::Polyhedron, delset, ::FourierMotzkin)     = error("Fourier-Motzkin elimination not implemented for $(typeof(p))")
supportselimination(p::Polyhedron, ::BlockElimination) = false
eliminate(p::Polyhedron, delset, ::BlockElimination)   = error("Block elimination not implemented for $(typeof(p))")

eliminate(p::Polyhedron, algo::EliminationAlgorithm) = eliminate(p, BitSet(fulldim(p)), algo)

# eliminate the last dimension by default
eliminate(p::Polyhedron, delset=BitSet(fulldim(p))) = eliminate(p, delset, DefaultElimination())

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
    eliminate(p, delset, algo)
end

function eliminate(p::Polyhedron, delset, ::ProjectGenerators)
    N = fulldim(p)
    Id = Matrix(1I, N, N)
    Iproj = Id[collect(setdiff(BitSet(1:N), delset)),:]
    Iproj * p
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
        throw(DimensionMismatch("The columns of P should have the same dimension than the polyhedron"))
    end
    m = size(P, 2)
    if m > N
        throw(DimensionMismatch("P should have more rows than columns"))
    end
    Q = Array{Float64}(P) # normalize will make it nonrational
    Proj = zeros(eltype(P), N, N)
    for i = 1:m
        Q[:,i] = normalize(Q[:,i] - Proj * Q[:,i])
        Proj += Q[:,i] * Q[:,i]'
    end
    if m == N
        basis = Q
    else
        # For the rest, we take the canonical basis and we look at
        # I - Proj * I
        Id = Matrix(1.0I, N, N)
        R = Id - Proj
        # We take the n-m that have highest norm
        order = sortperm([dot(R[:,i], R[:,i]) for i in 1:N])
        R = Id[:,order[m+1:N]]
        for i in 1:N-m
            R[:,i] = normalize(R[:,i] - Proj * R[:,i])
            Proj += R[:,i] * R[:,i]'
        end
        basis = [Q R]
    end
    eliminate(basis \ p, BitSet(m+1:N))
end

_fixelim(h::ElemT, I, J, v, dout::FullDim) where ElemT<:HRepElement = similar_type(ElemT, dout)(h.a[J], h.β - dot(h.a[I], v))

"""
    fixandeliminate(p::HRep{T}, I, v)

Fix the variables with indices in `I` to the corresponding value in `v`. This is equivalent to doing the following:
```julia
function ei(i)
    a = zeros(T, fulldim(p))
    a[i] = one(T)
    a
end
eliminate(p ∩ HyperPlane(ei(I[1], v[1]) ∩ ... ∩ HyperPlane(ei(I[1], v[1]))
```
but it is much more efficient. The code above does a polyhedral projection while this function simply replace
each halfspace `⟨a, x⟩ ≤ β` (resp. each hyperplane `⟨a, x⟩ = β`) by the halfspace `⟨a_J, x⟩ ≤ β - ⟨a_I, v⟩`
(resp. the hyperplane `⟨a_J, x⟩ = β - ⟨a_I, v⟩`) where `J = setdiff(1:fulldim(p), I)`.
"""
function fixandeliminate(p::HRep{S}, I, v) where {S}
    J = setdiff(1:fulldim(p), I)
    d = length(J)
    f = (i, h) -> _fixelim(h, I, J, v, d)
    T = promote_type(S, eltype(v))
    similar(p, d, T, hmap(f, d, T, p)...)
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

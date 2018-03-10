export implementseliminationmethod, eliminate, project, fixandeliminate

project(p::Polyhedron{N}, pset) where {N} = eliminate(p, setdiff(1:N, pset))
project(p::Polyhedron{N}, pset, method) where {N} = eliminate(p, setdiff(1:N, pset), method)

implementseliminationmethod(p::Polyhedron, ::Type{Val{:FourierMotzkin}})   = false
eliminate(p::Polyhedron, delset, ::Type{Val{:FourierMotzkin}})             = error("Fourier-Motzkin elimination not implemented for $(typeof(p))")
implementseliminationmethod(p::Polyhedron, ::Type{Val{:BlockElimination}}) = false
eliminate(p::Polyhedron, delset, ::Type{Val{:BlockElimination}})           = error("Block elimination not implemented for $(typeof(p))")

eliminate(p::Polyhedron, method::Symbol) = eliminate(p, Val{method})
eliminate(p::Polyhedron, delset, method::Symbol) = eliminate(p, delset, Val{method})

eliminate(p::Polyhedron{N}, method::Type{Val{:ProjectGenerators}}) where {N} = eliminate(p, IntSet(N), method)

# eliminate the last dimension by default
function eliminate(p::Polyhedron{N}, delset=IntSet(N)) where N
    fm = implementseliminationmethod(p, Val{:FourierMotzkin})
    be = implementseliminationmethod(p, Val{:BlockElimination})
    if (!fm && !be) || vrepiscomputed(p)
        method = :ProjectGenerators
    elseif fm && (!be || (length(delset) == 1 && N in delset))
        method = :FourierMotzkin
    else
        method = :BlockElimination
    end
    eliminate(p, delset, Val{method})
end

function eliminate(p::Polyhedron{N}, delset, ::Type{Val{:ProjectGenerators}}) where N
    I = eye(Int, N)
    Iproj = I[collect(setdiff(IntSet(1:N), delset)),:]
    Iproj * p
end

function project(p::Polyhedron{N,T}, P::AbstractMatrix) where {N,T}
    # Function to make x orthogonal to an orthonormal basis in Q
    # We first make the columns of P orthonormal
    if size(P, 1) != N
        error("The columns of P should have the same dimension than the polyhedron")
    end
    m = size(P, 2)
    if m > N
        error("P should have more columns than rows")
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
        I = eye(Float64, N)
        R = I - Proj
        # We take the n-m that have highest norm
        order = sortperm([dot(R[:,i], R[:,i]) for i in 1:N])
        R = I[:,order[m+1:N]]
        for i in 1:N-m
            R[:,i] = normalize(R[:,i] - Proj * R[:,i])
            Proj += R[:,i] * R[:,i]'
        end
        basis = [Q R]
    end
    eliminate(basis \ p, IntSet(m+1:N))
end

_fixelim(h::ElemT, I, J, v, dout::FullDim) where ElemT<:HRepElement = similar_type(ElemT, dout)(h.a[J], h.β - dot(h.a[I], v))

"""
    fixandeliminate(p::HRep{N, T}, I, v)

Fix the variables with indices in `I` to the corresponding value in `v`. This is equivalent to doing the following:
```julia
function ei(i)
    a = zeros(T, N)
    a[i] = one(T)
    a
end
eliminate(p ∩ HyperPlane(ei(I[1], v[1]) ∩ ... ∩ HyperPlane(ei(I[1], v[1]))
```
but it is much more efficient. The code above does a polyhedral projection while this function simply replace
each halfspace `⟨a, x⟩ ≤ β` (resp. each hyperplane `⟨a, x⟩ = β`) by the halfspace `⟨a_J, x⟩ ≤ β - ⟨a_I, v⟩`
(resp. the hyperplane `⟨a_J, x⟩ = β - ⟨a_I, v⟩`) where `J = setdiff(1:N, I)`.
"""
function fixandeliminate(p::HRep{N, T}, I, v) where {N, T}
    J = setdiff(1:N, I)
    dout = FullDim{length(J)}()
    f = (i, h) -> _fixelim(h, I, J, v, dout)
    Tout = promote_type(T, eltype(v))
    similar_type(typeof(p), dout, Tout)(hmap(f, dout, Tout, p)...)
end

# TODO rewrite, it is just cutting a cone with a half-space, nothing more
# export radialprojectoncut
# function radialprojectoncut{N}(p::Polyhedron{N}, cut::Vector, at)
#   if isapproxzero(at)
#     error("at is zero")
#   end
#   if length(cut) != N
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

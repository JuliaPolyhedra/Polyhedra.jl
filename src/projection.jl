export eliminate, project
implementseliminationmethod(p::Polyhedron, ::Type{Val{:FourierMotzkin}})             = false
eliminate(p::Polyhedron, delset, ::Type{Val{:FourierMotzkin}})               = error("Fourier-Motzkin elimination not implemented for $(typeof(p))")
implementseliminationmethod(p::Polyhedron, ::Type{Val{:BlockElimination}})           = false
eliminate(p::Polyhedron, delset, ::Type{Val{:BlockElimination}})             = error("Block elimination not implemented for $(typeof(p))")

eliminate(p::Polyhedron, method::Symbol) = eliminate(p, Val{method})
eliminate(p::Polyhedron, delset, method::Symbol) = eliminate(p, delset, Val{method})

eliminate{N}(p::Polyhedron{N}, method::Type{Val{:ProjectGenerators}}) = eliminate(p, IntSet(N), method)

function eliminate{N}(p::Polyhedron{N}, delset=IntSet(N))
    fm = implementseliminationmethod(p, Val{:FourierMotzkin})
    be = implementseliminationmethod(p, Val{:BlockElimination})
    if (!fm && !be) || generatorsarecomputed(p)
        method = :ProjectGenerators
    elseif fm && (!be || (length(delset) == 1 && N in delset))
        method = :FourierMotzkin
    else
        method = :BlockElimination
    end
    eliminate(p, delset, Val{method})
end

function eliminate{N}(p::Polyhedron{N}, delset, ::Type{Val{:ProjectGenerators}})
    ext = vrep(p)
    I = eye(Int, N)
    polyhedron(I[setdiff(IntSet(1:N), collect(delset)),:] * ext, getlibrary(p))
end

function Base.convert{N, S, T}(::Type{Polyhedron{N, S}}, p::Polyhedron{N, T})
    if !hrepiscomputed(p) && vrepiscomputed(p)
        f = (i, x) -> changeeltype(typeof(x), S)(x)
        if decomposedvfast(p)
            polyhedron(PointIterator(p, f), RayIterator(p, f), getlibraryfor(p, N, S))
        else
            polyhedron(VRepIterator(p, f), getlibraryfor(p, N, S))
        end
    else
        if decomposedvfast(p)
            polyhedron(IneqIterator(p, f), EqIterator(p, f), getlibraryfor(p, N, S))
        else
            polyhedron(HRepIterator(p, f), getlibraryfor(p, N, S))
        end
    end
end

# eliminate the last dimension by default
eliminate{N,T}(p::Polyhedron{N,T})  = eliminate(p::Polyhedron, IntSet([N]))

function project{N,T}(p::Polyhedron{N,T}, P::AbstractMatrix)
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
    eliminate(p * basis, IntSet(m+1:N))
end


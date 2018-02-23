myeq(x::Real, y::Real) = myeq(promote(x, y)...)
myeq(x::T, y::T) where {T<:Real} = x == y
myeq(x::T, y::T) where {T<:AbstractFloat} = y < x+1024*eps(T) && x < y+1024*eps(T)
myeq(x::Vector{S}, y::Vector{T}) where {S<:Real,T<:Real} = myeq(promote(x, y)...)
myeq(x::Vector{T}, y::Vector{T}) where {T<:Real} = x == y
myeq(x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat} = myeq(norm(x - y), zero(T))
myeqzero(x::T) where {T<:Real} = myeq(x, zero(T))

tomatrix(M::Matrix) = M
function tomatrix(v::Vector)
    M = Matrix{eltype(v)}(length(v), 1)
    M[:,1] = v
    M
end

function myparallel(x, y)
    @assert length(x) == length(y)
    cstchosen = false
    cst = zero(Base.promote_op(/, eltype(y), eltype(x)))
    for i in 1:length(x)
        if myeqzero(x[i])
            myeqzero(y[i]) || return false
        elseif myeqzero(y[i])
            myeqzero(x[i]) || return false
        else
            c = y[i] / x[i]
            if cstchosen
                myeq(cst, c) || return false
            else
                cst = c
            end
        end
    end
    return true
end

function inaffspace(x, y, L)
    for i in 1:size(L, 1)
        z = @view L[i, :]
        # remove component
        x = x * dot(z, z) - z * dot(z, x)
        y = y * dot(z, z) - z * dot(z, y)
    end
    myparallel(x, y)
end

function inequality_fulltest(ine::MixedMatHRep, A, b, linset, aff = ine[collect(linset)])
    @test size(ine.A) == size(A)
    @test length(ine.linset) == length(linset)

    affAb = [aff.b aff.A]
    inaff(x, y) = inaffspace(x, y, affAb)

    for i in 1:size(A, 1)
        found = false
        for j in 1:size(ine.A, 1)
            # vec for julia 0.4
            if !((i in linset) ⊻ (j in ine.linset)) && inaff([b[i]; A[i,:]], [ine.b[j]; ine.A[j,:]])
                found = true
                break
            end
        end
        @test found
    end
end
function inequality_fulltest(p::Polyhedron, A, b, linset)
    A = tomatrix(A)
    detecthlinearities!(p)
    removehredundancy!(p)
    inequality_fulltest(MixedMatHRep(hrep(p)), A, b, linset, MixedMatHRep(affinehull(p)))
end

function generator_fulltest(ext::MixedMatVRep, V, R=Matrix{eltype(V)}(0, size(V, 2)), Vlinset = IntSet(), Rlinset = IntSet())
    @test size(ext.V) == size(V)
    @test size(ext.R) == size(R)
    @test length(ext.Vlinset) == length(Vlinset)
    @test length(ext.Rlinset) == length(Rlinset)
    for i in 1:size(V, 1)
        found = false
        for j in 1:size(ext.V, 1)
            if myeq(V[i, :], ext.V[j, :])
                found = true
                break
            end
        end
        @test found
    end
    linspace = ext.R[collect(ext.Rlinset),:]
    inaff(x, y) = inaffspace(x, y, linspace)
    for i in 1:size(R, 1)
        found = false
        for j in 1:size(ext.R, 1)
            if !((i in Rlinset) ⊻ (j in ext.Rlinset)) && inaff(R[i,:], ext.R[j,:])
                #if parallel(vec(R[i, :]), vec(ext.R[j, :]), (i in Rlinset) || (j in ext.Rlinset))
                found = true
                break
            end
        end
        @test found
    end
end
function generator_fulltest(p::Polyhedron, V, R=Matrix{eltype(V)}(0, size(V, 2)), Vlinset = IntSet(), Rlinset = IntSet())
    V = tomatrix(V)
    R = tomatrix(R)
    detectvlinearities!(p)
    removevredundancy!(p)
    generator_fulltest(MixedMatVRep(p), V, R, Vlinset, Rlinset)
end
#generator_fulltest(p::Polyhedron, V) = generator_fulltest(p, V, Matrix{eltype(V)}(0, size(V, 2)))

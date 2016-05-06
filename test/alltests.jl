using Polyhedra

include("simplex.jl")
include("permutahedron.jl")
include("board.jl")

myeq(x::Real, y::Real) = myeq(promote(x, y)...)
myeq{T<:Real}(x::T, y::T) = x == y
myeq{T<:AbstractFloat}(x::T, y::T) = y < x+1024*eps(T) && x < y+1024*eps(T)
myeq{S<:Real,T<:Real}(x::Vector{S}, y::Vector{T}) = myeq(promote(x, y)...)
myeq{T<:Real}(x::Vector{T}, y::Vector{T}) = x == y
myeq{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}) = myeq(norm(x - y), zero(T))
myeqzero{T<:Real}(x::T) = myeq(x, zero(T))

tomatrix(M::Matrix) = M
function tomatrix(v::Vector)
  M = Matrix{eltype(v)}(length(v), 1)
  M[:,1] = v
  M
end

function inlinspace(x, L)
  for i in 1:size(L, 1)
    y = vec(L[i,:])
    # remove component
    x = x * dot(y, y) - y * dot(y, x)
  end
  myeqzero(norm(x))
end

function inequality_fulltest(p::Polyhedron, A, b, linset)
  A = tomatrix(A)
  detecthlinearities!(p)
  removeredundantinequalities!(p)
  ine = SimpleHRepresentation(getinequalities(p))
  @test size(ine.A) == size(A)
  @test length(ine.linset) == length(linset)

  aff = SimpleHRepresentation(getinequalities(affinehull(p)))
  affAb = [aff.b aff.A]
  inaff(x) = inlinspace(x, affAb)

  for i in 1:size(A, 1)
    found = false
    for j in 1:size(ine.A, 1)
      # vec for julia 0.4
      if !((i in linset) $ (j in ine.linset)) && inaff([b[i]-ine.b[j];vec(A[i,:]-ine.A[j,:])])
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
  removeredundantgenerators!(p)
  ext = SimpleVRepresentation(getgenerators(p))
  @test size(ext.V) == size(V)
  @test size(ext.R) == size(R)
  @test length(ext.Vlinset) == length(Vlinset)
  @test length(ext.Rlinset) == length(Rlinset)
  for i in 1:size(V, 1)
    found = false
    for j in 1:size(ext.V, 1)
      if myeq(vec(V[i, :]), vec(ext.V[j, :]))
        found = true
        break
      end
    end
    @test found
  end
  linspace = ext.R[collect(ext.Rlinset),:]
  inlin(x) = inlinspace(x, linspace)
  for i in 1:size(R, 1)
    found = false
    for j in 1:size(ext.R, 1)
      if !((i in Rlinset) $ (j in ext.Rlinset)) && inlin(R[i,:]-ext.R[j,:])
      #if parallel(vec(R[i, :]), vec(ext.R[j, :]), (i in Rlinset) || (j in ext.Rlinset))
        found = true
        break
      end
    end
    @test found
  end
end
#generator_fulltest(p::Polyhedron, V) = generator_fulltest(p, V, Matrix{eltype(V)}(0, size(V, 2)))

function alltests{Lib<:PolyhedraLibrary}(lib::Lib)
  simplextest(lib)
  permutahedrontest(lib)
  boardtest(lib)
end

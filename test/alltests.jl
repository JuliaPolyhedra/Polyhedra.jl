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

tomatrix(M::Matrix) = M
function tomatrix(v::Vector)
  M = Matrix{eltype(v)}(length(v), 1)
  M[:,1] = v
  M
end

function inequality_fulltest(p::Polyhedron, A, b, linset)
  A = tomatrix(A)
  removeredundantinequalities!(p)
  ine = getinequalities(p)
  @test size(ine.A) == size(A)
  @test length(ine.linset) == length(linset)

  aff = getinequalities(affinehull(p))
  affAb = [aff.b aff.A]
  function inaff(x)
    for i in 1:size(affAb, 1)
      y = vec(affAb[i,:])
      # remove component
      x = x * dot(y, y) - y * dot(y, x)
    end
    myeq(norm(x), zero(eltype(x)))
  end

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
function generator_fulltest(p::Polyhedron, V, R)
  V = tomatrix(V)
  R = tomatrix(R)
  removeredundantgenerators!(p)
  ext = getgenerators(p)
  Polyhedra.splitvertexrays!(ext)
  @test size(ext.V) == size(V)
  @test size(ext.R) == size(R)
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
  for i in 1:size(R, 1)
    found = false
    for j in 1:size(ext.R, 1)
      if myeq(vec(R[i, :]), vec(ext.R[j, :]))
        found = true
        break
      end
    end
    @test found
  end
end
generator_fulltest(p::Polyhedron, V) = generator_fulltest(p, V, Matrix{eltype(V)}(0, size(V, 2)))

function alltests{Lib<:PolyhedraLibrary}(lib::Lib)
  simplextest(lib)
  permutahedrontest(lib)
  boardtest(lib)
end

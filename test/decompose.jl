using Test
using StaticArrays
import GeometryBasics
using Polyhedra

struct Face
    points::Vector{Vector{Float64}}
end
Face(args...) = Face([args...])

function Base.isapprox(a::Face, b::Face; atol=1e-8, rtol=1e-8)
    for i in 1:3
        ok = false
        for j in 1:3
            if !isapprox(a.points[j], b.points[((j+i-1)%3)+1]; atol=atol, rtol=rtol)
                ok = true
            end
        end
        ok && return true
    end
    return false
end

nfaces(d::Dict{<:Any, Face}) = length(d)
nfaces(d::Dict{<:Any, <:Vector}) = sum(map(length, values(d)))
function test_decompose(p, d::Dict)
    N = ndims(p)
    P = GeometryBasics.Point{N, Float64}
    NR = GeometryBasics.Point{3, Float64}
    points = GeometryBasics.decompose(P, p)
    faces = GeometryBasics.decompose(GeometryBasics.GLTriangleFace, p)
    normals = GeometryBasics.decompose(GeometryBasics.Normal{NR}(), p)
    nf = nfaces(d)
    @test length(points) == 3nf
    @test length(faces) == nf
    @test length(normals) == 3nf
    for i in eachindex(faces)
        a, b, c = faces[i]
        @test normals[3i] == normals[3i-1] == normals[3i-2]
        normal = map(x -> x === -0.0 ? 0.0 : x, round.(normals[3i], digits=6))
        face = Face(points[3i-2],
                    points[3i-1],
                    points[3i])
        @test normal in keys(d)
        f = d[normal]
        if f isa Face
            @test f ≈ face
            delete!(d, normal)
        else
            j = findfirst(fa -> fa ≈ face, f)
            @test j isa Int && j > 0
            deleteat!(f, j)
            if isempty(f)
                delete!(d, normal)
            end
        end
    end
    @test isempty(d)
end
test_decompose(p::Polyhedron, d::Dict) = test_decompose(Polyhedra.Mesh(p), d)

function _quadrantdecomposetest(lib::Polyhedra.Library, V)
    v = vrep([Ray(convert(V, [1., 0])),
              Ray(convert(V, [0, 1.]))])
    p = polyhedron(v, lib)
    d = Dict([0.0,  0.0, 1.0] => Face([0.0, 0.0],
                                      [2.0, 0.0],
                                      [0.0, 2.0]))
    test_decompose(p, d)
end
function quadrantdecomposetest(lib::Polyhedra.Library)
    _quadrantdecomposetest(lib, Vector{Float64})
    _quadrantdecomposetest(lib, SVector{2,Float64})
end
function _orthantdecomposetest(lib::Polyhedra.Library, V)
    v = vrep([Ray(convert(V, [1., 0, 0])),
              Ray(convert(V, [0, 1., 0])),
              Ray(convert(V, [0, 0, 1.]))])
    p = polyhedron(v, lib)
    d = Dict([-1.0,  0.0,  0.0] => Face([0.0, 0.0, 0.0],
                                        [0.0, 2.0, 0.0],
                                        [0.0, 0.0, 2.0]),
             [ 0.0, -1.0,  0.0] => Face([0.0, 0.0, 0.0],
                                        [0.0, 0.0, 2.0],
                                        [2.0, 0.0, 0.0]),
             [ 0.0,  0.0, -1.0] => Face([0.0, 0.0, 0.0],
                                        [2.0, 0.0, 0.0],
                                        [0.0, 2.0, 0.0]))
    test_decompose(p, d)
end
function orthantdecomposetest(lib::Polyhedra.Library)
    _orthantdecomposetest(lib, Vector{Float64})
    _orthantdecomposetest(lib, SVector{3,Float64})
end

function _squaredecomposetest(lib::Polyhedra.Library, V)
    p = polyhedron(convexhull(convert(V, [ 1,  1]),
                              convert(V, [-1,  1]),
                              convert(V, [ 1, -1]),
                              convert(V, [-1, -1])), lib)
    d = Dict([0.0,  0.0,  1.0] => [Face([-1.0,  1.0], [ 1.0,  1.0], [ 1.0, -1.0]),
                                   Face([ 1.0, -1.0], [-1.0, -1.0], [-1.0,  1.0])])
    test_decompose(p, d)
end
function squaredecomposetest(lib::Polyhedra.Library)
    _squaredecomposetest(lib, Vector{Float64})
    _squaredecomposetest(lib, SVector{2,Float64})
end

function _cubedecomposetest(lib::Polyhedra.Library, V)
    p = polyhedron(convexhull(convert(V, [ 1.,  1,  1]),
                              convert(V, [ 1., -1,  1]),
                              convert(V, [ 1.,  1, -1]),
                              convert(V, [ 1., -1, -1]),
                              convert(V, [-1.,  1,  1]),
                              convert(V, [-1., -1,  1]),
                              convert(V, [-1.,  1, -1]),
                              convert(V, [-1., -1, -1])), lib)
    # FIXME failing
    #p = polyhedron(convexhull(SymPoint([1., 1, 1]),
    #                          SymPoint([1., -1, 1]),
    #                          SymPoint([1., 1, -1]),
    #                          SymPoint([1., -1, -1])))
    d = Dict([ 0.0,  0.0, -1.0] => [Face([-1.0, -1.0, -1.0], [ 1.0, -1.0, -1.0], [ 1.0,  1.0, -1.0]),
                                    Face([ 1.0,  1.0, -1.0], [-1.0,  1.0, -1.0], [-1.0, -1.0, -1.0])],
             [-1.0,  0.0,  0.0] => [Face([-1.0, -1.0, -1.0], [-1.0,  1.0, -1.0], [-1.0,  1.0,  1.0]),
                                    Face([-1.0,  1.0,  1.0], [-1.0, -1.0,  1.0], [-1.0, -1.0, -1.0])],
             [ 0.0, -1.0,  0.0] => [Face([-1.0, -1.0,  1.0], [ 1.0, -1.0,  1.0], [ 1.0, -1.0, -1.0]),
                                    Face([ 1.0, -1.0, -1.0], [-1.0, -1.0, -1.0], [-1.0, -1.0,  1.0])],
             [ 0.0,  1.0,  0.0] => [Face([-1.0,  1.0, -1.0], [ 1.0,  1.0, -1.0], [ 1.0,  1.0,  1.0]),
                                    Face([ 1.0,  1.0,  1.0], [-1.0,  1.0,  1.0], [-1.0,  1.0, -1.0])],
             [ 0.0,  0.0,  1.0] => [Face([-1.0,  1.0,  1.0], [ 1.0,  1.0,  1.0], [ 1.0, -1.0,  1.0]),
                                    Face([ 1.0, -1.0,  1.0], [-1.0, -1.0,  1.0], [-1.0,  1.0,  1.0])],
             [ 1.0,  0.0,  0.0] => [Face([ 1.0, -1.0,  1.0], [ 1.0,  1.0,  1.0], [ 1.0,  1.0, -1.0]),
                                    Face([ 1.0,  1.0, -1.0], [ 1.0, -1.0, -1.0], [ 1.0, -1.0,  1.0])])
    test_decompose(p, d)
end
function cubedecomposetest(lib::Polyhedra.Library)
    _cubedecomposetest(lib, Vector{Float64})
    _cubedecomposetest(lib, SVector{3,Float64})
end

function largedecomposetest(lib::Polyhedra.Library)
    V = [-1 -1  1;
         -1  1  1;
          1 -1  1;
          1  1  1;
          0  1 -1;
          1  0 -1;
         -1  0 -1;
          0 -1 -1]
    R = [0 0 1]
    v = vrep(V, R)
    p = polyhedron(v, lib)
    d13 = round(1/3, digits=6)
    d23 = round(2/3, digits=6)
    md13 = round(-1/3, digits=6)
    md23 = round(-2/3, digits=6)
    d = Dict([ 1.0, 0.0, 0.0] => [Face([1.0,1.0,1.0], [1.0,-1.0,1.0], [1.0,1.0,2.0]),
                                  Face([1.0,-1.0,1.0], [1.0,-1.0,2.0], [1.0,1.0,2.0]),
                                  Face([1.0,1.0,1.0], [1.0,0.0,-1.0], [1.0,-1.0,1.0])],
             [ 0.0, 1.0, 0.0] => [Face([-1.0,1.0,1.0], [1.0,1.0,1.0], [-1.0,1.0,2.0]),
                                  Face([1.0,1.0,1.0], [1.0,1.0,2.0], [-1.0,1.0,2.0]),
                                  Face([-1.0,1.0,1.0], [0.0,1.0,-1.0], [1.0,1.0,1.0])],
             [ d23, d23,md13] => [Face([1.0,1.0,1.0], [0.0,1.0,-1.0], [1.0,0.0,-1.0])],
             [ 0.0,-1.0, 0.0] => [Face([1.0,-1.0,1.0], [-1.0,-1.0,1.0], [1.0,-1.0,2.0]),
                                  Face([-1.0,-1.0,1.0], [-1.0,-1.0,2.0], [1.0,-1.0,2.0]),
                                  Face([1.0,-1.0,1.0], [0.0,-1.0,-1.0], [-1.0,-1.0,1.0])],
             [-1.0, 0.0, 0.0] => [Face([-1.0,-1.0,1.0], [-1.0,1.0,1.0], [-1.0,-1.0,2.0]),
                                  Face([-1.0,1.0,1.0], [-1.0,1.0,2.0], [-1.0,-1.0,2.0]),
                                  Face([-1.0,-1.0,1.0], [-1.0,0.0,-1.0], [-1.0,1.0,1.0])],
             [md23,md23,md13] => [Face([-1.0,-1.0,1.0], [0.0,-1.0,-1.0], [-1.0,0.0,-1.0])],
             [md23, d23,md13] => [Face([-1.0,0.0,-1.0], [0.0,1.0,-1.0], [-1.0,1.0,1.0])],
             [ d23,md23,md13] => [Face([1.0,0.0,-1.0], [0.0,-1.0,-1.0], [1.0,-1.0,1.0])],
             [ 0.0, 0.0,-1.0] => [Face([0.0,-1.0,-1.0], [1.0,0.0,-1.0], [0.0,1.0,-1.0]),
                                  Face([0.0,1.0,-1.0], [-1.0,0.0,-1.0], [0.0,-1.0,-1.0])])
    test_decompose(p, d)
end

decomposetests = Dict("orthantdecompose" => orthantdecomposetest,
                      "squaredecompose"  => squaredecomposetest,
                      "cubedecompose"    => cubedecomposetest,
                      "largedecompose"   => largedecomposetest)

@polytestset decompose

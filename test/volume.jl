using Test
using Polyhedra

# simple algorithm for calculating area of polygons
# requires vertices to be sorted (counter)clockwise
function shoelace(verts::AbstractMatrix{<:Real})
    @assert size(verts, 1) == 2  "shoelace only works for polygons"
    xs = verts[begin,:]
    ys = verts[end,:]
    A = (ys[end]+ys[begin])*(xs[end]-xs[begin])
    for i in axes(verts,2)[begin:end-1]
        A += (ys[i]+ys[i+1])*(xs[i]-xs[i+1])
    end
    A = abs(A)
    A isa AbstractFloat ? A/2 : A//2
end
isvertsapprox(verts, points) = all(any(isapprox(p, v) for v in verts) for p in points)

# issue 285: area of square [-0.5, -0.5] x [0.5, 0.5]
function check_vol_issue285(lib)
    ineq = [
        HalfSpace([-2.0, -2.0], 4.0),
        HalfSpace([-2.0, -1.0], 2.5),
        HalfSpace([-2.0, 0.0], 2.0),
        HalfSpace([-2.0, 1.0], 2.5),
        HalfSpace([-2.0, 2.0], 4.0),
        HalfSpace([-1.0, -2.0], 2.5),
        HalfSpace([-1.0, -1.0], 1.0),
        HalfSpace([-1.0, 0.0], 0.5),
        HalfSpace([-1.0, 1.0], 1.0),
        HalfSpace([-1.0, 2.0], 2.5),
        HalfSpace([0.0, -2.0], 2.0),
        HalfSpace([0.0, -1.0], 0.5),
        HalfSpace([0.0, 0.0], 0.0),
        HalfSpace([0.0, 1.0], 0.5),
        HalfSpace([0.0, 2.0], 2.0),
        HalfSpace([1.0, -2.0], 2.5),
        HalfSpace([1.0, -1.0], 1.0),
        HalfSpace([1.0, 0.0], 0.5),
        HalfSpace([1.0, 1.0], 1.0),
        HalfSpace([1.0, 2.0], 2.5),
        HalfSpace([2.0, -2.0], 4.0),
        HalfSpace([2.0, -1.0], 2.5),
        HalfSpace([2.0, 0.0], 2.0),
        HalfSpace([2.0, 1.0], 2.5),
        HalfSpace([2.0, 2.0], 4.0),
    ]
    square = polyhedron(reduce(intersect, ineq), lib)
    sqverts = [-1 1 1 -1; -1 -1 1 1]/2
    @assert isvertsapprox(eachcol(sqverts), points(square))
    return volume(square) ≈ shoelace(sqverts)
end

# issue 249
function check_vol_issue249_1(lib)
    poly = polyhedron(HalfSpace([-1.0, 0.0], 0.0) ∩
                    HalfSpace([0.0, -1.0], 0.0) ∩
                    HalfSpace([1.0, 0.0], 1.0) ∩
                    HalfSpace([0.0, 1.0], 1.0) ∩
                    HalfSpace([-0.2, -0.8], -0.0) ∩
                    HalfSpace([0.6, 0.4], 0.6), lib)

    verts = [1/3 1 0 0; 1 0 0 1]
    @assert isvertsapprox(eachcol(verts), points(poly))
    return volume(poly) ≈ shoelace(verts)
end

function check_vol_issue249_2(lib)
    poly2 = polyhedron(# HalfSpace([-1.0, 0.0], 0.0) ∩
                    HalfSpace([0.0, -1.0], 0.0) ∩ # Comment here to get correct area
                    HalfSpace([1.0, 0.0], 1.0) ∩
                    HalfSpace([0.0, 1.0], 1.0) ∩
                    HalfSpace([-0.6, -0.4], -0.6) ∩
                    HalfSpace([0.2, 0.8], 0.95), lib)
    verts2 = [1/3 3/4 1 1; 1 1 15/16 0]
    @assert isvertsapprox(eachcol(verts2), points(poly2))
    return volume(poly2) ≈ shoelace(verts2)
end

function check_vol_issue249_3(lib)
    L3 = polyhedron(vrep([
        0 0 0 0 0 0
        0 0 1 0 0 0
        0 1 0 0 0 0
        0 1 1 0 0 1
        1 0 0 0 0 0
        1 0 1 0 1 0
        1 1 0 1 0 0
        1 1 1 1 1 1
    ]),lib)
    # reportedly solution taken from https://doi.org/10.1016/j.dam.2018.10.038
    return volume(L3) == 1//180
end
# examples similar to Cohen-Hickey paper, table 1
# https://doi.org/10.1145/322139.322141
function check_vol_simplex(n, lib)
    s = polyhedron(vrep([i==j for i in 0:n, j in 1:n]), lib)
    return volume(s) == 1//factorial(n)
end


function check_vol_isocahedron(lib)
    ϕ = (1 + √5)/2
    isoc = polyhedron(vrep([
        0  1  ϕ
        0  1 -ϕ
        0 -1  ϕ
        0 -1 -ϕ
        1  ϕ  0
        1 -ϕ  0
        -1  ϕ  0
        -1 -ϕ  0
        ϕ  0  1
        ϕ  0 -1
        -ϕ  0  1
        -ϕ  0 -1
        ]), lib)
    # https://en.wikipedia.org/wiki/Regular_icosahedron
    return volume(isoc) ≈ (5/12)*(3+√5)*2^3
end

function check_vol_dodecahedron(lib)
    ϕ = (1 + √5)/2
    ϕ² = ϕ^2
    dodec = polyhedron(vrep([
        ϕ  ϕ  ϕ
        ϕ  ϕ -ϕ
        ϕ -ϕ  ϕ
        ϕ -ϕ -ϕ
        -ϕ  ϕ  ϕ
        -ϕ  ϕ -ϕ
        -ϕ -ϕ  ϕ
        -ϕ -ϕ -ϕ
        0  1  ϕ²
        0  1 -ϕ²
        0 -1  ϕ²
        0 -1 -ϕ²
        1  ϕ² 0
        1 -ϕ² 0
        -1  ϕ² 0
        -1 -ϕ² 0
        ϕ² 0  1
        ϕ² 0 -1
        -ϕ² 0  1
        -ϕ² 0 -1
    ]), lib)
    # https://en.wikipedia.org/wiki/Regular_dodecahedron
    return volume(dodec) ≈ (1/4)*(15+7*√5)*2^3
end

@testset "volumes" begin
    lib_float = DefaultLibrary{Float64}()
    lib_exact = DefaultLibrary{Int}()
    @testset "simplex" for n in 1:5
        @test check_vol_simplex(n, lib_exact)
    end
    @test check_vol_issue285(lib_float)
    @test check_vol_issue249_1(lib_float)
    @test check_vol_issue249_2(lib_float)
    @test check_vol_issue249_3(lib_exact)
    @test check_vol_isocahedron(lib_float)
    @test check_vol_dodecahedron(lib_float)
end

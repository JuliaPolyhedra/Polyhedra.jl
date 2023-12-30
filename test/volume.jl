using Test
using Polyhedra
using CDDLib

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
isvertsapprox = (verts, points) -> all(any(isapprox(p, v) for v in verts) for p in points)

# issue 285: area of square [-0.5, -0.5] x [0.5, 0.5]

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
square = polyhedron(reduce(intersect, ineq))
sqverts = [-1 1 1 -1; -1 -1 1 1]/2
@test isvertsapprox(eachcol(sqverts), points(square))
@test volume(square) ≈ shoelace(sqverts)

# issue 249

poly = polyhedron(HalfSpace([-1.0, 0.0], 0.0) ∩
                  HalfSpace([0.0, -1.0], 0.0) ∩
                  HalfSpace([1.0, 0.0], 1.0) ∩
                  HalfSpace([0.0, 1.0], 1.0) ∩
                  HalfSpace([-0.2, -0.8], -0.0) ∩
                  HalfSpace([0.6, 0.4], 0.6))

verts = [1/3 1 0 0; 1 0 0 1]
@test isvertsapprox(eachcol(verts), points(poly))
@test volume(poly) ≈ shoelace(verts)

poly2 = polyhedron(# HalfSpace([-1.0, 0.0], 0.0) ∩
                    HalfSpace([0.0, -1.0], 0.0) ∩ # Comment here to get correct area
                    HalfSpace([1.0, 0.0], 1.0) ∩
                    HalfSpace([0.0, 1.0], 1.0) ∩
                    HalfSpace([-0.6, -0.4], -0.6) ∩
                    HalfSpace([0.2, 0.8], 0.95))
verts2 = [1/3 3/4 1 1; 1 1 15/16 0]
@test isvertsapprox(eachcol(verts2), points(poly2))
@test volume(poly2) ≈ shoelace(verts2)

lib=CDDLib.Library(:exact)
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
@test volume(L3) == 1//180

# examples similar to Cohen-Hickey paper, table 1
# https://doi.org/10.1145/322139.322141
s5 = polyhedron(vrep([
    0 0 0 0 0
    0 0 0 0 1
    0 0 0 1 0
    0 0 1 0 0
    0 1 0 0 0
    1 0 0 0 0
    ]), lib)
@test volume(s5) == 1//120

flib = CDDLib.Library()

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
    ]), flib)
# https://en.wikipedia.org/wiki/Regular_icosahedron
@test volume(isoc) ≈ (5/12)*(3+√5)*2^3

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
    ]), flib)
# https://en.wikipedia.org/wiki/Regular_dodecahedron
@test volume(dodec) ≈ (1/4)*(15+7*√5)*2^3

using SparseArrays

function sparsetest(lib::Polyhedra.Library)
    # Sparse H-representation, fixandeliminate and ⊆
    # 3x_1 + 2x_2        ≤  2
    # 2x_1        + 3x_3 ≥  1
    #  x_1               ≥  0
    #         x_2        ≥ -1
    #                x_3 ≤  1
    h = HalfSpace(sparsevec([1, 2], [3, 2], 3), 2) ∩
        HalfSpace(sparsevec([1, 3], [-2, -3]), -1) ∩
        HalfSpace(sparsevec([1], [-1], 3),      0) ∩
        HalfSpace(sparsevec([2], [-1], 3),      1) ∩
        HalfSpace(sparsevec([3], [1]),          1)
    v = convexhull([0,     1,  1//3],
                   [4//3, -1, -5//9],
                   [0,    -1,  1//3],
                   [0,     1,  1],
                   [4//3, -1,  1],
                   [0,    -1,  1])
    p = polyhedron(h, lib)
    @test !vrepiscomputed(p)
    @test p ⊆ HalfSpace([3, 2, 0], 2)
    @test !vrepiscomputed(p)
    @test !(p ⊆ HalfSpace([3, 2, 0], 1))
    @test !vrepiscomputed(p)
    inequality_fulltest(p, h)
    generator_fulltest(p, v)
    q = polyhedron(h, lib)
    @test !vrepiscomputed(q)
    q = fixandeliminate(q, [1], [1])
    # If lib supported sparse, hvectortype(q) should be sparse hence hvectortype(p) should also be sparse
    @test Polyhedra.hvectortype(typeof(p)) == Polyhedra.hvectortype(typeof(q))
    @test !vrepiscomputed(q)
    @test q ⊆ HalfSpace([2, -3], 0) # This also checks that the solver was ot dropped during fixandeliminate
    @test !vrepiscomputed(q)
    @test !(q ⊆ HalfSpace([2, -3], -1))
    @test !vrepiscomputed(q)
    hfix = HalfSpace([ 2,  0], -1) ∩
           HalfSpace([ 0, -3],  1) ∩
           HalfSpace([-1,  0],  1) ∩
           HalfSpace([ 0,  1],  1)
    vfix = convexhull([-1//2, -1//3],
                      [-1,    -1//3],
                      [-1//2,  1],
                      [-1,     1])
    inequality_fulltest(q, hfix)
    generator_fulltest(q, vfix)
end

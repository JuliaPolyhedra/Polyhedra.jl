using SparseArrays

function sparserecttest(lib::Polyhedra.Library)
    # Sparse H-representation, cartesian product
    # -1 ≤ x_1 ≤ 1
    # 2 ≤ x_2 ≤ 3
    h1 = HalfSpace(sparsevec([1], [ 1]),  1) ∩
         HalfSpace(sparsevec([1], [-1]),  1)
    h2 = HalfSpace(sparsevec([1], [ 1]),  3) ∩
         HalfSpace(sparsevec([1], [-1]), -2)
    h = h1 * h2
    v = convexhull([-1, 2],
                   [ 1, 2],
                   [-1, 3],
                   [ 1, 3])
    p = polyhedron(h, lib)
    @test !vrepiscomputed(p)
    inequality_fulltest(p, h)
    generator_fulltest(p, v)
    p1 = polyhedron(h1, lib)
    @test !vrepiscomputed(p1)
    p2 = polyhedron(h2, lib)
    @test !vrepiscomputed(p2)
    q = p1 * p2
    @test !vrepiscomputed(q)
    inequality_fulltest(q, h)
    generator_fulltest(q, v)
end

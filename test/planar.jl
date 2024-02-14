using Test
using StaticArrays
using Polyhedra

function test_planar_square(VT, d)
    ei = ntuple(i -> Polyhedra.basis(VT, d, i), 2)
    so = -ei[1] + ei[2]
    no = -ei[1] - ei[2]
    ne = ei[1] - ei[2]
    se = ei[1] + ei[2]
    expected = [no, so, se, ne]
    for ps in [
        [no, so, ne, se],
        [no, se, so, ne],
        [se, so, no, ne],
        expected
    ]
        hull = Polyhedra.planar_hull(vrep(ps))
        @test collect(points(hull)) == expected
        @test !hasallrays(hull)
        v = vrep(ps) + conichull(se)
        hull = Polyhedra.planar_hull(v)
        @test collect(points(hull)) == [ne, no, so]
        @test collect(rays(hull)) == [Ray(se)]
        @test !haslines(hull)
        v = vrep(ps) + conichull(se, se)
        hull = Polyhedra.planar_hull(v)
        @test collect(points(hull)) == [ne, no, so]
        @test collect(rays(hull)) == [Ray(se)]
        @test !haslines(hull)
        ne = ei[1] - ei[2]
        for v in [vrep(ps) + conichull(se, ne),
                  vrep(ps) + conichull(se, ne, ei[1])]
            hull = Polyhedra.planar_hull(v)
            @test collect(points(hull)) == [no, so]
            @test collect(rays(hull)) == [Ray(ne), Ray(se)]
            @test !haslines(hull)
        end
        v = vrep(ps) + conichull(se, ne, -se)
        hull = Polyhedra.planar_hull(v)
        @test collect(points(hull)) == [ne]
        @test collect(rays(hull)) == [Ray(ne)]
        @test collect(lines(hull)) == [Line(se)]
        for v in [vrep(ps) + conichull(se, ne, -ne),
                  vrep(ps) + conichull(se, Line(ne)),
                  vrep(ps) + conichull(ne, se, -ne)]
            hull = Polyhedra.planar_hull(v)
            @test collect(points(hull)) == [no]
            @test collect(rays(hull)) == [Ray(se)]
            @test collect(lines(hull)) == [Line(ne)]
        end
        for v in [vrep(ps) + conichull(se, ne, -ne, -se),
                  vrep(ps) + conichull(se, ne, -ei[1])]
            hull = Polyhedra.planar_hull(v)
            @test collect(points(hull)) == [no]
            @test !hasrays(hull)
            @test collect(lines(hull)) == [Line(ne), Line(se)]
        end
    end
end

# Square + redundant points + [+-1.5, -+0.25]
function test_planar_square_with_more()
    expected = [[-1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [1.0, -1.0]]
    v0 = convexhull([-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [0.0, 1.0])
    v1 = Polyhedra.planar_hull(v0)
    @test collect(points(v1)) == expected
    v0 = convexhull([-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [0.0, 1.0 + 1e-10])
    v1 = Polyhedra.planar_hull(v0)
    @test collect(points(v1)) == expected
    v0 = convexhull([0.0, 1.0 + 1e-10], [-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0])
    v1 = Polyhedra.planar_hull(v0)
    @test collect(points(v1)) == expected
    expected = [[1.0, -1.0], [-1.0, -1.0], [-1.0, 1.0], [0.0, 1.0 + 1e-4], [1.0, 1.0]]
    v0 = convexhull([1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [0.0, 1.0 + 1e-4], [-1.0, -1.0])
    v1 = Polyhedra.planar_hull(v0)
    @test collect(points(v1)) == expected
    expected = [[-1.0, -1.0], [-1.5, 0.25], [-1.0, 1.0], [1.0, 1.0], [1.5, -0.25], [1.0, -1.0]]
    v0 = convexhull([-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [0.0, 1.0], [-1.5, 0.25], [1.5, -0.25])
    v1 = Polyhedra.planar_hull(v0)
    @test collect(points(v1)) == expected
    v0 = convexhull([-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [0.0, 1.0], [-0.5, -1.0], [-1.5, 0.25], [0.0, -0.75], [0.0, 0.75], [1.5, -0.25])
    v1 = Polyhedra.planar_hull(v0)
    @test collect(points(v1)) == expected
    v0 = convexhull([-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [0.0, -1.0], [0.0, 1.0], [-0.5, -1.0], [-1.5, 0.25], [0.0, -0.75], [0.0, 0.75], [1.5, -0.25])
    v1 = Polyhedra.planar_hull(v0)
    @test collect(points(v1)) == expected
    v0 = convexhull([-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [0.5, 1.0], [0.0, -1.0], [0.0, 1.0], [-0.5, -1.0], [-1.5, 0.25], [0.0, -0.75], [0.0, 0.75], [1.5, -0.25])
    v1 = Polyhedra.planar_hull(v0)
    @test collect(points(v1)) == expected
end

function test_issue_271()
    expected = [[-0.5, 0.0], [0.5, 0.5], [1.0, 0.0]]
    v = convexhull([0.0, 0.0], [1, 0], [0.5, 0.0], [-0.5, 0.0], [0.5, 0.5])
    p = Polyhedra.planar_hull(v)
    @test collect(points(p)) == expected
    v = convexhull([0.0, 0.0], [0.5, 0.0], [-0.5, 0.0], [0.5, 0.5], [1, 0])
    p = Polyhedra.planar_hull(v)
    @test collect(points(p)) == expected

    expected = [[1.0, 0.0], [-0.5, 0.0], [0.0, 0.5], [1.0, 0.5]]
    v = vrep([[1, 0], [0.0, 0.0], [0.5, 0.0], [-0.5, 0.0], [0.5, 0.5], [1.0, 0.5], [0.0, 0.5], [0.2, 0.5]])
    p = Polyhedra.planar_hull(v)
    @test collect(points(p)) == expected
    v = vrep([[1, 0], [0.0, 0.0], [0.5, 0.0], [-0.5, 0.0], [1.0, 0.5], [0.2, 0.5], [0.5, 0.5], [0.0, 0.5]])
    p = Polyhedra.planar_hull(v)
    @test collect(points(p)) == expected
end

function test_issue_333()
    vs = [
        0.5 -0.4
        -0.5  1
        0    1
        0.8  1
    ]

    h = Polyhedra.planar_hull(vrep(vs))
    @test all(points(h) .== [vs[i, :] for i in [1, 2, 4]])

    vs = [-0.96     -0.4;
        -0.5      -0.4;
        0.0      -0.4;
        0.5      -0.4;
        1.0      -0.36065655;
        -0.995    -0.05;
        1.0      -0.05;
        -1.0       0.3;
        0.97      0.3;
        -1.0       0.65;
        0.930928  0.65;
        -1.0       1.0;
        -0.5       1.0;
        0.0       1.0;
        0.5       1.0;
        0.865979 1.0;]

    h = Polyhedra.planar_hull(vrep(vs))
    @test all(points(h) .== [vs[i, :] for i in [1, 6, 8, 12, 16, 11, 9, 7, 5, 4]])
end

@testset "Planar" begin
    for T in [Float64, Int]
        for (VT, d) in [(Vector{T}, 2), (SVector{2, T}, StaticArrays.Size(2))]
            test_planar_square(VT, d)
        end
    end
    test_planar_square_with_more()
    test_issue_271()
    test_issue_333()
end

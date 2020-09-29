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

@testset "Planar" begin
    for T in [Float64, Int]
        for (VT, d) in [(Vector{T}, 2), (SVector{2, T}, StaticArrays.Size(2))]
            test_planar_square(VT, d)
        end
    end
end

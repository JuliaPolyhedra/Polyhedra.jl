@testset "Duplicate removal" begin
    h = removeduplicates(hrep([1 1; -1 -1; 1 0; 0 -1], [1, -1, 1, 0]))
    @test nhyperplanes(h) == 1
    @test first(hyperplanes(h)) == HyperPlane([1, 1], 1)
    @test first(hyperplanes(h)) == HyperPlane([-1, -1], -1)
    @test nhalfspaces(h) == 1
    h = removeduplicates(hrep([1 1; -1 -1], [2, -1]))
    @test !hashyperplanes(h)
    @test nhalfspaces(h) == 2
    v = removeduplicates(convexhull([1, 2], SymPoint([1, 2]), [-1, 1]))
    @test typeof(v) == Polyhedra.PointsHull{2,Int,Vector{Int}}
    @test collect(sympoints(v)) == [SymPoint([1, 2])]
    @test collect(points(v)) == [[-1, 1]]
    @test !hasallrays(v)
    for v in (removeduplicates(convexhull([1, 2], [2, 1], SymPoint([1, 2]), [-1, 1]) + Line([0, 1])),
              removeduplicates(Line([0, 1]) + convexhull([1, 2], [2, 1], SymPoint([1, 2]), [-1, 1])),
              removeduplicates(Line([0, 1]) + convexhull([-1, -2], [2, 1], [1, 2], [-1, 1], [-1, -2])))
        @test typeof(v) == Polyhedra.Hull{2,Int,Vector{Int}}
        @test collect(sympoints(v)) == [SymPoint([1, 2])]
        @test collect(points(v)) == [[2, 1]]
        @test !hasrays(v)
        @test collect(lines(v)) == [Line([0, 1])]
    end
end

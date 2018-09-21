@testset "Redundancy removal" begin
    @testset "isredundant" begin
        @test isredundant(HyperPlane([0, 1], 0) ∩ HalfSpace([1, 1], 1), [0, 0], nl=0)
        @test !isredundant(HalfSpace([1, 1], 1) ∩ HyperPlane([0, 1], 0) ∩ HalfSpace([-1, -1], 1), [1, 0], nl=0)
        #@test !ishredundant(convexhull(SymPoint([1, 0])), HyperPlane([0, 1], 0), d=1)
        #@test isredundant(convexhull(SymPoint([1, 0])), HalfSpace([0, 1], 1), d=1)
    end
#    @testset "Non-redundant SymPoint" begin
#        vr = convexhull(SymPoint([-1, 0]), SymPoint([0, 1]))
#        hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1) ∩ HalfSpace([1, 1], 1) ∩ HalfSpace([-1, 1], 1)
#        vr = @inferred Polyhedra.removevredundancy(vr, hr)
#        @test vr isa Polyhedra.PointsHull{2,Int,Vector{Int}}
#        @test collect(sympoints(vr)) == [SymPoint([-1, 0]), SymPoint([0, 1])]
#        @test !haspoints(vr)
#        @test !hasallrays(vr)
#    end
    @testset "Point" begin
        vr = convexhull([-1, 0], [0, -1]) + conichull([1, 1], [-1, 1])
        hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1)
        vr = Polyhedra.removevredundancy(vr, hr)
        @test collect(points(vr)) == [[0, -1]]
        @test !haslines(vr)
        @test collect(rays(vr)) == [Ray([1, 1]), Ray([-1, 1])]
    end
#    @testset "Split SymPoint" begin
#        vr = convexhull(SymPoint([-1, 0]), SymPoint([0, 1])) + conichull([1, 1], [-1, 1])
#        hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1)
#        vr = Polyhedra.removevredundancy(vr, hr)
#        @test !hassympoints(vr)
#        @test collect(points(vr)) == [[0, -1]]
#        @test !haslines(vr)
#        @test collect(rays(vr)) == [Ray([1, 1]), Ray([-1, 1])]
#
#        vr = convexhull(SymPoint([1, 0]), SymPoint([0, 1])) + conichull([-1, 1])
#        hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1) ∩ HalfSpace([1, 1], 1)
#        vr = Polyhedra.removevredundancy(vr, hr)
#        @test !hassympoints(vr)
#        @test collect(points(vr)) == [[1, 0], [0, -1]]
#        @test !haslines(vr)
#        @test collect(rays(vr)) == [Ray([-1, 1])]
#    end
end
@testset "Duplicate removal" begin
    h = removeduplicates(hrep([1 1; -1 -1; 1 0; 0 -1], [1, -1, 1, 0]))
    @test nhyperplanes(h) == 1
    @test first(hyperplanes(h)) == HyperPlane([1, 1], 1)
    @test first(hyperplanes(h)) == HyperPlane([-1, -1], -1)
    @test nhalfspaces(h) == 1
    h = removeduplicates(hrep([1 1; -1 -1], [2, -1]))
    @test !hashyperplanes(h)
    @test nhalfspaces(h) == 2
#    v = removeduplicates(convexhull([1, 2], SymPoint([1, 2]), [-1, 1]))
#    @test typeof(v) == Polyhedra.PointsHull{2,Int,Vector{Int}}
#    @test collect(sympoints(v)) == [SymPoint([1, 2])]
#    @test collect(points(v)) == [[-1, 1]]
#    @test !hasallrays(v)
    for v in (#removeduplicates(convexhull([1, 2], [2, 1], SymPoint([1, 2]), [-1, 1]) + Line([0, 1])),
              #removeduplicates(Line([0, 1]) + convexhull([1, 2], [2, 1], SymPoint([1, 2]), [-1, 1])),
              removeduplicates(Line([0, 1]) + convexhull([-1, -2], [2, 1], [1, 2], [-1, 1], [-1, -2])),)
        @test typeof(v) == Polyhedra.Hull{Int, Vector{Int}, Int}
        #@test collect(sympoints(v)) == [SymPoint([1, 2])]
        @test collect(points(v)) == [[-1, -2], [2, 1], [1, 2]]
        @test !hasrays(v)
        @test collect(lines(v)) == [Line([0, 1])]
    end
end

function convexhull_test()
    v = convexhull([1, 0, 0]) + conichull(Ray([0, 1, 0]), Line([0, 0, 1]))
    @test npoints(v) == 1
    @test nlines(v) == 1
    @test nrays(v) == 1
    ps = v.points
    _ps = ps.points
    rs = v.rays
    _ls = rs.lines
    __ls = _ls.lines
    _rs = rs.rays
    function _inplace()
        @test v.points === ps
        @test v.rays === rs
        @test v.points.points === _ps
        @test v.rays.lines === _ls
        @test v.rays.lines.lines === __ls
        @test v.rays.rays === _rs
    end
    convexhull!(v, [2, 0, 0])
    _inplace()
    @test npoints(v) == 2
    @test nlines(v) == 1
    @test nrays(v) == 1
    convexhull!(v, Ray([0, 2, 0]))
    _inplace()
    @test npoints(v) == 2
    @test nlines(v) == 1
    @test nrays(v) == 2
    convexhull!(v, Line([1, 0, 0]))
    _inplace()
    @test npoints(v) == 2
    @test nlines(v) == 2
    @test nrays(v) == 2
end

function intersect_test()
    h = HalfSpace([1, 0], 1) ∩ HyperPlane([0, 1], 1)
    @test nhyperplanes(h) == 1
    @test nhalfspaces(h) == 1
    hs = h.halfspaces
    hp = h.hyperplanes
    _hp = h.hyperplanes.hyperplanes
    intersect!(h, HalfSpace([1, 1], 2))
    @test h.halfspaces === hs
    @test h.hyperplanes === hp
    @test h.hyperplanes.hyperplanes === _hp
    @test nhyperplanes(h) == 1
    @test nhalfspaces(h) == 2
    intersect!(h, HyperPlane([1, 1], 3))
    @test h.halfspaces === hs
    @test h.hyperplanes === hp
    @test h.hyperplanes.hyperplanes === _hp
    @test nhyperplanes(h) == 2
    @test nhalfspaces(h) == 2
end

@testset "VecRep" begin
    @testset "Convex hull" begin
        convexhull_test()
    end
    @testset "Intersection" begin
        intersect_test()
    end
end

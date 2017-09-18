@testset "Element" begin
    @testset "RemProj" begin
        r = Ray([1//1, -5//1, 0//1])
        l = Line([-15//1, 75//1, -26//1])
        r = Polyhedra.remproj(r, l)
        @test r isa Ray
        @test coord(r) == [1//1, -5//1, -15//1]
    end
    @testset "Simplify" begin
        r = @inferred simplify(Ray([0, 0]))
        @test r isa Ray{2, Int}
        @test coord(r) == [0, 0]
        p = @inferred simplify([6, 4])
        @test p isa Vector{Int}
        @test p == [6, 4]
        l = @inferred simplify(Line([6, 4]))
        @test l isa Line{2, Int}
        @test coord(l) == [3, 2]
        hp = @inferred simplify(HyperPlane([6, 4], 2))
        @test hp isa HyperPlane{2, Int}
        @test hp.a == [3, 2]
        @test hp.β == 1
        hs = @inferred simplify(HalfSpace([0, 0], 0))
        @test hs isa HalfSpace{2, Int}
        @test hs.a == [0, 0]
        @test hs.β == 0
        hpr = @inferred simplify(HyperPlane([0//1, 0//1], 0//1))
        @test hpr isa HyperPlane{2, Rational{Int}}
        @test hpr.a == [0, 0]
        @test hpr.β == 0
    end
end

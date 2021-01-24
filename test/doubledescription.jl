using Test
using StaticArrays

function empty_test(h::HRepresentation{T}) where T
    #v = @inferred doubledescription(h)
    v = doubledescription(h)
    U = T == Float64 ? Float64 : Rational{BigInt}
    v = doubledescription(h)
    @test v isa Polyhedra.Hull{U,Vector{U}}
    @test !haspoints(v)
    @test !hasallrays(v)
end
empty_hyperplane_test(z) = empty_test(intersect(HyperPlane([z, z], 1)))
empty_halfspace_test(z) = empty_test(intersect(HalfSpace([z, z], -1)))
empty_space_test(z) = empty_test(HyperPlane([z, 1], 1) ∩ HyperPlane([z, 1], -1))
empty_range_test(z) = empty_test(HalfSpace([1, z], -1) ∩ HalfSpace([-1, z], -1))

@testset "Double Description" begin
    @test Polyhedra.polytypefor(Float32) == Float64
    @testset "H-representation -> V-representation" begin
        @testset "Intersection" begin
            @testset "Vector" begin
                @testset "Exact" begin
                    h = hrep([HalfSpace([ 1,  1], 1),
                              HalfSpace([ 1, -1], 0),
                              HalfSpace([-1,  0], 0)])
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Rational{BigInt},Vector{Rational{BigInt}}}
                    @test collect(points(v)) == [[0, 0], [0, 1], [1//2, 1//2]]
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = hrep([HalfSpace([ 1.,  1], 1),
                              HalfSpace([ 1., -1], 0),
                              HalfSpace([-1.,  0], 0)])
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Float64,Vector{Float64}}
                    @test collect(points(v)) == [[0.0, 0.0], [0.0, 1.0], [1/2, 1/2]]
                    @test !hasallrays(v)
                end
            end
            @testset "SVector" begin
                @testset "Exact" begin
                    h = hrep([HalfSpace((@SVector [ 1,  1]), 1),
                              HalfSpace((@SVector [ 1, -1]), 0),
                              HalfSpace((@SVector [-1,  0]), 0)])
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Rational{BigInt},SVector{2,Rational{BigInt}}}
                    @test collect(points(v)) == [(@SVector [0, 0]), (@SVector [0, 1]), (@SVector [1//2, 1//2])]
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = hrep([HalfSpace((@SVector [ 1.,  1]), 1),
                              HalfSpace((@SVector [ 1., -1]), 0),
                              HalfSpace((@SVector [-1.,  0]), 0)])
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Float64,SVector{2,Float64}}
                    @test collect(points(v)) == [(@SVector [0.0, 0.0]), (@SVector [0.0, 1.0]), (@SVector [1/2, 1/2])]
                    @test !hasallrays(v)
                end
            end
        end
        @testset "MixedMatHRep" begin
            @testset "Exact" begin
                h = hrep([ 1  1
                           1 -1
                          -1  0],
                         [1, 0, 0])
                #v = @inferred doubledescription(h)
                v = doubledescription(h)
                @test v isa Polyhedra.MixedMatVRep{Rational{BigInt}}
                @test v.V == [0//1 0//1; 0//1 1//1; 1//2 1//2]
                @test v.R == zeros(Rational{BigInt}, 0, 2)
                @test isempty(v.Rlinset)
            end
            @testset "Numerical" begin
                h = hrep([ 1  1
                           1 -1
                          -1  0],
                         [1., 0, 0])
                #v = @inferred doubledescription(h)
                v = doubledescription(h)
                @test v isa Polyhedra.MixedMatVRep{Float64}
                @test v.V == [0 0; 0 1; 1/2 1/2]
                @test v.R == zeros(0, 2)
                @test isempty(v.Rlinset)
            end
        end
        _str(::Float64) = "Numerical"
        _str(::Int) = "Exact"
        @testset "Empty Intersection $(_str(z))" for z in [0, 0.0]
            @testset "0x_1 + 0x_2 = 1" begin
                empty_hyperplane_test(z)
            end
            @testset "0x_1 + 0x_2 ≤ -1" begin
                empty_halfspace_test(z)
            end
            @testset "-1 = 0x_1 + x_2 = 1" begin
                empty_space_test(z)
            end
            @testset "1 ≤ x_1 + 0x_2 ≤ -1" begin
                empty_range_test(z)
            end
        end
        @testset "Simple hyperplane $(_str(o))" for o in [1, 1.0]
            h = intersect(HyperPlane([o, o], o))
            v = doubledescription(h)
            @test nlines(v) == 1
            @test first(lines(v)) == Line([-1, 1])
            @test !hasrays(v) == 1
            @test npoints(v) == 1
            @test first(points(v)) == [1, 0]
        end
        @testset "Quadrilateral $(_str(z))" for z in [0, 0.0]
            h = HalfSpace([2, 1], 4) ∩ HalfSpace([1, 2], 4) ∩
                HalfSpace([-1, z], z) ∩ HalfSpace([z, -1], z)
            v = doubledescription(h)
            @test !haslines(v)
            @test !hasrays(v)
            @test npoints(v) == 4
            @test all(points(v) .≈ [[0, 0], [2, 0], [0, 2], [4//3, 4//3]])
        end
    end
    @testset "V-representation -> H-representation" begin
        @testset "Hull" begin
            @testset "Vector" begin
                @testset "Exact" begin
                    v = conichull([1, 0],
                                  [0, 1])
                    h = doubledescription(v)
                    #h = @inferred doubledescription(v)
                    @test !hashyperplanes(h)
                    @test collect(halfspaces(h)) == [HalfSpace([0, 0], 1), HalfSpace([-1, 0], 0), HalfSpace([0, -1], 0)] # FIXME get rid of (0, 0) 1
                end
                @testset "Numerical" begin
                    v = conichull([1., 0.],
                                  [0., 1.])
                    h = doubledescription(v)
                    #h = @inferred doubledescription(v)
                    @test !hashyperplanes(h)
                    @test collect(halfspaces(h)) == [HalfSpace([0, 0], 1), HalfSpace([-1, 0], 0), HalfSpace([0, -1], 0)] # FIXME get rid of (0, 0) 1
                end
            end
        end
    end
end

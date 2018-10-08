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
                    @test collect(points(v)) == [[1//2, 1//2], [0, 0], [0, 1]]
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = hrep([HalfSpace([ 1.,  1], 1),
                              HalfSpace([ 1., -1], 0),
                              HalfSpace([-1.,  0], 0)])
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Float64,Vector{Float64}}
                    @test collect(points(v)) == [[1/2, 1/2], [0.0, 0.0], [0.0, 1.0]]
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
                    @test collect(points(v)) == [(@SVector [1//2, 1//2]), (@SVector [0, 0]), (@SVector [0, 1])]
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = hrep([HalfSpace((@SVector [ 1.,  1]), 1),
                              HalfSpace((@SVector [ 1., -1]), 0),
                              HalfSpace((@SVector [-1.,  0]), 0)])
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Float64,SVector{2,Float64}}
                    @test collect(points(v)) == [(@SVector [1/2, 1/2]), (@SVector [0.0, 0.0]), (@SVector [0.0, 1.0])]
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
                @test v.V == [1//2 1//2; 0//1 0//1; 0//1 1//1]
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
                @test v.V == [1/2 1/2; 0 0; 0 1]
                @test v.R == zeros(0, 2)
                @test isempty(v.Rlinset)
            end
        end
        @testset "Empty Intersection" begin
            h = HyperPlane([0, 1], 1) ∩ HyperPlane([0, 1], -1)
            @testset "0x_1 + 0x_2 = 1" begin
                @testset "Exact" begin
                    h = intersect(HyperPlane([0, 0], 1))
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Rational{BigInt},Vector{Rational{BigInt}}}
                    @test !haspoints(v)
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = intersect(HyperPlane([0., 0], 1))
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Float64,Vector{Float64}}
                    @test !haspoints(v)
                    @test !hasallrays(v)
                end
            end
            @testset "-1 = 0x_1 + x_2 = 1" begin
                @testset "Exact" begin
                    h = HyperPlane([0, 1], 1) ∩ HyperPlane([0, 1], -1)
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Rational{BigInt},Vector{Rational{BigInt}}}
                    @test !haspoints(v)
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = HyperPlane([0, 1.], 1) ∩ HyperPlane([0, 1.], -1)
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Float64,Vector{Float64}}
                    @test !haspoints(v)
                    @test !hasallrays(v)
                end
            end
            @testset "1 ≤ x_1 + 0x_2 ≤ -1" begin
                @testset "Exact" begin
                    h = HalfSpace([1, 0], -1) ∩ HalfSpace([-1, 0], -1)
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Rational{BigInt},Vector{Rational{BigInt}}}
                    @test !haspoints(v)
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = HalfSpace([1, 0.], -1) ∩ HalfSpace([-1, 0.], -1)
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Float64,Vector{Float64}}
                    @test !haspoints(v)
                    @test !hasallrays(v)
                end
            end
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
                    @test collect(halfspaces(h)) == [HalfSpace([0, -1], 0), HalfSpace([-1, 0], 0), HalfSpace([0, 0], 1)] # FIXME get rid of (0, 0) 1
                end
                @testset "Numerical" begin
                    v = conichull([1., 0.],
                                  [0., 1.])
                    h = doubledescription(v)
                    #h = @inferred doubledescription(v)
                    @test !hashyperplanes(h)
                    @test collect(halfspaces(h)) == [HalfSpace([0, -1], 0), HalfSpace([-1, 0], 0), HalfSpace([0, 0], 1)] # FIXME get rid of (0, 0) 1
                end
            end
        end
    end
end

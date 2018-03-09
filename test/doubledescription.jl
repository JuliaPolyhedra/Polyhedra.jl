@testset "Double Description" begin
    @testset "H-representation -> V-representation" begin
        @testset "Intersection" begin
            @testset "Vector" begin
                @testset "Exact" begin
                    h = hrep([HalfSpace([ 1,  1], 1),
                              HalfSpace([ 1, -1], 0),
                              HalfSpace([-1,  0], 0)])
                    v = @inferred doubledescription(h)
                    @test v isa Polyhedra.Hull{2,Rational{BigInt},Vector{Rational{BigInt}}}
                    @test collect(points(v)) == [[1//2, 1//2], [0, 0], [0, 1]]
                    @test !hassympoints(v)
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = hrep([HalfSpace([ 1.,  1], 1),
                              HalfSpace([ 1., -1], 0),
                              HalfSpace([-1.,  0], 0)])
                    v = @inferred doubledescription(h)
                    @test v isa Polyhedra.Hull{2,Float64,Vector{Float64}}
                    @test collect(points(v)) == [[1/2, 1/2], [0.0, 0.0], [0.0, 1.0]]
                    @test !hassympoints(v)
                    @test !hasallrays(v)
                end
            end
            @testset "SVector" begin
                @testset "Exact" begin
                    h = hrep([HalfSpace((@SVector [ 1,  1]), 1),
                              HalfSpace((@SVector [ 1, -1]), 0),
                              HalfSpace((@SVector [-1,  0]), 0)])
                    v = @inferred doubledescription(h)
                    @test v isa Polyhedra.Hull{2,Rational{BigInt},SVector{2,Rational{BigInt}}}
                    @test collect(points(v)) == [(@SVector [1//2, 1//2]), (@SVector [0, 0]), (@SVector [0, 1])]
                    @test !hassympoints(v)
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = hrep([HalfSpace((@SVector [ 1.,  1]), 1),
                              HalfSpace((@SVector [ 1., -1]), 0),
                              HalfSpace((@SVector [-1.,  0]), 0)])
                    v = @inferred doubledescription(h)
                    @test v isa Polyhedra.Hull{2,Float64,SVector{2,Float64}}
                    @test collect(points(v)) == [(@SVector [1/2, 1/2]), (@SVector [0.0, 0.0]), (@SVector [0.0, 1.0])]
                    @test !hassympoints(v)
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
                v = @inferred doubledescription(h)
                @test v isa Polyhedra.MixedMatVRep{2,Rational{BigInt}}
                @test v.V == [1//2 1//2; 0//1 0//1; 0//1 1//1]
                @test isempty(v.Vlinset)
                @test v.R == zeros(Rational{BigInt}, 0, 2)
                @test isempty(v.Rlinset)
            end
            @testset "Numerical" begin
                h = hrep([ 1  1
                           1 -1
                          -1  0],
                         [1., 0, 0])
                v = @inferred doubledescription(h)
                @test v isa Polyhedra.MixedMatVRep{2,Float64}
                @test v.V == [1/2 1/2; 0 0; 0 1]
                @test isempty(v.Vlinset)
                @test v.R == zeros(0, 2)
                @test isempty(v.Rlinset)
            end
        end
    end
    @testset "V-representation -> H-representation" begin
    end
end

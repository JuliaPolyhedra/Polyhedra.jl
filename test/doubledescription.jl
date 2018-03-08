@testset "Double Description" begin
    @testset "Intersection" begin
        @testset "Vector" begin
            h = hrep([HalfSpace([ 1.,  1], 1),
                      HalfSpace([ 1., -1], 0),
                      HalfSpace([-1.,  0], 0)])
            #v = @inferred doubledescription(h) # FIXME
            v = doubledescription(h)
            @test v isa Polyhedra.Hull{2,Float64,Array{Float64,1}}
            @test collect(points(v)) == [[1/2, 1/2], [0.0, 0.0], [0.0, 1.0]]
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

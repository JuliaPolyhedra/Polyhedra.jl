@testset "LPHRepresentation" begin
    @testset "with bad arguments" begin
        @test_throws DimensionMismatch LPHRepresentation(ones(2, 2), [1], [1], [1, 2], [1, 2])
        @test_throws DimensionMismatch LPHRepresentation(ones(2, 2), [1, 2], [1, 2], [1], [1])
    end
end

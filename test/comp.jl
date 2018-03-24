@testset "Comparison" begin
    #@test SymPoint([1, 2]) ≈ SymPoint([-1, -2])
    #@test !(SymPoint([1, 2]) ≈ SymPoint([-2, -1]))
    @test Polyhedra._lt(1, 1.5)
    @test !Polyhedra._lt(2, 1.5)
end

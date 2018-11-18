@testset "Lifted Representation" begin
    change_fulldim_test(LiftedHRepresentation)
    change_fulldim_test(LiftedVRepresentation)

    isempty_vtest(LiftedVRepresentation(zeros(0, 2)), 0, 0)
    isempty_vtest(LiftedVRepresentation([1 0; 1 1]), 0, 2)
    isempty_vtest(LiftedVRepresentation([0 1; 0 2]), 2, 0)
    isempty_vtest(LiftedVRepresentation([0 0; 1 1]), 1, 1)

    isempty_htest(LiftedHRepresentation(ones(0, 2)), 0, 0)
    isempty_htest(LiftedHRepresentation([0 1; 0 2]), 0, 2)
    isempty_htest(LiftedHRepresentation([0 1; 0 2], BitSet(1:2)), 2, 0)
    isempty_htest(LiftedHRepresentation([0 1; 0 2], BitSet([2])), 1, 1)

    @testset "with bad arguments" begin
        A = [1 -1 -1; 0 1 0; 0 0 1]
        ls = BitSet([1])

        @test_throws ErrorException LiftedHRepresentation(A, BitSet([4]))
        ine = copy(LiftedHRepresentation(A, ls))
        @test ine.A == A
        @test ine.A !== A
        #@test linset(ine) == ls
        @test ine.linset !== ls
        @test Polyhedra.similar_type(LiftedHRepresentation{Int, Matrix{Int}}, Float64) == LiftedHRepresentation{Float64, Matrix{Float64}}
        @test Polyhedra.similar_type(LiftedHRepresentation{Int, SparseMatrixCSC{Int, Int}}, 3, Float64) == LiftedHRepresentation{Float64, SparseMatrixCSC{Float64, Int}}

        A2 = [1 1; -1 0; 0 -1]
        b2 = [1, 0, 0]
        linset2 = BitSet([1])
        ine2 = hrep(A2, b2, linset2)

        ine = LiftedHRepresentation(ine2)
        @test ine.A == A
        @test ine.linset == ls

        V = [1 0 1; 1 1 0]
        Vlinset = BitSet(2)
        @test_throws ErrorException LiftedVRepresentation(V, BitSet([4]))
        ext = copy(LiftedVRepresentation(V, Vlinset))
        @test ext.R == V
        @test ext.R !== V
        #@test linset(ext) == Vlinset
        @test ext.linset !== Vlinset

        @test Polyhedra.similar_type(LiftedVRepresentation{Int, SparseMatrixCSC{Int, Int}}, Float64) == LiftedVRepresentation{Float64, SparseMatrixCSC{Float64, Int}}
        @test Polyhedra.similar_type(LiftedVRepresentation{Int, Matrix{Int}}, 3, Float64) == LiftedVRepresentation{Float64, Matrix{Float64}}
    end
end

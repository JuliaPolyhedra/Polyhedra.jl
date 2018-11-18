@testset "MixMatRep" begin
    change_fulldim_test(MixedMatHRep)
    change_fulldim_test(MixedMatVRep)

    isempty_vtest(vrep(zeros(0, 3)), 0, 0)
    isempty_vtest(vrep(zeros(1, 2)), 0, 1)
    isempty_vtest(vrep(zeros(1, 4), ones(2, 4)), 2, 1)
    isempty_vtest(vrep(zeros(2, 1), ones(1, 1)), 1, 2)

    isempty_htest(hrep(ones(0, 3), zeros(0)), 0, 0)
    isempty_htest(hrep(ones(1, 2), zeros(1)), 0, 1)
    isempty_htest(hrep(ones(2, 4), zeros(2), BitSet(1:2)), 2, 0)
    isempty_htest(hrep(ones(3, 1), zeros(3), BitSet([2])), 1, 2)

    @testset "with bad arguments" begin
        A = [1 1; -1 0; 0 -1]
        b = [1, 0, 0]
        linset = BitSet([1])

        @test_throws ErrorException hrep(A, [0, 0], linset)
        @test_throws ErrorException hrep(A, b, BitSet([4]))
        ine = hrep(A, b, linset)
        @test fulldim(ine) == 2
        @test (@inferred Polyhedra.FullDim(ine)) == 2
        @test Polyhedra.coefficient_type(ine) == Int
        @test translate(ine, [1, 0]).b == [2, -1, 0]

        V = [0 1; 1 0]
        @test_throws ErrorException vrep(zeros(0, 2), [1 0]) # V-consistency
        @test_throws ErrorException vrep(V, [1 0 0], BitSet())
        @test_throws ErrorException vrep(V, [1 1], BitSet([2]))
        ext = vrep(V)
        @test fulldim(ext) == 2
        @test (@inferred Polyhedra.FullDim(ine)) == 2
        @test Polyhedra.coefficient_type(ext) == Int
        @test translate(ext, [1, 0]).V == [1 1; 2 0]
    end

    @testset "Cartesian product" begin
        A = [1 2; 3 4; 5 6]
        a = [7, 8, 9]
        B = [10 11 12; 13 14 15]
        b = [16, 17]
        p1 = hrep(A, a, BitSet([2]))
        p2 = hrep(B, b, BitSet([1]))
        p = p1 * p2
        @test p.A == [A[2,:]' zeros(1, 3)
                      zeros(1, 2) B[1, :]'
                      A[[1,3],:] zeros(2, 3)
                      zeros(1, 2) B[2, :]']
        @test p.b == [a[2]; b[1]; a[[1,3]]; b[2]]
        @test p.linset == BitSet([1, 2])
    end
end

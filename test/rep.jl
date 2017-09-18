@testset "Representation tests" begin
    @testset "Simple Representation with bad arguments" begin
        A = [1 1; -1 0; 0 -1]
        b = [1, 0, 0]
        linset = IntSet([1])

        @test_throws ErrorException SimpleHRepresentation(A, [0, 0], linset)
        @test_throws ErrorException SimpleHRepresentation(A, b, IntSet([4]))
        @test_throws ErrorException SimpleHRepresentation{3, Int}(A, b)
        ine = SimpleHRepresentation(A, b, linset)
        @test fulldim(ine) == 2
        @test eltype(ine) == Int
        @test translate(ine, [1, 0]).b == [2, -1, 0]

        V = [0 1; 1 0]
        @test_throws ErrorException SimpleVRepresentation{3, Int}(V, [1 0], IntSet(), IntSet())
        @test_throws ErrorException SimpleVRepresentation(V, [1 0 0], IntSet(), IntSet())
        @test_throws ErrorException SimpleVRepresentation(V, [1 1], IntSet(), IntSet([2]))
        @test_throws ErrorException SimpleVRepresentation(V, [1 1], IntSet([4]), IntSet())
        @test_throws ErrorException SimpleVRepresentation(V, IntSet([4]))
        ext = SimpleVRepresentation(V)
        @test fulldim(ext) == 2
        @test eltype(ext) == Int
        @test translate(ext, [1, 0]).V == [1 1; 2 0]
    end

    @testset "Lifted Representation with bad arguments" begin
        A = [1 -1 -1; 0 1 0; 0 0 1]
        ls = IntSet([1])

        @test_throws ErrorException LiftedHRepresentation(A, IntSet([4]))
        @test_throws ErrorException LiftedHRepresentation{3, Int}(A, ls)
        ine = copy(LiftedHRepresentation(A, ls))
        @test ine.A == A
        @test ine.A !== A
        @test linset(ine) == ls
        @test ine.linset !== ls
        @test changeeltype(LiftedHRepresentation{2, Int}, Float64) == LiftedHRepresentation{2, Float64}
        @test changeboth(LiftedHRepresentation{2, Int}, 3, Float64) == LiftedHRepresentation{3, Float64}

        A2 = [1 1; -1 0; 0 -1]
        b2 = [1, 0, 0]
        linset2 = IntSet([1])
        ine2 = SimpleHRepresentation(A2, b2, linset2)

        ine = LiftedHRepresentation(ine2)
        @test ine.A == A
        @test ine.linset == ls

        V = [1 0 1; 1 1 0]
        Vlinset = IntSet(2)
        @test_throws ErrorException LiftedVRepresentation{3, Int}(V, Vlinset)
        @test_throws ErrorException LiftedVRepresentation(V, IntSet([4]))
        ext = copy(LiftedVRepresentation(V, Vlinset))
        @test ext.R == V
        @test ext.R !== V
        @test linset(ext) == Vlinset
        @test ext.linset !== Vlinset

        @test changeeltype(LiftedVRepresentation{2, Int}, Float64) == LiftedVRepresentation{2, Float64}
        @test changeboth(LiftedVRepresentation{2, Int}, 3, Float64) == LiftedVRepresentation{3, Float64}
    end

    @testset "eltype for some iterators is incorrect #7" begin
        function collecttest(it, exp_type)
            @test eltype(it) == exp_type
            a = collect(it)
            @test typeof(a) = Vector{exp_type}
        end
        hr = SimpleHRepresentation([1 2 3; 4 5 6], [7., 8], IntSet([1]))
        @test eltype(hr) == Float64
        @test eltype(hreps(hr)) == HRepElement{3, Float64}
        @test eltype(ineqs(hr)) == HalfSpace{3, Float64}
        @test eltype(eqs(hr)) == HyperPlane{3, Float64}
        vr = SimpleVRepresentation([1 2; 3 4])
        @test eltype(vreps(vr)) == VRepElement{2, Int}
        @test eltype(points(vr)) == AbstractPoint{2, Int}
        @test eltype(rays(vr)) == AbstractRay{2, Int}
    end

    @testset "Iterating over ineqs of a SimpleHRepresentation broken #9" begin
        A = [1 2; 3 4; 5 6]
        b = [1, 2, 3]
        ineq = [1, 3]
        eq = [2]
        linset = IntSet(2)
        hrep = SimpleHRepresentation(A, b, linset)
        for (i, h) in enumerate(hreps(hrep))
            @test h.a == A[i, :]
            @test h.β == b[i]
            @test isa(h, i in linset ? HyperPlane{2, Int} : HalfSpace{2, Int})
        end
        for (i, h) in enumerate(ineqs(hrep))
            @test h.a == A[ineq[i], :]
            @test h.β == b[ineq[i]]
            @test isa(h, HalfSpace{2, Int})
        end
        for (i, h) in enumerate(eqs(hrep))
            @test h.a == A[eq[i], :]
            @test h.β == b[eq[i]]
            @test isa(h, HyperPlane{2, Int})
        end
    end

    @testset "changefulldim" begin
        N = 5
        M = 10
        T = Int64
        reps = [SimpleHRepresentation{N, T}, SimpleVRepresentation{N, T}, LiftedHRepresentation{N, T}, LiftedVRepresentation{N, T}]
        for rep in reps
            changedrep = changefulldim(rep, M)
            @test fulldim(changedrep) == M
            @test eltype(changedrep) == T
        end
    end

    @testset "Cartesian product" begin
        A = [1 2; 3 4; 5 6]
        a = [7, 8, 9]
        B = [10 11 12; 13 14 15]
        b = [16, 17]
        p1 = SimpleHRepresentation(A, a, IntSet([2]))
        p2 = SimpleHRepresentation(B, b, IntSet([1]))
        p = p1 * p2
        @test p.A == [A zeros(3, 3); zeros(2, 2) B]
        @test p.b == [a; b]
        @test p.linset == IntSet([2, 4])
    end

    @testset "isempty not working correctly for iterators #17" begin
        function vtest(vr, nr, np)
            hasr = nr > 0
            hasp = np > 0
            @test nvreps(vr) == length(vreps(vr)) == nr + np
            @test nrays(vr) == length(rays(vr)) == nr
            @test npoints(vr) == length(points(vr)) == np
            @test hasvreps(vr) == !isempty(vreps(vr)) == hasr || hasp
            @test hasrays(vr) == !isempty(rays(vr)) == hasr
            @test haspoints(vr) == !isempty(points(vr)) == hasp
        end
        vtest(SimpleVRepresentation(zeros(0, 3)), 0, 0)
        vtest(SimpleVRepresentation(zeros(1, 2)), 0, 1)
        vtest(SimpleVRepresentation(zeros(0, 4), ones(2, 4)), 2, 0)
        vtest(SimpleVRepresentation(zeros(2, 1), ones(1, 1)), 1, 2)
        vtest(LiftedVRepresentation(zeros(0, 2)), 0, 0)
        vtest(LiftedVRepresentation([1 0; 1 1]), 0, 2)
        vtest(LiftedVRepresentation([0 1; 0 2]), 2, 0)
        vtest(LiftedVRepresentation([0 0; 1 1]), 1, 1)
        function htest(hr, ne, ni)
            hase = ne > 0
            hasi = ni > 0
            @test nhreps(hr) == length(hreps(hr)) == ne + ni
            @test neqs(hr) == length(eqs(hr)) == ne
            @test nineqs(hr) == length(ineqs(hr)) == ni
            @test hashreps(hr) == !isempty(hreps(hr)) == hase || hasi
            @test haseqs(hr) == !isempty(eqs(hr)) == hase
            @test hasineqs(hr) == !isempty(ineqs(hr)) == hasi
        end
        htest(SimpleHRepresentation(ones(0, 3), zeros(0)), 0, 0)
        htest(SimpleHRepresentation(ones(1, 2), zeros(1)), 0, 1)
        htest(SimpleHRepresentation(ones(2, 4), zeros(2), IntSet(1:2)), 2, 0)
        htest(SimpleHRepresentation(ones(3, 1), zeros(3), IntSet([2])), 1, 2)
        htest(LiftedHRepresentation(ones(0, 2)), 0, 0)
        htest(LiftedHRepresentation([0 1; 0 2]), 0, 2)
        htest(LiftedHRepresentation([0 1; 0 2], IntSet(1:2)), 2, 0)
        htest(LiftedHRepresentation([0 1; 0 2], IntSet([2])), 1, 1)
    end

    @testset "Building rep with different type" begin
        @test eltype(SimpleHRepresentation{2, Float64}([1 2; 3 4], [1, 2])) == Float64
        @test eltype(SimpleVRepresentation{2, Float64}([1 2; 3 4], [1 2; 3 4])) == Float64
        @test eltype(LiftedHRepresentation{1, Float64}([1 2; 3 4])) == Float64
        @test eltype(LiftedVRepresentation{1, Float64}([1 2; 3 4])) == Float64
    end

    @testset "Chebyshev center" begin
        p = SimpleHRepresentation(eye(2), zeros(2))
        @test_throws ErrorException chebyshevcenter(p, lpsolver) # unbounded

        p = SimpleHRepresentation([1 1; -1 -1], [0, -1])
        @test_throws ErrorException chebyshevcenter(p, lpsolver) # empty

        # examples/chebyshevcenter.ipynb
        A = [ 2  1
              2 -1
             -1  2
             -1 -2]
        b = ones(4)
        p = SimpleHRepresentation(A, b)
        c, r = chebyshevcenter(p, lpsolver)
        @test c ≈ [0, 0] atol=1e-6
        @test r ≈ 0.4472135955 atol=1e-6
    end
end

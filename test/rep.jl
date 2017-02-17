@testset "Simple Representation with bad arguments" begin
    A = [1 1; -1 0; 0 -1]
    b = [1, 0, 0]
    linset = IntSet([1])

    @test_throws ErrorException SimpleHRepresentation(A, [0, 0], linset)
    @test_throws ErrorException SimpleHRepresentation(A, [0, 0], IntSet([4]))
    ine = SimpleHRepresentation(A, b, linset)
    @test fulldim(ine) == 2
    @test eltype(ine) == Int

    V = [0 1; 1 0]
    @test_throws ErrorException SimpleVRepresentation(V, [1 0 0], IntSet(), IntSet())
    @test_throws ErrorException SimpleVRepresentation(V, [1 1], IntSet(), IntSet([2]))
    @test_throws ErrorException SimpleVRepresentation(V, [1 1], IntSet([4]), IntSet())
    @test_throws ErrorException SimpleVRepresentation(V, IntSet([4]))
    ext = SimpleVRepresentation(V)
    @test fulldim(ext) == 2
    @test eltype(ext) == Int
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

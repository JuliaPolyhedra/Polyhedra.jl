struct InconsistentVRep{N, T, AT} <: VRepresentation{N, T}
    points::Polyhedra.PointsHull{N, T, AT}
    rays::Polyhedra.RaysHull{N, T, AT}
    function InconsistentVRep{N, T, AT}(sympoints, points, lines, rays) where {N, T, AT}
        new{N, T, AT}(Polyhedra.PointsHull(sympoints, points), Polyhedra.RaysHull(lines, rays))
    end
end
Polyhedra.dualtype(::Type{InconsistentVRep{N,T,AT}}, ::Type{AT}) where {N, T, AT} = Polyhedra.Intersection{N, T, AT}
Polyhedra.arraytype(::Union{InconsistentVRep{N, T, AT}, Type{InconsistentVRep{N, T, AT}}}) where {N, T, AT} = AT
Polyhedra.similar_type(PT::Type{<:InconsistentVRep}, d::FullDim{N}, ::Type{T}) where {N, T} = InconsistentVRep{N, T, Polyhedra.similar_type(Polyhedra.arraytype(PT), d, T)}
Polyhedra.fulltype(::Type{InconsistentVRep{N, T, AT}}) where {N, T, AT} = InconsistentVRep{N, T, AT}
Polyhedra.@subrepelem InconsistentVRep SymPoint points
Polyhedra.@subrepelem InconsistentVRep Point points
Polyhedra.@subrepelem InconsistentVRep Line rays
Polyhedra.@subrepelem InconsistentVRep Ray rays
@testset "Representation tests" begin
    @testset "MixMatRep with bad arguments" begin
        A = [1 1; -1 0; 0 -1]
        b = [1, 0, 0]
        linset = IntSet([1])

        @test_throws ErrorException hrep(A, [0, 0], linset)
        @test_throws ErrorException hrep(A, b, IntSet([4]))
        @test_throws ErrorException MixedMatHRep{3, Int}(A, b, IntSet())
        ine = hrep(A, b, linset)
        @test fulldim(ine) == 2
        @test (@inferred FullDim(ine)) == FullDim{2}()
        @test coefficienttype(ine) == Int
        @test translate(ine, [1, 0]).b == [2, -1, 0]

        V = [0 1; 1 0]
        @test_throws ErrorException MixedMatVRep{3, Int}(V, [1 0], IntSet(), IntSet())
        @test_throws ErrorException vrep(zeros(0, 2), [1 0]) # V-consistency
        @test_throws ErrorException vrep(V, [1 0 0], IntSet(), IntSet())
        @test_throws ErrorException vrep(V, [1 1], IntSet(), IntSet([2]))
        @test_throws ErrorException vrep(V, [1 1], IntSet([4]), IntSet())
        @test_throws ErrorException vrep(V, IntSet([4]))
        ext = vrep(V)
        @test fulldim(ext) == 2
        @test (@inferred FullDim(ine)) == FullDim{2}()
        @test MP.coefficienttype(ext) == Int
        @test translate(ext, [1, 0]).V == [1 1; 2 0]
    end

    @testset "LPHRepresentation with bad arguments" begin
        @test_throws DimensionMismatch LPHRepresentation(ones(2, 2), [1], [1], [1, 2], [1, 2])
        @test_throws DimensionMismatch LPHRepresentation(ones(2, 2), [1, 2], [1, 2], [1], [1])
        @test_throws DimensionMismatch LPHRepresentation{2, Int, Matrix{Int}}(ones(Int, 1, 1), [1], [1], [1], [1])
    end

    @testset "Lifted Representation with bad arguments" begin
        A = [1 -1 -1; 0 1 0; 0 0 1]
        ls = IntSet([1])

        @test_throws ErrorException LiftedHRepresentation(A, IntSet([4]))
        @test_throws ErrorException LiftedHRepresentation{3, Int}(A, ls)
        ine = copy(LiftedHRepresentation(A, ls))
        @test ine.A == A
        @test ine.A !== A
        #@test linset(ine) == ls
        @test ine.linset !== ls
        @test Polyhedra.similar_type(LiftedHRepresentation{2, Int}, Float64) == LiftedHRepresentation{2, Float64}
        @test Polyhedra.similar_type(LiftedHRepresentation{2, Int}, FullDim{3}(), Float64) == LiftedHRepresentation{3, Float64}

        A2 = [1 1; -1 0; 0 -1]
        b2 = [1, 0, 0]
        linset2 = IntSet([1])
        ine2 = hrep(A2, b2, linset2)

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
        #@test linset(ext) == Vlinset
        @test ext.linset !== Vlinset

        @test Polyhedra.similar_type(LiftedVRepresentation{2, Int}, Float64) == LiftedVRepresentation{2, Float64}
        @test Polyhedra.similar_type(LiftedVRepresentation{2, Int}, FullDim{3}(), Float64) == LiftedVRepresentation{3, Float64}
    end

    @testset "eltype for some iterators is incorrect #7" begin
        function collecttest(it, exp_type)
            @test MP.coefficienttype(it) == exp_type
            a = collect(it)
            @test typeof(a) = Vector{exp_type}
        end
        hps = [HyperPlane([1, 2, 3], 7.)]
        shps = [@inferred HyperPlane((@SVector [1, 2, 3]), 7.)]
        @test eltype(shps) == HyperPlane{3, Float64, SVector{3, Float64}}
        hss = [HalfSpace([4, 5, 6.], 8)]
        shss = [@inferred HalfSpace((@SVector [4., 5., 6.]), 8)]
        @test eltype(shss) == HalfSpace{3, Float64, SVector{3, Float64}}
        for (hr, AT) in (((@inferred hrep(hps)), Vector{Float64}),
                         ((@inferred hrep(shps)), SVector{3, Float64}),
                         ((@inferred hrep(hss)), Vector{Float64}),
                         ((@inferred hrep(shss)), SVector{3, Float64}),
                         ((@inferred hrep(hps, hss)), Vector{Float64}),
                         ((@inferred hrep(shps, shss)), SVector{3, Float64}),
                         (hrep([1 2 3; 4 5 6], [7., 8], IntSet([1])), Vector{Float64}),
                         (SimpleHRepresentation([1 2 3; 4 5 6], [7., 8], IntSet([1])), Vector{Float64}))
            @test (@inferred coefficienttype(hr)) == Float64
            @test                                               (@inferred eltype(allhalfspaces(hr)))  == HalfSpace{3, Float64, AT}
            @test                                               (@inferred collect(allhalfspaces(hr))) isa Vector{HalfSpace{3, Float64, AT}}
            @test isempty(allhalfspaces(hr)) == iszero(nallhalfspaces(hr))
            @test (@inferred Polyhedra.halfspacetype(hr))    == (@inferred eltype(halfspaces(hr)))     == HalfSpace{3, Float64, AT}
            @test                                               (@inferred collect(halfspaces(hr)))    isa Vector{HalfSpace{3, Float64, AT}}
            @test isempty(halfspaces(hr)) == iszero(nhalfspaces(hr))
            @test (@inferred Polyhedra.hyperplanetype(hr))   == (@inferred eltype(hyperplanes(hr)))    == HyperPlane{3, Float64, AT}
            @test                                               (@inferred collect(hyperplanes(hr)))   isa Vector{HyperPlane{3, Float64, AT}}
            @test isempty(hyperplanes(hr)) == iszero(nhyperplanes(hr))
            @test_throws DimensionMismatch hr * ones(2, 3) # TODO replace by ones(2, 3) \ hr
        end
        symps = [SymPoint([0, 1])]
        ssymps = [SymPoint(@SVector [0, 1])]
        ps = [[1, 2], [3, 4]]
        sps = [(@SVector [1, 2]), (@SVector [3, 4])]
        rs = [Ray([0, 1])]
        srs = [Ray(@SVector [0, 1])]
        ls = [Line([0, 1])]
        sls = [Line(@SVector [0, 1])]
        @test eltype(sps) == SVector{2, Int}
        for (vr, AT) in (((@inferred vrep(symps)), Vector{Int}),
                         ((@inferred vrep(ssymps)), SVector{2, Int}),
                         (vrep(ps), Vector{Int}),
                         ((@inferred vrep(sps)), SVector{2, Int}),
                         ((@inferred vrep(symps, ps)), Vector{Int}),
                         (convexhull(symps..., ps...), Vector{Int}),
                         (convexhull(ps..., symps...), Vector{Int}),
                         ((@inferred vrep(ssymps, sps)), SVector{2, Int}),
                         ((@inferred convexhull(ssymps..., sps...)), SVector{2, Int}),
                         ((@inferred convexhull(sps..., ssymps...)), SVector{2, Int}),
                         ((@inferred vrep(ls)), Vector{Int}),
                         ((@inferred vrep(sls)), SVector{2, Int}),
                         ((@inferred vrep(rs)), Vector{Int}),
                         ((@inferred vrep(srs)), SVector{2, Int}),
                         ((@inferred vrep(ls, rs)), Vector{Int}),
                         ((@inferred vrep(sls, srs)), SVector{2, Int}),
                         ((@inferred vrep(symps, ps, ls, rs)), Vector{Int}),
                         ((@inferred vrep(ssymps, sps, sls, srs)), SVector{2, Int}),
                         (vrep([1 2; 3 4]), Vector{Int}),
                         (SimpleVRepresentation([1 2; 3 4], zeros(Int, 0, 0), IntSet(), IntSet()), Vector{Int}))
            @test (@inferred coefficienttype(vr)) == Int
            @test                                           (@inferred eltype(allpoints(vr)))  == AT
            @test                                           (@inferred collect(allpoints(vr))) isa Vector{AT}
            @test isempty(allpoints(vr)) == iszero(nallpoints(vr))
            @test (@inferred Polyhedra.sympointtype(vr)) == (@inferred eltype(sympoints(vr)))  == SymPoint{2, Int, AT}
            @test                                           (@inferred collect(sympoints(vr))) isa Vector{SymPoint{2, Int, AT}}
            @test isempty(sympoints(vr)) == iszero(nsympoints(vr))
            @test (@inferred Polyhedra.pointtype(vr))    == (@inferred eltype(points(vr)))     == AT
            @test                                           (@inferred collect(points(vr)))    isa Vector{AT}
            @test isempty(points(vr)) == iszero(npoints(vr))
            @test                                           (@inferred eltype(allrays(vr)))    == Ray{2, Int, AT}
            @test                                           (@inferred collect(allrays(vr)))   isa Vector{Ray{2, Int, AT}}
            @test isempty(allrays(vr)) == iszero(nallrays(vr))
            @test (@inferred Polyhedra.linetype(vr))     == (@inferred eltype(lines(vr)))      == Line{2, Int, AT}
            @test                                           (@inferred collect(lines(vr)))     isa Vector{Line{2, Int, AT}}
            @test isempty(lines(vr)) == iszero(nlines(vr))
            @test (@inferred Polyhedra.raytype(vr))      == (@inferred eltype(rays(vr)))       == Ray{2, Int, AT}
            @test                                           (@inferred collect(rays(vr)))      isa Vector{Ray{2, Int, AT}}
            @test isempty(rays(vr)) == iszero(nrays(vr))
            @test_throws DimensionMismatch ones(2, 1) * vr
        end
    end

    @testset "Iterating over halfspaces of a MixedMatHRep broken #9" begin
        A = [1 2; 3 4; 5 6]
        b = [1, 2, 3]
        halfspace = [1, 3]
        hyperplane = [2]
        linset = IntSet(2)
        hr = hrep(A, b, linset)
        Aall = [3 4; -3 -4; 1 2; 5 6]
        ball = [2, -2, 1, 3]
        for (i, h) in enumerate(allhalfspaces(hr))
            @test h.a == Aall[i, :]
            @test h.β == ball[i]
            @test isa(h, HalfSpace{2, Int})
        end
        for (i, h) in enumerate(halfspaces(hr))
            @test h.a == A[halfspace[i], :]
            @test h.β == b[halfspace[i]]
            @test isa(h, HalfSpace{2, Int})
        end
        for (i, h) in enumerate(hyperplanes(hr))
            @test h.a == A[hyperplane[i], :]
            @test h.β == b[hyperplane[i]]
            @test isa(h, HyperPlane{2, Int})
        end
    end

    @testset "Change FullDim" begin
        N = 5
        M = 10
        T = Int64
        reps = [MixedMatHRep{N, T}, MixedMatVRep{N, T}, LiftedHRepresentation{N, T}, LiftedVRepresentation{N, T}]
        for rep in reps
            changedrep = Polyhedra.similar_type(rep, FullDim{M}())
            @test fulldim(changedrep) == M
            @test (@inferred FullDim(changedrep)) == FullDim{M}()
            @test MP.coefficienttype(changedrep) == T
        end
    end

    @testset "Cartesian product" begin
        A = [1 2; 3 4; 5 6]
        a = [7, 8, 9]
        B = [10 11 12; 13 14 15]
        b = [16, 17]
        p1 = hrep(A, a, IntSet([2]))
        p2 = hrep(B, b, IntSet([1]))
        p = p1 * p2
        @test p.A == [A[2,:]' zeros(1, 3)
                      zeros(1, 2) B[1, :]'
                      A[[1,3],:] zeros(2, 3)
                      zeros(1, 2) B[2, :]']
        @test p.b == [a[2]; b[1]; a[[1,3]]; b[2]]
        @test p.linset == IntSet([1, 2])
    end

    @testset "isempty not working correctly for iterators #17" begin
        function vtest(vr, nr, np)
            hasr = nr > 0
            hasp = np > 0
            @test nallpoints(vr) == length(allpoints(vr)) == np
            @test nsympoints(vr) == length(sympoints(vr)) == 0
            @test npoints(vr) == length(points(vr)) == np
            @test nallrays(vr) == length(allrays(vr)) == nr
            @test nlines(vr) == length(lines(vr)) == 0
            @test nrays(vr) == length(rays(vr)) == nr
            @test hasallpoints(vr) == !isempty(allpoints(vr)) == hasp
            @test hassympoints(vr) == !isempty(sympoints(vr)) == false
            @test haspoints(vr) == !isempty(points(vr)) == hasp
            @test hasallrays(vr) == !isempty(allrays(vr)) == hasr
            @test haslines(vr) == !isempty(lines(vr)) == false
            @test hasrays(vr) == !isempty(rays(vr)) == hasr
        end
        vtest(vrep(zeros(0, 3)), 0, 0)
        vtest(vrep(zeros(1, 2)), 0, 1)
        vtest(vrep(zeros(1, 4), ones(2, 4)), 2, 1)
        vtest(vrep(zeros(2, 1), ones(1, 1)), 1, 2)
        vtest(LiftedVRepresentation(zeros(0, 2)), 0, 0)
        vtest(LiftedVRepresentation([1 0; 1 1]), 0, 2)
        vtest(LiftedVRepresentation([0 1; 0 2]), 2, 0)
        vtest(LiftedVRepresentation([0 0; 1 1]), 1, 1)
        function htest(hr, ne, ni)
            hase = ne > 0
            hasi = ni > 0
            @test nallhalfspaces(hr) == length(allhalfspaces(hr)) == 2ne + ni
            @test nhyperplanes(hr) == length(hyperplanes(hr)) == ne
            @test nhalfspaces(hr) == length(halfspaces(hr)) == ni
            @test hasallhalfspaces(hr) == !isempty(allhalfspaces(hr)) == hase || hasi
            @test hashyperplanes(hr) == !isempty(hyperplanes(hr)) == hase
            @test hashalfspaces(hr) == !isempty(halfspaces(hr)) == hasi
        end
        htest(hrep(ones(0, 3), zeros(0)), 0, 0)
        htest(hrep(ones(1, 2), zeros(1)), 0, 1)
        htest(hrep(ones(2, 4), zeros(2), IntSet(1:2)), 2, 0)
        htest(hrep(ones(3, 1), zeros(3), IntSet([2])), 1, 2)
        htest(LiftedHRepresentation(ones(0, 2)), 0, 0)
        htest(LiftedHRepresentation([0 1; 0 2]), 0, 2)
        htest(LiftedHRepresentation([0 1; 0 2], IntSet(1:2)), 2, 0)
        htest(LiftedHRepresentation([0 1; 0 2], IntSet([2])), 1, 1)
    end

    @testset "Building rep with different type" begin
        @test coefficienttype(MixedMatHRep{2, Float64}([1 2; 3 4], [1, 2], IntSet())) == Float64
        @test coefficienttype(MixedMatVRep{2, Float64}([1 2; 3 4], [1 2; 3 4], IntSet(), IntSet())) == Float64
        @test coefficienttype(LiftedHRepresentation{1, Float64}([1 2; 3 4])) == Float64
        @test coefficienttype(LiftedVRepresentation{1, Float64}([1 2; 3 4])) == Float64
    end

    @testset "Chebyshev center" begin
        p = hrep(eye(2), zeros(2))
        @test_throws ErrorException chebyshevcenter(p, lpsolver) # unbounded

        p = hrep([1 1; -1 -1], [0, -1])
        @test_throws ErrorException chebyshevcenter(p, lpsolver) # empty

        # examples/chebyshevcenter.ipynb
        A = [ 2  1
              2 -1
             -1  2
             -1 -2]
        b = ones(4)
        p = hrep(A, b)
        c, r = chebyshevcenter(p, lpsolver)
        @test c ≈ [0, 0] atol=1e-6
        @test r ≈ 0.4472135955 atol=1e-6

        p = convexhull([0, 0], [0, 1], [1, 0])
        @test_throws ErrorException chebyshevcenter(p) # Not yet implemented
    end

    @testset "V-consistency with iterator constructor" begin
        T = Int
        AT = Vector{Int}
        for VRepType in (Polyhedra.LiftedVRepresentation{2, T},
                         Polyhedra.MixedMatVRep{2, T},
                         Polyhedra.Hull{2, T, AT})
            @test_throws ErrorException VRepType(SymPoint{2, T, AT}[], AT[], [Line([1, 2])])
            @test_throws ErrorException VRepType(SymPoint{2, T, AT}[], AT[], Line{2, T, AT}[], [Ray([1, 2])])
            @test_throws ErrorException VRepType(SymPoint{2, T, AT}[], AT[], [Line([1, 2])], [Ray([1, 2])])
            @test isempty(VRepType(SymPoint{2, T, AT}[], AT[], Line{2, T, AT}[], Ray{2, T, AT}[]))
            @test isempty(VRepType(Line{2, T, AT}[], Ray{2, T, AT}[]))
            v = VRepType([Line([1, 2])])
            @test !hassympoints(v)
            @test collect(points(v)) == [[0, 0]]
            @test collect(lines(v)) == [Line([1, 2])]
            @test !hasrays(v)
        end
        for vinc in (InconsistentVRep{2, T, AT}(SymPoint{2, T, AT}[], AT[], Line{2, T, AT}[], [Ray([1, 2])]),
                     InconsistentVRep{2, T, AT}(SymPoint{2, T, AT}[], AT[], [Line([1, 2])], Ray{2, T, AT}[]),
                     InconsistentVRep{2, T, AT}(SymPoint{2, T, AT}[], AT[], [Line([1, 2])], [Ray([1, 2])]))
            @test_throws ErrorException Polyhedra.checkvconsistency(vinc)
            pinc = polyhedron(vinc)
            @test_throws ErrorException Polyhedra.checkvconsistency(pinc)
        end
    end

    @testset "Conversion with different array type" begin
        @testset "V-representation" begin
            vv = convexhull(@SVector [0, 1]) + conichull(Ray(@SVector [1, 1]), Line(@SVector [1, 0]))
            mv = vrep([0 1], [1 1; 1 0], IntSet(), IntSet([2]))
            generator_fulltest(vv, mv)
            mvv = @inferred typeof(mv)(vv)
            vmv = @inferred typeof(vv)(mv)
            generator_fulltest(vv, mvv)
            generator_fulltest(vmv, vv)
        end
        @testset "H-representation" begin
            vh = HalfSpace((@SVector [1, 0]), 0) ∩ HyperPlane((@SVector [0, 1]), 1)
            mh = hrep([0 1; 1 0], [1, 0], IntSet(1))
            inequality_fulltest(vh, mh)
            mvh = @inferred typeof(mh)(vh)
            vmh = @inferred typeof(vh)(mh)
            inequality_fulltest(vh, mvh)
            inequality_fulltest(vmh, vh)
        end
    end
    @testset "Copy test" begin
        @testset "V-representation" begin
            vr = convexhull(@SVector [0, 1])
            vc = copy(vr)
            @test vc !== vr
            @test vc.points !== vr.points
            generator_fulltest(vc, vr)
            vr = conichull(Line(@SVector [1, 1]), Line(@SVector [1, 0]))
            vc = copy(vr)
            @test vc !== vr
            @test vc.lines !== vr.lines
            generator_fulltest(vc, vr)
            vr = conichull(Line(@SVector [1, 1]), Line(@SVector [1, 0]), Ray(@SVector [1, 1]))
            vc = copy(vr)
            @test vc !== vr
            @test vc.rays !== vr.rays
            @test vc.lines !== vr.lines
            @test vc.lines.lines !== vr.lines.lines
            generator_fulltest(vc, vr)
            vr = convexhull(@SVector [0, 1]) + conichull(Line(@SVector [1, 1]), Ray(@SVector [1, 1]), Line(@SVector [1, 0]))
            vc = copy(vr)
            @test vc !== vr
            @test vc.points !== vr.points
            @test vc.points.points !== vr.points.points
            @test vc.rays !== vr.rays
            @test vc.rays.lines !== vr.rays.lines
            @test vc.rays.lines.lines !== vr.rays.lines.lines
            @test vc.rays.rays !== vr.rays.rays
            generator_fulltest(vc, vr)
        end
        @testset "H-representation" begin
            hr = HyperPlane([1, 0], 1) ∩ HyperPlane([0, 1], -1)
            hc = copy(hr)
            @test hc !== hr
            @test hc.hyperplanes !== hr.hyperplanes
            inequality_fulltest(hc, hr)
            hr = HyperPlane([1, 0], 1) ∩ HyperPlane([1, 1], 0) ∩ HalfSpace([0, 1], -1)
            hc = copy(hr)
            @test hc !== hr
            @test hc.hyperplanes !== hr.hyperplanes
            @test hc.hyperplanes.hyperplanes !== hr.hyperplanes.hyperplanes
            @test hc.halfspaces !== hr.halfspaces
            inequality_fulltest(hc, hr)
        end
    end
end

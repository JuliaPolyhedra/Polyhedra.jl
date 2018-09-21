struct InconsistentVRep{T, AT} <: VRepresentation{T}
    points::Polyhedra.PointsHull{T, AT}
    rays::Polyhedra.RaysHull{T, AT}
    function InconsistentVRep{T, AT}(points, lines, rays) where {T, AT}
        new{T, AT}(Polyhedra.PointsHull{T, AT}(points), Polyhedra.RaysHull(lines, rays))
    end
end
Polyhedra.FullDim(rep::InconsistentVRep{T, AT}) where {T, AT <: StaticArrays.SVector} = FullDim(AT)
Polyhedra.FullDim(rep::InconsistentVRep) = Polyhedra.FullDim_rec(rep.points, rep.rays)
Polyhedra.dualtype(::Type{InconsistentVRep{T,AT}}, ::Type{AT}) where {T, AT} = Polyhedra.Intersection{T, AT}
Polyhedra.hvectortype(::Type{InconsistentVRep{T, AT}}) where {T, AT} = AT
Polyhedra.vvectortype(::Type{InconsistentVRep{T, AT}}) where {T, AT} = AT
Polyhedra.similar_type(PT::Type{<:InconsistentVRep}, d::Polyhedra.FullDim, ::Type{T}) where {T} = InconsistentVRep{T, Polyhedra.similar_type(Polyhedra.hvectortype(PT), d, T)}
Polyhedra.fulltype(::Type{InconsistentVRep{T, AT}}) where {T, AT} = InconsistentVRep{T, AT}
#Polyhedra.@subrepelem InconsistentVRep SymPoint points
Polyhedra.@subrepelem InconsistentVRep Point points
Polyhedra.@subrepelem InconsistentVRep Line rays
Polyhedra.@subrepelem InconsistentVRep Ray rays
@testset "Representation tests" begin
    @testset "MixMatRep with bad arguments" begin
        A = [1 1; -1 0; 0 -1]
        b = [1, 0, 0]
        linset = BitSet([1])

        @test_throws ErrorException hrep(A, [0, 0], linset)
        @test_throws ErrorException hrep(A, b, BitSet([4]))
        ine = hrep(A, b, linset)
        @test fulldim(ine) == 2
        @test (@inferred Polyhedra.FullDim(ine)) == 2
        @test coefficienttype(ine) == Int
        @test translate(ine, [1, 0]).b == [2, -1, 0]

        V = [0 1; 1 0]
        @test_throws ErrorException vrep(zeros(0, 2), [1 0]) # V-consistency
        @test_throws ErrorException vrep(V, [1 0 0], BitSet())
        @test_throws ErrorException vrep(V, [1 1], BitSet([2]))
        ext = vrep(V)
        @test fulldim(ext) == 2
        @test (@inferred Polyhedra.FullDim(ine)) == 2
        @test MP.coefficienttype(ext) == Int
        @test translate(ext, [1, 0]).V == [1 1; 2 0]
    end

    @testset "LPHRepresentation with bad arguments" begin
        @test_throws DimensionMismatch LPHRepresentation(ones(2, 2), [1], [1], [1, 2], [1, 2])
        @test_throws DimensionMismatch LPHRepresentation(ones(2, 2), [1, 2], [1, 2], [1], [1])
    end

    @testset "Lifted Representation with bad arguments" begin
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

    @testset "eltype for some iterators is incorrect #7" begin
        function collecttest(it, exp_type)
            @test MP.coefficienttype(it) == exp_type
            a = collect(it)
            @test typeof(a) = Vector{exp_type}
        end
        hps = [HyperPlane([1, 2, 3], 7.)]
        shps = [@inferred HyperPlane((@SVector [1, 2, 3]), 7.)]
        @test eltype(shps) == HyperPlane{Float64, SVector{3, Float64}}
        hss = [HalfSpace([4, 5, 6.], 8)]
        shss = [@inferred HalfSpace((@SVector [4., 5., 6.]), 8)]
        @test eltype(shss) == HalfSpace{Float64, SVector{3, Float64}}
        function htest(hr::Polyhedra.HRepresentation, AT::Type{<:AbstractVector})
            @test (@inferred coefficienttype(hr)) == Float64
            @test                                               (@inferred eltype(allhalfspaces(hr)))  == HalfSpace{Float64, AT}
            @test                                               (@inferred collect(allhalfspaces(hr))) isa Vector{HalfSpace{Float64, AT}}
            @test isempty(allhalfspaces(hr)) == iszero(nallhalfspaces(hr))
            @test (@inferred Polyhedra.halfspacetype(hr))    == (@inferred eltype(halfspaces(hr)))     == HalfSpace{Float64, AT}
            @test                                               (@inferred collect(halfspaces(hr)))    isa Vector{HalfSpace{Float64, AT}}
            @test isempty(halfspaces(hr)) == iszero(nhalfspaces(hr))
            @test (@inferred Polyhedra.hyperplanetype(hr))   == (@inferred eltype(hyperplanes(hr)))    == HyperPlane{Float64, AT}
            @test                                               (@inferred collect(hyperplanes(hr)))   isa Vector{HyperPlane{Float64, AT}}
            @test isempty(hyperplanes(hr)) == iszero(nhyperplanes(hr))
            @test_throws DimensionMismatch ones(2, 3) \ hr
        end
        htest((@inferred hrep(hps)), Vector{Float64})
        htest((@inferred hrep(shps)), SVector{3, Float64})
        htest((@inferred hrep(hss)), Vector{Float64})
        htest((@inferred hrep(shss)), SVector{3, Float64})
        htest((@inferred hrep(hps, hss)), Vector{Float64})
        htest((@inferred hrep(shps, shss)), SVector{3, Float64})
        htest(hrep([1 2 3; 4 5 6], [7., 8], BitSet([1])), Vector{Float64})
        htest(hrep(spzeros(2, 3), [7., 8], BitSet([1])), SparseVector{Float64, Int})
        htest(hrep([1 2 3; 4 5 6], [7., 8]), Vector{Float64})
        ps = [[1, 2], [3, 4]]
        sps = [(@SVector [1, 2]), (@SVector [3, 4])]
        rs = [Ray([0, 1])]
        srs = [Ray(@SVector [0, 1])]
        ls = [Line([0, 1])]
        sls = [Line(@SVector [0, 1])]
        @test eltype(sps) == SVector{2, Int}
        function vtest(vr::VRepresentation, AT::Type{<:AbstractVector})
            @test (@inferred coefficienttype(vr)) == Int
            @test (@inferred Polyhedra.pointtype(vr))    == (@inferred eltype(points(vr)))     == AT
            @test                                           (@inferred collect(points(vr)))    isa Vector{AT}
            @test isempty(points(vr)) == iszero(npoints(vr))
            @test                                           (@inferred eltype(allrays(vr)))    == Ray{Int, AT}
            @test                                           (@inferred collect(allrays(vr)))   isa Vector{Ray{Int, AT}}
            @test isempty(allrays(vr)) == iszero(nallrays(vr))
            @test (@inferred Polyhedra.linetype(vr))     == (@inferred eltype(lines(vr)))      == Line{Int, AT}
            @test                                           (@inferred collect(lines(vr)))     isa Vector{Line{Int, AT}}
            @test isempty(lines(vr)) == iszero(nlines(vr))
            @test (@inferred Polyhedra.raytype(vr))      == (@inferred eltype(rays(vr)))       == Ray{Int, AT}
            @test                                           (@inferred collect(rays(vr)))      isa Vector{Ray{Int, AT}}
            @test isempty(rays(vr)) == iszero(nrays(vr))
            @test_throws DimensionMismatch ones(2, 1) * vr
        end
        vtest(vrep(ps), Vector{Int})
        vtest((@inferred vrep(sps)), SVector{2, Int})
        vtest(convexhull(ps...), Vector{Int})
        vtest((@inferred convexhull(sps...)), SVector{2, Int})
        vtest((@inferred vrep(ls)), Vector{Int})
        vtest((@inferred vrep(sls)), SVector{2, Int})
        vtest((@inferred vrep(rs)), Vector{Int})
        vtest((@inferred vrep(srs)), SVector{2, Int})
        vtest((@inferred vrep(ls, rs)), Vector{Int})
        vtest((@inferred vrep(sls, srs)), SVector{2, Int})
        vtest((@inferred vrep(ps, ls, rs)), Vector{Int})
        vtest((@inferred vrep(sps, sls, srs)), SVector{2, Int})
        vtest(vrep([1 2; 3 4]), Vector{Int})
        vtest(vrep(spzeros(Int, 2, 2)), SparseVector{Int, Int})
        vtest(vrep([1 2; 3 4], zeros(Int, 0, 0), BitSet()), Vector{Int})
        vtest(vrep([1 2; 3 4], zeros(Int, 0, 0)), Vector{Int})
        vtest(vrep([1 2; 3 4]), Vector{Int})
    end

    @testset "Iterating over halfspaces of a MixedMatHRep broken #9" begin
        A = [1 2; 3 4; 5 6]
        b = [1, 2, 3]
        halfspace = [1, 3]
        hyperplane = [2]
        linset = BitSet(2)
        hr = hrep(A, b, linset)
        Aall = [3 4; -3 -4; 1 2; 5 6]
        ball = [2, -2, 1, 3]
        for (i, h) in enumerate(allhalfspaces(hr))
            @test h.a == Aall[i, :]
            @test h.β == ball[i]
            @test isa(h, HalfSpace{Int})
        end
        for (i, h) in enumerate(halfspaces(hr))
            @test h.a == A[halfspace[i], :]
            @test h.β == b[halfspace[i]]
            @test isa(h, HalfSpace{Int})
        end
        for (i, h) in enumerate(hyperplanes(hr))
            @test h.a == A[hyperplane[i], :]
            @test h.β == b[hyperplane[i]]
            @test isa(h, HyperPlane{Int})
        end
    end

    @testset "Change Polyhedra.FullDim" begin
        N = 5
        M = 10
        T = Int64
        reps = [MixedMatHRep{T, Matrix{T}}, MixedMatVRep{T, Matrix{T}}, LiftedHRepresentation{T, Matrix{T}}, LiftedVRepresentation{T, Matrix{T}}]
        for rep in reps
            @test fulldim(rep) == -1
            @test MP.coefficienttype(rep) == T
            changedrep = Polyhedra.similar_type(rep, M)
            @test fulldim(changedrep) == -1
            @test Polyhedra.FullDim(changedrep) == -1
            #@test (@inferred Polyhedra.FullDim(changedrep)) == M
            @test MP.coefficienttype(changedrep) == T
        end
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

    @testset "isempty not working correctly for iterators #17" begin
        function vtest(vr, nr, np)
            hasr = nr > 0
            hasp = np > 0
            @test npoints(vr) == length(points(vr)) == np
            @test nallrays(vr) == length(allrays(vr)) == nr
            @test nlines(vr) == length(lines(vr)) == 0
            @test nrays(vr) == length(rays(vr)) == nr
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
        htest(hrep(ones(2, 4), zeros(2), BitSet(1:2)), 2, 0)
        htest(hrep(ones(3, 1), zeros(3), BitSet([2])), 1, 2)
        htest(LiftedHRepresentation(ones(0, 2)), 0, 0)
        htest(LiftedHRepresentation([0 1; 0 2]), 0, 2)
        htest(LiftedHRepresentation([0 1; 0 2], BitSet(1:2)), 2, 0)
        htest(LiftedHRepresentation([0 1; 0 2], BitSet([2])), 1, 1)
    end

    @testset "Building rep with different type" begin
        @test coefficienttype(MixedMatHRep{Float64}([1 2; 3 4], [1, 2], BitSet())) == Float64
        @test coefficienttype(MixedMatVRep{Float64}([1 2; 3 4], [1 2; 3 4], BitSet())) == Float64
        @test coefficienttype(LiftedHRepresentation{Float64}([1 2; 3 4])) == Float64
        @test coefficienttype(LiftedVRepresentation{Float64}([1 2; 3 4])) == Float64
    end

    @testset "Chebyshev center" begin
        p = hrep(eye(2), zeros(2))
        @test_throws ErrorException chebyshevcenter(p, lp_solver) # unbounded

        p = hrep([1 1; -1 -1], [0, -1])
        @test_throws ErrorException chebyshevcenter(p, lp_solver) # empty

        # examples/chebyshevcenter.ipynb
        A = [ 2  1
              2 -1
             -1  2
             -1 -2]
        b = ones(4)
        p = hrep(A, b)
        c, r = chebyshevcenter(p, lp_solver)
        @test c ≈ [0, 0] atol=1e-6
        @test r ≈ 0.4472135955 atol=1e-6

        p = convexhull([0, 0], [0, 1], [1, 0])
        @test_throws ErrorException chebyshevcenter(p) # Not yet implemented
    end

    @testset "V-consistency with iterator constructor" begin
        T = Int
        AT = Vector{Int}
        for VRepType in (Polyhedra.LiftedVRepresentation{T, Matrix{T}},
                         Polyhedra.MixedMatVRep{T, Matrix{T}},
                         Polyhedra.Hull{T, AT})
            @test_throws ErrorException VRepType(AT[], [Line([1, 2])])
            @test_throws ErrorException VRepType(AT[], Line{T, AT}[], [Ray([1, 2])])
            @test_throws ErrorException VRepType(AT[], [Line([1, 2])], [Ray([1, 2])])
            @test isempty(VRepType(AT[], Line{T, AT}[], Ray{T, AT}[]))
            @test isempty(VRepType(Line{T, AT}[], Ray{T, AT}[]))
            v = VRepType([Line([1, 2])])
            @test collect(points(v)) == [[0, 0]]
            @test collect(lines(v)) == [Line([1, 2])]
            @test !hasrays(v)
        end
        for vinc in (InconsistentVRep{T, AT}(AT[], Line{T, AT}[], [Ray([1, 2])]),
                     InconsistentVRep{T, AT}(AT[], [Line([1, 2])], Ray{T, AT}[]),
                     InconsistentVRep{T, AT}(AT[], [Line([1, 2])], [Ray([1, 2])]))
            @test_throws ErrorException Polyhedra.checkvconsistency(vinc)
            pinc = polyhedron(vinc)
            @test_throws ErrorException Polyhedra.checkvconsistency(pinc)
        end
    end
    @testset "Combination of different coefficient type" begin
        @testset "V-representation" begin
            generator_fulltest(convexhull([1, 0], Line([0, 1.])), convexhull(Line([0, 1]), [1, 0.]))
            @test conichull(convexhull([1, 0.]), conichull([0, 1])) isa Polyhedra.RaysHull{Float64}
            @test convexhull(convexhull([1, 0]), Line([0, 1.])) isa Polyhedra.Hull{Float64}
            @test convexhull(Line([0, 1.]), convexhull([1, 0])) isa Polyhedra.Hull{Float64}
            @test convexhull(convexhull(Line([0, 1.])), [1, 0]) isa Polyhedra.Hull{Float64}
            @test convexhull(convexhull(Line([1, 0])), Line([0, 1.])) isa Polyhedra.LinesHull{Float64}
            @test convexhull(conichull([1, 0.]), Line([0, 1])) isa Polyhedra.RaysHull{Float64}
            @test convexhull(conichull([1, 0.]), [0, 1]) isa Polyhedra.Hull{Float64}
        end
        @testset "H-representation" begin
            #FIXME inference fails @test (@inferred (HyperPlane([1, 1], 0) ∩ HyperPlane([1, 0], 1)) ∩ (HyperPlane([1, 1], 0) ∩ HyperPlane([1., 0.], 1))) isa Polyhedra.HyperPlanesIntersection{Float64}
            @test ((HyperPlane([1, 1], 0) ∩ HyperPlane([1, 0], 1)) ∩ (HyperPlane([1, 1], 0) ∩ HyperPlane([1., 0.], 1))) isa Polyhedra.HyperPlanesIntersection{Float64}
            @test HyperPlane([1, 1], 0) ∩ HalfSpace([1., 0.], 1.) isa Polyhedra.Intersection{Float64}
            #FIXME inference fails @test (@inferred (HalfSpace([1., 0.], 1.) ∩ (HyperPlane([1, 1], 0) ∩ HyperPlane([1, 0], 1)))) isa Polyhedra.Intersection{Float64}
            @test ((HalfSpace([1., 0.], 1.) ∩ (HyperPlane([1, 1], 0) ∩ HyperPlane([1, 0], 1)))) isa Polyhedra.Intersection{Float64}
        end
    end
    @testset "Conversion with different array type" begin
        @testset "V-representation" begin
            vv = convexhull(@SVector [0, 1]) + conichull((@SVector [1, 1]), Line(@SVector [1, 0]))
            mv = vrep([0 1], [1 1; 1 0], BitSet([2]))
            generator_fulltest(vv, mv)
            mvv = @inferred typeof(mv)(vv)
            vmv = @inferred typeof(vv)(mv)
            generator_fulltest(vv, mvv)
            generator_fulltest(vmv, vv)
        end
        @testset "H-representation" begin
            vh = HalfSpace((@SVector [1, 0]), 0) ∩ HyperPlane((@SVector [0, 1]), 1)
            mh = hrep([0 1; 1 0], [1, 0], BitSet(1))
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
            vr = conichull(Line(@SVector [1, 1]), Line(@SVector [1, 0]), @SVector [1, 1])
            vc = copy(vr)
            @test vc !== vr
            @test vc.rays !== vr.rays
            @test vc.lines !== vr.lines
            @test vc.lines.lines !== vr.lines.lines
            generator_fulltest(vc, vr)
            vr = convexhull(@SVector [0, 1]) + conichull(Line(@SVector [1, 1]), (@SVector [1, 1]), Line(@SVector [1, 0]))
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
    @testset "Preserving sparsity" begin
        for h in (HalfSpace(sparsevec([1], [1], 2), 1) ∩ HyperPlane(sparsevec([2], [-1]), 3), hrep(sparse([1, 2], [1, 2], [1, -1]), [1, 3]))
            @test Polyhedra.Intersection(h) isa Polyhedra.Intersection{Int, SparseVector{Int, Int}}
            @test LPHRepresentation(h) isa LPHRepresentation{Int, SparseMatrixCSC{Int, Int}}
            @test MixedMatHRep(h) isa MixedMatHRep{Int, SparseMatrixCSC{Int, Int}}
            @test LiftedHRepresentation(h) isa LiftedHRepresentation{Int, SparseMatrixCSC{Int, Int}}
            @test Polyhedra.Intersection{Float64}(h) isa Polyhedra.Intersection{Float64, SparseVector{Float64, Int}}
            @test LPHRepresentation{Float64}(h) isa LPHRepresentation{Float64, SparseMatrixCSC{Float64, Int}}
            @test MixedMatHRep{Float64}(h) isa MixedMatHRep{Float64, SparseMatrixCSC{Float64, Int}}
            @test LiftedHRepresentation{Float64}(h) isa LiftedHRepresentation{Float64, SparseMatrixCSC{Float64, Int}}
        end
        for v in (convexhull(sparsevec([1], [1], 2), sparsevec([2], [-1])), vrep(sparse([1, 2], [1, 2], [1, -1])))
            @test Polyhedra.Hull(v) isa Polyhedra.Hull{Int, SparseVector{Int, Int}}
            @test MixedMatVRep(v) isa MixedMatVRep{Int, SparseMatrixCSC{Int, Int}}
            @test LiftedVRepresentation(v) isa LiftedVRepresentation{Int, SparseMatrixCSC{Int, Int}}
            @test Polyhedra.Hull{Float64}(v) isa Polyhedra.Hull{Float64, SparseVector{Float64, Int}}
            @test MixedMatVRep{Float64}(v) isa MixedMatVRep{Float64, SparseMatrixCSC{Float64, Int}}
            @test LiftedVRepresentation{Float64}(v) isa LiftedVRepresentation{Float64, SparseMatrixCSC{Float64, Int}}
        end
    end
end

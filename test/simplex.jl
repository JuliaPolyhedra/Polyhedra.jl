import JuMP
function simplextest(lib::Polyhedra.Library)
    hsim = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HyperPlane([1, 1], 1)
    vsim = convexhull([0, 1], [1, 0])
    poly1 = polyhedron(hsim, lib)
    @test !isempty(poly1)
    center, radius = chebyshevcenter(poly1)
    @test center ≈ [1/2, 1/2]
    @test radius ≈ 1/2
    @test dim(poly1) == 1 # FIXME doing dim earlier makes chebyshevcenter fail
    inequality_fulltest(poly1, hsim)
    generator_fulltest(poly1, vsim)

    # Test incidence
    hpidx = first(eachindex(hyperplanes(poly1)))
    for ps in (incidentpoints(poly1, hpidx), get.(poly1, incidentpointindices(poly1, hpidx)))
        @test (ps == [[0, 1], [1, 0]] || ps == [[1, 0], [0, 1]])
    end
    for hsidx in eachindex(halfspaces(poly1))
        h = get(poly1, hsidx)
        if dot(h.a, [1, 0]) ≈ h.β
            expps = [[1, 0]]
        else
            @assert dot(h.a, [0, 1]) ≈ h.β
            expps = [[0, 1]]
        end
        for ps in (incidentpoints(poly1, hsidx), get.(poly1, incidentpointindices(poly1, hsidx)))
            @test ps == expps
        end
    end
    for pidx in eachindex(points(poly1))
        for hps in (incidenthyperplanes(poly1, pidx), get.(poly1, incidenthyperplaneindices(poly1, pidx)))
            @test hps == [HyperPlane([1, 1], 1)]
        end
        for hss in (incidenthalfspaces(poly1, pidx), get.(poly1, incidenthalfspaceindices(poly1, pidx)))
            h = hss[1]
            @test dot(h.a, get(poly1, pidx)) ≈ h.β
        end
    end

    @testset "Optimize with objective : max 2x_1" begin
        model, T = Polyhedra.layered_optimizer(Polyhedra.linear_objective_solver(poly1))
        x = MOI.add_variables(model, fulldim(poly1))
        MOI.add_constraint(model, MOI.VectorOfVariables(x),
                           Polyhedra.PolyhedraOptSet(poly1))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                MOI.ScalarAffineFunction(
                    [MOI.ScalarAffineTerm{T}(2, x[1])], zero(T)))
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) == 2
        @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get(model, MOI.VariablePrimal(), x[1]) == 1
        @test MOI.get(model, MOI.VariablePrimal(), x[2]) == 0
    end
    @testset "Optimize with objective : max x_1 + 3x_2" begin
        model, T = Polyhedra.layered_optimizer(Polyhedra.linear_objective_solver(poly1))
        x = MOI.add_variables(model, fulldim(poly1))
        MOI.add_constraint(model, MOI.VectorOfVariables(x),
                           Polyhedra.PolyhedraOptSet(poly1))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                MOI.ScalarAffineFunction(
                    [MOI.ScalarAffineTerm{T}(1, x[1]),
                     MOI.ScalarAffineTerm{T}(3, x[2])], zero(T)))
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) == 3
        @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get(model, MOI.VariablePrimal(), x[1]) == 0
        @test MOI.get(model, MOI.VariablePrimal(), x[2]) == 1
    end

    poly2 = polyhedron(vsim, lib)
    @test dim(poly2) == 1
    @test !isempty(poly2)
    inequality_fulltest(poly2, hsim)
    generator_fulltest(poly2, vsim)

    @testset "x_1 cannot be 2" begin
        hempty = hsim ∩ HyperPlane([1, 0], 2)
        poly = polyhedron(hempty, lib)
        @test isempty(poly)
    end

    # We now add the vertex (0, 0)
    ext0 = convexhull([0, 0])
    @test collect(points(translate(ext0, [1, 0]))) == [[1, 0]]

    htri = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HalfSpace([1, 1], 1)
    vtri = convexhull(vsim, ext0)
    convexhull!(poly1, ext0)
    inequality_fulltest(poly1, htri)
    generator_fulltest(poly1, vtri)
    convexhull!(poly2, ext0)
    inequality_fulltest(poly2, htri)
    generator_fulltest(poly2, vtri)

    # nonnegative orthant cut by x_1 + x_2 = 1
    vray = conichull([1, 0], [0, 1])
    poly3 = polyhedron(vray, lib)
    @test_throws ErrorException chebyshevcenter(poly3)
    @test dim(poly3) == 2

    @testset "Optimize with objective : min x_1 + x_2" begin
        model, T = Polyhedra.layered_optimizer(Polyhedra.linear_objective_solver(poly3))
        x = MOI.add_variables(model, fulldim(poly3))
        MOI.add_constraint(model, MOI.VectorOfVariables(x),
                           Polyhedra.PolyhedraOptSet(poly3))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                MOI.ScalarAffineFunction(
                    [MOI.ScalarAffineTerm{T}(1, x[1]),
                     MOI.ScalarAffineTerm{T}(1, x[2])], zero(T)))
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) == 0
        @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get(model, MOI.VariablePrimal(), x[1]) == 0
        @test MOI.get(model, MOI.VariablePrimal(), x[2]) == 0
    end

    @testset "Optimize with objective : max x_2" begin
        model, T = Polyhedra.layered_optimizer(Polyhedra.linear_objective_solver(poly3))
        x = MOI.add_variables(model, fulldim(poly3))
        MOI.add_constraint(model, MOI.VectorOfVariables(x),
                           Polyhedra.PolyhedraOptSet(poly3))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                MOI.ScalarAffineFunction(
                    [MOI.ScalarAffineTerm{T}(1, x[2])], zero(T)))
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
        @test MOI.get(model, MOI.ObjectiveValue()) == 1
        @test MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
        @test MOI.get(model, MOI.VariablePrimal(), x[1]) == 0
        @test MOI.get(model, MOI.VariablePrimal(), x[2]) == 1
    end

    hcutel = HyperPlane([1, 1], 1)
    hcut = intersect(hcutel)
    vcut = convexhull([1, 0]) + conichull(Line([1, -1]))
    @test !ininterior([1/2, 1/2], hcut)
    @test inrelativeinterior([1/2, 1/2], hcut)

    poly4 = copy(poly3)

    polycut3 = poly3 ∩ hcutel
    @test dim(polycut3) == 1
    inequality_fulltest(polycut3, hsim)
    generator_fulltest(polycut3, vsim)
    intersect!(poly3, hcutel)
    @test dim(poly3) == 1
    inequality_fulltest(poly3, hsim)
    generator_fulltest(poly3, vsim)

    # It should not have been cut as it is a copy of poly3
    @test dim(poly4) == 2

    polycut4 = poly4 ∩ hcut
    @test dim(polycut4) == 1
    inequality_fulltest(polycut4, hsim)
    generator_fulltest(polycut4, vsim)
    intersect!(poly4, hcut)
    @test dim(poly4) == 1
    inequality_fulltest(poly4, hsim)
    generator_fulltest(poly4, vsim)

    # FIXME needs float currently but should be updated
    # poly4 = project(poly1, [1; 0])
    # inequality_fulltest(poly4, [-1; 1], [0, 1], BitSet())
    # generator_fulltest(poly4, [0; 1], [])

    #\
    # \
    # |\
    # |_\
    #    \
    hlin = HalfSpace([1, 1], 1) ∩ HalfSpace([-1, -1], -1)
    plin = polyhedron(hlin, lib)
    @test dim(plin) == 1
    inequality_fulltest(plin, hcut)
    generator_fulltest(plin, vcut)
    #ineout = hrep(plin)
    #@test linset(ineout) == BitSet(1)
    vlin = convexhull([1, 0]) + conichull([1, -1], [-1, 1])
    plin = polyhedron(vlin, lib)
    inequality_fulltest(plin, hcut)
    generator_fulltest(plin, vcut)
    #extout = vrep(plin)
    #@test linset(extout) == BitSet(1)
end

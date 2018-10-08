using LinearAlgebra # for I
function boardtest(lib::Polyhedra.Library)
    A1 = -Matrix(1I, 9, 9) # x >= 0
    b1 = zeros(Int, 9)
    A2 = Matrix(1I, 9, 9) # x <= 1
    b2 = ones(Int, 9)
    A3 = zeros(Int, 9, 9)
    b3 = 3 * ones(Int, 9)
    i = 1
    for a = 1:3
        for b = (a+1):3
            for c = 1:3
                for d = (c+1):3
                    ac = a + (c-1) * 3
                    ad = a + (d-1) * 3
                    bc = b + (c-1) * 3
                    bd = b + (d-1) * 3
                    A3[i, ac] = 1
                    A3[i, ad] = 1
                    A3[i, bc] = 1
                    A3[i, bd] = 1
                    i += 1
                end
            end
        end
    end
    A = [A1; A2; A3]
    b = [b1; b2; b3]
    ine = hrep(A, b)
    poly = polyhedron(ine, lib)
    @test nrays(poly) == 0
    @test npoints(poly) == 614
    @test !isempty(poly)
    ext  = MixedMatVRep(vrep(poly))
    target = ones(Int, 9) * (3 // 4)
    ok = false
    for i = 1:size(ext.V, 1)
        # In julia v0.4 [i,:] returns a row matrix and in v0.5 it is
        # a 1D vector hence the use of vec
        if vec(ext.V[i,:]) == target
            ok = true
        end
    end
    @test ok

    cutA = ones(Int, 1, 9)
    cutb = 6
    Acut = [cutA; A]
    bcut = [cutb; b]
    inecut = hrep(Acut, bcut)
    polycut = polyhedron(inecut, lib)
    @test !isempty(polycut)
    #(isredundant, certificate) = isredundantinequality(polycut, 1)
    #@test !isredundant
    #@test certificate == target
    @test !isredundant(polycut, first(eachindex(halfspaces(polycut))))
    # @test BitSet() == gethredundantindices(polycut) # TODO reactivate it when I figure out why it makes LRS tests fail
    #(issredundant, scertificate) = isstronglyredundantinequality(polycut, 1)
    #@test !issredundant
    #@test scertificate == target
    #@test BitSet([]) == getstronglyredundantinequalities(polycut)
end

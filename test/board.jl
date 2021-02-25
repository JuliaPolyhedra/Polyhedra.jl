using LinearAlgebra # for I
function boardtest(lib::Polyhedra.Library, n)
    N = n^2
    A1 = -Matrix(1I, N, N) # x >= 0
    b1 = zeros(Int, N)
    A2 = Matrix(1I, N, N) # x <= 1
    b2 = ones(Int, N)
    A3 = zeros(Int, N, N)
    b3 = n * ones(Int, N)
    i = 1
    for a = 1:n
        for b = (a + 1):n
            for c = 1:n
                for d = (c + 1):n
                    ac = a + (c - 1) * n
                    ad = a + (d - 1) * n
                    bc = b + (c - 1) * n
                    bd = b + (d - 1) * n
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
    @test !hasrays(poly)
    @test nrays(poly) == 0
    expected = [2, 11, 614]
    if n in eachindex(expected)
        @test npoints(poly) == expected[n]
    end
    @test !isempty(poly)
    if n == 3
        ext  = MixedMatVRep(vrep(poly))
        target = ones(Int, N) * (3 // 4)
        @test any(i -> ext.V[i,:] == target, 1:size(ext.V, 1))
    end

    if n == 3
        cutA = ones(Int, 1, N)
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
end

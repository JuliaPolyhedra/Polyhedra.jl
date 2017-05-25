function doctest(lib::PolyhedraLibrary)
    A = [1 1;1 -1;-1 0]
    b = [1,0,0]
    V = [1//2 1//2; 0 1; 0 0]
    hrep = SimpleHRepresentation(A, b)
    p = polyhedron(hrep, lib)
    inequality_fulltest(p, A, b, IntSet())
    generator_fulltest(p, V)
    Ap = [1; -1]
    bp = [1, 0]
    inequality_fulltest(eliminate(p, [1]), Ap, bp, IntSet())
    inequality_fulltest(project(p, [2]), Ap, bp, IntSet())
    inequality_fulltest(eliminate(p, [1], Val{:ProjectGenerators}), Ap, bp, IntSet())
    inequality_fulltest(project(p, [2], Val{:ProjectGenerators}), Ap, bp, IntSet())
    if implementseliminationmethod(p, Val{:BlockElimination})
        inequality_fulltest(eliminate(p, [1], Val{:BlockElimination}), Ap, bp, IntSet())
        inequality_fulltest(project(p, [2], Val{:BlockElimination}), Ap, bp, IntSet())
    end
    inequality_fulltest(project(p, [0; 1]), Ap, bp, IntSet())
end

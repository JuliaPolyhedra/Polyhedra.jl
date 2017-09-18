function issue48test(lib::PolyhedraLibrary)
    v = [-5 -15; 2 3; -4 2; -2 5]
    V = SimpleVRepresentation(v)
    p = polyhedron(V, lib)
    b = [70, 16, 8, 15]
    A = [-17  1
          -3  2
           1  2
          18 -7]
    inequality_fulltest(p, A, b, IntSet([]))
    generator_fulltest(p, v)
end

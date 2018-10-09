function issue48test(lib::Polyhedra.Library)
    v = [-5 -15; 2 3; -4 2; -2 5]
    V = vrep(v)
    p = polyhedron(V, lib)
    b = [70, 16, 8, 15]
    A = [-17  1
          -3  2
           1  2
          18 -7]
    inequality_fulltest(p, A, b, BitSet([]))
    generator_fulltest(p, v)

    ind = [3, 4, 6, 1, 7, 8, 9, 2, 5]
    v = [3 13; -1 0; 16 9; 21 1; -5 -15; 9 0; 10 4; 5 -4; 0 3]
    A = [-13   4
           4  13
           8   5
         -15   4
           8 -13]
    b = [13, 181, 173, 15, 155]

    V1 = vrep(v)
    p = polyhedron(V1, lib)
    inequality_fulltest(p, A, b, BitSet([]))

    V2 = vrep(v[ind, :])
    p = polyhedron(V2, lib)
    inequality_fulltest(p, A, b, BitSet([]))
end

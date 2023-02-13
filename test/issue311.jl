function issue311test(lib::Polyhedra.Library)
    W = [
        [4, 8],
        [3, 9],
        [9, 3],
        [8, 4],
        [9, 8],
        [8, 9],
        [4, 8],
        [8, 4],
        [8, 8],
        [3, 8],
        [8, 3],
        [8, 8],
    ]
    V = [
        [3, 9],
        [3, 8],
        [9, 8],
        [8, 9],
        [8, 3],
        [9, 3],
    ]
    p = polyhedron(vrep(W), lib)
    generator_fulltest(p, V)
end

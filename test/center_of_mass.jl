using Polyhedra

"""
Test `center_of_mass` in a case where the center of mass is different from the
centroid (namely, a square pyramid).
"""
function comsquarepyramidtest(lib::Polyhedra.Library)
    poly = polyhedron(convexhull([0, 0, 0],
                                 [1, 0, 0],
                                 [0, 1, 0],
                                 [1, 1, 0],
                                 [1/2, 1/2, 1]),
                      lib)
    # Centroid is [1/2, 1/2, 1/5].
    # Center of mass is [1/2, 1/2, 1/4].
    @test center_of_mass(poly) ≈ [1/2, 1/2, 1/4]
end

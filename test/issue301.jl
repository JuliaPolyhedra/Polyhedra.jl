"""
Test for https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/301
"""
function issue301test(lib::Polyhedra.Library)
    p = polyhedron(vrep([[3, 0], [0, 3], [0, 0], [1, 1]]), lib)
    removevredundancy!(p; verbose=0, ztol=1e-7)
    @test npoints(p) == 3
    removehredundancy!(p; verbose=0, ztol=1e-7)
    @test nhalfspaces(p) == 3
end

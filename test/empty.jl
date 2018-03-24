using StaticArrays

function emptytest(lib::PolyhedraLibrary)
    vr = vrep(SVector{2, Float64}[])
    p = polyhedron(vr, lib)
    generator_fulltest(p, vr)
end

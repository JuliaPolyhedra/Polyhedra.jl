using StaticArrays

function emptytest(lib::Polyhedra.Library)
    vr = vrep(SVector{2,Float64}[])
    p = polyhedron(vr, lib)
    return generator_fulltest(p, vr)
end

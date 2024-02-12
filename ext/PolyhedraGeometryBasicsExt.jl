module PolyhedraGeometryBasicsExt

using LinearAlgebra
import GeometryBasics
using Polyhedra
using Polyhedra: FullDim, typed_fulldim, isapproxzero, _planar_hull, counterclockwise, rotate
using StaticArrays

"""
    struct Mesh{N, T, PT <: Polyhedron{T}} <: GeometryBasics.GeometryPrimitive{N, T}
        polyhedron::PT
        coordinates::Vector{GeometryBasics.Point{3, T}}
        faces::Vector{GeometryBasics.TriangleFace{Int}}
        normals::Vector{GeometryBasics.Point{3, T}}
    end

Mesh wrapper type that inherits from `GeometryPrimitive` to be used for plotting
a polyhedron. Note that `Mesh(p)` is type unstable but one can use `Mesh{3}(p)`
instead if it is known that `p` is defined in a 3-dimensional space.
"""
mutable struct Mesh{N, T, PT <: Polyhedron{T}} <: GeometryBasics.GeometryPrimitive{N, T}
    polyhedron::PT
    coordinates::Vector{GeometryBasics.Point{N, T}}
    faces::Vector{GeometryBasics.TriangleFace{Int}}
    normals::Vector{GeometryBasics.Point{3, T}}
end
function Mesh{N}(polyhedron::Polyhedron{T}) where {N, T}
    return Mesh{N, T, typeof(polyhedron)}(polyhedron, [], [], [])
end
function Mesh(polyhedron::Polyhedron, ::StaticArrays.Size{N}) where N
    return Mesh{N[1]}(polyhedron)
end
function Mesh(polyhedron::Polyhedron, N::Int)
    # This is type unstable but there is no way around that,
    # use polyhedron built from StaticArrays vector to avoid that.
    return Mesh{N}(polyhedron)
end
function Polyhedra.Mesh(polyhedron::Polyhedron)
    return Mesh(polyhedron, typed_fulldim(polyhedron))
end

function fulldecompose!(mesh::Mesh)
    if isempty(mesh.coordinates)
        mesh.coordinates, mesh.faces, mesh.normals = fulldecompose(mesh)
    end
    return
end

# Creates a scene for the vizualisation to be used to truncate the lines and rays
function scene(vr::Mesh{N}, ::Type{T}) where {N,T}
    # First compute the smallest rectangle containing the P-representation (i.e. the points).
    ps = points(vr.polyhedron)
    coord_min = ntuple(i -> minimum(Base.Fix2(getindex, i), ps), Val(N))
    coord_max = ntuple(i -> maximum(Base.Fix2(getindex, i), ps), Val(N))
    width = maximum(coord_max .- coord_min)
    if iszero(width)
        width = 2
    end
    lower = StaticArrays.SVector((coord_min .+ coord_max) ./ 2 .- width)
    width_vector = StaticArrays.SVector(ntuple(_ -> convert(T, 2width), Val(N)))
    scene = GeometryBasics.HyperRectangle{N, T}(lower, width_vector)
    # Intersection of rays with the limits of the scene
    (start, r) -> begin
        ray = coord(r)
        λ = nothing
        min_scene = minimum(scene)
        max_scene = maximum(scene)
        for i in 1:N
            r = ray[i]
            if !iszero(r)
                cur = max((min_scene[i] - start[i]) / r, (max_scene[i] - start[i]) / r)
                if λ === nothing || cur < λ
                    λ = cur
                end
            end
        end
        start + λ * ray
    end
end

function _isdup(zray, triangles)
    for tri in triangles
        normal = tri[2]
        if isapproxzero(cross(zray, normal)) && dot(zray, normal) > 0 # If A[j,:] is almost 0, it is always true...
            # parallel and equality or inequality and same sense
            return true
        end
    end
    false
end
_isdup(poly, hidx, triangles) = _isdup(get(poly, hidx).a, triangles)

function decompose_plane!(triangles::Vector, d::FullDim, zray, incident_points, incident_lines, incident_rays, exit_point::Function, counterclockwise::Function, rotate::Function)
    # xray should be the rightmost ray
    xray = nothing
    # xray should be the leftmost ray
    yray = nothing
    isapproxzero(zray) && return

    # Checking rays
    hull, lines, rays = _planar_hull(d, incident_points, incident_lines, incident_rays, counterclockwise, rotate)
    if isempty(lines)
        if length(hull) + length(rays) < 3
            return
        end
        @assert length(rays) <= 2
        if !isempty(rays)
            if length(rays) + length(hull) >= 2
                push!(hull, exit_point(last(hull), last(rays)))
            end
            push!(hull, exit_point(first(hull), first(rays)))
        end
    else
        if length(hull) == 2
            @assert length(lines) == 1 && isempty(rays)
            a, b = hull
            line = first(lines)
            empty!(hull)
            push!(hull, exit_point(a, line))
            push!(hull, exit_point(a, -line))
            push!(hull, exit_point(b, -line))
            push!(hull, exit_point(b, line))
        else
            @assert length(hull) == 1 && isempty(rays)
            center = first(hull)
            empty!(hull)
            a = first(lines)
            b = nothing
            if length(lines) == 2
                @assert isempty(rays)
                b = last(lines)
            elseif !isempty(rays)
                @assert length(lines) == 1
                @assert length(rays) == 1
                b = linearize(first(rays))
            end
            push!(hull, exit_point(center, a))
            if b !== nothing
                push!(hull, exit_point(center, b))
            end
            push!(hull, exit_point(center, -a))
            if b !== nothing && length(lines) == 2 || length(rays) >= 2
                @assert length(rays) == 2
                push!(hull, exit_point(center, -b))
            end
        end
    end

    if length(hull) >= 3
        a = pop!(hull)
        b = pop!(hull)
        while !isempty(hull)
            c = pop!(hull)
            push!(triangles, ((a, b, c), zray))
            b = c
        end
    end
end

const _Tri{N,T} = Tuple{Tuple{StaticArrays.SVector{N,T},StaticArrays.SVector{N,T},StaticArrays.SVector{N,T}},StaticArrays.SVector{3,T}}

function fulldecompose(triangles::Vector{_Tri{N,T}}) where {N,T}
    ntri = length(triangles)
    pts = Vector{GeometryBasics.Point{N, T}}(undef, 3ntri)
    faces = Vector{GeometryBasics.TriangleFace{Int}}(undef, ntri)
    ns = Vector{GeometryBasics.Point{3, T}}(undef, 3ntri)
    for i in 1:ntri
        tri = pop!(triangles)
        normal = tri[2]
        for j = 1:3
            idx = 3*(i-1)+j
            #ns[idx] = -normal
            ns[idx] = normal
        end
        faces[i] = collect(3*(i-1) .+ (1:3))
        k = 1
        for k = 1:3
            # reverse order of the 3 vertices so that if I compute the
            # normals with the `normals` function, they are in the good
            # sense.
            # I know I don't use the `normals` function but I don't know
            # what is the OpenGL convention so I don't know if it cares
            # about the order of the vertices.
            pts[3*i-k+1] = tri[1][k]
        end
    end
    # If the type of ns is Rational, it also works.
    # The normalized array in in float but then it it recast into Rational
    map!(normalize, ns, ns)
    return pts, faces, ns
end

function fulldecompose(poly_geom::Mesh, ::Type{T}) where T
    poly = poly_geom.polyhedron
    exit_point = scene(poly_geom, T)
    triangles = _Tri{2,T}[]
    z = StaticArrays.SVector(zero(T), zero(T), one(T))
    decompose_plane!(triangles, typed_fulldim(poly), z, collect(points(poly)), lines(poly), rays(poly), exit_point, counterclockwise, rotate)
    return fulldecompose(triangles)
end

function fulldecompose(poly_geom::Mesh{3}, ::Type{T}) where T
    poly = poly_geom.polyhedron
    exit_point = scene(poly_geom, T)
    triangles = _Tri{3,T}[]
    function decompose_plane(hidx)
        zray = get(poly, hidx).a
        counterclockwise(a, b) = dot(cross(a, b), zray)
        rotate(r) = cross(zray, r)
        decompose_plane!(triangles, typed_fulldim(poly), zray, incidentpoints(poly, hidx), incidentlines(poly, hidx), incidentrays(poly, hidx), exit_point, counterclockwise, rotate)
    end
    for hidx in eachindex(hyperplanes(poly))
        decompose_plane(hidx)
    end
    # If there is already a triangle, its normal is a hyperplane and it is the only face
    if isempty(triangles)
        for hidx in eachindex(halfspaces(poly))
            if !_isdup(poly, hidx, triangles)
                decompose_plane(hidx)
            end
        end
    end
    return fulldecompose(triangles)
end

fulldecompose(poly::Mesh{N, T}) where {N, T} = fulldecompose(poly, typeof(one(T)/2))

GeometryBasics.coordinates(poly::Mesh) = (fulldecompose!(poly); poly.coordinates)
GeometryBasics.faces(poly::Mesh) = (fulldecompose!(poly); poly.faces)
GeometryBasics.texturecoordinates(poly::Mesh) = nothing
GeometryBasics.normals(poly::Mesh) = (fulldecompose!(poly); poly.normals)

end

import GeometryBasics

"""
    struct Mesh{N, T, PT <: Polyhedron{T}} <: GeometryBasics.GeometryPrimitive{N, T}
        polyhedron::PT
        coordinates::Union{Nothing, Vector{GeometryBasics.Point{3, T}}}
        faces::Union{Nothing, Vector{GeometryBasics.TriangleFace{Int}}}
        normals::Union{Nothing, Vector{GeometryBasics.Point{3, T}}}
    end

Mesh wrapper type that inherits from `GeometryPrimitive` to be used for plotting
a polyhedron. Note that `Mesh(p)` is type unstable but one can use `Mesh{3}(p)`
instead if it is known that `p` is defined in a 3-dimensional space.
"""
mutable struct Mesh{N, T, PT <: Polyhedron{T}} <: GeometryBasics.GeometryPrimitive{N, T}
    polyhedron::PT
    coordinates::Union{Nothing, Vector{GeometryBasics.Point{N, T}}}
    faces::Union{Nothing, Vector{GeometryBasics.TriangleFace{Int}}}
    normals::Union{Nothing, Vector{GeometryBasics.Point{N, T}}}
end
function Mesh{N}(polyhedron::Polyhedron{T}) where {N, T}
    return Mesh{N, T, typeof(polyhedron)}(polyhedron, nothing, nothing, nothing)
end
function Mesh(polyhedron::Polyhedron, ::StaticArrays.Size{N}) where N
    return Mesh{N[1]}(polyhedron)
end
function Mesh(polyhedron::Polyhedron, N::Int)
    # This is type unstable but there is no way around that,
    # use polyhedron built from StaticArrays vector to avoid that.
    return Mesh{N}(polyhedron)
end
function Mesh(polyhedron::Polyhedron)
    return Mesh(polyhedron, FullDim(polyhedron))
end

function fulldecompose!(mesh::Mesh)
    if mesh.coordinates === nothing
        mesh.coordinates, mesh.faces, mesh.normals = fulldecompose(mesh)
    end
    return
end

# Creates a scene for the vizualisation to be used to truncate the lines and rays
function scene(vr::VRep, ::Type{T}) where T
    # First compute the smallest rectangle containing the P-representation (i.e. the points).
    (xmin, xmax) = extrema(map((x)->x[1], points(vr)))
    (ymin, ymax) = extrema(map((x)->x[2], points(vr)))
    (zmin, zmax) = extrema(map((x)->x[3], points(vr)))
    width = max(xmax-xmin, ymax-ymin, zmax-zmin)
    if width == zero(T)
        width = 2
    end
    scene = GeometryBasics.HyperRectangle{3, T}([(xmin + xmax) / 2 - width,
                                                 (ymin + ymax) / 2 - width,
                                                 (zmin + zmax) / 2 - width],
                                                2 * width * ones(T, 3))
    # Intersection of rays with the limits of the scene
    (start, r) -> begin
        ray = coord(r)
        times = max.((Vector(minimum(scene))-start) ./ ray, (Vector(maximum(scene))-start) ./ ray)
        times[ray .== 0] .= Inf # To avoid -Inf with .../(-0)
        time = minimum(times)
        start + time * ray
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

function fulldecompose(poly_geom::Mesh{3}, ::Type{T}) where T
    poly = poly_geom.polyhedron
    exit_point = scene(poly, T)

    triangles = Tuple{Tuple{Vector{T},Vector{T},Vector{T}}, Vector{T}}[]

    function decomposeplane(hidx)
        h = get(poly, hidx)
        # xray should be the rightmost ray
        xray = nothing
        # xray should be the leftmost ray
        yray = nothing
        zray = h.a
        isapproxzero(zray) && return

        # Checking rays
        counterclockwise(a, b) = dot(cross(a, b), zray)
        face_vert = pointtype(poly)[]
        for x in points(poly)
            if _isapprox(dot(x, zray), h.β)
                push!(face_vert, x)
            end
        end
        hull, lines, rays = _planar_hull(3, face_vert, incidentlines(poly, hidx), incidentrays(poly, hidx), counterclockwise, r -> cross(zray, r))
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

    for hidx in eachindex(hyperplanes(poly))
        decomposeplane(hidx)
    end
    # If there is already a triangle, his normal is an hyperplane and it is the only face
    if isempty(triangles)
        for hidx in eachindex(halfspaces(poly))
            if !_isdup(poly, hidx, triangles)
                decomposeplane(hidx)
            end
        end
    end

    ntri = length(triangles)
    pts = Vector{GeometryBasics.Point{3, T}}(undef, 3ntri)
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
    (pts, faces, ns)
end

fulldecompose(poly::Mesh{N, T}) where {N, T} = fulldecompose(poly, typeof(one(T)/2))

GeometryBasics.coordinates(poly::Mesh) = (fulldecompose!(poly); poly.coordinates)
GeometryBasics.faces(poly::Mesh) = (fulldecompose!(poly); poly.faces)
GeometryBasics.texturecoordinates(poly::Mesh) = nothing
GeometryBasics.normals(poly::Mesh) = (fulldecompose!(poly); poly.normals)

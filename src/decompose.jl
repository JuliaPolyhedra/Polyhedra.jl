import GeometryTypes

"""
    struct Mesh{N, T, PT <: Polyhedron{T}} <: GeometryTypes.GeometryPrimitive{N, T}
        polyhedron::PT
    end

Mesh wrapper type that inherits from `GeometryPrimitive` to be used for plotting
a polyhedron. Note that `Mesh(p)` is type unstable but one can use `Mesh{3}(p)`
instead if it is known that `p` is defined in a 3-dimensional space.
"""
struct Mesh{N, T, PT <: Polyhedron{T}} <: GeometryTypes.GeometryPrimitive{N, T}
    polyhedron::PT
end
function Mesh{N}(polyhedron::Polyhedron{T}) where {N, T}
    return Mesh{N, T, typeof(polyhedron)}(polyhedron)
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
    scene = GeometryTypes.HyperRectangle{3, T}([(xmin + xmax) / 2 - width,
                                                (ymin + ymax) / 2 - width,
                                                (zmin + zmax) / 2 - width],
                                               2 * width * ones(T, 3))
    # Intersection of rays with the limits of the scene
    (start, ray) -> begin
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
        line = nothing
        lineleft = false
        lineright = false
        function checkleftright(r::Union{Ray, Line})
            cc = counterclockwise(r, line)
            if !isapproxzero(cc)
                if cc < 0 || islin(r)
                    lineleft = true
                end
                if cc > 0 || islin(r)
                    lineright = true
                end
            end
        end
        for l in incidentlines(poly, hidx)
            if !isapproxzero(l)
                if line === nothing
                    line = l
                else
                    checkleftright(l)
                end
            end
        end
        for r in incidentrays(poly, hidx)
            if !isapproxzero(r)
                if line === nothing
                    if xray === nothing || counterclockwise(r, xray) > 0
                        xray = coord(r) # r is more right than xray
                    end
                    if yray === nothing || counterclockwise(r, yray) < 0
                        yray = coord(r) # r is more left than xray
                    end
                else
                    checkleftright(r)
                end
            end
        end

        # Checking vertices
        face_vert = pointtype(poly)[]
        for x in points(poly)
            if _isapprox(dot(x, zray), h.β)
                push!(face_vert, x)
            end
        end

        if line !== nothing
            if isempty(face_vert)
                center = origin(pointtype(poly), 3)
            else
                center = first(face_vert)
            end
            hull = pointtype(poly)[]
            push!(hull, exit_point(center, line))
            if lineleft
                push!(hull, exit_point(center, cross(zray, line)))
            end
            push!(hull, exit_point(center, -line))
            if lineright
                push!(hull, exit_point(center, cross(line, zray)))
            end
            hulls = (hull,)
        else
            #if length(face_vert) < 3 # Wrong, they are also the rays
            #  error("Not enough vertices and rays to form a face, it may be because of numerical rounding. Otherwise, please report this bug.")
            #end
            if length(face_vert) < 3 && (xray == nothing || (length(face_vert) < 2 && (yray == xray || length(face_vert) < 1)))
                return
            end
            if xray == nothing
                sweep_norm = cross(zray, [1,0,0])
                if sum(abs, sweep_norm) == 0
                    sweep_norm = cross(zray, [0,1,0])
                end
            else
                sweep_norm = cross(zray, xray)
            end
            sort!(face_vert, by = x -> dot(x, sweep_norm))
            xtoy_hull = getsemihull(face_vert, 1, counterclockwise, yray)
            if yray == nothing
                ytox_hull = getsemihull(face_vert, -1, counterclockwise, yray)
            else
                ytox_hull = Any[]
                push!(ytox_hull, face_vert[1])
                if last(xtoy_hull) != face_vert[1]
                    push!(ytox_hull, last(xtoy_hull))
                end
                push!(ytox_hull, exit_point(last(xtoy_hull), yray))
                push!(ytox_hull, exit_point(face_vert[1], xray))
            end
            hulls = (xtoy_hull, ytox_hull)
        end
        for hull in hulls
            if length(hull) >= 3
                a = pop!(hull)
                b = pop!(hull)
                while !isempty(hull)
                    c = pop!(hull)
                    push!(triangles, ((a,b,c), zray))
                    b = c
                end
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
    pts = Vector{GeometryTypes.Point{3, T}}(undef, 3ntri)
    faces = Vector{GeometryTypes.Face{3, Int}}(undef, ntri)
    ns = Vector{GeometryTypes.Normal{3, T}}(undef, 3ntri)
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

GeometryTypes.isdecomposable(::Type{<:GeometryTypes.Point{3}},  ::Type{<:Mesh{3}}) = true
GeometryTypes.isdecomposable(::Type{<:GeometryTypes.Face{3}},   ::Type{<:Mesh{3}}) = true
GeometryTypes.isdecomposable(::Type{<:GeometryTypes.Normal{3}}, ::Type{<:Mesh{3}}) = true
function GeometryTypes.decompose(PT::Type{<:GeometryTypes.Point}, poly::Mesh)
    points = fulldecompose(poly)[1]
    map(PT, points)
end
function GeometryTypes.decompose(FT::Type{<:GeometryTypes.Face}, poly::Mesh)
    faces = fulldecompose(poly)[2]
    GeometryTypes.decompose(FT, faces)
end
function GeometryTypes.decompose(NT::Type{<:GeometryTypes.Normal}, poly::Mesh)
    ns = fulldecompose(poly)[3]
    map(NT, ns)
end

# In AbstractPlotting, when asking to plot an object, it calls this constructor
# which is only defined for `GeometryTypes.GeometryPrimitive` which is a
# supertype of `Polyhedra.Mesh`
GeometryTypes.GLNormalMesh(p::Polyhedron) = GeometryTypes.GLNormalMesh(Mesh(p))

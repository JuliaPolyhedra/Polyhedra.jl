# I only import it and do not use "using" so that Datastructures.status does not conflict with MathProgBase.status
import DataStructures
import GeometryTypes.decompose, GeometryTypes.isdecomposable

function fulldecompose{T}(poly::Polyhedron{3,T})
    ine = SimpleHRepresentation(poly)

    # I need to do division so if T is e.g. Integer, I need to use another type
    RT = typeof(one(T)/2)

    A = ine.A
    #rayinface{T<:Real}(r::Vector{T}, i::Integer) = myeqzero(dot(r, A[i,:])) && !myeqzero(r)
    #vertinface{T<:Real}(r::Vector{T}, i::Integer) = myeqzero(dot(r, A[i,:])) && !myeqzero(r)

    ps = haspoints(poly) ? points(poly) : [zeros(T, 3)]

    # Intersection of rays with the limits of the scene
    (xmin, xmax) = extrema(map((x)->x[1], ps))
    (ymin, ymax) = extrema(map((x)->x[2], ps))
    (zmin, zmax) = extrema(map((x)->x[3], ps))
    width = max(xmax-xmin, ymax-ymin, zmax-zmin)
    if width == zero(T)
        width = 2
    end
    scene = HyperRectangle{3,RT}([(xmin+xmax)/2-width, (ymin+ymax)/2-width, (zmin+zmax)/2-width], 2*width*ones(RT,3))
    function exit_point(start, ray)
        times = max((Vector(minimum(scene))-start) ./ ray, (Vector(maximum(scene))-start) ./ ray)
        times[ray .== 0] = Inf # To avoid -Inf with .../(-0)
        time = minimum(times)
        start + time * ray
    end

    triangles = DataStructures.Stack(Tuple{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}},Int64})
    for i in 1:size(A, 1)
        xray = nothing
        yray = nothing
        zray = A[i,:]
        if myeqzero(zray)
            continue
        end
        newface = true
        for j in 1:i-1
            if myeqzero(cross(zray, A[j,:])) && (i in ine.linset || dot(zray, A[j,:]) > 0) # If A[j,:] is almost 0, it is always true...
                # parallel and equality or inequality and same sense
                # TODO is it possible that A[i,:] is stronger than A[j,:] ?
                newface = false
            end
        end
        if !newface
            continue
        end
        if i > 1 && A[i,:] == A[i-1,:] # Same row, only need to check i-1 since the rows are sorted
            continue
        end

        # Checking rays
        counterclockwise(a, b) = dot(cross(a, b), zray)
        line = nothing
        lineleft = false
        lineright = false
        function checkleftright(r::Vector, lin::Bool)
            cc = counterclockwise(r, line)
            if !myeqzero(cc)
                if cc < 0 || lin
                    lineleft = true
                end
                if cc > 0 || lin
                    lineright = true
                end
            end
        end
        for r in rays(poly)
            if myeqzero(dot(r, zray)) && !myeqzero(r)
                if line != nothing
                    checkleftright(r, islin(r))
                else
                    if islin(r)
                        line = vec(r)
                        if xray != nothing
                            checkleftright(xray, false) # false otherwise line wouldn't be nothing
                        end
                        if yray != nothing
                            checkleftright(yray, false)
                        end
                    end
                    if xray == nothing || counterclockwise(r, xray) > 0
                        xray = vec(r)
                    end
                    if yray == nothing || counterclockwise(r, yray)  < 0
                        yray = vec(r)
                    end
                end
            end
        end

        # Checking vertices
        face_vert = []
        for x in ps
            if myeq(dot(x, zray), ine.b[i])
                push!(face_vert, x)
            end
        end

        if line != nothing
            if isempty(face_vert)
                center = zeros(T, 3)
            else
                center = vec(first(face_vert))
            end
            hull = DataStructures.Stack(Any)
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
                continue
            end
            if xray == nothing
                sweep_norm = cross(zray, [1,0,0])
                if sum(abs(sweep_norm)) == 0
                    sweep_norm = cross(zray, [0,1,0])
                end
            else
                sweep_norm = cross(zray, xray)
            end
            sort!(face_vert, by = x -> dot(x, sweep_norm))
            function getsemihull(sign_sense)
                hull = DataStructures.Stack(Vector{T})
                prev = sign_sense == 1 ? face_vert[1] : face_vert[length(face_vert)]
                cur = prev
                for j in (sign_sense == 1 ? (2:length(face_vert)) : ((length(face_vert)-1):-1:1))
                    while prev != cur && counterclockwise(cur - prev, face_vert[j] - prev) >= 0
                        cur = prev
                        pop!(hull)
                        if !isempty(hull)
                            prev = DataStructures.top(hull)
                        end
                    end
                    if yray != nothing && counterclockwise(face_vert[j] - cur, yray) >= 0
                        break
                    else
                        push!(hull, cur)
                        prev = cur
                        cur = face_vert[j]
                    end
                end
                push!(hull, cur)
                hull
            end
            xtoy_hull = getsemihull(1)
            if yray == nothing
                ytox_hull = getsemihull(-1)
            else
                ytox_hull = DataStructures.Stack(Any)
                push!(ytox_hull, face_vert[1])
                if DataStructures.top(xtoy_hull) != face_vert[1]
                    push!(ytox_hull, DataStructures.top(xtoy_hull))
                end
                push!(ytox_hull, exit_point(DataStructures.top(xtoy_hull), yray))
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
                    push!(triangles, ((a,b,c), i))
                    b = c
                end
            end
        end

    end
    ntri = length(triangles)
    pts  = Vector{GeometryTypes.Point{3,RT}}(3*ntri)
    faces   = Vector{GeometryTypes.Face{3,Int,0}}(ntri)
    ns = Vector{GeometryTypes.Normal{3,RT}}(3*ntri)
    for i in 1:ntri
        tri = pop!(triangles)
        normal = vec(A[tri[2],:])
        for j = 1:3
            idx = 3*(i-1)+j
            #ns[idx] = -normal
            ns[idx] = normal
        end
        faces[i] = Array(3*(i-1)+(1:3))
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
    map!(normalize, ns)
    (pts, faces, ns)
end


isdecomposable{T<:Point, S<:Polyhedron}(::Type{T}, ::Type{S}) = true
isdecomposable{T<:Face, S<:Polyhedron}(::Type{T}, ::Type{S}) = true
isdecomposable{T<:Normal, S<:Polyhedron}(::Type{T}, ::Type{S}) = true
function decompose{N, T1, T2}(PT::Type{Point{N, T1}}, poly::Polyhedron{N, T2})
    points = fulldecompose(poly)[1]
    map(PT, points)
end
function decompose{N, T, T2}(FT::Type{Face{N, T}}, poly::Polyhedron{3, T2})
    faces = fulldecompose(poly)[2]
    decompose(FT, faces)
end
function decompose{NT<:Normal, T}(::Type{NT}, poly::Polyhedron{3,T})
    ns = fulldecompose(poly)[3]
    map(NT, ns)
end

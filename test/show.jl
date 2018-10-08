# Taken from Base: test/show.jl
replstr(x, lim=true) = sprint((io,x) -> show(IOContext(io, :limit => lim, :displaysize => (24, 80)), MIME("text/plain"), x), x)
showstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), x), x)

macro test_show(expr)
    esc(quote
            @test string($expr) == "$($expr)"
    end)
end

@testset "Show" begin
    @testset "Simple show" begin
        @test_show convexhull([1, 0]) + convexhull(Line([1, 1]))
        @test_show HyperPlane([1, 0], 0) ∩ HyperPlane([0, 1], 0) ∩ HalfSpace([1, 1], 1)
    end
    vshort = vrep([[i, -i] for i in 1:3])
    vlong = vrep([[i, -i] for i in 1:30])
    @testset "Iterator" begin
        @test replstr(points(vshort))        == "3-element iterator of Array{$Int,1}:\n [1, -1]\n [2, -2]\n [3, -3]"
        @test replstr(points(vshort), false) == "3-element iterator of Array{$Int,1}:\n [1, -1]\n [2, -2]\n [3, -3]"
        @test replstr(points(vlong)) == "30-element iterator of Array{$Int,1}:\n [1, -1]\n [2, -2]\n [3, -3]\n [4, -4]\n [5, -5]\n [6, -6]\n [7, -7]\n [8, -8]\n [9, -9]\n [10, -10]\n [11, -11]\n [12, -12]\n [13, -13]\n [14, -14]\n [15, -15]\n [16, -16]\n [17, -17]\n [18, -18]\n [19, -19]\n [20, -20]\n  ⋮"
        @test replstr(points(vlong), false) == "30-element iterator of Array{$Int,1}:\n [1, -1]\n [2, -2]\n [3, -3]\n [4, -4]\n [5, -5]\n [6, -6]\n [7, -7]\n [8, -8]\n [9, -9]\n [10, -10]\n [11, -11]\n [12, -12]\n [13, -13]\n [14, -14]\n [15, -15]\n [16, -16]\n [17, -17]\n [18, -18]\n [19, -19]\n [20, -20]\n [21, -21]\n [22, -22]\n [23, -23]\n [24, -24]\n [25, -25]\n [26, -26]\n [27, -27]\n [28, -28]\n [29, -29]\n [30, -30]"
        @test showstr(points(vshort)) == "Array{$Int,1}[[1, -1], [2, -2], [3, -3]]"
        @test showstr(points(vlong)) == "Array{$Int,1}[[1, -1], [2, -2], [3, -3], [4, -4], [5, -5], [6, -6], [7, -7], [8, -8], [9, -9], [10, -10], [11, -11], [12, -12], [13, -13], [14, -14], [15, -15], [16, -16], [17, -17], [18, -18], [19, -19], [20, -20] … ]"
    end
    v = vshort + Ray([1, 0])
    p = polyhedron(v)
    @testset "Polyhedron without H-representation" begin
        @test replstr(p)        == "Polyhedron DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},$Int},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},$Int}}:\n3-element iterator of Array{Rational{BigInt},1}:\n Rational{BigInt}[1//1, -1//1]\n Rational{BigInt}[2//1, -2//1]\n Rational{BigInt}[3//1, -3//1],\n1-element iterator of Ray{Rational{BigInt},Array{Rational{BigInt},1}}:\n Ray(Rational{BigInt}[1//1, 0//1])"
        @test replstr(p, false) == "Polyhedron DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},$Int},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},$Int}}:\n3-element iterator of Array{Rational{BigInt},1}:\n Rational{BigInt}[1//1, -1//1]\n Rational{BigInt}[2//1, -2//1]\n Rational{BigInt}[3//1, -3//1],\n1-element iterator of Ray{Rational{BigInt},Array{Rational{BigInt},1}}:\n Ray(Rational{BigInt}[1//1, 0//1])"
        @test showstr(p) == "convexhull([1//1, -1//1], [2//1, -2//1], [3//1, -3//1]) + convexhull(Ray(Rational{BigInt}[1//1, 0//1]))"
    end
    h = hrep(p)
    @testset "V-Representation" begin
        @test replstr(v)        == "V-representation Polyhedra.Hull{$Int,Array{$Int,1},$Int}:\n3-element iterator of Array{$Int,1}:\n [1, -1]\n [2, -2]\n [3, -3],\n1-element iterator of Ray{$Int,Array{$Int,1}}:\n Ray([1, 0])"
        @test replstr(v, false) == "V-representation Polyhedra.Hull{$Int,Array{$Int,1},$Int}:\n3-element iterator of Array{$Int,1}:\n [1, -1]\n [2, -2]\n [3, -3],\n1-element iterator of Ray{$Int,Array{$Int,1}}:\n Ray([1, 0])"
        @test showstr(v) == "convexhull([1, -1], [2, -2], [3, -3]) + convexhull(Ray([1, 0]))"
    end
    @testset "H-Representation" begin
        @test replstr(h)        == "H-representation Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},$Int}:\n3-element iterator of HalfSpace{Rational{BigInt},Array{Rational{BigInt},1}}:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
        @test replstr(h, false) == "H-representation Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},$Int}:\n3-element iterator of HalfSpace{Rational{BigInt},Array{Rational{BigInt},1}}:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
        @test showstr(h) == "HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1) ∩ HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4) ∩ HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
    end
    @testset "Polyhedron" begin
        @test replstr(p)        == "Polyhedron DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},$Int},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},$Int}}:\n3-element iterator of HalfSpace{Rational{BigInt},Array{Rational{BigInt},1}}:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2):\n3-element iterator of Array{Rational{BigInt},1}:\n Rational{BigInt}[1//1, -1//1]\n Rational{BigInt}[2//1, -2//1]\n Rational{BigInt}[3//1, -3//1],\n1-element iterator of Ray{Rational{BigInt},Array{Rational{BigInt},1}}:\n Ray(Rational{BigInt}[1//1, 0//1])"
        @test replstr(p, false) == "Polyhedron DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},$Int},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},$Int}}:\n3-element iterator of HalfSpace{Rational{BigInt},Array{Rational{BigInt},1}}:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2):\n3-element iterator of Array{Rational{BigInt},1}:\n Rational{BigInt}[1//1, -1//1]\n Rational{BigInt}[2//1, -2//1]\n Rational{BigInt}[3//1, -3//1],\n1-element iterator of Ray{Rational{BigInt},Array{Rational{BigInt},1}}:\n Ray(Rational{BigInt}[1//1, 0//1])"
        @test showstr(p) == "HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1) ∩ HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4) ∩ HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2) : convexhull([1//1, -1//1], [2//1, -2//1], [3//1, -3//1]) + convexhull(Ray(Rational{BigInt}[1//1, 0//1]))"
    end
    p = polyhedron(h)
    @testset "Polyhedron without V-representation" begin
        @test replstr(p)        == "Polyhedron DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},$Int},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},$Int}}:\n3-element iterator of HalfSpace{Rational{BigInt},Array{Rational{BigInt},1}}:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
        @test replstr(p, false) == "Polyhedron DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},$Int},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},$Int}}:\n3-element iterator of HalfSpace{Rational{BigInt},Array{Rational{BigInt},1}}:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
        @test showstr(p) == "HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1) ∩ HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4) ∩ HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
    end
end

using Test, Polyhedra

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
    PT = Polyhedra.pointtype(vshort)
    RT = Polyhedra.raytype(vshort)
    @testset "Iterator" begin
        @test replstr(points(vshort))        == "3-element iterator of $PT:\n [1, -1]\n [2, -2]\n [3, -3]"
        @test replstr(points(vshort), false) == "3-element iterator of $PT:\n [1, -1]\n [2, -2]\n [3, -3]"
        @test replstr(points(vlong)) == "30-element iterator of $PT:\n [1, -1]\n [2, -2]\n [3, -3]\n [4, -4]\n [5, -5]\n [6, -6]\n [7, -7]\n [8, -8]\n [9, -9]\n [10, -10]\n [11, -11]\n [12, -12]\n [13, -13]\n [14, -14]\n [15, -15]\n [16, -16]\n [17, -17]\n [18, -18]\n [19, -19]\n [20, -20]\n  ⋮"
        @test replstr(points(vlong), false) == "30-element iterator of $PT:\n [1, -1]\n [2, -2]\n [3, -3]\n [4, -4]\n [5, -5]\n [6, -6]\n [7, -7]\n [8, -8]\n [9, -9]\n [10, -10]\n [11, -11]\n [12, -12]\n [13, -13]\n [14, -14]\n [15, -15]\n [16, -16]\n [17, -17]\n [18, -18]\n [19, -19]\n [20, -20]\n [21, -21]\n [22, -22]\n [23, -23]\n [24, -24]\n [25, -25]\n [26, -26]\n [27, -27]\n [28, -28]\n [29, -29]\n [30, -30]"
        @test showstr(points(vshort)) == string([[1, -1], [2, -2], [3, -3]])
        @test showstr(points(vlong)) == string([[i, -i] for i in 1:20])[1:end-1] * " … ]"
    end
    v = vshort + Ray([1, 0])
    p = polyhedron(v)
    PPT = Polyhedra.pointtype(p)
    PRT = Polyhedra.raytype(p)
    @testset "Polyhedron without H-representation" begin
        @test replstr(p)        == "Polyhedron $(typeof(p)):\n3-element iterator of $PPT:\n Rational{BigInt}[1//1, -1//1]\n Rational{BigInt}[2//1, -2//1]\n Rational{BigInt}[3//1, -3//1],\n1-element iterator of $PRT:\n Ray(Rational{BigInt}[1//1, 0//1])"
        @test replstr(p, false) == "Polyhedron $(typeof(p)):\n3-element iterator of $PPT:\n Rational{BigInt}[1//1, -1//1]\n Rational{BigInt}[2//1, -2//1]\n Rational{BigInt}[3//1, -3//1],\n1-element iterator of $PRT:\n Ray(Rational{BigInt}[1//1, 0//1])"
        @test showstr(p) == "convexhull([1//1, -1//1], [2//1, -2//1], [3//1, -3//1]) + convexhull(Ray(Rational{BigInt}[1//1, 0//1]))"
    end
    h = hrep(p)
    HS = Polyhedra.halfspacetype(p)
    @testset "V-Representation" begin
        @test replstr(v)        == "V-representation $(typeof(v)):\n3-element iterator of $PT:\n [1, -1]\n [2, -2]\n [3, -3],\n1-element iterator of $RT:\n Ray([1, 0])"
        @test replstr(v, false) == "V-representation $(typeof(v)):\n3-element iterator of $PT:\n [1, -1]\n [2, -2]\n [3, -3],\n1-element iterator of $RT:\n Ray([1, 0])"
        @test showstr(v) == "convexhull([1, -1], [2, -2], [3, -3]) + convexhull(Ray([1, 0]))"
    end
    @testset "H-Representation" begin
        @test replstr(h)        == "H-representation $(typeof(h)):\n3-element iterator of $HS:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
        @test replstr(h, false) == "H-representation $(typeof(h)):\n3-element iterator of $HS:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
        @test showstr(h) == "HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1) ∩ HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4) ∩ HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
    end
    @testset "Polyhedron" begin
        @test replstr(p)        == "Polyhedron $(typeof(p)):\n3-element iterator of $HS:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2):\n3-element iterator of $PPT:\n Rational{BigInt}[1//1, -1//1]\n Rational{BigInt}[2//1, -2//1]\n Rational{BigInt}[3//1, -3//1],\n1-element iterator of $PRT:\n Ray(Rational{BigInt}[1//1, 0//1])"
        @test replstr(p, false) == "Polyhedron $(typeof(p)):\n3-element iterator of $HS:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2):\n3-element iterator of $PPT:\n Rational{BigInt}[1//1, -1//1]\n Rational{BigInt}[2//1, -2//1]\n Rational{BigInt}[3//1, -3//1],\n1-element iterator of $PRT:\n Ray(Rational{BigInt}[1//1, 0//1])"
        @test showstr(p) == "HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1) ∩ HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4) ∩ HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2) : convexhull([1//1, -1//1], [2//1, -2//1], [3//1, -3//1]) + convexhull(Ray(Rational{BigInt}[1//1, 0//1]))"
    end
    p = polyhedron(h)
    @testset "Polyhedron without V-representation" begin
        @test replstr(p)        == "Polyhedron $(typeof(p)):\n3-element iterator of $HS:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
        @test replstr(p, false) == "Polyhedron $(typeof(p)):\n3-element iterator of $HS:\n HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4)\n HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
        @test showstr(p) == "HalfSpace(Rational{BigInt}[-1//2, -1//2], 0//1) ∩ HalfSpace(Rational{BigInt}[0//1, -1//4], 3//4) ∩ HalfSpace(Rational{BigInt}[0//1, 1//2], -1//2)"
    end
end

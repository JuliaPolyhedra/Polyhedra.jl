# Taken from Base: test/show.jl
replstr(x, lim=true) = sprint((io,x) -> show(IOContext(io, :limit => lim, :displaysize => (24, 80)), MIME("text/plain"), x), x)
showstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), x), x)

@testset "Show" begin
    vshort = vrep([[i, -i] for i in 1:3])
    vlong = vrep([[i, -i] for i in 1:30])
    @testset "Iterator" begin
        @test replstr(points(vshort)) == "3-element iterator of Array{Int64,1}:\n [1, -1]\n [2, -2]\n [3, -3]"
        @test replstr(points(vshort), false) == "3-element iterator of Array{Int64,1}:\n [1, -1]\n [2, -2]\n [3, -3]"
        @test replstr(points(vlong)) == "30-element iterator of Array{Int64,1}:\n [1, -1]\n [2, -2]\n [3, -3]\n [4, -4]\n [5, -5]\n [6, -6]\n [7, -7]\n [8, -8]\n [9, -9]\n [10, -10]\n [11, -11]\n [12, -12]\n [13, -13]\n [14, -14]\n [15, -15]\n [16, -16]\n [17, -17]\n [18, -18]\n [19, -19]\n [20, -20]\n  ⋮"
        @test replstr(points(vlong), false) == "30-element iterator of Array{Int64,1}:\n [1, -1]\n [2, -2]\n [3, -3]\n [4, -4]\n [5, -5]\n [6, -6]\n [7, -7]\n [8, -8]\n [9, -9]\n [10, -10]\n [11, -11]\n [12, -12]\n [13, -13]\n [14, -14]\n [15, -15]\n [16, -16]\n [17, -17]\n [18, -18]\n [19, -19]\n [20, -20]\n [21, -21]\n [22, -22]\n [23, -23]\n [24, -24]\n [25, -25]\n [26, -26]\n [27, -27]\n [28, -28]\n [29, -29]\n [30, -30]"
        @test showstr(points(vshort)) == "Array{Int64,1}[[1, -1], [2, -2], [3, -3]]"
        @test showstr(points(vlong)) == "Array{Int64,1}[[1, -1], [2, -2], [3, -3], [4, -4], [5, -5], [6, -6], [7, -7], [8, -8], [9, -9], [10, -10], [11, -11], [12, -12], [13, -13], [14, -14], [15, -15], [16, -16], [17, -17], [18, -18], [19, -19], [20, -20] … ]"
    end
end

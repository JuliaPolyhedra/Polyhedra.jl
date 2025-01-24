using Test

import MutableArithmetics
const MA = MutableArithmetics

using Polyhedra

function translate_test(T::Type)
    o = one(T)
    v = [o, -o]
    # We need to create independent copies of `o` in case these are `BigInt`.
    pr = convexhull([1o, 2o])
    @test MA.promote_operation(translate, typeof(pr), typeof(v)) == typeof(pr)
    @test MA.mutability(typeof(pr)) isa MA.IsMutable
    @test MA.mutability(typeof(pr), translate, typeof(pr), typeof(v)) isa MA.IsMutable
    rr = conichull([3o, -o])
    @test MA.promote_operation(translate, typeof(rr), typeof(v)) == typeof(rr)
    @test MA.mutability(typeof(rr)) isa MA.IsMutable
    @test MA.mutability(typeof(rr), translate, typeof(rr), typeof(v)) isa MA.IsMutable
    vr = pr + rr
    @test MA.promote_operation(translate, typeof(vr), typeof(v)) == typeof(vr)
    alloc_test(() -> MA.promote_operation(translate, typeof(vr), typeof(v)), 0)
    @test MA.mutability(typeof(vr)) isa MA.IsMutable
    @test MA.mutability(typeof(vr), translate, typeof(vr), typeof(v)) isa MA.IsMutable
    alloc_test(() -> MA.mutability(typeof(vr), translate, typeof(vr), typeof(v)), 0)
    vr2 = MA.operate!(translate, vr, v)
    rec_identical(vr, vr2)
    alloc_test(() -> MA.mutable_operate!(+, a, b), 0)
    alloc_test(() -> MA.mutable_operate!(translate, a, b), 0)
    alloc_test(() -> MA.mutable_operate!(translate, vr.points.points[1], v), 0)
    alloc_test(() -> MA.mutable_operate!(translate, vr.points, v), 0)
    alloc_test(() -> MA.mutable_operate!(translate, vr.rays, v), 0)
    alloc_test(() -> MA.mutable_operate!(translate, vr, v), 0)
    alloc_test(() -> MA.operate!(translate, vr, v), 0)

#    hr = HalfSpace([1o, 1], 2o) ∩ HyperPlane([-o, 2o, 0o])
#    @test MA.promote_operation(translate, typeof(hr), typeof(v)) == typeof(vr)
#    @test MA.mutability(translate, typeof(hr), typeof(v)) isa MA.IsMutable
end

@testset "Translate" begin
    translate_test(Float64)
end

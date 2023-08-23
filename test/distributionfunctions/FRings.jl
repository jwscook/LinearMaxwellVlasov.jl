using Dates
println("$(now()) $(@__FILE__)")

using Test, QuadGK
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


@testset "FRings" begin
  for _ in 1:10
    mi = 1836*LMV.mₑ
    vth = thermalspeed(1 + rand(), mi)
    vd = thermalspeed(10 + 10 * rand(), mi)

    ring = FRing(vth, vd)
    @test QuadGK.quadgk(v->2π*v*ring(v), 0, 12 * (vth + vd))[1] ≈ 1
    @test LMV.is_normalised(ring)
    deriv = rand(Bool)
    kernel = LMV.PerpendicularKernel(rand() / vth,
      Pair(rand(-3:3), rand(-3:3)), rand(0:5))
    result = LMV.integrate(ring, kernel, deriv)
    expected = QuadGK.quadgk(v->kernel(v)*ring(v, deriv),
      0, Inf, rtol=10eps())[1]
    @test result ≈ expected
    result0 = LMV.integrate(FRing(vth, 0), kernel, deriv)
    expected0 = LMV.integrate(FPerpendicularMaxwellian(vth), kernel, deriv)
    @test result0 ≈ expected0
  end
end

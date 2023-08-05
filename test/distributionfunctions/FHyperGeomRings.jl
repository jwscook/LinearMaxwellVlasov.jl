using Dates
println("$(now()) $(@__FILE__)")

using Test, QuadGK
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


@testset "FHyperGeomRings" begin
  for _ in 1:10
    vth = 1 + rand()
    vd = vth / 10 #10 + 10 * rand()

    normalring = LMV.FRing(vth, vd)
    ring = LMV.FHyperGeomRing(vth, vd)

    @test QuadGK.quadgk(v->2π*v*ring(v), 0, 12 * (vth + vd))[1] ≈ 1
    @test LMV.is_normalised(ring)
    deriv = rand(Bool)
    kernel = LMV.PerpendicularKernel(rand() * vth,
      Pair(rand(-3:3), rand(-3:3)), rand(0:5))
    result = LMV.integrate(ring, kernel, deriv)
    resultnormal = LMV.integrate(normalring, kernel, deriv)
    expected = QuadGK.quadgk(v->kernel(v)*ring(v, deriv), 0, Inf, rtol=10eps())[1]
    expectednormal = QuadGK.quadgk(v->kernel(v)*normalring(v, deriv), 0, Inf, rtol=10eps())[1]
    @test result ≈ expected
    @show resultnormal, result
    @show expectednormal, expected
  end
 
end

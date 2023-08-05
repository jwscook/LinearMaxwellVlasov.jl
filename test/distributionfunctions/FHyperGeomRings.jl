using Dates
println("$(now()) $(@__FILE__)")

using Test, QuadGK
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


@testset "FHyperGeomRings" begin
  for z in 0:2:26

  @testset "Ratio $z" begin
    for _ in 1:10
    vth = rand()
    vd = z * vth# + 2 * rand()

    normalring = LMV.FRing(vth, vd)
    ring = LMV.FHyperGeomRing(vth, vd)

    @test QuadGK.quadgk(v->2π*v*ring(v), 0, 12 * (vth + vd))[1] ≈ 1
    @test LMV.is_normalised(ring)
    deriv = rand(Bool)
    n = rand(-10:10)
    m = n + rand(-1:1)
    kernel = LMV.PerpendicularKernel(rand() * vth,
      Pair(n, m), rand(0:3))
    result = LMV.integrate(ring, kernel, deriv)
    expected = QuadGK.quadgk(v->kernel(v)*ring(v, deriv), 0, Inf, rtol=eps())[1]
    @test result ≈ expected
    #resultnormal = LMV.integrate(normalring, kernel, deriv)
    #expectednormal = QuadGK.quadgk(v->kernel(v)*normalring(v, deriv), 0, Inf, rtol=10eps())[1]
    #@show resultnormal, result
    #@show expectednormal, expected
  end
  end
  end
 
end

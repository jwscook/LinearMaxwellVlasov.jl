using Dates
println("$(now()) $(@__FILE__)")

using Test, QuadGK
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


@testset "FHyperGeomRings" begin
  for z in 0:20#:200

  vth = 1.0#rand()
  for n in 1:5:50
    vd = z * vth# + 2 * rand()

    normalring = LMV.FRing(vth, vd)
    ring = LMV.FHyperGeomRing(vth, vd)

    @test QuadGK.quadgk(v->2π*v*ring(v), 0, 12 * (vth + vd))[1] ≈ 1
    @test LMV.is_normalised(ring)
  for deriv in (true, false)
  for m in n-1:n+1
  @testset "Ratio $z, bessel index $n, $m, and deriv $deriv" begin
    kernel = LMV.PerpendicularKernel(randn() * vth,
      Pair(n, m), rand(0:3))
    result = LMV.integrate(ring, kernel, deriv)
    lo = max(0.0, vd - 12 * vth)
    hi = min(Inf, vd + 12 * vth)
    expected = QuadGK.quadgk(v->kernel(v)*ring(v, deriv), lo, hi,
                             rtol=2eps(), order=57)[1]
    @test result ≈ expected rtol=1e-4
    #if !isapprox(result, expected, rtol=1e-4)
    #  resultnormal = LMV.integrate(normalring, kernel, deriv)
    #  @show result, resultnormal
    #end
    #expectednormal = QuadGK.quadgk(v->kernel(v)*normalring(v, deriv), 0, Inf, rtol=10eps())[1]
    #@show resultnormal, result
    #@show expectednormal, expected
  end
  end
  end
  end
  end
 
end

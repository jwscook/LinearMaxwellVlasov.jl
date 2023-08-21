using Dates
println("$(now()) $(@__FILE__)")

using Test, QuadGK, Random
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

Random.seed!(0)

@testset "FWideRings" begin
  # vth/vd ≈ sqrt(3.6e6) / sqrt(1e4)
  # vth ≈ vd sqrt(3.6e2)
  # vd ≈ vth / sqrt(3.6e2) = vth / 20
  for z in vcat(1:5, 6:2:10, 15:5:20)#, 30:10:50)
    vth = 1.0
    vd = z * vth# + 2 * rand()
    #normalring = LMV.FRing(vth, vd) # similar to the wide ring
    ring = LMV.FWideRing(vth, vd)
    #for n in (-10, -5, 0, 5, 10, 40)
    n = rand((-50:50))
    m = rand(-1:1) + n #for m in (n+1,)#n-1:n+1
    @test QuadGK.quadgk(v->2π*v*ring(v), 0, 12 * (vth + vd))[1] ≈ 1
    @test LMV.is_normalised(ring)
    for deriv in (true, false)
      @testset "Ratio $z, bessel index $n, $m, and deriv $deriv" begin
        kernel = LMV.PerpendicularKernel(randn() * vth,
          Pair(n, m), rand(0:3))
        t1 = @elapsed result = LMV.integrate(ring, kernel, deriv)
        lo = max(0.0, vd - 32 * vth)
        hi = min(Inf, vd + 32 * vth)
        expected = QuadGK.quadgk(v->kernel(v)*ring(v, deriv), lo, hi,
                                 rtol=2eps(), order=57)[1]
        @test result ≈ expected rtol=1e-2 atol=1e-64
        #resultnormal = LMV.integrate(normalring, kernel, deriv)
        #@show result, resultnormal
        #expectednormal = QuadGK.quadgk(v->kernel(v)*normalring(v, deriv), 0, Inf, rtol=10eps())[1]
        #@show resultnormal, result
        #@show expectednormal, expected
      end
      #end # for m in ...
    end
    #end # for n in ...
  end
  @test_throws ArgumentError FWideRing(1.0, 0.0)
  @test_warn "unreliable" FWideRing(1.0, 20.0 * (1 + eps()))
end

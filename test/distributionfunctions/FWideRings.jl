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
    ring = LMV.FWideRing(BigFloat(vth, precision=128), BigFloat(vd, precision=128))
    ringfloat = LMV.FWideRing(vth, vd)
    for n in (-10, -5, 0, 5, 10, 40)
    #n = rand((-50:50))
    for m in (n-1:n+1)
    #m = rand(-1:1) + n 
    result = QuadGK.quadgk(v->2π*v*ring(v), 0, 12 * (vth + vd))[1]
    @test result ≈ 1 rtol=sqrt(eps())
    @test LMV.is_normalised(ring)
    for deriv in (true, false)
      for k⊥_Ω in 0.25:0.25:4
        #@show "Ratio $z, bessel index $n, $m, deriv $deriv, vth k⊥_Ω $(vth*k⊥_Ω)"
        @testset "Ratio $z, bessel index $n, $m, deriv $deriv, vth k⊥_Ω $(vth*k⊥_Ω)" begin
          kernel = LMV.PerpendicularKernel(k⊥_Ω, Pair(n, m), rand(0:3))
          tr = @elapsed result = LMV.integrate(ring, kernel, deriv)
          lo = LMV.lower(ring)
          hi = LMV.upper(ring)
          te = @elapsed expected = QuadGK.quadgk(v->kernel(v)*ringfloat(v, deriv),
           lo, 2hi, rtol=2eps(), order=57)[1]
          #@show te / tr
          @test result ≈ expected rtol=1e-2 atol=1e-64
          #resultnormal = LMV.integrate(normalring, kernel, deriv)
          #@show result, resultnormal
          #expectednormal = QuadGK.quadgk(v->kernel(v)*normalring(v, deriv), 0, Inf, rtol=10eps())[1]
          #@show resultnormal, result
          #@show expectednormal, expected
        end
      end
      end # for m in ...
    end
    end # for n in ...
  end

  @test LMV.is_normalised(FWideRing(1.0, 1.0))
  @test_throws ArgumentError FWideRing(1.0, 0.0)
  @test_warn "unreliable" FWideRing(1.0, 20.0 * (1 + eps()))
end

using Dates
println("$(now()) $(@__FILE__)")

using Test, SpecialFunctions, QuadGK
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

verbose = false

@testset "Plasma dispersion functions vs Quadrature" begin
  @inline function test_maxwellian(Z, pow)
    for (i, z) in enumerate(Z)
      principle(x) = x.^pow .* exp(-x.^2)/sqrt(pi)
      integrand(x) = principle(x) / (x-z)
      ap = 0.0
      if iszero(imag(z))
        folded_integrand(v) = (principle(v + z) - principle(-v + z)) ./ v
        ap = QuadGK.quadgk(folded_integrand, eps(), 10*abs(z), rtol=eps())[1]
        @assert iszero(imag(ap))
      else
        ap += QuadGK.quadgk(integrand, -10*abs(z), 10*abs(z), rtol=eps())[1]
      end
      ar = LMV.residue(principle, LMV.Pole(z, 1, 0.0)) # no deformation for this test
      a = ap + ar
      b = LMV.plasma_dispersion_function(z, pow)
      outcome = isapprox(a, b, rtol=sqrt(eps()), atol=eps())
      @test real(a) ≈ real(b) rtol=sqrt(eps()) atol=eps()
      @test imag(a) ≈ imag(b) rtol=sqrt(eps()) atol=eps()
    end
  end
  Zs = []
  push!(Zs,  1.0 + im*0)
  push!(Zs, -1.0 + im*0)
  push!(Zs,  1.0 - im/4)
  push!(Zs, -1.0 - im/4)
  push!(Zs,  1.0 + im/4)
  push!(Zs, -1.0 + im/4)
  for z in Zs
    @testset "QuadGK of maxwellian when z = $z" begin
      @testset "QuadGK of maxwellian with 0th moment" begin
        test_maxwellian(z, 0)
      end
      @testset "QuadGK of maxwellian with 1st moment" begin
        test_maxwellian(z, 1)
      end
      @testset "QuadGK of maxwellian with 2nd moment" begin
        test_maxwellian(z, 2)
      end
    end
  end
end

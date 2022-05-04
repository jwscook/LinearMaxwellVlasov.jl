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
      ar = LMV.residue(principle, z)
      #ar = ar # 2 wrong
      #ar = sign(real(z)) < 0 ? conj(ar) : ar # 4 wrong
      #ar = sign(real(z)) < 0 ? -ar : ar # 5 wrong
      #ar = sign(real(z)) < 0 ? -conj(ar) : ar # 3 wrong
      b = LMV.plasma_dispersion_function(z, pow)
      a = ap + ar
      #iszero(imag(z)) && (a += ar)
      outcome = isapprox(a, b, rtol=sqrt(eps()), atol=eps())
      @test real(a) ≈ real(b) rtol=sqrt(eps()) atol=eps()
      @test imag(a) ≈ imag(b) rtol=sqrt(eps()) atol=eps()
      #verbose && outcome || @show pow, i, z
      #verbose && outcome || @show b, a, ap, ar
    end
  end
  Zs = []
  push!(Zs,  1.0 + im*0)
  push!(Zs, -1.0 + im*0)
  push!(Zs,  1.0 - im/4)
  push!(Zs, -1.0 - im/4)
  push!(Zs,  1.0 + im/4)
  push!(Zs, -1.0 + im/4)
  @testset "QuadGK of maxwellian with 0th moment" begin
    test_maxwellian(Zs, 0)
  end
  @testset "QuadGK of maxwellian with 1st moment" begin
    test_maxwellian(Zs, 1)
  end
  @testset "QuadGK of maxwellian with 2nd moment" begin
    test_maxwellian(Zs, 2)
  end
end

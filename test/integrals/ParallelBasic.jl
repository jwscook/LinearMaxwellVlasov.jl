using Dates
println("$(now()) $(@__FILE__)")

using Test, QuadGK, SpecialFunctions
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

verbose = false

@testset "Parallel: basic integral logic" begin

  function singular_integral_basic(num, den, pole)
    @assert imag(pole) == 0.0
    a, b = pole*(1.0 - 10eps())-eps(), pole*(1.0 + 10eps())+eps()
    a, b = min(a, b), max(a, b)
    A = QuadGK.quadgk(x->num(x)/den(x), -Inf,a)[1]
    B = QuadGK.quadgk(x->num(x)/den(x), b, Inf)[1]
    return A+B
  end
  function singular_integral_advanced(num, pole)
    @assert imag(pole) == 0.0
    return QuadGK.quadgk(x->(num(x+pole)-num(-x+pole))/x, eps(),Inf)[1]
  end

  num(x, n=0) = x.^n.*exp.(-x.^2)/sqrt(pi)
  den(x, pole) = x-pole
  for n in [0, 1, 2]
    for p in range(-6.0, stop=6.0, length=10)
      a = singular_integral_basic(x->num(x, n), x->den(x, p), p)
      b = singular_integral_advanced(x->num(x, n), p)
      try
        @test a ≈ b rtol=1.0e-1
      catch
        println(n, " ", p, " ", a, " ", b)
      end
    end
  end

end

@testset "Plasma dispersion functions vs Quadrature" begin
  @inline function test_maxwellian_quad(Z, pow)
    for z in Z
      pole = LMV.Pole(z, 1)#sign(real(z)))
      principal(x) = x.^pow .* exp(-x.^2)/sqrt(pi)
      integrand(x) = principal(x) / (x-z)
      ap = QuadGK.quadgk(integrand, -12 + im * pole.deformation, 12 + im * pole.deformation, rtol=eps())[1]
      ar = LMV.residue(principal, pole)
      b = LMV.plasma_dispersion_function(z, pow)
      a = ap + ar
      @test real(a) ≈ real(b) rtol=1.0e-8
      @test imag(a) ≈ imag(b) rtol=1.0e-8 atol=1.0e-8
      #verbose && isapprox(a, b, rtol=1.0e-8) || @show z, pow, ap, ar, a, b
    end
  end

  for z in (1.0 + im*0, -1.0 + im*0,  1.0 + im/4, -1.0 + im/4,  1.0 - im/4, -1.0 - im/4)
    @testset "z = $z" begin
      test_maxwellian_quad(z, 0)
      test_maxwellian_quad(z, 1)
      test_maxwellian_quad(z, 2)
      test_maxwellian_quad(z, 10)
    end
  end
end

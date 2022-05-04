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
      principal(x) = x.^pow .* exp(-x.^2)/sqrt(pi)
      integrand(x) = principal(x) / (x-z)
      ap = 0.0
      if iszero(imag(z))
        folded = LMV.foldnumeratoraboutpole(principal, real(z))
        ap += QuadGK.quadgk(folded, eps(), 6 + 6*abs(z), rtol=eps())[1]
      else
        ap += QuadGK.quadgk(integrand, -10*abs(z), 10*abs(z), rtol=eps())[1]
      end
      ar = LMV.residue(principal, z)
      b = LMV.plasma_dispersion_function(z, pow)
      a = ap + ar
      @test real(a) ≈ real(b) rtol=1.0e-8
      @test imag(a) ≈ imag(b) rtol=1.0e-8 atol=1.0e-8
      #verbose && isapprox(a, b, rtol=1.0e-8) || @show z, pow, ap, ar, a, b
    end
  end
  @inline function test_maxwellian_log(Z, pow)
    for z in Z
      principal(x) = x.^pow .* exp(-x.^2)/sqrt(pi)
      denom(x) = x - z
      integrand(x) = principal(x) / denom(x)
      h = abs(z) * sqrt(eps()) * 10
      dprincipaldx(x) = (principal(x + h) - principal(x - h)) / 2h
      ddenomdx(x) = 1.0 # do derivative by hand because of complex z
      function logdenom(x)
        output = log(Complex(denom(x)))
        return isfinite(output) ? output : zero(output)
      end
      target(x) = - logdenom(x) * ddenomdx(x) * dprincipaldx(x)
      a = QuadGK.quadgk(target, -6, 6, rtol=eps())[1]
      b = LMV.plasma_dispersion_function(z, pow)
      if imag(z) < 0
        # the log method doesn't work for negative imaginary part
        @test !isapprox(real(a), real(b), rtol=1.0e-8)
        @test !isapprox(imag(a), imag(b), rtol=1.0e-8, atol=1.0e-8)
      else
        @test real(a) ≈ real(b) rtol=1.0e-8
        @test imag(a) ≈ imag(b) rtol=1.0e-8 atol=1.0e-8
      end
      #verbose && isapprox(a, b, rtol=1.0e-8) || @show z, pow, ap, ar, a, b
    end
  end

  Zs = []
  push!(Zs,  1.0 + im*0)
  push!(Zs, -1.0 + im*0)
  push!(Zs,  1.0 + im/4)
  push!(Zs, -1.0 + im/4)
  push!(Zs,  1.0 - im/4) # log fails
  push!(Zs, -1.0 - im/4) # log fails
  @testset "QuadGK of maxwellian with 0th moment" begin
    test_maxwellian_quad(Zs, 0)
  end
  @testset "QuadGK of maxwellian with 1st moment" begin
    test_maxwellian_quad(Zs, 1)
  end
  @testset "QuadGK of maxwellian with 2nd moment" begin
    test_maxwellian_quad(Zs, 2)
  end
  @testset "QuadGK of maxwellian with 10th moment" begin
    test_maxwellian_quad(Zs, 10)
  end
  @testset "Log integral of maxwellian with 0th moment" begin
    test_maxwellian_log(Zs, 0)
  end
  @testset "Log integral of maxwellian with 1st moment" begin
    test_maxwellian_log(Zs, 1)
  end
  @testset "Log integral of maxwellian with 2nd moment" begin
    test_maxwellian_log(Zs, 2)
  end
end

using Dates
println("$(now()) $(@__FILE__)")

using Random, Test, InteractiveUtils, QuadGK

using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov
using PlasmaDispersionFunctions

@testset "Poles" begin

@testset "CauchyResidues" begin
  for σ ∈ (1, -1)
    @testset "principalpart, wavenumber sign is $σ" begin
      for i ∈ 1:10
        a = 10.0^rand(-5:5) * (rand() - 0.5) + im * 10.0^rand(-5:5) * (rand() - 0.5)
        z = 10.0^rand(-5:5) * (rand() - 0.5) + im * 10.0^rand(-5:5) * (rand() - 0.5)
        pole = LMV.Pole(z, σ)
        f(x) = a / (x - pole)
        RP1 = LMV.residuepart(f, pole)
        RP2 = LMV.residuepartadaptive(f, pole, 8)
        @test RP1 ≈ a atol=0 rtol=sqrt(eps())
        @test RP2 ≈ a atol=0 rtol=sqrt(eps())
      end
    end
  end
end
@testset "isreal isfinite" begin
  for r in (1.0, Inf), s in (-1, 0, 1)
    @test !isreal(LMV.Pole(r + 1im, s))
    @test isreal(LMV.Pole(r + 0im, s))
  end
  @test !isfinite(LMV.Pole(Inf + 0im, 1))
  @test !isfinite(LMV.Pole(0 + Inf*im, 1))
  @test !isfinite(LMV.Pole(Inf + Inf*im, 1))
end

@testset "integral contour deformation logic" begin
  numerator(x) = exp(-x*x) / sqrt(π)
  foobles(x, z) = numerator(x) ./ (x - z)
  for r in (1.0, ), i in (0.1,-0.1, 0.0), deformation in (0.0, 0.01, -0.01, -0.2, 0.2)
    d = deformation
    z = r + im * i
    iszero(d - i) && continue # otherwise test logic is broken ...
    pole = LMV.Pole(z, 1)

    expected = plasma_dispersion_function(z)

    σ = LMV.residuesigma(pole - im * d)
    rs = im * π * σ * numerator(z)
    ab = QuadGK.quadgk(x->foobles(x, z), -12 + im * d, 12 + im * d)[1] # ... here
    manualresult = ab + rs

    deformedpole = LMV.Pole(z, 1, d)
    result = QuadGK.quadgk(x->foobles(x, z), -12 + im * d, 12 + im * d)[1] + LMV.residue(numerator, deformedpole)

    @testset "r=$r, i=$i, d=$d" begin
      @test expected ≈ manualresult
      @test expected ≈ result
    end
  end
end

@testset "integral contour deformation" begin
  numerator(x) = exp(-x*x) / sqrt(π)
  foobles(x, z) = numerator(x) ./ (x - z)
  for r in (1.0,), i in (0.1, -0.1, 0.0, 1e-10, -1e-10)
    z = r + im * i

    expected = plasma_dispersion_function(z)

    d = LMV.imagcontourdeformation(z)
    deformedpole = LMV.Pole(z, 1, d)
    result = QuadGK.quadgk(x->foobles(x, z), -12 + im * d, 12 + im * d)[1] + LMV.residue(numerator, deformedpole)

    @testset "r=$r, i=$i, d=$d" begin
      @test expected ≈ result
    end
  end
end

end

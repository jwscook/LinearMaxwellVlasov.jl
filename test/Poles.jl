using Dates
println("$(now()) $(@__FILE__)")

using Random, Test, InteractiveUtils

using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


@testset "Poles" begin
@testset "Pole fix" begin
  pole⁺ = LMV.Pole(1.0, 1)
  z = 1.0 + im
  for op ∈ (pole⁺, LMV.wavedirectionalityhandler(pole⁺))
    @test op(1) == 1
    @test op(z) == z
    @test op(conj(z)) == conj(z)
  end
  pole⁻ = LMV.Pole(1.0, -1)
  for op ∈ (pole⁻, LMV.wavedirectionalityhandler(pole⁻))
    @test op(1) == 1
    @test op(z) == conj(z)
    @test op(conj(z)) == z
  end
end
@testset "Pole fix code_warntype" begin
  Random.seed!(0) # seed rand
  verbose = false

  pole = LMV.Pole(1.0, 1)
  z = 1.0 + im
  inferredworks = true
  try
    @inferred pole(z)
  catch
    inferredworks = false
  end
  @test inferredworks
end

@testset "CauchyResidues" begin
  for σ ∈ (1, -1)
    @testset "principalpart, wavenumber sign is $σ" begin
      for i ∈ 1:10
        a = 10.0^rand(-5:5) * (rand() - 0.5) + im * 10.0^rand(-5:5) * (rand() - 0.5)
        z = 10.0^rand(-5:5) * (rand() - 0.5) + im * 10.0^rand(-5:5) * (rand() - 0.5)
        pole = LMV.Pole(z, σ)
        f(x) = a / (x - pole)
        PP1 = LMV.principalpart(f, pole)
        PP2 = LMV.principalpartadaptive(f, pole, 8)
        @test PP1 ≈ a atol=0 rtol=sqrt(eps())
        @test PP2 ≈ a atol=0 rtol=sqrt(eps())
      end
    end
  end
end
end

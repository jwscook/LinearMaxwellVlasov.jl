using Dates
println("$(now()) $(@__FILE__)")

using Test, SpecialFunctions, QuadGK
using LinearMaxwellVlasov


@testset "Parallel: sign of real component of parallel wavenumber" begin

unitG = SeparableVelocitySpecies(1.0, 1.0,
  FParallelNumerical(1.0, -1.0),
  FPerpendicularNumerical(1.0))
unitR = RingBeamSpecies(1.0, 1.0, 1.0, 1.0, -1.0)
F = ComplexF64(0.5, 0.5)
δs = [-1.0, -0.1, -sqrt(eps()), -eps(), 0.0, eps(), sqrt(eps()), 0.1, 1.0]
δs = [-0.1, 0.0, 0.1]
θs = [-π, -3π/4, -π/2, -π/4, 0.0, π/4, π/2, 3π/4, π]
for θ ∈ θs, p ∈ Unsigned.(0:2)
  @testset "θ = $θ, pow = $p, diffbool = false" begin
    for δ ∈ δs
      K = Wavenumber(wavenumber=1.0 + δ*im, propagationangle=θ)
      C = Configuration(F, K)
      a = parallel(unitR, C, 0, p, false)
      b = parallel(unitG, C, 0, p, false)
      @test a ≈ b
    end
  end

  p > 1 && continue
  @testset "θ = $θ, pow = $p, diffbool = true" begin
    for δ ∈ δs
      K = Wavenumber(wavenumber=1.0 + δ*im, propagationangle=π/8)
      C = Configuration(F, K)
      a = parallel(unitR, C, 0, p, true)
      b = parallel(unitG, C, 0, p, true)
      @test a ≈ b
    end
  end
end

end
#

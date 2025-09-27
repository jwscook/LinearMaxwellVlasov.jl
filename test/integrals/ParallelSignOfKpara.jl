using Dates
println("$(now()) $(@__FILE__)")

using Test, SpecialFunctions, QuadGK
using LinearMaxwellVlasov


@testset "Parallel: sign of real component of parallel wavenumber" begin

unitG = SeparableVelocitySpecies(1.0, 1.0,
  FParallelNumerical(1.0, -1.0),
  FPerpendicularNumerical(1.0))
unitR = RingBeamSpecies(1.0, 1.0, 1.0, 1.0, -1.0)
F = ComplexF64(0.5, 0.05)
δs = [-1.0, -0.1, -sqrt(eps()), -eps(), 0.0, eps(), sqrt(eps()), 0.1, 1.0]
#δs = [-0.1, 0.0, 0.1]
θs = [-π, -3π/4, -π/2, -π/4, 0.0, π/4, π/2, 3π/4, π]
for θ ∈ θs, p ∈ Unsigned.(0:2), diffbool = (true, false)
  diffbool && p > 1 && continue
  @testset "θ = $θ, pow = $p, diffbool=$diffbool" begin
    for δ ∈ δs
      K = Wavenumber(wavenumber=1.0 + δ*im, propagationangle=θ)
      C = Configuration(F, K)
      a = parallel(unitR, C, 0, p, diffbool)
      b = parallel(unitG, C, 0, p, diffbool)
      @test isapprox(a, b, rtol=1e-5, atol=1e-8)
    end
  end
end

end
#

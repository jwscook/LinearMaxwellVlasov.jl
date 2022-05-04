using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


@testset "Parameters" begin
  Ωi = 100 # low frequency limit
  K⊥ = Wavenumber(parallel=0.0, perpendicular=1.0)
  kz = Wavenumber(parallel=1.0, perpendicular=0.0)

  Va⁺ = magnetoacousticfrequency(1.0, 0.0, K⊥, 1)
  @test Va⁺ == 1.0
  Va⁻ = magnetoacousticfrequency(1.0, 0.0, K⊥, -1)
  @test Va⁻ == 0.0
  Va⁺ = magnetoacousticfrequency(1.0, 1.0, K⊥, 1)
  Va⁻ = magnetoacousticfrequency(1.0, 1.0, K⊥, -1)
  @test Va⁺ > Va⁻
  Vaslow = slowmagnetoacousticfrequency(1.0, 1.0, K⊥)
  @test Vaslow == Va⁻
  Vafast = fastmagnetoacousticfrequency(1.0, 1.0, K⊥)
  @test Vafast == Va⁺
  Va = shearfrequency(1.0, kz)
  Va⁺para = fastmagnetoacousticfrequency(1.0, 1.0, kz)

  ωᵥ = LMV.zerobetamagnetoacousticfrequency(1.0, K⊥, Ωi, 1)
  @test ωᵥ == 1.0
  ωᵥ⁻ = slowzerobetamagnetoacousticfrequency(1.0, K⊥, Ωi)
  ωᵥ⁺ = fastzerobetamagnetoacousticfrequency(1.0, K⊥, Ωi)
  @test ωᵥ⁺ > ωᵥ⁻

  Va = LMV.shearspeed(1.0, kz)
  Va⁺para = LMV.fastmagnetoacousticspeed(1.0, 1.0, kz)
  @test Va == Va⁺para
end

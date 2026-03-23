using Dates
println("$(now()) $(@__FILE__)")

using Test, Random, DualNumbers, SpecialFunctions, LinearAlgebra
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov
using Statistics
using ForwardDiff

Random.seed!(0)

@testset "Separable vs Coupled velocity tensors" begin
  mₑ = LMV.mₑ
  M = 1
  Z = -1

  rtol = 1e-5
  B0 = 3.0
  n0 = 1e20
  Π = 4.242118193111994e11
  Ω = -3.344812404889745e11
  ϵV = 10.0e3
  vth = thermalspeed(ϵV, mₑ)
  maxwellian = CoupledRelativisticSpecies(Π, Ω, mₑ, mₑ * vth)
  kz = 0.0
  k⊥ = 2830.7429179878372#ForwardDiff.Dual(2830.7429179878372,0.015707317311820675,0.9998766324816606)
  ω = 5.3545952051433167e11
  F = ComplexF64(ω, 0.0)
  K = Wavenumber(kz=kz, k⊥=k⊥)
  options = LMV.Options(rtols=1e-5, cubature_maxevals=1_000_000,
                        erroruponcubaturenonconvergence=true,)

  config = Configuration(F, K, options)
  output = LMV.contribution(maxwellian, config)
  @test all(isfinite, output)
end

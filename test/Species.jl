using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


include("species/NumericalSpecies.jl")

@testset "Species" begin
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  B0 = 1.0
  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωi = cyclotronfrequency(B0, mi, 1)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πi = plasmafrequency(n0, mi, 1)
  ϵV = 1.0e3
  vth = thermalspeed(ϵV, mi)
  lri = vth / Ωi
  Va = sqrt(B0^2/LMV.μ₀/n0/mi)
  small = 100*sqrt(eps())
  vth = Va * small
  cold = ColdSpecies(Πi, Ωi)
  
  @test plasmafrequency(cold) == cold.Π
  @test cyclotronfrequency(cold) == cold.Ω
  @test LMV.is_normalised(cold)
  
  warm = WarmSpecies(Πi, Ωi, vth)
  
  @test ColdSpecies(warm) == cold
  
  maxwellian = MaxwellianSpecies(Πi, Ωi, vth, vth)
  
  @test LMV.is_normalised(maxwellian)
  @test ColdSpecies(maxwellian) == cold
  @test WarmSpecies(maxwellian, 1) == warm
  
  separable = SeparableVelocitySpecies(Πi, Ωi,
      FParallelNumerical(vth),
      FPerpendicularNumerical(vth))
  
  @test maxwellian(vth, vth) == maxwellian([vth, vth])
  @test separable(vth, vth) ≈ maxwellian(vth, vth) rtol=1e-4
  
  @test ColdSpecies(separable) == cold
  
  ringbeam = RingBeamSpecies(Πi, Ωi, vth, vth)
  
  @test LMV.is_normalised(ringbeam)
  
  coupled = CoupledVelocitySpecies(Πi, Ωi, vth, vth)
  
  @test LMV.is_normalised(coupled)
  
  coupled = CoupledVelocitySpecies(Πi, Ωi, vth)
  @test coupled(vth, vth) ≈ maxwellian(vth, vth)
  relativistic = CoupledRelativisticSpecies(Πi, Ωi, mi, vth * mi)
  @test maxwellian(vth, vth) ≈ relativistic(mi * vth, mi * vth) * mi^3
end

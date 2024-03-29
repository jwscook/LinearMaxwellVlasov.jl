using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


@testset "Cold Limit conductivities" begin

mₑ = LMV.mₑ
mi = 1836*mₑ
n0 = 1.0e19
B0 = 1.0
Ωe = cyclotronfrequency(B0, mₑ, -1)
Ωi = cyclotronfrequency(B0, mi, 1)
Πe = plasmafrequency(n0, mₑ, -1)
Πi = plasmafrequency(n0, mi, 1)
ϵV = 1.0
vth = thermalspeed(ϵV, mi)
Va = sqrt(B0^2/LMV.μ₀/n0/mi)
small = 100*sqrt(eps())
lri = vth / Ωi
cold = ColdSpecies(Πi, Ωi)
warm = WarmSpecies(Πi, Ωi, vth)
delta = SeparableVelocitySpecies(Πi, Ωi,
  FParallelDiracDelta(0.0),
  FPerpendicularDiracDelta(eps()))
maxwellian = MaxwellianSpecies(Πi, Ωi, vth, vth)
numerical = SeparableVelocitySpecies(Πi, Ωi,
    FParallelNumerical(vth),
    FPerpendicularNumerical(vth))
coupled = CoupledVelocitySpecies(Πi, Ωi, vth, vth)
relativistic = CoupledRelativisticSpecies(Πi, Ωi, mi, vth * mi)

k = Ωi / vth
ω = k * vth * 20
γ = ω / 100
K = Wavenumber(k=k, θ=π/5)
F = ComplexF64(ω, γ)
O = Options(rtols=100eps())
C = Configuration(F, K, O)

tol = 1.0e-3

outputCold = conductivity(cold, C)
outputDelta = conductivity(delta, C)
outputWarm = conductivity(warm, C)
outputMaxwellian = conductivity(maxwellian, C)
outputNumerical = conductivity(numerical, C)
O = Options(rtols=1e-8)
C = Configuration(F, K, O)
outputCoupled = conductivity(coupled, C)
outputRelativistic = conductivity(relativistic, C)
for i ∈ 1:3, j ∈ 1:3
  expected = outputCold[i, j]
  @test expected ≈ outputDelta[i, j] rtol=tol atol=tol
  @test expected ≈ outputWarm[i, j] rtol=tol atol=tol
  @test expected ≈ outputMaxwellian[i, j] rtol=tol atol=tol
  @test expected ≈ outputNumerical[i, j] rtol=tol atol=tol
  @test expected ≈ outputCoupled[i, j] rtol=tol atol=tol
  @test expected ≈ outputRelativistic[i, j] rtol=tol atol=tol
end
# TODO this is a rubbish test, improve
for s ∈ (maxwellian, numerical, coupled, relativistic)
  @test LMV.lowerintegralbounds(s)[1] < 0
  @test LMV.lowerintegralbounds(s)[2] >= 0
  @test LMV.upperintegralbounds(s)[1] > 0
  @test LMV.upperintegralbounds(s)[2] >= 0
end

end

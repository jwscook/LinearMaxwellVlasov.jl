using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


@testset "Cold Limit Dielectrics" begin

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

k = Ωi / vth # tests for num and gyro fail if no multiply by 10
ω = k * vth * 100
γ = ω / 100
K = Wavenumber(k=k, θ=π/5)
F = ComplexF64(ω, γ)
O = Options(rtols=100eps())
C = Configuration(F, K, O)

tol = 1.0e-3

outputCold = dielectric(cold, C)
outputDelta = dielectric(delta, C)
outputWarm = dielectric(warm, C)
outputMaxwellian = dielectric(maxwellian, C)
outputNumerical = dielectric(numerical, C)
outputCoupled = dielectric(coupled, C)
outputRelativistic = dielectric(relativistic, C)
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

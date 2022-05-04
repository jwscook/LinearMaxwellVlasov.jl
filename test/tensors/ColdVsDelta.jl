using Dates
println("$(now()) $(@__FILE__)")
using Test
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

@testset "Tensors: Cold vs Delta" begin
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  B0 = 1.0
  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωi = cyclotronfrequency(B0, mi, 1)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πi = plasmafrequency(n0, mi, 1)
  Va = sqrt(B0^2/LMV.μ₀/n0/mi)
  protonD = SeparableVelocitySpecies(Πi, Ωi,
    FParallelDiracDelta(0.0),
    FPerpendicularDiracDelta(eps()))
  protonC = ColdSpecies(Πi, Ωi)
  
  k = Ωi/Va*0.75
  K = Wavenumber(k, π/4)
  F = ComplexF64(0.25*Ωi, 0.0)
  C = Configuration(F, K)
  output20 = zeros(ComplexF64, 3, 3)
  function sum_tensors(num)
    outputC = zeros(ComplexF64, 3, 3)
    outputD = zeros(ComplexF64, 3, 3)
    for n in -num:num
      outputD += LMV.contribution(protonD, C, n)
    end
    outputC = LMV.contribution(protonC, C) # not a function of n
    return outputC, outputD
  end
  outputC40, outputD40 = sum_tensors(40)
  outputC = LMV.contribution(protonC, C)
  outputD = LMV.contribution(protonD, C)
  for i in 1:3, j in 1:3
    @test real(outputD[i,j])≈real(outputD40[i,j]) rtol=100sqrt(eps()) atol=10eps()
    @test imag(outputD[i,j])≈imag(outputD40[i,j]) rtol=100sqrt(eps()) atol=10eps()
    @test real(outputD[i,j])≈real(outputC40[i,j]) rtol=100sqrt(eps()) atol=10eps()
    @test imag(outputD[i,j])≈imag(outputC40[i,j]) rtol=100sqrt(eps()) atol=10eps()
    @test real(outputD[i,j])≈real(outputC[i,j]) rtol=100sqrt(eps()) atol=10eps()
    @test imag(outputD[i,j])≈imag(outputC[i,j]) rtol=100sqrt(eps()) atol=10eps()
  end
end

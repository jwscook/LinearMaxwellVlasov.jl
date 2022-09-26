using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

@testset "Maxwellian vs Numerical tensors" begin
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  B0 = 1.0
  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωi = cyclotronfrequency(B0, mi, 1)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πi = plasmafrequency(n0, mi, 1)
  ϵV = 1.0e3
  vthe = thermalspeed(ϵV, mₑ)
  vthi = thermalspeed(ϵV, mi)
  λD = vthe / Πe
  lri = vthi / Ωi
  Va = sqrt(B0^2/LMV.μ₀/n0/mi)
  #electronN = SeparableVelocitySpecies(Πe, Ωe,
  #    FParallelNumerical(vthe),
  #    FPerpendicularNumerical(vthe))
  protonN = SeparableVelocitySpecies(Πi, Ωi,
      FParallelNumerical(vthi),
      FPerpendicularNumerical(vthi))
  #N = [electronN, protonN]

  electronA = MaxwellianSpecies(Πe, Ωe, vthe, vthe)
  protonA = MaxwellianSpecies(Πi, Ωi, vthi, vthi)

  k = 2π/lri*0.5
  K = Wavenumber(k=k, θ=π/4)
  F = ComplexF64(0.5*Ωi, 0.5*Ωi)
  C = Configuration(F, K)
  output20 = zeros(ComplexF64, 3, 3)
  function sum_tensors(num)
    outputA = zeros(ComplexF64, 3, 3)
    outputN = zeros(ComplexF64, 3, 3)
    for n ∈ -num:num
      outputN += LMV.contribution(protonN, C, n)
      outputA += LMV.contribution(protonA, C, n)
    end
    return outputA, outputN
  end
  num = 100
  outputASUM, outputNSUM = sum_tensors(num)

  outputA = LMV.contribution(protonA, C)
  outputN = LMV.contribution(protonN, C)

  @testset "Inferred" begin
    for species ∈ (protonA, protonN)
      try
        @inferred LMV.contribution(species, C, 0)
        @test true
      catch
        @warn "contribution not inferred for $(nameof(typeof(species)))"
        @test_broken false
      end
      try
        @inferred LMV.contribution(species, C)
        @test true
      catch
        @warn "contribution not inferred for $(nameof(typeof(species)))"
        @test_broken false
      end
    end
  end

  for i ∈ 1:3, j ∈ 1:3
    outputNSUM[i,j] ≈ outputA[i, j]
    outputN[i,j] ≈ outputA[i, j]
    outputNSUM[i,j] ≈ outputN[i, j]
  end

  @testset "Maxwellian Converged vs Maxwellian summed to $num" begin
    @test real(outputA[1,1])≈real(outputASUM[1,1]) rtol=100sqrt(eps()) atol=0.0
    @test real(outputA[1,2])≈real(outputASUM[1,2]) rtol=100sqrt(eps()) atol=0.0
    @test real(outputA[1,3])≈real(outputASUM[1,3]) rtol=100sqrt(eps()) atol=0.0
    @test real(outputA[2,1])≈real(outputASUM[2,1]) rtol=100sqrt(eps()) atol=0.0
    @test real(outputA[2,2])≈real(outputASUM[2,2]) rtol=100sqrt(eps()) atol=0.0
    @test real(outputA[2,3])≈real(outputASUM[2,3]) rtol=100sqrt(eps()) atol=0.0
    @test real(outputA[3,1])≈real(outputASUM[3,1]) rtol=100sqrt(eps()) atol=0.0
    @test real(outputA[3,2])≈real(outputASUM[3,2]) rtol=100sqrt(eps()) atol=0.0
    @test real(outputA[3,3])≈real(outputASUM[3,3]) rtol=100sqrt(eps()) atol=0.0
  end
  @testset "Numerical vs numerical summed to $num" begin
    @testset "real" begin
      @test real(outputN[1,1])≈real(outputNSUM[1,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[1,2])≈real(outputNSUM[1,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[1,3])≈real(outputNSUM[1,3]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,1])≈real(outputNSUM[2,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,2])≈real(outputNSUM[2,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,3])≈real(outputNSUM[2,3]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,1])≈real(outputNSUM[3,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,2])≈real(outputNSUM[3,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,3])≈real(outputNSUM[3,3]) rtol=100sqrt(eps()) atol=0.0
    end
    @testset "imag" begin
      @test imag(outputN[1,1])≈imag(outputNSUM[1,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[1,2])≈imag(outputNSUM[1,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[1,3])≈imag(outputNSUM[1,3]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,1])≈imag(outputNSUM[2,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,2])≈imag(outputNSUM[2,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,3])≈imag(outputNSUM[2,3]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,1])≈imag(outputNSUM[3,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,2])≈imag(outputNSUM[3,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,3])≈imag(outputNSUM[3,3]) rtol=100sqrt(eps()) atol=0.0
    end
  end

  @testset "Numerical vs maxwellian summed to $num" begin
    @testset "real" begin
      @test real(outputN[1,1])≈real(outputASUM[1,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[1,2])≈real(outputASUM[1,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[1,3])≈real(outputASUM[1,3]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,1])≈real(outputASUM[2,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,2])≈real(outputASUM[2,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,3])≈real(outputASUM[2,3]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,1])≈real(outputASUM[3,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,2])≈real(outputASUM[3,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,3])≈real(outputASUM[3,3]) rtol=100sqrt(eps()) atol=0.0
    end
    @testset "imag" begin
      @test imag(outputN[1,1])≈imag(outputASUM[1,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[1,2])≈imag(outputASUM[1,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[1,3])≈imag(outputASUM[1,3]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,1])≈imag(outputASUM[2,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,2])≈imag(outputASUM[2,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,3])≈imag(outputASUM[2,3]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,1])≈imag(outputASUM[3,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,2])≈imag(outputASUM[3,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,3])≈imag(outputASUM[3,3]) rtol=100sqrt(eps()) atol=0.0
    end

  end
  @testset "Numerical and analytical tensors" begin
    @testset "real" begin
      @test real(outputN[1,1])≈real(outputA[1,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[1,2])≈real(outputA[1,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[1,3])≈real(outputA[1,3]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,1])≈real(outputA[2,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,2])≈real(outputA[2,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[2,3])≈real(outputA[2,3]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,1])≈real(outputA[3,1]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,2])≈real(outputA[3,2]) rtol=100sqrt(eps()) atol=0.0
      @test real(outputN[3,3])≈real(outputA[3,3]) rtol=100sqrt(eps()) atol=0.0
    end
    @testset "imag" begin
      @test imag(outputN[1,1])≈imag(outputA[1,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[1,2])≈imag(outputA[1,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[1,3])≈imag(outputA[1,3]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,1])≈imag(outputA[2,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,2])≈imag(outputA[2,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[2,3])≈imag(outputA[2,3]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,1])≈imag(outputA[3,1]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,2])≈imag(outputA[3,2]) rtol=100sqrt(eps()) atol=0.0
      @test imag(outputN[3,3])≈imag(outputA[3,3]) rtol=100sqrt(eps()) atol=0.0
    end
  end
end

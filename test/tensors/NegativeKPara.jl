using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

@testset "Tensors: negative parallel wavenumber are dealt with correctly" begin
  verbose = false
  
  unitM = MaxwellianSpecies(1.0, 1.0, 1.0, 1.0, -1.0)
  unitR = RingBeamSpecies(1.0, 1.0, 1.0, 1.0, -1.0)
  unitN = SeparableVelocitySpecies(1.0, 1.0,
      FParallelNumerical(1.0, -1.0),
      FPerpendicularNumerical(1.0))
  
  k = 2*π * 10
  K = Wavenumber(k=-k, θ=π/4)
  F = ComplexF64(0.5, 0.5)
  C = Configuration(F, K)

  @testset "contribution" begin
    for i ∈ 1:20
      n = rand(-10:10)
      outputN = LMV.contribution(unitN, C, n)
      outputR = LMV.contribution(unitR, C, n)
      for i ∈ 1:9
        @test isapprox(real(outputN[i]), real(outputR[i]),
                       rtol=1.0e-6, atol=eps())
        @test isapprox(imag(outputN[i]), imag(outputR[i]),
                       rtol=1.0e-6, atol=eps())
        v1, v2 = real(outputN[i]), real(outputR[i])
        w1, w2 = imag(outputN[i]), imag(outputR[i])
        iszero(v1) && iszero(v2) && iszero(w1) && iszero(w2) && continue
        #verbose && a || @show n, l, m, i, v1, v2, v1 / v2
        #verbose && b || @show n, l, m, i, w1, w2, w1 / w2
        #@test a
        #@test b
      end
    end
  end

  @testset "tensors" begin
    outputM = conductivity(unitM, C)
    outputN = conductivity(unitN, C)
    outputR = conductivity(unitR, C)
    for i ∈ 1:9
      @test real(outputN[i]) ≈ real(outputM[i]) rtol=1.0e-6
      @test imag(outputN[i]) ≈ imag(outputM[i]) rtol=1.0e-6
      @test real(outputR[i]) ≈ real(outputM[i]) rtol=1.0e-6
      @test imag(outputR[i]) ≈ imag(outputM[i]) rtol=1.0e-6
    end
  end
end

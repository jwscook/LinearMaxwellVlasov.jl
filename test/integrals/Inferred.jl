using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov
const DPR = LinearMaxwellVlasov


@testset "Inferred integrals" begin
  unitN = SeparableVelocitySpecies(1.0, 1.0,
    FParallelNumerical(1.0, -1.0),
    FPerpendicularNumerical(1.0))
  unitR = RingBeamSpecies(1.0, 1.0, 1.0, 1.0, -1.0)
  unitD = SeparableVelocitySpecies(1.0, 1.0,
    FParallelDiracDelta(-1.0),
    FPerpendicularDiracDelta(1.0))

  K = Wavenumber(wavenumber=1.0, propagationangle=π/4)
  F = ComplexF64(0.5, 0.5)
  C = Configuration(F, K)

  @testset "Parallel" begin
    for species ∈ (unitR, unitN, unitD)
      try
        @inferred LMV.parallel(species, C, 0, Unsigned(1), false)
        @test true
      catch
        @test false
      end
    end
  end

  @testset "Perpendicular" begin
    for species ∈ (unitR, unitN, unitD)
      try
        @inferred LMV.perpendicular(species, C, 0, 0, Unsigned(1), false)
        @test true
      catch
        @test false
      end
    end
  end

end


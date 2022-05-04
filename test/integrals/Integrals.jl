using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


@testset "Integrals memoisation" begin
  species = RingBeamSpecies(rand(6)...)
  @testset "parallel" begin
    @testset "do memoise" begin
      options = Options(memoiseparallel=true)
      config = Configuration(options)
      ∫ = LMV.parallel_integral(species, config)
      @test !(∫ === LMV.parallel)
    end
    @testset "dont memoise" begin
      options = Options(memoiseparallel=false)
      config = Configuration(options)
      ∫ = LMV.parallel_integral(species, config)
      @test ∫ === LMV.parallel
    end
  end

  @testset "perpendicular" begin
    @testset "do memoise" begin
      options = Options(memoiseperpendicular=true)
      config = Configuration(options)
      ∫ = LMV.perpendicular_integral(species, config)
      @test !(∫ === LMV.perpendicular)
    end
    @testset "dont memoise" begin
      options = Options(memoiseperpendicular=false)
      config = Configuration(options)
      ∫ = LMV.perpendicular_integral(species, config)
      @test ∫ === LMV.perpendicular
    end
  end
end


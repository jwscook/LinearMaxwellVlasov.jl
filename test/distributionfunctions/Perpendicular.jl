using Dates
println("$(now()) $(@__FILE__)")

using Test, Random
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

Random.seed!(0)

@testset "FPerpendicularNumerical derivative and norm" begin
  a = rand()
  b = rand()
  this = FPerpendicularNumerical(a, b)
  compare = FPerpendicularNumerical(this.F, this.lower, this.upper)
  @test this.F(a) ≈ compare.F(a) rtol=100eps()
  @test this.F(b) ≈ compare.F(b) rtol=100eps()
  @test this.dFdv(a) ≈ compare.dFdv(a)
  @test this.dFdv(b) ≈ compare.dFdv(b)
end
@testset "FRing vs FPerpendicularNumerical" begin
  a = rand()
  b = rand()
  
  ring = FRing(a, b)
  @test LMV.is_normalised(ring)
  num = FPerpendicularNumerical(a, b)

  for i ∈ 1:10
    v = rand()*4
    @test ring(v) ≈ num(v)
    @test ring(v, true) ≈ num(v, true)
    @test ring(v, false) ≈ num(v, false)
  end
end

@testset "FPerpendicularMaxwellian" begin
  fmx = FPerpendicularMaxwellian(1.0)
  @test LMV.is_normalised(fmx)
end

using Dates
println("$(now()) $(@__FILE__)")

using Test, Random
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

include("../species/NumericalSpecies.jl")

Random.seed!(0)
@testset "Perpendicular integrals" begin
  reltol = 1.0e-3
  tol = Tolerance(eps(), eps())

  @testset "Basic" begin
    for i in 1:5
      Π, Ω, vthb, vth⊥, vdb, vd⊥ = 2 .* rand(6)
      k⊥ = ComplexF64(rand() * 10)
      Ω *= 2 * rand() * sign(rand() - 0.5)
      vdb *= 2(rand() - 0.5)
      diffbool = rand() > 0.5
      N = NumericalSpecies(Π, Ω, vthb, vth⊥, vdb, vd⊥)
      R = RingBeamSpecies(Π, Ω, vthb, vth⊥, vdb, vd⊥)
      N0 = NumericalSpecies(Π, Ω, vthb, vth⊥, vdb, 0.0)
      R0 = MaxwellianSpecies(Π, Ω, vthb, vth⊥, vdb)
      pow = Unsigned(rand(0:3))
      n, l = rand(-5:5), rand(-5:5)
      nl = Pair(n, l)
      resultN = LMV.perpendicular(N.F⊥, N.Ω, nl, k⊥, pow, diffbool)[1]
      resultR = LMV.perpendicular(R.F⊥, R.Ω, nl, k⊥, pow, diffbool)
      resultN0 = LMV.perpendicular(N0.F⊥, N0.Ω, nl, k⊥, pow, diffbool)[1]
      resultR0 = LMV.perpendicular(R0.F⊥, R0.Ω, nl, k⊥, pow, diffbool)
      @test resultN ≈ resultR rtol=reltol
      @test resultN0 ≈ resultR0 rtol=reltol
      @test LMV.makepositive(-1.0) > 0
      @test all(reim(LMV.makepositive(-1.0 - im)) .> 0)
    end
  end

  @testset "Perpendicular short circuit key and factor" begin
    for i in 1:10
      Π, Ω, vth, vd = 1.0, 1.0, 1.0, 1.0
      # the number ξ = 2 * (k⊥ *vth / Ω)^2
      # integral gets difficult when ξ > 10
      k⊥ = 10 / sqrt(2) * Ω / vth
      k⊥ *= rand() < 0.5 ? -1 : 1
      A = FRing(vth, 0.0)
      B = FRing(vth, vd)
      pow = Unsigned(rand(0:3))
      for pow ∈ 0:3
        μνvector = sort([rand(-20:20), rand(-20:20)])
        μν = Pair(μνvector[1], μνvector[2])
        ∂F∂v = rand(Bool)
        aR = LMV.perpendicular(A, Ω, μν, k⊥, Unsigned(pow), ∂F∂v, tol)
        bR = LMV.perpendicular(B, Ω, μν, k⊥, Unsigned(pow), ∂F∂v, tol)
        kernel = LMV.PerpendicularKernel(k⊥ / Ω, μν, Unsigned(pow))
        aE = LMV.integrate(A, kernel, ∂F∂v, tol)
        bE = LMV.integrate(B, kernel, ∂F∂v, tol)
        @test isapprox(aE, aR)
        @test isapprox(bE, bR)
      end
    end
  end
  @testset "Sign K perp, even and odd μ ν" begin
    Π, Ω, vth, vd = 1.0, 1.0, 1.0, 1.0
    k⊥ = 10 / sqrt(2) * Ω / vth
    kz = k⊥
    A = FPerpendicularMaxwellian(vth) # FRing(vth, 0.0)
    B = FRing(vth, vd)
    # FBeam isn't whats being tested
    C = SeparableVelocitySpecies(Π, Ω, FBeam(1.0, 1.0), B)
    F = ComplexF64(1, 1)# doesn't matter
    # only do a subset of 0:1
    for ∂F∂v ∈ (false, true), pow ∈ 0:1, μ ∈ -2:2, ν ∈ -2:2, σ ∈ (-1, 1)
      K = Wavenumber(parallel=kz, perpendicular=σ * k⊥)
      config = Configuration(F, K)
      μν = Pair(μ, ν)
      aR = LMV.perpendicular(A, Ω, μν, σ * k⊥, Unsigned(pow), ∂F∂v, tol)
      bR = LMV.perpendicular(B, Ω, μν, σ * k⊥, Unsigned(pow), ∂F∂v, tol)
      kernel = LMV.PerpendicularKernel(σ * k⊥ / Ω, μν, Unsigned(pow))
      aE = LMV.integrate(A, kernel, ∂F∂v, tol)
      bE = LMV.integrate(B, kernel, ∂F∂v, tol)
      cN = LMV.perpendicular(C, config, μ, ν, Unsigned(pow), ∂F∂v)
      @test isapprox(aE, aR)
      @test isapprox(bE, bR)
      @test isapprox(cN, bE)
    end
    for n ∈ (-3:3), σ ∈ (-1, 1)
      K = Wavenumber(parallel=kz, perpendicular=σ * k⊥)
      config = Configuration(F, K)
      cN = LMV.perpendicular(C, config, n)
      memoisedperpendicular = LMV.perpendicular_integral(C, config)
      cM = memoisedperpendicular(C, config, n)
      @test all(isapprox.(cN, cM))
    end
  end
end


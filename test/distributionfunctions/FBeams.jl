using Dates
println("$(now()) $(@__FILE__)")

using Base.Threads, Random, LinearAlgebra, Test, StaticArrays, InteractiveUtils
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


@testset "RingBeam vs Numerical ring beams and pancakes" begin
  Random.seed!(0) # seed rand
  verbose = false
  
  mi = 1836*LMV.mₑ
  n0 = 2.0e20
  B0 = 2.1
  Va = sqrt(B0^2/LMV.μ₀/n0/mi)
  Ω = LMV.q₀ * B0 / mi

  vth = thermalspeed(1.0e4, mi)
  vd = thermalspeed(1.0e5, mi)

  numerical = FParallelNumerical(vth, vd)
  beam = FBeam(vth, vd)
  @test LMV.is_normalised(beam)
  powers = (0, 1, 2)
  @testset "Compare parallel integral of FBeam vs Numerical quadrature" begin
    for i ∈ 1:1, dfdv ∈ (true, false), p ∈ powers
      n = rand(-5:5)
      ω = Ω / Va * rand(ComplexF64)
      ω = rand() > 0.5 ? conj(ω) : ω
      KPara = 10 * (rand()-0.5) * Ω / Va
      dfdv && (p >= 2) && continue
      beamanswer = LMV.parallel(beam, ω, KPara, n, Ω, Unsigned(p), dfdv)
      numeanswer = LMV.parallel(numerical, ω, KPara, n, Ω, Unsigned(p), dfdv)
      @test beamanswer ≈ numeanswer rtol=sqrt(eps()) atol=eps()
    end
  end

  @testset "make sure @inferred passes for FBeam" begin
    ω = Ω / Va * rand(ComplexF64)
    KPara = Ω / Va
    inferred_false = true
    inferred_true = true
    for power ∈ (UInt64(0), UInt64(1))
      try
        @inferred LMV.parallel(beam, ω, KPara, 3, Ω, power, false)
      catch
        inferred_false = false
      end
      try
        @inferred LMV.parallel(beam, ω, KPara, 3, Ω, power, true)
      catch
        inferred_true = false
      end
      @test inferred_false
      @test inferred_true
    end
  end

end

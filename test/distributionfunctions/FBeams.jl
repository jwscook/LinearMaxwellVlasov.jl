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

  powers = (0, 1, 2)
  fails = []
  @testset "Compare parallel integral of FBeam vs Numerical quadrature" begin
    for ksign = (-1, 1), γsign = (-1, 1), vdsign = (-1, 1)
      vth = thermalspeed(1.0e4, mi)
      vd = rand() * 4 * vth * vdsign
      numerical = FParallelNumerical(vth, vd)
      beam = FBeam(vth, vd)
      for dfdv ∈ (true, false), p ∈ powers
        @testset "ksign=$ksign, γsign=$γsign, vdsign=$vdsign, dfdv=$dfdv, p=$p" begin
          for i ∈ 1:100
            n = rand(-5:5)
            ω = Ω / Va * Complex(rand(), γsign * rand() / 100)
            kpara = rand() * 10 * Ω / Va * ksign
            wavenumber = Wavenumber(kpara, 0.0)
            dfdv && (p >= 2) && continue
            beamanswertuple = LMV.parallel(beam, ω, wavenumber, n * Ω)
            ind = findfirst(i == (Unsigned(p), dfdv) for i in LMV.PARALLEL_TUPLE_ORDER)
            beamanswer = beamanswertuple[ind]
            numeanswer = LMV.parallel(numerical, ω, wavenumber, n * Ω, Unsigned(p), dfdv)
            diff = (beamanswer - numeanswer) ./ norm(beamanswer)
            if !isapprox(beamanswer, numeanswer, rtol=1000sqrt(eps()), atol=100eps())
              push!(fails, (ksign=ksign, γsign=γsign, vdsign=vdsign, ω=ω,
                            wavenumber=wavenumber, n=n, p=p, dfdv=dfdv, diff=diff))
            end
            @test beamanswer ≈ numeanswer rtol=100sqrt(eps()) atol=1000eps()
          end
        end
      end
    end
  end
  for f in fails
  @show f
  end

  @testset "make sure @inferred passes for FBeam" begin
    ω = Ω / Va * rand(ComplexF64)
    kpara = Ω / Va
    inferred = true
    vth = thermalspeed(1.0e4, mi)
    vd = (rand() - 0.5) * 8 * vth
    wavenumber = Wavenumber(kpara, 0.0)
    numerical = FParallelNumerical(vth, vd)
    beam = FBeam(vth, vd)
    @test LMV.is_normalised(beam)
    for power ∈ (UInt64(0), UInt64(1))
      try
        @inferred LMV.parallel(beam, ω, wavenumber, 3 * Ω)
      catch err
        @info err
        rethrow(err)
        inferred = false
      end
      @test inferred
    end
  end

end

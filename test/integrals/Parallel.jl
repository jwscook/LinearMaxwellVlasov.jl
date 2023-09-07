using Dates
println("$(now()) $(@__FILE__)")

using Base, SpecialFunctions, QuadGK, Test, UnitTestDesign
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

include("../species/NumericalSpecies.jl")

@testset "Parallel" begin
  verbose = false
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  B0 = 1.0
  Va = B0 / sqrt(n0*mi*LMV.μ₀)
  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωi = cyclotronfrequency(B0, mi, 1)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πi = plasmafrequency(n0, mi, 1)

  ϵV = 1.0e2
  vthe = thermalspeed(ϵV, mₑ)
  vthi = thermalspeed(ϵV, mi)
  λD = vthe / Πe
  ρi = vthi / Ωi
  Cs = sqrt(LMV.q₀*ϵV/mi)
  electron = NumericalSpecies(Πe, Ωe, vthe)
  proton = NumericalSpecies(Πi, Ωi, vthi)
  Ω1, Ω2 = 1.0, 2.0
  unit = NumericalSpecies(1.0, Ω1, 1.0)
  twos = NumericalSpecies(2.0, Ω2, 2.0)

  @assert LMV.is_normalised(electron)
  @assert LMV.is_normalised(proton)
  @assert LMV.is_normalised(unit)
  # Gradteyn and Rizhik 3.356
  @testset "QuadGK integrates normalised maxwellian to unity" begin
    @test QuadGK.quadgk(x->exp(-x.^2)/sqrt(pi), -Inf, Inf, rtol=eps())[1] ≈ 1.0
  end

  @assert LMV.is_normalised(unit)
  @assert LMV.is_normalised(twos)
  function para_integral(S, ω, kz, n, pow, diffbool, vth)
    @assert vth > 0.0
    M = RingBeamSpecies(S.Π, S.Ω, vth, vth, 0.0)
    return LMV.parallel(M.Fz, ω, kz, n, S.Ω, UInt64(pow), diffbool)
  end

  @testset "Test parallel integral numerical vs ring beam" begin
    a =      para_integral(unit,    1.0 + 0.0*im, 0.05 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, 0.05 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=1.0e-2           
                                      
    a =      para_integral(unit,    1.0 + 0.0*im, -0.05 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, -0.05 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=1.0e-2           
                                      
    a =      para_integral(unit,    1.0 + 0.0*im, 0.05 + 0.0*im, 0,     1        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, 0.05 + 0.0*im, 0, Ω1, UInt64(1), false)
    @test a ≈ b rtol=1.0e-2           
                                      
    a =      para_integral(unit,    1.0 + 0.0*im, -0.05 + 0.0*im, 0,     1        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, -0.05 + 0.0*im, 0, Ω1, UInt64(1), false)
    @test a ≈ b rtol=1.0e-2           
                                      
    a =      para_integral(unit,    1.0 + 0.0*im, 1.0 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, 1.0 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=1.0e-2                                       
                                                                  
    a =      para_integral(unit,    1.0 + 0.0*im, 1.0 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, 1.0 + 0.0*im, 0, Ω1, UInt64(0), false     )
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(unit,    1.0 + 0.0*im, 5.0 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, 5.0 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(unit,    5.0 + 0.0*im, 1.0 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 5.0 + 0.0*im, 1.0 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(unit,    5.0 + 0.0*im, 3.0 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 5.0 + 0.0*im, 3.0 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(unit,    1.0 + 0.0*im, 1.0 + 0.0*im, 0,     1        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, 1.0 + 0.0*im, 0, Ω1, UInt64(1), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(twos,    1.0 + 0.0*im, 1.0 + 0.0*im, 0,     0        , false, 2.0)
    b = parallel(twos.Fz, 1.0 + 0.0*im, 1.0 + 0.0*im, 0, Ω2, UInt64(0), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(twos,    2.0 + 0.0*im, 3.0 + 0.0*im, 1,     0        , false, 2.0)
    b = parallel(twos.Fz, 2.0 + 0.0*im, 3.0 + 0.0*im, 1, Ω2, UInt64(0), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(twos,    4.0 + 0.0*im, 7.0 + 0.0*im, 1,     1        , false, 2.0)
    b = parallel(twos.Fz, 4.0 + 0.0*im, 7.0 + 0.0*im, 1, Ω2, UInt64(1), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(twos,    4.0 - 0.0*im, 7.0 + 0.0*im, 1,     1        , false, 2.0)
    b = parallel(twos.Fz, 4.0 - 0.0*im, 7.0 + 0.0*im, 1, Ω2, UInt64(1), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(twos,    4.0 - 0.0*im, 7.0 + 0.0*im, 2,     1        , false, 2.0)
    b = parallel(twos.Fz, 4.0 - 0.0*im, 7.0 + 0.0*im, 2, Ω2, UInt64(1), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(unit,    1.0 + 0.1*im, 1.0 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.1*im, 1.0 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=sqrt(eps())                                  
    a =      para_integral(unit,    1.0 - 0.1*im, 1.0 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 1.0 - 0.1*im, 1.0 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=sqrt(eps())                                  
    a =      para_integral(unit,    1.0 + 0.0*im, 1.0 + 0.0*im, 0,     0        , false, 1.0)
    b = parallel(unit.Fz, 1.0 + 0.0*im, 1.0 + 0.0*im, 0, Ω1, UInt64(0), false)
    @test a ≈ b rtol=sqrt(eps())                                  
                                                                  
    a =      para_integral(twos,    1.0 + 0.1*im, 1.0 + 0.0*im, 2,     2        , false, 2.0)
    b = parallel(twos.Fz, 1.0 + 0.1*im, 1.0 + 0.0*im, 2, Ω2, UInt64(2), false)
    @test a ≈ b rtol=sqrt(eps())                                  
    a =      para_integral(twos,    1.0 - 0.1*im, 1.0 + 0.0*im, 2,     2        , false, 2.0)
    b = parallel(twos.Fz, 1.0 - 0.1*im, 1.0 + 0.0*im, 2, Ω2, UInt64(2), false)
    @test a ≈ b rtol=sqrt(eps())                                  
    a =      para_integral(twos,    1.0 + 0.0*im, 1.0 + 0.0*im, 2,     2        , false, 2.0)
    b = parallel(twos.Fz, 1.0 + 0.0*im, 1.0 + 0.0*im, 2, Ω2, UInt64(2), false)
    @test a ≈ b rtol=sqrt(eps())
  end

  function bigtest(speciesk0vth)
    (s, k0, vth) = speciesk0vth

    ktmp = 10.0.^range(-3, stop=3, length=7)

    ωs = range(-3, stop=3, length=2) * ComplexF64(s.Π)
    ks = [-ktmp; 0.0; ktmp] * ComplexF64(k0)

    for params ∈ all_pairs(ks, ωs, -1:1, 0:2, (true, false))
      (kz, ω, n, pow, diffbool) = params
      diffbool && abs(pow) == 2 && continue
      t1 = @elapsed output = LMV.parallel(s.Fz, ω, kz, n, s.Ω, UInt64(pow), diffbool)
      analytical = para_integral(s, ω, kz, n, UInt64(pow), diffbool, vth)
      Base.isnan(analytical) && continue
      @test output ≈ analytical  atol=1.0e-2 rtol=1.0e-2
    end
  end

  @testset "Numerical vs Analytical solution" begin
    ntest = 0
    for speciesk0vth in [(electron, 1.0/λD, vthe), (proton, 1.0/ρi, vthi)]
      ntest += 1
      @testset "test number $(ntest)" begin
        bigtest(speciesk0vth)
      end
    end
  end

end



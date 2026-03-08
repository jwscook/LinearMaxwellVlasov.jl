using Dates
println("$(now()) $(@__FILE__)")

using Base, SpecialFunctions, QuadGK, Test, UnitTestDesign
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

include("../species/NumericalSpecies.jl")


@testset "Parallel unnittest 4" begin
  verbose = false
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  B0 = 1.0
  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Πe = plasmafrequency(n0, mₑ, -1)
  Ωi = cyclotronfrequency(B0, mi, 1)
  Πi = plasmafrequency(n0, mi, 1)

  ϵV = 1.0e2
  vthe = thermalspeed(ϵV, mₑ)
  vthi = thermalspeed(ϵV, mi)
  electron = NumericalSpecies(Πe, Ωe, vthe)
  proton = NumericalSpecies(Πi, Ωi, vthi)

  @assert LMV.is_normalised(electron)
  @assert LMV.is_normalised(proton)

  function para_integral(S, ω, kz, n, pow, diffbool, vth)
    @assert vth > 0.0
    M = RingBeamSpecies(S.Π, S.Ω, vth, vth, 0.0)
    return LMV.parallel(M.Fz, ω, kz, n, S.Ω, UInt64(pow), diffbool)
  end

  for speciesk0vth in [(electron, vthe), (proton, vthi)]
    (s, vth) = (electron, vthe)
    (kz, ω, n, pow, diffbool) = (-69208.49198219362 - 0.0im, 1.2490402992112974e10 + 0.0im, 1, 2, false)
    k = Wavenumber(kz, 0.0)
    t1 = @elapsed output = LMV.parallel(s.Fz, ω, k, n * s.Ω, UInt64(pow), diffbool)
    analytical = para_integral(s, ω, k, n, UInt64(pow), diffbool, vth)
    @testset "test: kz=$(kz), ω=$(ω), n=$(n), pow=$(pow), diffbool=$(diffbool)" begin
      @test output ≈ analytical  atol=1.0e-2 rtol=1.0e-2
    end
  end

end



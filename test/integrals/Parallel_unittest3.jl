using Dates
println("$(now()) $(@__FILE__)")

using Test, Statistics
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

include("../species/NumericalSpecies.jl")

@testset "Parallel unittest 3" begin

  verbose = false
  number = exp(1) * π
  Ω = number
  vth = number
  twosnum = NumericalSpecies(number, Ω, vth)
  twosrb = RingBeamSpecies(number, Ω, vth)

  @assert LMV.is_normalised(twosnum)

  imagfails = []
  realfails = []
  imagratios = 0.0
  realratios = 0.0

  @testset "Test numerical parallel integral with against analytical result" begin
    p = 0
    for n ∈ (0,)#-20:20
      for ω ∈ (1.0 + 0*im, 1.0 + 0.1*im, 1.0 - 0.1im), kz in (1.0 + 0*im, 1.0 + 0.1*im, 1.0 - 0.1im)
        a = LMV.parallel(twosrb.Fz,  ω, kz, n, Ω, Unsigned(p), false)
        b = LMV.parallel(twosnum.Fz, ω, kz, n, Ω, Unsigned(p), false)
        @assert twosrb.Fz((ω - n* Ω) / kz) ≈ twosnum.Fz((ω - n* Ω) / kz)
        realratios += abs((real(a) - real(b)) / real(a))
        imagratios += abs((imag(a) - imag(b)) / imag(a))
        if !isapprox(real(a), real(b), rtol=1.0e-3)
          push!(realfails, n)
          verbose && @show "n = $n, ωandkz = $ωandkz"
        end
        if !isapprox(imag(a), imag(b), rtol=1.0e-3)
          push!(imagfails, n)
        end
        @testset "n = $n, ω = $ω, kz = $kz" begin
          expected = a
          result = b
          @test result ≈ expected rtol=1.0e-3
        end
      end
    end
    verbose && @show imagratios
    verbose && @show realratios
    verbose && @show unique(realfails)
    verbose && @show unique(imagfails)
  end

end






#

using Dates
println("$(now()) $(@__FILE__)")

using Test, Statistics
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

include("../species/NumericalSpecies.jl")

@testset "Parallel unittest 2" begin

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
    for n ∈ -20:20
      for ωandkz ∈ (1.0 + 0*im, 1.0 + 0.1*im, 1.0 - 0.1im)
        a = LMV.parallel(twosrb.Fz,  ωandkz, ωandkz, n, Ω, Unsigned(p), false)
        b = LMV.parallel(twosnum.Fz, ωandkz, ωandkz, n, Ω, Unsigned(p), false)
        realratios += abs((real(a) - real(b)) / real(a))
        imagratios += abs((imag(a) - imag(b)) / imag(a))
        if !isapprox(real(a), real(b), rtol=1.0e-3)
          r = (real(b) - real(a))/max(abs(real(a)), abs(real(a)))
          push!(realfails, n)
          verbose && @show "n = $n, ωandkz = $ωandkz"
        end
        if !isapprox(imag(a), imag(b), rtol=1.0e-3)
          r = (imag(b) - imag(a))/max(abs(imag(a)), abs(imag(a)))
          push!(imagfails, n)
        end
        @test a ≈ b rtol=1.0e-3
        #verbose && isapprox(a, b, rtol=1.0e-3) || @show n, ωandkz
      end
    end
    verbose && @show imagratios
    verbose && @show realratios
    verbose && @show unique(realfails)
    verbose && @show unique(imagfails)
  end

end






#

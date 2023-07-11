using Dates
println("$(now()) $(@__FILE__)")

using Base, SpecialFunctions, QuadGK, Test, UnitTestDesign, StaticArrays
using LinearAlgebra, Statistics
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

function exp_z_besseli_nz(n::Int, z::T) where {T<:Complex}
  return SpecialFunctions.besselix(n, z) * exp(-im * imag(z))
end
exp_z_besseli_nz(n::Int, z::T) where {T<:Real} = SpecialFunctions.besselix(n, z)

function brambilladielectriccontribution(species, config, cyclotronharmonic::Int)
  @assert !iszero(para(config.wavenumber)) # causes Inf
  @assert !iszero(species.Fz.vth) # causes Inf
  @assert species.F⊥.vth == species.Fz.vth
  @assert species.Fz.vd == 0

  ω = config.frequency
  kz, k⊥ = para(config.wavenumber), perp(config.wavenumber)
 
  n = cyclotronharmonic
  Ω = species.Ω
  nΩ = n * species.Ω
  vth = species.F⊥.vth
  μ = k⊥ * vth / Ω
  λ = μ^2 / 2
  Inexp = besseli(n, λ) * exp(-λ)
  Idnexp = (besseli(n - 1, λ) + besseli(n + 1, λ) ) / 2 * exp(-λ)
  xn = (ω - nΩ) / kz / vth 
  x0 = ω / kz / vth
  Zxn = LMV.plasma_dispersion_function(xn, 0)
  Zdxn = -2 * (1 + xn * Zxn)
  x0Zxn = x0 * Zxn
  nz = LMV.c₀ * kz / ω
  n⊥ = LMV.c₀ * k⊥ / ω

  output = @MArray zeros(ComplexF64, 3, 3)

  output[1, 1] = - n^2 / λ * Inexp * (-x0Zxn)
  output[1, 2] = - im * n * (Idnexp - Inexp) * (-x0Zxn)
  output[1, 3] = - 0.5 * n⊥ * nz * ω / Ω * vth^2 / LMV.c₀^2 * n / λ * Inexp * x0^2 * Zdxn
  output[2, 1] = -output[1, 2]
  output[2, 2] = - (n^2 / λ * Inexp - 2λ * (Idnexp - Inexp)) * (-x0Zxn)
  output[2, 3] = im /2 * n⊥ * nz * ω / Ω * vth^2 / LMV.c₀^2 * (Idnexp - Inexp) * x0^2 * Zdxn
  output[3, 1] = output[1, 3]
  output[3, 2] = -output[2, 3]
  output[3, 3] = -Inexp * x0 * xn * Zdxn

  return output
end

@testset "Compare against Brambilla" begin
  # Brambilla Eq 14.13 page 106
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
  electronmax = RingBeamSpecies(Πe, Ωe, vthe)
  protonmax = RingBeamSpecies(Πi, Ωi, vthi)
  electronrb = RingBeamSpecies(Πe, Ωe, vthe, vthe, 0.0, 0.0)
  protonrb = RingBeamSpecies(Πi, Ωi, vthi, vthi, 0.0, 0.0)

  λD = vthe / Πe
  lri = vthi / Ωi
  Va = sqrt(B0^2/LMV.μ₀/n0/mi)

  for (k, ω, species) in ((2π/λD, abs(Ωe), electronmax),
                          (2π/λD, abs(Ωe), electronrb),
                          (2π / lri, Ωi, protonmax),
                          (2π / lri, Ωi, protonrb))
    wavenumber = Wavenumber(k=k, θ=π/4)
    frequency = ComplexF64(0.5*ω, 0.5*ω)
    config = Configuration(frequency, wavenumber)
    for n in -5:5
      result = LMV.contribution(species, config, n)
      expected = brambilladielectriccontribution(species, config, n)
      normvals = max.(abs.(expected), abs.(result))
      normvals[normvals .== 0] .= 1
      for i in eachindex(result, expected)
        divvalue = iszero(expected[i]) ? one(eltype(expected)) : expected[i]
        normaliseddiff = (result[i] - expected[i]) / divvalue
        @test abs(normaliseddiff) < sqrt(eps())
      end
    end
  end

end

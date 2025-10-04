using Dates
println("$(now()) $(@__FILE__)")

using Test, Random, DualNumbers, SpecialFunctions, LinearAlgebra
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov
using Statistics

Random.seed!(0)

@testset "Separable vs Coupled velocity tensors" begin
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  verbose = false

  rtol = 1e-5
  for (M, Z) ∈ ((1836, 1),)#(100, -1))#, _ ∈ 1:2
    B0 = 3.0
    n0 = 1e20
    m = M * mₑ
    Ω = cyclotronfrequency(B0, m, Z)
    Π = plasmafrequency(n0, m, Z)
    Va = sqrt(B0^2 / LMV.μ₀ / n0 / mi)
    ϵV = 1.0e3
    vth = thermalspeed(ϵV, m)
    λD = vth / Π

    argsM = ones(3) * vth
    coupledMaxwellian = CoupledVelocitySpecies(Π, Ω, argsM...)
    separableMaxwellian = MaxwellianSpecies(Π, Ω, argsM...)
    argsR = ones(4) * vth
    coupledRingBeam = CoupledVelocitySpecies(Π, Ω, argsR...)
    separableRingBeam = RingBeamSpecies(Π, Ω, argsR...)
    norms = []

    for (coupled, separable) ∈ (
                                (coupledMaxwellian, separableMaxwellian),
                                (coupledRingBeam, separableRingBeam),
                               )
      k = abs(Ω / Va / 2)
      ωrs = (real(abs(vth * abs(k))), 0.9 * Ω, 2.5 * Ω, 5.5Ω, 10.5Ω) # real ωr must be > 0
      σs = (0, -1e-1, 1e-1, -1e-2, 1e-2, -1e-5, 1e-5, -1e-8, 1e-8, -1e-10, 1e-10)
      kzs = (2k, k/2, 0, -k/2, -2k)
      k⊥s = (k/2, 2k)
      for kz in kzs, σ ∈ σs, k⊥ in k⊥s, ωr in ωrs
        # this clearly isn't great, but how often is σ zero
        F = ComplexF64(ωr, σ * ωr)
        K = Wavenumber(kz=kz, k⊥=k⊥)
        iszero(K) && continue # kparallel and kperp cannot both be zero
        config = Configuration(F, K)
        config.options = Options(quadrature_rtol=1e-6, summation_rtol=1e-6)
        outputS = LMV.contribution(separable, config)
        config.options = Options(quadrature_rtol=1e-6, summation_rtol=1e-6, cubature_rtol=1e-6)
        outputC = LMV.contribution(coupled, config)

        @test separable(0.0, 0.0) ≈ coupled(0.0, 0.0)
        push!(norms, (norm(outputC - outputS)) / norm(outputS))

        atol=10eps() * norm(outputS)
        @testset "$kz, $σ, $k⊥, $M, $(ωr/Ω)" begin
          if verbose
            @show Float16.(abs.(outputC ./ outputS .- 1)), ωr, sign(kz)
          else
            @test outputC[1,1]≈outputS[1,1] rtol=rtol atol=atol
            @test outputC[1,2]≈outputS[1,2] rtol=rtol atol=atol
            @test outputC[1,3]≈outputS[1,3] rtol=rtol atol=atol
            @test outputC[2,1]≈outputS[2,1] rtol=rtol atol=atol
            @test outputC[2,2]≈outputS[2,2] rtol=rtol atol=atol
            @test outputC[2,3]≈outputS[2,3] rtol=rtol atol=atol
            @test outputC[3,1]≈outputS[3,1] rtol=rtol atol=atol
            @test outputC[3,2]≈outputS[3,2] rtol=rtol atol=atol
            @test outputC[3,3]≈outputS[3,3] rtol=rtol atol=atol
          end
        end
      end
    end
    verbose && @show mean(norms)
    verbose && @show std(norms)
  end
end

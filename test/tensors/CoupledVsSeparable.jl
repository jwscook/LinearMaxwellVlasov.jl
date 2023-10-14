using Dates
println("$(now()) $(@__FILE__)")

using Test, Random, DualNumbers, SpecialFunctions, LinearAlgebra
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

Random.seed!(0)

@testset "Separable vs Coupled velocity tensors" begin
  mₑ = LMV.mₑ
  mi = 1836*mₑ
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

    for (coupled, separable) ∈ (
                                (coupledMaxwellian, separableMaxwellian),
                                #(coupledRingBeam, separableRingBeam),
                               )
      k = abs(Ω / Va / 2)
      ωr = real(abs(vth * abs(k))) # real ωr must be > 0
      σs = (0, -1, 1)
      kzs = (2k, k/2, 0, -k/2, -2k)
      k⊥s = (k/2, 2k)
      for σ ∈ σs, kz in kzs, k⊥ in k⊥s
        # this clearly isn't great, but how often is σ zero
        rtol = iszero(σ) ? 1e-2 : 1e-5
        F = ComplexF64(ωr, σ * ωr / 100)
        K = Wavenumber(kz=kz, k⊥=k⊥)
        iszero(K) && continue
        config = Configuration(F, K)
        config.options = Options(quadrature_rtol=1.0e-15, summation_rtol=4eps())
        outputS = LMV.contribution(separable, config)
        config.options = Options(quadrature_rtol=1.0e-6, summation_rtol=1e-7)
        outputC = LMV.contribution(coupled, config)

        @test separable(0.0, 0.0) ≈ coupled(0.0, 0.0)

        atol=10eps() * norm(outputS)
        for op in (identity, )#real, imag, abs)
          @testset "$M, $σ, $kz, $k⊥, $op" begin
            @test op(outputC[1,1])≈op(outputS[1,1]) rtol=rtol atol=atol
            @test op(outputC[1,2])≈op(outputS[1,2]) rtol=rtol atol=atol
            @test op(outputC[1,3])≈op(outputS[1,3]) rtol=rtol atol=atol
            @test op(outputC[2,1])≈op(outputS[2,1]) rtol=rtol atol=atol
            @test op(outputC[2,2])≈op(outputS[2,2]) rtol=rtol atol=atol
            @test op(outputC[2,3])≈op(outputS[2,3]) rtol=rtol atol=atol
            @test op(outputC[3,1])≈op(outputS[3,1]) rtol=rtol atol=atol
            @test op(outputC[3,2])≈op(outputS[3,2]) rtol=rtol atol=atol
            @test op(outputC[3,3])≈op(outputS[3,3]) rtol=rtol atol=atol
          end
        end
      end
    end
  end
end

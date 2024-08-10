using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


@testset "Maxwellian: Classical vs Relativistic tensors" begin
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  B0 = 1.0
  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Πe = plasmafrequency(n0, mₑ, -1)
  ϵV = 1.0 # 1 eV is definitely a classical temperature
  vthe = thermalspeed(ϵV, mₑ)
  pthe = thermalmomentum(ϵV, mₑ)
  λD = vthe / Πe

  relativistic = CoupledRelativisticSpecies(Πe, Ωe, mₑ, pthe)
  classical = MaxwellianSpecies(Πe, Ωe, vthe, vthe)

  k = 2π/λD / 2
  ω = abs(vthe * k) # + Ωe
  for σ ∈ (0, -1, 1), ϕ ∈ (-1, 1)
    F = ComplexF64(ω, σ * ω / 100)
    K = Wavenumber(k=ϕ * k, θ=π/4)
    options = Options(summation_rtol=1e-3, quadrature_rtol=1e-6)
    config = Configuration(F, K, options)
    outputC = LMV.contribution(classical, config)
    outputR = LMV.contribution(relativistic, config)

    @testset "Inferred" begin
      loop = zip((classical, relativistic), ("classical", "relativistic"))
      for (species, name) ∈ loop
        try
          @inferred LMV.contribution(species, config, 0)
          @test true
        catch
          @warn "contribution not inferred for $(nameof(typeof(species)))"
          @test_broken false # return type of HCubature not inferrable
        end
        try
          @inferred LMV.contribution(species, config)
          @test true
        catch
          @warn "contribution not inferred for $(nameof(typeof(species)))"
          @test_broken false # return type of HCubature not inferrable
        end
      end
    end

    rtol=1.0e-3
    atol=1.0e-5

    @testset "real sign growth rate $σ, sign wavenumber $ϕ" begin
      @test real(outputC[1,1])≈real(outputR[1,1]) rtol=rtol atol=atol
      @test real(outputC[1,2])≈real(outputR[1,2]) rtol=rtol atol=atol
      @test real(outputC[1,3])≈real(outputR[1,3]) rtol=rtol atol=atol
      @test real(outputC[2,1])≈real(outputR[2,1]) rtol=rtol atol=atol
      @test real(outputC[2,2])≈real(outputR[2,2]) rtol=rtol atol=atol
      @test real(outputC[2,3])≈real(outputR[2,3]) rtol=rtol atol=atol
      @test real(outputC[3,1])≈real(outputR[3,1]) rtol=rtol atol=atol
      @test real(outputC[3,2])≈real(outputR[3,2]) rtol=rtol atol=atol
      @test real(outputC[3,3])≈real(outputR[3,3]) rtol=rtol atol=atol
    end
    @testset "imag sign growth rate $σ, sign wavenumber $ϕ" begin
      @test imag(outputC[1,1])≈imag(outputR[1,1]) rtol=rtol atol=atol
      @test imag(outputC[1,2])≈imag(outputR[1,2]) rtol=rtol atol=atol
      @test imag(outputC[1,3])≈imag(outputR[1,3]) rtol=rtol atol=atol
      @test imag(outputC[2,1])≈imag(outputR[2,1]) rtol=rtol atol=atol
      @test imag(outputC[2,2])≈imag(outputR[2,2]) rtol=rtol atol=atol
      @test imag(outputC[2,3])≈imag(outputR[2,3]) rtol=rtol atol=atol
      @test imag(outputC[3,1])≈imag(outputR[3,1]) rtol=rtol atol=atol
      @test imag(outputC[3,2])≈imag(outputR[3,2]) rtol=rtol atol=atol
      @test imag(outputC[3,3])≈imag(outputR[3,3]) rtol=rtol atol=atol
    end
  end
end

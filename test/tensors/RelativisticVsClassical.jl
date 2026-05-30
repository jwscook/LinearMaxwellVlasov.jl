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
  first = true
  for Z in (1, -1)
    Ωe = cyclotronfrequency(B0, mₑ, Z)
    Πe = plasmafrequency(n0, mₑ, Z)
    for ϵV in (1e0, 1e3)
      @testset "charge sign $Z, temparature $ϵV eV" begin
        vthe = thermalspeed(ϵV, mₑ)
        pthe = thermalmomentum(ϵV, mₑ)
        λD = vthe / Πe

        coupledrelativistic = CoupledRelativisticSpecies(Πe, Ωe, mₑ, pthe)
        maxwellainrelativistic = MaxwellJuttnerSpecies(Πe, Ωe, mₑ, ϵV)
        classical = MaxwellianSpecies(Πe, Ωe, vthe, vthe)

        k = 2π/λD / 2 / 10
        ω = abs(vthe * k) # + Ωe
        for σ ∈ (0, 1, -1), ϕ ∈ (1, -1, 0)
          iszero(σ) && iszero(ϕ) && continue
          @testset "sign growth rate $σ, sign wavenumber $ϕ" begin
            F = ComplexF64(ω, σ * ω / 100)
            K = Wavenumber(ϕ * k, k)
            config = Configuration(F, K)
            config.options = Options(summation_rtol=1e-8, quadrature_rtol=1e-8)
            tc = @elapsed outputC = LMV.contribution(classical, config)
            tm = @elapsed outputM = if σ < 0 && iszero(ϕ)
              @test_throws ArgumentError LMV.contribution(maxwellainrelativistic, config)
              nothing
            else
              LMV.contribution(maxwellainrelativistic, config)
            end
            config.options = Options(quadrature_rtol=1e-4, summation_rtol=1e-4, cubature_rtol=1e-5)
            tp = @elapsed outputP = LMV.contribution(coupledrelativistic, config)
            #@show tm / tc, tp / tc
            if first
              first = false
              @testset "Inferred" begin
                loop = zip((classical, coupledrelativistic, maxwellainrelativistic),
                           ("classical", "coupledrelativistic", "maxwellianrelativistic"))
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
            end

            rtol=1.0e-3
            atol=1.0e-5
            if ϵV == 1
              for (A, B, str) in ((outputC, outputP, "classical vs coupledrel"),
                                  (outputC, outputM, "classical vs maxwellianrel"),
                                  (outputP, outputM, "coupledrl vs maxwellianrel"))
                if isnothing(A) || isnothing(B)
                  continue
                end
                @testset "$str" begin
                  for op in (real, imag), i in 1:3, j in 1:3
                    @testset "$op $i $j" begin
                      @test op(A[i, j]) ≈ op(B[i, j]) rtol=rtol atol=atol
                    end
                  end
                end
              end
            else
              for (A, B, str) in ((outputP, outputM, "coupledrl vs maxwellianrel"),)
                if isnothing(A) || isnothing(B)
                  continue
                end
                @testset "$str" begin
                  for op in (real, imag), i in 1:3, j in 1:3
                    @testset "$op $i $j" begin
                      @test op(A[i, j]) ≈ op(B[i, j]) rtol=rtol atol=atol
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end

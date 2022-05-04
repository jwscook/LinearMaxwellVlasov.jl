using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


@testset "Plasmas" begin
  n0e = 1.0
  Ze = -1
  mₑ = LMV.mₑ
  for i in 1:1
    B = 10*(rand() - 0.5)
    mi = 1836 * mₑ
    Zi = rand(1:5)
    n0i = - Ze * n0e / Zi
    Va_ish = abs(B) / sqrt(LMV.μ₀ * mi * n0i)
    Πe = plasmafrequency(n0e, mₑ, Ze)
    Πi = plasmafrequency(n0i, mi, Zi)
    Ωe = cyclotronfrequency(B, mₑ, Ze)
    Ωi = cyclotronfrequency(B, mi, Zi)
    Se = ColdSpecies(Πe, Ωe)
    Si = ColdSpecies(Πi, Ωi)
    plasma = Plasma([Se, Si])
    @test length(plasma) == 2
    for (s, i) in enumerate(plasma)
      i == 1 && @test s == Se
      i == 2 && @test s == Si
    end

    @test isneutral(plasma)
    @test LMV.alfvenspeed(plasma) ≈ Va_ish rtol=0.01
    
      NeutralPlasma([Se, Si])
    try
      NeutralPlasma([Se, Si])
      @test true
    catch
      @test false
    end

    Se1 = ColdSpecies(Πe * (1 + 1e-8), Ωe)
    @test !isneutral(Plasma([Se1, Si]))
    Se2 = ColdSpecies(Πe, Ωe * (1 + 1e-8))
    @test !isneutral(Plasma([Se2, Si]))

    try
      NeutralPlasma([Se1, Si])
      @test false
    catch
      @test true
    end

    Si1 = ColdSpecies(Πi * (1 + 1e-8), Ωi)
    @test !isneutral(Plasma([Se, Si1]))
    Si2 = ColdSpecies(Πi, Ωi * (1 + 1e-8))
    @test !isneutral(Plasma([Se, Si2]))

  end 
end

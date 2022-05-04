using Dates
println("$(now()) $(@__FILE__)")

path = "../"
using Test, LinearAlgebra
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


verbose = false

function performtest()
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  B0 = 1.0
  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωi = cyclotronfrequency(B0, mi, 1)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πi = plasmafrequency(n0, mi, 1)
  ϵV = 1.0e3
  vthe = thermalspeed(ϵV, mₑ)
  vthi = thermalspeed(ϵV, mi)
  λD = vthe / Πe
  lri = vthi / Ωi
  Va = sqrt(B0^2/LMV.μ₀/n0/mi)
  electron_num = SeparableVelocitySpecies(Πe, Ωe,
    FParallelNumerical(vthe),
    FPerpendicularNumerical(vthe))
  proton_num = SeparableVelocitySpecies(Πi, Ωi,
    FParallelNumerical(vthi),
    FPerpendicularNumerical(vthi))
  electron_max = MaxwellianSpecies(Πe, Ωe, vthe, vthe)
  proton_max = MaxwellianSpecies(Πi, Ωi, vthi, vthi)
  electron_rb = RingBeamSpecies(Πe, Ωe, vthe, vthe)
  proton_rb = RingBeamSpecies(Πi, Ωi, vthi, vthi)
  S_num = Plasma([electron_num, proton_num])
  S_max = Plasma([electron_max, proton_max])
  S_rb = Plasma([electron_rb, proton_rb])

  @testset "Tensors: massive test over different configurations" begin
    count = 0
    freq₀ = abs(Ωi)
    ωs = Vector(range(0.01, stop=3.0, length=3)) * freq₀
    γs = Vector(range(-3.0, stop=3.0, length=3)) * freq₀ / 100
    for k_norm ∈ [-5.0, -1.0, -0.2, 0.2, 1.0, 5.0]
      k = 2π/lri*k_norm
      K = Wavenumber(k=k, θ=π/4)
      for (i, ω) ∈ enumerate(ωs), (j, γ) ∈ enumerate(γs)
        F = ComplexF64(ω, γ)
        C = Configuration(F, K)
        count += 1
        verbose && @show count, ω, γ, k
        output_num = det(tensor(S_num, C))
        output_max = det(tensor(S_max, C))
        output_rb = det(tensor(S_rb, C))
        abs(output_max) > 1.0e14 && continue
        @test output_num ≈ output_max rtol=0.01 atol=1.0e-4
        @test output_rb ≈ output_max rtol=0.1 atol=1.0e-4
      end
    end
  end

  @testset "electrostatic" begin
    k = 2π/lri
    K = Wavenumber(k=k, θ=0)
    C = Configuration((1.1 + 0.1*im)* Πe, K)
    output_tensor = tensor(S_num, C)[3,3]
    output = electrostaticdielectric(S_num, C)
    @test output_tensor ≈ output
  end

  @testset "two stream" begin
    Πe = plasmafrequency(1e19, mₑ, -1)
    vbeam = thermalspeed(1e3, mₑ)

    left = SeparableVelocitySpecies(Πe, eps(),
      FParallelDiracDelta(-vbeam), FPerpendicularDiracDelta(eps()))
    right = SeparableVelocitySpecies(Πe, eps(),
      FParallelDiracDelta(vbeam), FPerpendicularDiracDelta(eps()))

    config = Configuration(0 + im * Πe / 2,
                           Wavenumber(k=√3/2 * Πe / vbeam, θ=0.0))
    @test tensor(Plasma([left, right]), config)[3,3] ≈ 0 atol=eps()
    @test electrostaticdielectric(Plasma([left, right]), config) ≈ 0 atol=eps()
  end
end
performtest()

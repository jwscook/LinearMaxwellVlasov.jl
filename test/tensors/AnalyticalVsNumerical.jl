using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


@testset "Tensors: Analytical vs Numerical" begin
  mₑ = LMV.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  B0 = 1.0
  Ωi = cyclotronfrequency(B0, mi, 1)
  Πi = plasmafrequency(n0, mi, 1)
  ϵV = 1.0e3
  vthi = thermalspeed(ϵV, mi)
  lri = vthi / Ωi
  proton = SeparableVelocitySpecies(Πi, Ωi,
    FParallelNumerical(
      x->exp(-x^2/vthi^2), -12vthi, 12vthi),
    FPerpendicularNumerical(
      x->exp(-x^2/vthi^2), 0.0, 12vthi))
  S = proton
  # construct the same one for analytical purposes
  M = MaxwellianSpecies(Πi, Ωi, vthi, vthi)
  tols = Tolerance(rtol=1.0e-3, atol=0.0)
  O = Options(tols=tols)

  f0 = proton.Ω
  k0 = 1.0/lri
  N = 2
  rtol = 1.0e-3
  atol = 1.0e-3
  time_r, time_a, time_n = 0.0, 0.0, 0.0
  @testset "looping over k parallel, k perp, omega and gamma" begin
    ks = range(-3.0, stop=3.0, length=N) * k0
    θs = (π, 0.0, 3π/4, π/2, π/4, -π/4, -π/2, -3π/4)#range(-π, stop=π, length=5)
    ωs = range(0.1, stop=3.0, length=N) * f0
    γs = range(-0.5, stop=0.5, length=N) # not larger - could get nans
    for θ in θs, k in ks, ω in ωs, γ1 in γs
      iszero(k) && continue
      γ = γ1*ω
      K = Wavenumber(k=k, θ=θ)
      k⊥ = perp(K)
      kz = para(K)
#      z = k⊥ * vthi / Ωi
#      ζ = kz * vthi / ω
      F = ComplexF64(ω, γ)
      C = Configuration(F, K, O)
      time_a += @elapsed (ta = conductivity(M, C))
      time_n += @elapsed (tn = conductivity(S, C))
      @test all(isapprox.(tn, ta, rtol=rtol, atol=atol))
      for i in 1:3, j in 1:3
        if !isapprox(tn[i, j], ta[i, j], rtol=rtol, atol=atol)
          a, b = ta[i, j],  tn[i, j]
          @show i, j, θ / π, sign(kz), sign(k⊥), a, b
        end
      end
    end
  end

  @testset "Inferred" begin
    K = Wavenumber(k=50.0, θ=π/4)
    F = 1.0 + im
    C = Configuration(F, K)
    for Q ∈ (M, S)
      try
        @inferred conductivity(Q, C)
        @test true
      catch
        @warn "conductivity not inferred for $(nameof(typeof(Q)))"
        @test_broken false
      end
    end
  end
end

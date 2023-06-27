using Dates
println("$(now()) $(@__FILE__)")

using Base.Threads, Random, LinearAlgebra, Test
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

@testset "RingBeam vs Numerical ring beams and pancakes" begin
  Random.seed!(0) # seed rand
  verbose = false

  Ω    = rand()
  Π    = rand()
  v    = rand()
  vthb = rand()
  vth⊥ = rand()
  vb   = rand()
  v⊥   = rand()

  Va = 1.0
  f0 = abs(Ω)
  k0 = f0 / abs(Va)

  numerical_rb = SeparableVelocitySpecies(Π, Ω,
    FParallelNumerical(vthb, vb),
    FPerpendicularNumerical(vth⊥, v⊥))
  ringbeam_rb = SeparableVelocitySpecies(Π, Ω,
    FBeam(vthb, vb),
    FRing(vth⊥, v⊥))

  # beam pancake
  maxwellian_bp = MaxwellianSpecies(Π, Ω, vthb, vth⊥, vb)
  numerical_bp = SeparableVelocitySpecies(Π, Ω,
    FParallelNumerical(vthb, vb),
    FPerpendicularNumerical(vth⊥, 0.0))
  ringbeam_bp = SeparableVelocitySpecies(Π, Ω,
    FBeam(vthb, vb),
    FRing(vth⊥, 0.0))

  O = Options()

  ks_positive = range(0.1, stop=100.0, length=2^3) * k0
  ks = sort(vcat(-ks_positive, ks_positive))

  vals_numerical_rb = []
  vals_ringbeam_rb = []
  vals_maxwellian_bp = []
  vals_numerical_bp = []
  vals_ringbeam_bp = []
  configurations = []
  for k ∈ ks, m ∈ 0.2:0.4:2.0
    verbose && @show k, m
    K = Wavenumber(k=k, θ=π/4)
    ωi = abs(k * Va) * m
    γi = (rand() * 2 - 1) * abs(k * Va) / 100
    F = ComplexF64(ωi, γi)
    C = Configuration(F, K, O)
    push!(configurations, C)

    push!(vals_numerical_rb, conductivity(numerical_rb, C))
    push!(vals_ringbeam_rb, conductivity(ringbeam_rb, C))

    push!(vals_maxwellian_bp, conductivity(maxwellian_bp, C))
    push!(vals_numerical_bp, conductivity(numerical_bp, C))
    push!(vals_ringbeam_bp, conductivity(ringbeam_bp, C))
  end

  # these ones are different
  for k ∈ [-50, 50] .* k0, m ∈ 0.01:0.1:0.31
    K = Wavenumber(k=k, θ=π/4)
    ωi = abs(k * Va) * m
    γi = (rand() * 2 - 1) * abs(k * Va) / 100
    F = ComplexF64(ωi, γi)
    C = Configuration(F, K, O)
    push!(configurations, C)

    push!(vals_numerical_rb, conductivity(numerical_rb, C))
    push!(vals_ringbeam_rb, conductivity(ringbeam_rb, C))

    push!(vals_maxwellian_bp, conductivity(maxwellian_bp, C))
    push!(vals_numerical_bp, conductivity(numerical_bp, C))
    push!(vals_ringbeam_bp, conductivity(ringbeam_bp, C))
  end

  @testset "Inferred" begin
    K = Wavenumber(k=50.0, θ=π/4)
    F = ComplexF64(1.0 + im)
    C = Configuration(F, K, O)
    try
      @inferred conductivity(numerical_rb, C)
      @test true
    catch
      @warn "conductivity not inferred for $(nameof(typeof(numerical_rb)))"
      @test_broken false
    end
    try
      @inferred conductivity(ringbeam_rb, C)
      @test true
    catch
      @warn "conductivity not inferred for $(nameof(typeof(ringbeam_rb)))"
      @test_broken false
    end
  end

  verbose && @show length(vals_numerical_bp)

  rtol = sqrt(eps())
  atol = eps()
  function dotest(i, j, k, a, b)
    outcome = isapprox(a[k][i,j], b[k][i,j], rtol=rtol, atol=atol)
    @test isapprox(a[k][i,j], b[k][i,j], rtol=rtol, atol=atol)
    K = configurations[k].wavenumber
    if !outcome
      kz = para(K)
      k⊥ = perp(K)
      ω = configurations[k].frequency
      verbose && @show a[k][i,j], b[k][i,j]
      verbose && @show k, i, j, K.k, K.θ, kz, k⊥, ω
    end
  end

  @testset "RingBeam vs Numerical, ring beam" begin
    for i ∈ 1:3, j ∈ 1:3
      for k ∈ eachindex(vals_ringbeam_rb)
        dotest(i, j, k, vals_ringbeam_rb, vals_numerical_rb)
      end
    end
  end

  @testset "RingBeam vs Numerical, beam pancake" begin
    for i ∈ 1:3, j ∈ 1:3
      for k ∈ eachindex(vals_ringbeam_bp)
        dotest(i, j, k, vals_ringbeam_bp, vals_numerical_bp)
      end
    end
  end

  @testset "RingBeam vs Maxwellian, beam pancake" begin
    for i ∈ 1:3, j ∈ 1:3
      for k ∈ eachindex(vals_ringbeam_bp)
        dotest(i, j, k, vals_ringbeam_bp, vals_maxwellian_bp)
      end
    end
  end

end

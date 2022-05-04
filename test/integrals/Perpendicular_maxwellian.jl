using Dates
println("$(now()) $(@__FILE__)")

using Test, Random, QuadGK, SpecialFunctions
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


Random.seed!(0)
@testset "Perpendicular maxwellian integrals" begin
  # The missing factor of 2π in the expected values is hidden in F⊥
  J = besselj
  for i ∈ 1:10
    n = rand(-3:3)
    vth = rand()
    k⊥_Ω = rand()
    F = FRing(vth)
    mi⊥ = LMV.MaxwellianIntegralsPerpendicular(F.vth, k⊥_Ω, n)
    @testset "∫Jn²F⊥2πv⊥" begin
      result = LMV.∫Jn²F⊥2πv⊥(mi⊥)
      expected = quadgk(v -> 2π * v * F(v) * J(n, v * k⊥_Ω)^2, 0.0, 8 * vth,
        rtol=100eps(), atol=0.0)[1]
      @test result ≈ expected rtol=sqrt(eps()) atol=0.0
    end
    @testset "∫Jn²∂F⊥2π" begin
      result = LMV.∫Jn²∂F⊥2π(mi⊥)
      expected = quadgk(v -> 2π * F(v, true) * J(n, v * k⊥_Ω)^2, 0.0, 8 * vth,
        rtol=100eps(), atol=0.0)[1]
      @test result ≈ expected rtol=sqrt(eps()) atol=0.0
    end
    @testset "∫Jn∂JnF⊥2πv⊥²" begin
      result = LMV.∫Jn∂JnF⊥2πv⊥²(mi⊥)
      J0(v) = J(n, v * k⊥_Ω)
      J1(v) = J(n+1, v * k⊥_Ω)
      J_1(v) = J(n-1, v * k⊥_Ω)
      expected = quadgk(
        v -> 2π * v^2 * F(v, false) * J0(v) * (J_1(v) - J1(v)) / 2,
        0.0, 8 * vth, rtol=100eps(), atol=0.0)[1]
      @test result ≈ expected rtol=sqrt(eps()) atol=0.0
    end
    @testset "∫Jn∂Jn∂F⊥2πv⊥" begin
      result = LMV.∫Jn∂Jn∂F⊥2πv⊥(mi⊥)
      J0(v) = J(n, v * k⊥_Ω)
      J1(v) = J(n+1, v * k⊥_Ω)
      J_1(v) = J(n-1, v * k⊥_Ω)
      expected = quadgk(
        v -> 2π * v * F(v, true) * J0(v) * (J_1(v) - J1(v)) / 2,
        0.0, 8 * vth, rtol=100eps(), atol=0.0)[1]
      @test result ≈ expected rtol=sqrt(eps()) atol=0.0
    end
    @testset "∫∂Jn²F⊥2πv⊥³" begin
      result = LMV.∫∂Jn²F⊥2πv⊥³(mi⊥)
      J0(v) = J(n, v * k⊥_Ω)
      J1(v) = J(n+1, v * k⊥_Ω)
      J_1(v) = J(n-1, v * k⊥_Ω)
      expected = quadgk(
        v -> 2π * v^3 * F(v) * (J_1(v) - J1(v))^2 / 4,
        0.0, 8 * vth, rtol=100eps(), atol=0.0)[1]
      @test result ≈ expected rtol=sqrt(eps()) atol=0.0
    end
    @testset "∫∂Jn²∂F⊥2πv⊥²" begin
      result = LMV.∫∂Jn²∂F⊥2πv⊥²(mi⊥)
      J0(v) = J(n, v * k⊥_Ω)
      J1(v) = J(n+1, v * k⊥_Ω)
      J_1(v) = J(n-1, v * k⊥_Ω)
      expected = quadgk(
        v -> 2π * v^2 * F(v, true) * (J_1(v) - J1(v))^2 / 4,
        0.0, 8 * vth, rtol=100eps(), atol=0.0)[1]
      @test result ≈ expected rtol=sqrt(eps()) atol=0.0
    end

  end
end

@testset "Perpendicular Maxwellian Integrals in limit" begin
  vth = rand()*1e5
  Ω = rand()*1e9
  k⊥ = eps()
  k⊥_Ω = k⊥ / Ω
  for n in -2:2
    nΩ = n * Ω
    mi⊥ = LMV.MaxwellianIntegralsPerpendicular(vth, k⊥_Ω, n)

    ⊥2T11_2⊥2T1_1⊥2T_1_1 = 4 * LMV.∫∂Jn²∂F⊥2πv⊥²(mi⊥)
    @test isapprox(⊥2T11_2⊥2T1_1⊥2T_1_1, -2 * isone(abs(n)), atol=eps(), rtol=1e-5)

    ⊥3F11_2⊥3F1_1⊥3F_1_1 = 4 * LMV.∫∂Jn²F⊥2πv⊥³(mi⊥)
    @test isapprox(⊥3F11_2⊥3F1_1⊥3F_1_1, mi⊥.vth²₂ * 2 * isone(abs(n)),
      atol=eps(), rtol=1e-5)

    ⊥2T_1_1_⊥2T11 = 4nΩ / k⊥ * LMV.∫Jn∂Jn∂F⊥2πv⊥(mi⊥)
    @test isapprox(⊥2T_1_1_⊥2T11, (abs(n)==1) * sign(n) * (-2), atol=eps(), rtol=1e-5)

    ⊥3F_1_1_⊥3F11 = 4nΩ / k⊥ * LMV.∫Jn∂JnF⊥2πv⊥²(mi⊥)
    @test isapprox(⊥3F_1_1_⊥3F11, (abs(n)==1) * sign(n) * mi⊥.vth²₂ * 2, atol=eps(), rtol=1e-5)

    ⊥1T0_1⊥1T10 = 2nΩ / k⊥ * LMV.∫Jn²∂F⊥2π(mi⊥)
    @test isapprox(⊥1T0_1⊥1T10, 0.0, atol=eps(), rtol=1e-5)

    ⊥2F0_1⊥2F10 = 2nΩ / k⊥ * LMV.∫Jn²F⊥2πv⊥(mi⊥)
    @test isapprox(⊥2F0_1⊥2F10, 0.0, atol=100eps(), rtol=1e-5)

    ⊥2T112⊥2T1_1⊥2T_1_1 = 4nΩ^2 / k⊥^2 * LMV.∫Jn²∂F⊥2π(mi⊥)
    @test isapprox(⊥2T112⊥2T1_1⊥2T_1_1, isone(abs(n)) * (-2), atol=eps(), rtol=1e-5)

    ⊥3F112⊥3F1_1⊥3F_1_1 = 4nΩ^2 / k⊥^2 * LMV.∫Jn²F⊥2πv⊥(mi⊥)
    @test isapprox(⊥3F112⊥3F1_1⊥3F_1_1, isone(abs(n)) * mi⊥.vth²₂ * 2, atol=eps(), rtol=1e-5)
  end
end


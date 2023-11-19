using Dates
println("$(now()) $(@__FILE__)")

using Test, LinearAlgebra
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

const k0 = -1.234
@testset "Wavenumbers" begin
  @testset "parallel" begin
    for θ in (eps(), -eps())
      K = Wavenumber(wavenumber=k0, propagationangle=θ)
      @test para(K) ≈ k0
      @test perp(K) ≈ 0.0 atol=10eps()
    end
  end

  @testset "anti parallel" begin
    for θ in (π, -π)
      K = Wavenumber(wavenumber=k0, propagationangle=θ)
      @test para(K) ≈ -k0
      @test perp(K) ≈ 0.0 atol=10eps()
    end
  end

  @testset "perpendicular" begin
    for θ in (π/2, -π/2)
      K = Wavenumber(wavenumber=k0, propagationangle=θ)
      #@test para(K) ≈ 0.0
      @test perp(K) ≈ sign(θ) * k0
    end
  end

  @testset "errors" begin
    try
      Wavenumber(wavenumber=1, parallel=1)
      @test false
    catch
      @test true
    end
    try
      Wavenumber(asdfgsdga=1, apsngaegaerigbou=1)
      @test false
    catch
      @test true
    end
  end

  @testset "curl curl operator" begin
    for k1 in (k0, k0 + im)
      K = Wavenumber(wavenumber=k0, propagationangle=π/4)
      ∇x∇x = LMV.curlcurl(K)
      c = LMV.cartesian_vector(K)

      curl = im * [0.0 -c[3] c[2]; c[3] 0.0 -c[1]; [-c[2] c[1] 0.0]]
      @test all(curl * curl .== ∇x∇x)
    end
  end

  @testset "Costructor wavenumber propagationangle" begin
    K1 = Wavenumber(wavenumber=1.0, propagationangle=0.0)
    K2 = Wavenumber(wavenumber=1.0, propagationangle=0)
    K3 = Wavenumber(wavenumber=1.0+im, propagationangle=0.0)
    K4 = Wavenumber(wavenumber=1.0+im, propagationangle=0)
    @test true
  end

  @testset "Costructor propagationangle π/2" begin
    K1 = Wavenumber(wavenumber=1.0, propagationangle=π/2)
    K2 = Wavenumber(parallel=0.0, perpendicular=1.0)
    @test para(K1) == para(K2)
    @test perp(K1) == perp(K2)
  end

  @testset "angle" begin
    K = Wavenumber(parallel=1, perpendicular=0)
    @test angle(K) == 0
    K = Wavenumber(parallel=-1, perpendicular=0)
    @test angle(K) == Float64(π)
    K = Wavenumber(parallel=0, perpendicular=1)
    @test angle(K) == π/2
    K = Wavenumber(parallel=0, perpendicular=-1)
    @test angle(K) == -π/2
  end

  @testset "Minus" begin
    K12 = Wavenumber(parallel=1.0, perpendicular=2.0)
    K36 = Wavenumber(parallel=3.0, perpendicular=6.0)
    K = K36 - K12
    @test K == Wavenumber(parallel=2.0, perpendicular=4.0)
  end

  @testset "scalar value" begin
    K = Wavenumber(parallel=3.0, perpendicular=4.0)
    @test abs(K) ≈ 5.0
  end

  @testset "iszero" begin
    K0 = Wavenumber(parallel=0.0, perpendicular=0.0)
    @test iszero(K0)
    K1 = Wavenumber(parallel=1.0, perpendicular=1.0)
    @test !iszero(K1)
  end
end

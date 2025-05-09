using Dates
println("$(now()) $(@__FILE__)")

using LinearMaxwellVlasov
using Test, SpecialFunctions, QuadGK, HCubature, StaticArrays, Random
using DualNumbers, ForwardDiff

const LMV = LinearMaxwellVlasov

@testset "Maths core" begin
  Random.seed!(0)

  verbose = false

  tolerance = Tolerance()
  @testset "Plasma dispersion function" begin
    @test LMV.plasma_dispersion_function(0.0, 0) ≈ im*sqrt(pi) rtol=0.001
    @test LMV.plasma_dispersion_function(im, 0) ≈ im*0.757872156141312 rtol=0.001
    @test LMV.plasma_dispersion_function(ComplexF64(-1.52, 0.47), 0) ≈ ComplexF64(0.6088888957234254, 0.33494583882874024) rtol=0.001
    Z0 = LMV.plasma_dispersion_function(0.0, 0)
    Z1 = LMV.plasma_dispersion_function(0.0, 1)
    @test Z1 == LMV.plasma_dispersion_function(0.0, 1, Z0)
  end

  @testset "BesselJ generation function: used in derivation" begin
    for i ∈ 1:10, s ∈ (-1, 1)
      ϕ, ϕ₀, z, m = 2π .* (rand(3) .- 0.5)..., Int(round(20*(rand()-0.5)))

      cosexp = cos(m * ϕ) * exp(s * im * z * sin(ϕ))
      besselcosexpequiv = mapreduce(
        n -> 1/2 * (besselj(n - m, z) + besselj(n + m, z)) * exp(s * im * n * ϕ),
        +, -40:40)
      @test cosexp ≈ besselcosexpequiv rtol=sqrt(eps())

      sinexp = sin(m * ϕ) * exp(s * im * z * sin(ϕ))
      besselsinexpequiv = mapreduce(
        n -> -s*im/2 * (besselj(n - m, z) - besselj(n + m, z)) * exp(s * im * n * ϕ),
        +, -40:40)
      @test sinexp ≈ besselsinexpequiv rtol=sqrt(eps())
    end
  end

  @testset "transformaboutroots" begin
    vth = rand()*100
    bnd = 12*vth
    f(x) = exp(-x^2/2/vth^2) / sqrt(2 * π) / vth
    a, b = sort(2 .* (rand(2) .- 0.5) * vth)
    o = LMV.transformaboutroots(f, a)
    @test QuadGK.quadgk(o, -1, 1)[1] ≈ 1
    p = LMV.transformaboutroots(f, b)
    @test QuadGK.quadgk(p, -1, 1)[1] ≈ 1
    q = LMV.transformaboutroots(f, a, b)
    @test QuadGK.quadgk(q, -1, 1)[1] ≈ 1
  end

  dynamicrange = (1e-40, 1e-20, 1e-10, 1.0, 1e10, 1e20, 1e40)

  @testset "transformfrominfinity" begin
    for vth ∈ dynamicrange
      bnd = 12*vth
      f(x) = exp(-x^2/2/vth^2) / sqrt(2 * π) / vth
      @test QuadGK.quadgk(f, -bnd, bnd)[1] ≈ 1
      g = LMV.TransformFromInfinity(f, vth)
      @test QuadGK.quadgk(g, -1 + eps(), 1 - eps())[1] ≈ 1
    end
  end

  @testset "slow fourier transform" begin
    for n ∈ -5:5
      c = rand()
      f(x) = c * sin(n * x)
      fn = LMV.discretefouriertransform(f, n)
      @test c * im * (n != 0) ≈ fn atol=sqrt(eps()) rtol=sqrt(eps())
      g(x) = c * cos(n * x)
      gn = LMV.discretefouriertransform(g, n)
      @test c ≈ gn atol=sqrt(eps()) rtol=sqrt(eps())
    end
  end

  @testset "fold numerator about pole" begin
    vth = rand()*100
    bnd = 12*vth
    z = 3 * rand() * vth
    numerator(x) = exp.(-x.^2/2/vth^2) / sqrt(2 * π) / vth
    f(x) = numerator(x) ./ (x - z)
    g = LMV.foldnumeratoraboutpole(numerator, z)
    expected = QuadGK.quadgk(g, 2eps(), 12*vth)[1]
  end

  @testset "transform to polar" begin
    f(x) = exp(-sum(x.^2)) / π
    p = LMV.transformtopolar(f)
    n1 = HCubature.hcubature(f, [-6, -6], [6, 6], initdiv=3)[1]
    n2 = HCubature.hcubature(rθ -> rθ[1] * p(rθ), [0, -π], [6, π])[1]
    @test n1 ≈ 1
    @test n2 ≈ 1

    for i ∈ 1:10
      vb = rand() * 2 - 1
      v⊥ = rand()
      vbv⊥ = LMV.parallelperpfrompolar(LMV.polarfromparallelperp([vb, v⊥]))
      @test vb ≈ vbv⊥[1]
      @test v⊥ ≈ vbv⊥[2]
    end

    for i ∈ 1:10
      r = rand()
      θ = rand() * 2 * π - π
      rθ = LMV.parallelperpfrompolar(LMV.polarfromparallelperp([r, θ]))
      @test r ≈ rθ[1]
      @test θ ≈ rθ[2]
    end
  end

  @testset "unitsemicircleintegrandtransform" begin
    pth = 1e8 * 9.11e-31
    f(x) = 2π * x[2] * exp(-sum(x.^2) / pth^2) / sqrt(π)^3 / pth^3
    expected = HCubature.hcubature(f, pth .* [-12.0, 0.0], pth .* [12.0, 12.0])[1]
    g = LMV.UnitSemicircleIntegrandTransform(f, pth)
    result = HCubature.hcubature(g, [0, -π/2], [1, π/2])[1]
    @test expected ≈ result
    @test all(LMV.coordinates(g, [0.0, 0.0]) .== [0.0, 0.0])
    @test all(LMV.coordinates(g, [(sqrt(5)-1)/2, π/7]) .≈ [pth, π/7])
  end

  @testset "Ensure besselj composes with complex indices and Duals" begin
    @test DualNumbers.dualpart(besselj(1.0 + im, Dual(1.0, 1))) ≈
      (besselj(0+im, 1.0) - besselj(2.0+im, 1.0)) / 2
  end

  @testset "Duals vs ForwardDiff" begin
    for n in (3, 4, -3, -4), x in (-5.0, 5.0)
      @test DualNumbers.dualpart(besselix(n, Dual(x, 1))) ≈
        ForwardDiff.derivative(z->besselix(n, z), x)
      @test DualNumbers.dualpart(besselj(n, Dual(x, 1))) ≈
        ForwardDiff.derivative(z->besselj(n, z), x)
    end
    for n in (3, 4), x in (-5.0, 5.0)
      @assert ForwardDiff.derivative(z->besseli(n, z), x) ≈
        (besseli(n-1, x) + besseli(n+1, x)) / 2
      @test DualNumbers.dualpart(n^Dual(x, 1)) ≈
        ForwardDiff.derivative(z->n^z, x)
    end
  end

  @testset "Dual ^ Complex" begin
    a, b = rand(Float64), rand(ComplexF64)
    @test DualNumbers.realpart(DualNumbers.Dual(a, 1)^b) ≈ a^b
    @test DualNumbers.dualpart(DualNumbers.Dual(a, 1)^b) ≈ b * a^(b - 1)
    a, b = rand(ComplexF64), rand(ComplexF64)
    @test DualNumbers.realpart(DualNumbers.Dual(a, 1)^b) ≈ a^b
    @test DualNumbers.dualpart(DualNumbers.Dual(a, 1)^b) ≈ b * a^(b - 1)
  end

  @testset "isapproxinteger" begin
    @test LMV.isapproxinteger(1.0 + 0im, 0.0)
    @test LMV.isapproxinteger(1.0 + 0im, eps())
    @test LMV.isapproxinteger(1.0 + im * eps(), eps())
    @test LMV.isapproxinteger(0.0 + im * eps(), eps())
    @test !LMV.isapproxinteger(1.0 + 1im, eps())
    @test !LMV.isapproxinteger(0.0 + im * 2eps(), eps())
    @test !LMV.isapproxinteger(1.0 + im * 2eps(), eps())
    @test !LMV.isapproxinteger(2*(1 + 2eps()) + 0im, eps())
    @test LMV.isapproxinteger(2*(1 + eps()) + 0im, eps())
  end

  @testset "Concertinas.jl" begin
    Δ = 1e-14
    @testset "sin(x)/sin(x), Δ = $Δ" begin
      for i in 1:100
        L⊥ = rand() # doesnt matter!
        @test 2i ≈ HCubature.hcubature(LMV.ConcertinaSinpi((t, v⊥)->sinpi(t) / L⊥, (-i, i)), (Δ, 0.0), (1 - Δ, L⊥), rtol=1e-6)[1]
      end
    end

    for i in -4:-1:-9
      Δ = 10.0^i
      expected = (log(tan((1 - Δ) * π / 2)) - log(tan(Δ * π / 2))) / 2π
      L⊥ = rand() # doesnt matter!
      result = HCubature.hcubature(LMV.ConcertinaSinpi((t, v⊥) -> identity(t) / L⊥, (0, 1)), (Δ, 0.0), (1 - Δ, L⊥), rtol=1e-6)[1]
      @testset "1 / sin(x), Δ = $Δ" begin
        @test expected ≈ result
      end

      @testset "cos(x)^2 / sin(x), Δ = $Δ" begin
        fcos²_sin(x) = (cos(x) + log(sin(x/2)) - log(cos(x/2))) / π
        expected = fcos²_sin((1 - Δ) * π) - fcos²_sin(Δ * π)
        L⊥ = rand() # doesnt matter!
        result = HCubature.hcubature(LMV.ConcertinaSinpi((t, v⊥)->cospi(t)^2 / L⊥, (0, 1)), (Δ, 0.0), (1 - Δ, L⊥), rtol=1e-6)[1]
        @test expected ≈ result
      end
    end
  end
end

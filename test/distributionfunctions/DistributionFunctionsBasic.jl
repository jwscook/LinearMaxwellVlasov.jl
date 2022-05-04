using Dates
println("$(now()) $(@__FILE__)")

using SpecialFunctions, Test, ForwardDiff, QuadGK, Random
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

Random.seed!(0)

@testset "DistributionFunction" begin

@testset "shifted maxwellian parallel is normalised" begin
  for i ∈ 1:100
    vth, vd = rand()*10, (rand() - 0.5)*10
    result = QuadGK.quadgk(
      LMV.ShiftedMaxwellianParallel(vth, vd),
      -Inf, Inf, rtol=eps(), atol=0.0)[1]
    @test result ≈ 1.0 rtol=100*eps() atol=0.0
  end
end

@testset "shifted maxwellian perp is normalised" begin
  for i ∈ 1:100
    vth, vd = rand()*10, rand()*10
    i < 10 && (vd = 37 * vth + rand()*10)
    sm⊥ = LMV.ShiftedMaxwellianPerpendicular(vth, vd)
    result = QuadGK.quadgk(x -> 2*pi*x*sm⊥(x), max(0.0, vd - 12vth), vd + 12vth,
                           rtol=eps(), atol=0.0)[1]
    @test result ≈ 1.0 rtol=sqrt(eps()) atol=0.0
  end
end

@testset "shifted maxwellian derivatives" begin
  for c ∈ 1:2, i ∈ 1:10
    SMP = c == 1 ? LMV.ShiftedMaxwellianParallel :
      LMV.ShiftedMaxwellianPerpendicular
    vth, vd = rand()*10, (rand() - 0.5)*10
    c == 2 && (vd = abs(vd))
    sm = SMP(vth, vd)
    g(x::Number) = ForwardDiff.derivative(sm, x)
    for j ∈ 1:10
      v = rand()*6*vth - 3*vth - vd
      c == 2 && (v = abs(v))
      result = sm(v, true)
      @test result ≈ g(v) rtol=2*eps() atol=2*eps()
    end
  end
end

Fvbe = FParallelNumerical(Float64(π), -exp(1.0))
Fvpe = FPerpendicularNumerical(Float64(π), exp(1.0))

#Profile.clear()
#Profile.init(n=10^7)
#@profile for i in 1:1000 integrate(Fvbe) end
#Profile.print(format=:flat, combine=true, sortedby=:count)

@testset "distribution function para is normalised" begin
  @test LMV.is_normalised(Fvbe)
end
@testset "distribution function perp is normalised" begin
  @test LMV.is_normalised(Fvpe)
end

@testset "auto derivative against finite difference derivative" begin
  reltol = abstol = 1.0e-4
  @test LMV.is_normalised(Fvbe)
  @test LMV.is_normalised(Fvpe)

  F1 = FParallelNumerical(1.0, 1.0)
  gF1(x) = ForwardDiff.derivative(F1.F, x)
  tmp = 10.0.^range(-3, stop=3, length=1000)
  for v in sort(vcat(-tmp, [0.0], tmp))*10.0
    @test F1.dFdv(v) ≈ gF1(v) rtol=reltol atol=abstol
  end
  F2 = FParallelNumerical(x->exp.(-(x+2).^2), -10.0, 10.0)
  gF2(x) = ForwardDiff.derivative(F2.F, x)
  for v in sort(vcat(-tmp, [0.0], tmp))*10.0
    @test F2.dFdv(v) ≈ gF2(v) rtol=reltol atol=abstol
  end
end

@testset "integration by parts" begin
  reltol = abstol = 1.0e-4
  F3 = FParallelNumerical(Float64(π), exp(1.0))
  integrand(x) = x
  dintegranddx(x) = 1.0
  a = QuadGK.quadgk(x->F3.F(x)*dintegranddx(x), F3.lower, F3.upper, rtol=eps())
  b = QuadGK.quadgk(x->-F3.dFdv(x)*integrand(x), F3.lower, F3.upper, rtol=eps())
  @test a[1] ≈ b[1] rtol=reltol atol=abstol
end

end

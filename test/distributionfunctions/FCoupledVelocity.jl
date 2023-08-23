using Dates
println("$(now()) $(@__FILE__)")

using Test, HCubature, SpecialFunctions, Random
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

Random.seed!(0)

@testset "Test FCoupledVelocityNumerical basic interface" begin
  smc = LMV.FCoupledVelocityNumerical(1.0, 2.0).F
  @test smc((1.0, 2.0)) == smc(1.0, 2.0)
end
@testset "Test CoupledVelocityNumerical Shell:" begin


@testset "is normalised" begin
  for i ∈ 1:10
    vth = 2 * rand()
    vshell = 2 * rand()
    Lshell = vshell + 6*vth
    shell = FShell(vth, vshell)
    vbeam = 10 * rand() + 1
    vcrit = rand()
    vcutoffwidth = rand()
    Lsd = (vbeam + vcutoffwidth * 6)
    slowingdown = FSlowingDown(vbeam, vcrit, vcutoffwidth)
    for (S, L) in ((shell, Lshell), (slowingdown, Lsd))
      integral = HCubature.hcubature(x->2π * x[2] * S(x), [-L, 0], [L, L],
                                     rtol=1.0e-8, initdiv=3)[1]
      @test integral ≈ 1 rtol=1.0e-3

      S1 = LinearMaxwellVlasov.FCoupledVelocityNumerical(S.F, S.normalisation;
                                                         autonormalise=true)
      integral = HCubature.hcubature(x->2π * x[2] * S1(x), [-L, 0], [L, L],
                                     rtol=1.0e-8, initdiv=3)[1]
      @test integral ≈ 1 rtol=1.0e-6
    end
  end
end


@testset "integrals in polar coords" begin
  cartesiancounter = 0
  polarcounter = 0
  for i ∈ 1:10
    vth = 2 * rand()
    vshell = 2 * rand()
    Lshell = vshell + 6*vth
    shell = FShell(vth, vshell)
    vbeam = 10 * rand() + 1
    vcrit = rand()
    vcutoffwidth = rand()
    Lsd = (vbeam + vcutoffwidth * 6)
    slowingdown =FSlowingDown(vbeam, vcrit, vcutoffwidth)
    for (S, L) in ((shell, Lshell), (slowingdown, Lsd))
      n, l = rand(-4:4), rand(-4:4)
      f(x) = 2π * x[2] * x[1]^2 * besselj(n, x[2]) * besselj(l, x[2]) * S(x)
      fcounter = 0
      function fclosed(x)
        fcounter += 1
        return f(x)
      end
      a = HCubature.hcubature(fclosed, [-L, 0], [L, L],
                              rtol=1.0e-5, initdiv=3)[1]
      g(vrθ) = vrθ[1] * LMV.transformtopolar(f)(vrθ)
      gcounter = 0
      function gclosed(x)
        gcounter += 1
        return g(x)
      end
      b = HCubature.hcubature(gclosed, [0, -π/2], [L, π/2],
                              rtol=1.0e-5, initdiv=3)[1]
      h(vrθ) = vrθ[1] * f(LMV.parallelperpfrompolar(vrθ))
      c = HCubature.hcubature(h, [0, -π/2], [L, π/2],
                              rtol=1.0e-5, initdiv=3)[1]
      @test a ≈ b rtol=1.0e-3
      @test b ≈ c rtol=1.0e-3
      cartesiancounter += fcounter
      polarcounter += gcounter
    end
  end
  @test polarcounter < cartesiancounter
end

end

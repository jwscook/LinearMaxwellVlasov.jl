using Dates, Test
println("$(now()) $(@__FILE__)")

using Base.Threads, Random, LinearAlgebra
using LinearMaxwellVlasov

const DPR = LinearMaxwellVlasov

Random.seed!(0) # seed rand

function run()
  mα = 4*1836*LMV.mₑ
  vα = thermalspeed(3.52e6, mα)
  θpseudopitch = -135.0 * π/180
  vα⊥ = vα * abs(sin(θpseudopitch))
  vαb = vα * cos(θpseudopitch)
  vαth = thermalspeed(1.0e5, mα)
  Πα = 1.0 # doesn't matter
  (t, Ωα, n, ω, KPara, ∂F∂v, tolrel, tolabs) = (5.584536409, 1.0058612261143655e8, 21, 2.0596206058532245e9 - 1.85666570697272im, 5.17591148221635, false, 1.0e-14, 2.220446049250313e-16)
  power = Unsigned.([0, 1, 2])

  tol = LMV.Tolerance(tolrel, tolabs)
  alpha_cold = ColdSpecies(Πα, Ωα)
  alpha_hot = SeparableVelocitySpecies(Πα, Ωα,
    FParallelNumerical(vαth, vαb),
    FPerpendicularNumerical(vαth, vα⊥))

  output = [LMV.parallel(alpha_hot.Fz, ω, KPara, n, Ωα, p, ∂F∂v, tol)
             for p in power]
  @test all(isfinite.(output))
end
@testset "Isfinite" begin
run()
end

using Dates
println("$(now()) $(@__FILE__)")

using Test, SpecialFunctions, QuadGK
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov


@testset "Parallel unittest 1" begin
  mₑ = LMV.mₑ
  md = 2*1836*mₑ
  mα = 2*md
  n0 = 1.0e20
  B0 = 1.0
  ξ = 1.0e-3
  nα = ξ*n0 / (1.0 + 2*ξ)
  θ = 85.0/90.0 * π/2 # 80.0/90 * π/2
  
  Ωα = cyclotronfrequency(B0, mα, 2)
  Πα = plasmafrequency(nα, mα, 2)
  ϵV = 1.0e3
  vα = thermalspeed(3.52e6, mα)
  vα⊥ = vα / sqrt(2.0)
  vαb = - vα / sqrt(2.0)
  
  fpara = FParallelNumerical(vα/100, vαb)
  fperp = FPerpendicularNumerical(vα/100, vα⊥)
  alpha_hot = SeparableVelocitySpecies(Πα, Ωα, fpara, fperp)
  
  (Π, h, ω, KPara, pow, diffbool) = (4.15931043118463e8, 2, 6.635776485087542e7 + 1.5664971352130033e7im, 3.086205861011212 + 0.0im, UInt64[0x0000000000000000, 0x0000000000000001, 0x0000000000000002], false)
  
  Fz = FParallelNumerical(vα/100, vαb)
  pole = -9.538748782126304e6 + 5.075802476441841e6im
  
  function drifting_maxwellian(v::T, vth::Real, vd::Real) where {T}
    big = exp.(BigFloat(real(-0.5*((v-vd)/vth).^2)))
    imaginary = exp.(im*imag(-0.5*((v-vd)/vth).^2))
    return (Float64(big) + imaginary) / (vth * sqrt(2*pi))
    #return exp.(-0.5*((v-vd)/vth).^2 - log(vth * sqrt(2*pi)))
  end
  vIb0 = parallel(Fz, ω, KPara, h, Ωα, Unsigned.(0), false)
  vIb1 = parallel(Fz, ω, KPara, h, Ωα, Unsigned.(1), false)
  vIb2 = parallel(Fz, ω, KPara, h, Ωα, Unsigned.(2), false)
  @test isfinite(vIb0)
  @test isfinite(vIb1)
  @test isfinite(vIb2)
end

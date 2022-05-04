using Dates
println("$(now()) $(@__FILE__)")

using Test
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


@testset "FRings" begin
  mi = 1836*LMV.mₑ
  vth = thermalspeed(1.0e4, mi)
  vd = thermalspeed(1.0e5, mi)

  ring = FRing(vth, vd)
  @test LMV.is_normalised(ring)
end

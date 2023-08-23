using Dates
println("$(now()) $(@__FILE__)")

using Test, HCubature
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

@testset "Relativity tests" begin

# commented out while testing other things
@testset "Relativistic Maxwellian is normalised" begin
  m0 = LMV.mₑ
  ϵV = 1e3
  pth = LMV.thermalmomentum(ϵV, m0)
  f = LMV.RelativisticMaxwellian(pth)
  limit = 1 .- 100eps()
  normalisation = HCubature.hcubature(p -> 2π * p[2] * f(p),
    [-limit, eps()] .* 12 * pth, [limit, limit] .* 12 * pth,
    rtol=1e4*eps(), atol=0.0)[1]
  @test normalisation ≈ 1 rtol=1.0e-3

  f1 = LMV.TransformFromInfinity(p -> 2π * p[2] * f(p), [pth, pth])
  normalisation = HCubature.hcubature(f1,
    [-limit, eps()], [limit, limit],
    rtol=1e4*eps(), atol=0.0)[1]
  @test normalisation ≈ 1 # rtol=1.0e-6
end

end

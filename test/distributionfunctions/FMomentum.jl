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
  f = LMV.MaxwellJuttner(m0, pth)
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

@testset "Relativistic Maxwellian approaches classical" begin
  m0 = LMV.mₑ
  ϵV = 1.0
  pth = LMV.thermalmomentum(ϵV, m0)
  f = LMV.MaxwellJuttner(m0, pth)
  fcz = LMV.FBeam(pth / m0)
  fc⊥ = LMV.FRing(pth / m0)
  for i in 1:100
    pz, p⊥ = randn(2) .* pth
    result = f((pz, p⊥))
    expected = fc⊥(p⊥ / m0) * fcz(pz / m0) / m0^3
    @test result ≈ expected rtol=1e-4
  end
end


end

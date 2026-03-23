using Dates
println("$(now()) $(@__FILE__)")

using SpecialFunctions, Test, Random, QuadGK
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov


Random.seed!(0)

@testset "DistributionFunction FParallel" begin

@testset "FParallelNumerical" begin
  vf = 5.360866737263309e6
  lv = 18vf
  fvz⁴(vz) = (vz/vf)^4 * exp(-vz^2/vf^2) # obviously even
  nvz⁴ = QuadGK.quadgk(fvz⁴, -lv, lv)[1]
  fvz = FParallelNumerical(vz->fvz⁴(vz) / nvz⁴, -lv, lv)

  for power in 0:2, ∂F∂v in (false, true)
    kernel = LMV.ParallelKernelNumerator(Unsigned(power))
    result = LMV.integrate(fvz, kernel, ∂F∂v)
    if (power == 0) && !∂F∂v
      @test result ≈ -1
    elseif isodd(power + ∂F∂v)
      @test result < 1e-6
    end
  end

end

end



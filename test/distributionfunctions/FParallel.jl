using Dates
println("$(now()) $(@__FILE__)")

using SpecialFunctions, Test, Random
using LinearMaxwellVlasov


Random.seed!(0)

@testset "DistributionFunction FParallel" begin

@testset "shifted maxwellian implementations are consistent" begin
  for i ∈ 1:100
    vth, vd = rand()*10, (rand() - 0.5)*10
    x = (rand() - 0.5) * 2 * 3 * vth
    a = shifted_maxwellian(x, vth, vd)
    b = shifted_maxwellian(vth, vd)(x)
    @test a ≈ b rtol=100*eps() atol=0.0
  end
end


@testset "shifted maxwellian derivative implementations are consistent" begin
  for i ∈ 1:100
    vth, vd = rand()*10, (rand() - 0.5)*10
    x = (rand() - 0.5) * 2 * 3 * vth
    a = shifted_maxwellian_derivative(x, vth, vd)
    b = shifted_maxwellian_derivative(vth, vd)(x)
    @test a ≈ b rtol=100*eps() atol=0.0
  end
end

end



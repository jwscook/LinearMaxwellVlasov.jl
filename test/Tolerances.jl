using Dates
println("$(now()) $(@__FILE__)")

using LinearMaxwellVlasov
using Test


@testset "Tolerances" begin
  @test Tolerance(Float64) == Tolerance(0.0, eps()^(3//4))
  tol = Tolerance(Float64)
  @test LinearMaxwellVlasov.uniqueid(tol) == tol._uniqueid
end

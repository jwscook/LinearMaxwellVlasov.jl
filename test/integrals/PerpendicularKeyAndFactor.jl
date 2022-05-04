using Dates
println("$(now()) $(@__FILE__)")

using Test, Random, SpecialFunctions

Random.seed!(0)

@testset "Testing perp memoisation logic" begin
  for i in -3:3, j in -3:3, x in (-1.0, 1.0)
    expected = besselj(i, x) * besselj(j, x)
    factor = sign(i)^i * sign(j)^j
    a, b = sort([abs(i), abs(j)])
    factor *= x < 0 && isodd(a + b) ? -1 : 1
    y = abs(x)
    result = factor * besselj(a, y) * besselj(b, y)
    @test expected == result
  end
end

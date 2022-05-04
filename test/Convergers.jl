using Dates
println("$(now()) $(@__FILE__)")

using LinearMaxwellVlasov: converge
using Test

@testset "Convergers" begin
  function f(n::Int)
    k = Float64(abs(n))
    return iszero(n) ? 2.0 : Float64(k)^-k
  end
  expected = sum([f(i) for i ∈ -13:13])
  result = converge(f)
  @test expected ≈ result rtol=sqrt(eps())
  result = converge(x->[f(x)])[1]
  @test expected ≈ result rtol=sqrt(eps())

  function g(n::Int)
    k = Float64(abs(n))
    f = n < 0 ? 2 : 1
    return iszero(n) ? 2.0 : Float64(k)^-(f * k)
  end
  expected = sum([g(i) for i ∈ -20:20])
  result = converge(g)
  @test expected ≈ result rtol=sqrt(eps())
  result = converge(x->[g(x)])[1]
  @test expected ≈ result rtol=sqrt(eps())
end



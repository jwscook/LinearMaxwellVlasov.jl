using Dates
println("$(now()) $(@__FILE__)")

using LinearMaxwellVlasov: converge, fastisapprox
using Test, LinearAlgebra

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

  @testset "fastisapprox" begin
    yeses = 0
    count = 0
    for log10rtol in 0:5, log10atol in 0:5
      for _ in 1:100
        rtol = 10.0.^log10rtol * eps()
        atol = 10.0.^log10atol * eps()
        a = rand(ComplexF64, 3, 3)
        b = deepcopy(a)
        # the following means that approx half instances are approximately equal
        b .+= rtol .* randn(ComplexF64, 3, 3) * norm(a) / 2
        b .+= randn(ComplexF64, 3, 3) .* atol / 7
        expected = isapprox(a, b; atol=atol, rtol=rtol, nans=true)
        result = fastisapprox(a, b; atol=atol, rtol=rtol, nans=true)
        @test result == expected
        yeses += expected
        count += 1
      end
    end
    @assert 0.25 <= yeses / count <= 0.75
  end
end



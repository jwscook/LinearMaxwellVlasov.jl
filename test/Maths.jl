using Dates
println("$(now()) $(@__FILE__)")

using Test
@testset "Maths tests" begin
  include("maths/Maths.jl")
  include("maths/Maths_unittest0.jl")
end

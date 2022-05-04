using Dates
println("$(now()) $(@__FILE__)")

using Test
@testset "tensors/" begin
  include("tensors/ColdLimit.jl")
  include("tensors/ColdVsDelta.jl")
  include("tensors/NumericalVsMaxwellian.jl")
  include("tensors/AnalyticalVsNumerical.jl")
  include("tensors/NegativeKPara.jl")
  include("tensors/Basic.jl")
  include("tensors/RingBeamsAndPancakes.jl")
  include("tensors/CoupledVsSeparable.jl")
  include("tensors/RelativisticVsClassical.jl")
end

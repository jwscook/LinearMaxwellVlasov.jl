using Dates
println("$(now()) $(@__FILE__)")

using Test
@testset "Integrals tests" begin
  include("integrals/DiracDelta.jl")
  include("integrals/Inferred.jl")
  include("integrals/Integrals.jl")
  include("integrals/Parallel.jl")
  include("integrals/ParallelBasic.jl")
  include("integrals/ParallelSignOfKpara.jl")
  include("integrals/Parallel_unittest1.jl")
  include("integrals/Parallel_unittest2.jl")
  include("integrals/Parallel_unittest3.jl")
  include("integrals/Perpendicular.jl")
  include("integrals/PerpendicularKeyAndFactor.jl")
  include("integrals/Perpendicular_maxwellian.jl")
  include("integrals/Memoisation.jl")
end

using LinearMaxwellVlasov
using Test
using Dates
println("$(now()) $(@__FILE__)")

# in order of inter-dependency
@testset "LinearMaxwellVlasov" begin
  include("./Tolerances.jl")
  include("./Convergers.jl")
  include("./Maths.jl")
  include("./Wavenumbers.jl")
  include("./Poles.jl")
  include("./Options.jl")
  include("./Parameters.jl")
  include("./Configurations.jl")
  include("./DistributionFunctions.jl")
  include("./Integrals.jl")
  include("./Species.jl")
  include("./Plasmas.jl")
  include("./Tensors.jl")
end



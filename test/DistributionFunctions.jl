using Dates
println("$(now()) $(@__FILE__)")

include("distributionfunctions/DistributionFunctionsBasic.jl")
include("distributionfunctions/Perpendicular.jl")
include("distributionfunctions/FMomentum.jl")
include("distributionfunctions/FBeams.jl")
include("distributionfunctions/FRings.jl")
include("distributionfunctions/FWideRings.jl")
include("distributionfunctions/FCoupledVelocity.jl")

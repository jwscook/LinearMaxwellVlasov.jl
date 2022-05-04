module LinearMaxwellVlasov

const q₀ = 1.602176487e-19 # absolute value of electron change, SI units [C]
const ϵ₀ = 8.854187817e-12 # electric constant / permitivitty, SI units [F/m]
const mₑ = 9.10938188e-31 # electron mass, SI units [kg]
const c₀ = 2.99792458e8 # speed of light, SI units [m/s]
const μ₀ = 1.25663706144e-6 # magnetic constant / permeability, SI units [N/A²]

# in order of inter-dependency
include("./Tolerances.jl")
include("./Convergers.jl")
include("./Maths.jl")
include("./Wavenumbers.jl")
include("./Poles.jl")
include("./Options.jl")
include("./Parameters.jl")
include("./Configurations.jl")
include("./DistributionFunctions.jl")
include("./Species.jl")
include("./Integrals.jl")
include("./Plasmas.jl")
include("./Tensors.jl")

# Tolerances
export Tolerance
# Convergers
# Maths
# Wavenumbers
export Wavenumber, parallel, perpendicular, para, perp
# Poles
export Options
# Parameters
export cyclotronfrequency, plasmafrequency, thermalspeed, thermalmomentum
export magnetoacousticfrequency, shearfrequency
export slowzerobetamagnetoacousticfrequency, fastzerobetamagnetoacousticfrequency
export slowmagnetoacousticfrequency, fastmagnetoacousticfrequency
# Configurations
export Configuration
# DistributionFunctions
export FRing, FBeam, FParallelNumerical, FPerpendicularNumerical
export FPerpendicularMaxwellian, FRelativisticNumerical
export FCoupledVelocityNumerical, FShell
export FParallelDiracDelta, FPerpendicularDiracDelta
# Species
export ColdSpecies, WarmSpecies, MaxwellianSpecies, RingBeamSpecies
export SeparableVelocitySpecies, CoupledVelocitySpecies, CoupledRelativisticSpecies
# Integrals
export Cache
# Plasmas
export Plasma, NeutralPlasma, isneutral
# Tensors
export tensor, dielectric, electrostaticdielectric

__precompile__(true) # precompile dependencies
# precompile for species types
for s in (
    ColdSpecies{Float64,Float64},
    SeparableVelocitySpecies{Float64,Float64,FBeam,FRing},
    CoupledVelocitySpecies{Float64,Float64,FCoupledVelocityNumerical{ShiftedMaxwellianCoupled{Float64,Float64,Float64,Float64},Float64}},
    CoupledRelativisticSpecies{Float64,Float64,Float64,FRelativisticNumerical{RelativisticMaxwellian,Float64}})
  precompile(dielectric, (s, Configuration{ComplexF64, Float64, Float64, Float64}))
end

end # module

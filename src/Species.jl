abstract type AbstractSpecies end
abstract type AbstractFluidSpecies <: AbstractSpecies end
abstract type AbstractKineticSpecies <: AbstractSpecies end

abstract type AbstractClassicalSpecies <: AbstractKineticSpecies end
abstract type AbstractRelativisticSpecies <: AbstractKineticSpecies end

abstract type AbstractCoupledVelocitySpecies{Tz⊥} <:
  AbstractClassicalSpecies end

abstract type AbstractCoupledRelativisticSpecies{Tz⊥} <:
  AbstractRelativisticSpecies end

abstract type AbstractSeparableVelocitySpecies{Tz, T⊥} <:
  AbstractClassicalSpecies end

(S::AbstractKineticSpecies)(v) = S(v...)

"""
    ColdSpecies <: AbstractFluidSpecies
Cold plasma species.
**Fields:**
- `Π` Plasma frequency [rad / s]
- `Ω` Cyclotron Frequency [rad / s]
"""
struct ColdSpecies{TΠ<:Number, TΩ<:Number} <: AbstractFluidSpecies
  Π::TΠ # Plasma frequency in SI
  Ω::TΩ # Plasma frequency in SI
end

"""
Warm plasma species, with speeds of sound parallel and
perpendicular to the magnetic field.
**Fields:**
- `Π` Plasma frequency [rad / s]
- `Ω` Cyclotron Frequency [rad / s]
- `soundspeed` Sound speed [m/s]
"""
struct WarmSpecies{TΠ<:Number, TΩ<:Number, V<:Number} <: AbstractFluidSpecies
  Π::TΠ
  Ω::TΩ
  soundspeed::V
end

"""
    WarmSpecies(Π,Ω,thermalspeed,adiabiaticindex::Number)

Warm plasma species - accept thermalspeed and ratio of specific heats to get
sound speed

...
# Arguments
- `Π`: Plasma frequency [rad / s]
- `Ω`: Cyclotron Frequency [rad / s]
- `thermalspeed`: Thermal speed of Maxwellian distribution
- `adiabiaticindex`: Equation of state Gruneisen gamma (ratio of specific heats)
...
"""
function WarmSpecies(Π, Ω, thermalspeed, adiabiaticindex)
  0 < adiabiaticindex < 10 || @warn "The adiabiatic index is $adiabiaticindex"
  return WarmSpecies(Π, Ω, thermalspeed * sqrt(adiabiaticindex))
end

"""
Kinetic plasma species with separable distribution functions parallel and
perpendicular to the magnetic field.
**Fields:**
- `Π` Plasma frequency [rad / s]
- `Ω` Cyclotron Frequency [rad / s]
- `Fz :: AbstractFParallel`
  Distribution function parallel to magnetic field (normalised)
- `F⊥ :: AbstractFPerpendicular`
  Distribution function perpendicular to magnetic field (normalised)
"""
struct SeparableVelocitySpecies{
    TΠ<:Number, TΩ<:Number,
    Tz<:AbstractFParallel,
    T⊥<:AbstractFPerpendicular
    } <: AbstractSeparableVelocitySpecies{Tz, T⊥}
  Π::TΠ # plasma frequency
  Ω::TΩ # cyclotron frequency
  Fz::Tz
  F⊥::T⊥
  minharmonics::Int
end
function SeparableVelocitySpecies(Π, Ω, Fz, F⊥; minharmonics=DEFAULT_MIN_HARMONICS)
  return SeparableVelocitySpecies(Π, Ω, Fz, F⊥, minharmonics)
end
(S::SeparableVelocitySpecies)(vz, v⊥) = S.Fz(vz) * S.F⊥(v⊥)

"""
Kinetic plasma species defined by one coupled distribution function in velocity
space, parallel and perpendicular to the background magnetic field
**Fields:**
- `Π` Plasma frequency [rad / s]
- `Ω` Cyclotron Frequency [rad / s]
- `F :: AbstractCoupledVelocity`
  Distribution function in velocity space parallel and perpendicular to the
  background magnetic field (normalised)
"""
struct CoupledVelocitySpecies{
    TΠ<:Number, TΩ<:Number,
    TF<:AbstractCoupledVelocity
    } <: AbstractCoupledVelocitySpecies{TF}
  Π::TΠ # plasma frequency with rest mass
  Ω::TΩ # cyclotron frequency with rest mass
  F::TF
  minharmonics::Int
end
function CoupledVelocitySpecies(Π, Ω, F::AbstractCoupledVelocity; minharmonics=DEFAULT_MIN_HARMONICS)
  return CoupledVelocitySpecies(Π, Ω, F, minharmonics)
end
(S::CoupledVelocitySpecies)(vz, v⊥) = S.F(vz, v⊥)
"""
    CoupledVelocitySpecies(Π,Ω,vthz,vth⊥=vthz,vzdrift=0.0,v⊥drift=0.0)

...
# Arguments
- `Π`: plasma frequency [rad/s]
- `Ω`: cyclotron frequency [rad/s]
- `vthz`: parallel thermal speed [m/s]
- `vth⊥=vthz`: perpendicular thermal speed [m/s]
- `vzdrift=0.0`: parallel bulk speed [m/s]
- `v⊥drift=0.0`: perpendicular bulk speed [m/s]
...

"""
function CoupledVelocitySpecies(Π, Ω, vthz, vth⊥=vthz, vzdrift=0.0, v⊥drift=0.0;
    minharmonics=DEFAULT_MIN_HARMONICS)
  return CoupledVelocitySpecies(Π, Ω,
    FCoupledVelocityNumerical(vthz, vth⊥, vzdrift, v⊥drift); minharmonics)
end

"""
Kinetic plasma species defined by one coupled distribution function in momentum
space such that the relativistic dielectric tensor can be calculated.
**Fields:**
- `Π` Plasma frequency [rad / s]
- `Ω` Cyclotron Frequency [rad / s]
- `mass` Species particle mass [kg]
- `F :: AbstractFRelativisticMomentum`
  Distribution function in momentum space parallel and perpendicular to the
  background magnetic field (normalised)
"""
struct CoupledRelativisticSpecies{
    TΠ<:Number, TΩ<:Number, Tm<:Number,
    TF<:AbstractFRelativisticMomentum
    } <: AbstractCoupledRelativisticSpecies{TF}
  Π::TΠ # plasma frequency with rest mass
  Ω::TΩ # cyclotron frequency with rest mass
  m::Tm # rest mass of single particle
  F::TF # the distribution function
  minharmonics::Int
  function CoupledRelativisticSpecies(Π::TΠ, Ω::TΩ, m::Tm, F::TF, h::Int
      ) where {TΠ, TΩ, Tm, TF}
    @warn "CoupledRelativisticSpecies not stress tested"
    return new{TΠ,TΩ,Tm,TF}(Π, Ω, m, F, h)
  end
end
function CoupledRelativisticSpecies(Π, Ω, m, F; minharmonics=DEFAULT_MIN_HARMONICS)
  return CoupledRelativisticSpecies(Π, Ω, m, F, minharmonics)
end
(S::CoupledRelativisticSpecies)(pz, p⊥) = S.F(pz, p⊥)

"""
    CoupledRelativisticSpecies(Π,Ω,m,pthz::Number,pth⊥=pthz,pzdrift=0)

...
# Arguments
- `Π`: classical plasma frequency [rad/s]
- `Ω`: classical cyclotron frequency [rad/s]
- `m`: mass [kg]
- `pthz::Number`: parallel thermal momentum [kg m/s]
- `pth⊥=pthz`: perpendicular thermal moment [kg m/s]
- `pzdrift=0`: parallel bulk momentum [kg m/s]
...

"""
function CoupledRelativisticSpecies(Π, Ω, m, pthz::Number, pth⊥=pthz, pzdrift=zero(pthz);
    minharmonics=DEFAULT_MIN_HARMONICS)
  return CoupledRelativisticSpecies(Π, Ω, m,
    FRelativisticNumerical(m, pthz, pth⊥, pzdrift); minharmonics)
end

"""
    MaxwellJuttnerSpecies(Π,Ω,m,TeV)

MaxwellJuttnerSpecies implements a species with a relativistic Maxwell-Juttner
distribution function.

...
# Arguments
- `Π`: classical plasma frequency [rad/s]
- `Ω`: classical cyclotron frequency [rad/s]
- `m`: mass [kg]
- `TeV`: temperature [eV]
...
"""
struct MaxwellJuttnerSpecies{
    TΠ<:Number, TΩ<:Number, Tm<:Number, T<:Number,
    } <: AbstractRelativisticSpecies
  Π::TΠ # plasma frequency with rest mass
  Ω::TΩ # cyclotron frequency with rest mass
  m::Tm # rest mass of single particle in kg
  TeV::T # temperature in eV
  function MaxwellJuttnerSpecies(Π::TΠ, Ω::TΩ, m::Tm, TeV::T
      ) where {TΠ, TΩ, Tm, T}
    @warn "MaxwellJuttnerSpecies not stress tested"
    return new{TΠ,TΩ,Tm,T}(Π, Ω, m, TeV)
  end
end


"""
    MaxwellianSpecies(Π,Ω,vthb,vth⊥=vthb,vdb=0.0)

Kinetic Maxwellian Plasma species that can optionally have a drift along the
magnetic field

...
# Arguments
- `Π`: Plasma frequency [rad / s]
- `Ω`: Cyclotron Frequency [rad / s]
- `vthb`: parallel thermal speed [m/s]
- `vth⊥=vthb`: perpendicular thermal speed [m/s]
- `vdb=0.0`: parallel beam speed [m/s]
# Returns
- `SeparableVelocitySpecies(Π, Ω, FBeam(vthb, vdb), FPerpendicularMaxwellian(vth⊥))`
...

# Example
```julia
```
"""
function MaxwellianSpecies(Π, Ω, vthb, vth⊥=vthb, vdb=0.0;
    minharmonics=DEFAULT_MIN_HARMONICS)
  @assert vthb > 0.0 && vth⊥ > 0.0
  Fz = FBeam(vthb, vdb)
  F⊥ = FPerpendicularMaxwellian(vth⊥)
  return SeparableVelocitySpecies(Π, Ω, Fz, F⊥, minharmonics)
end

"""
Create a kinetic plasma species with separable distribution functions
parallel ``f(v_\\parallel)`` and perpendicular ``f(v_\\perpendicular)``
to the magnetic field, which are defined as a drifting beam and
a ring respectively.
...
# Arguments
- `Π`: Plasma frequency [rad / s]
- `Ω`: Cyclotron Frequency [rad / s]
- `vthb`: parallel thermal speed [m/s]
- `vth⊥=vthb`: perpendicular thermal speed [m/s]
- `vdb=0.0`: parallel beam speed [m/s]
- `vd⊥=0.0`: perpendicular ring speed [m/s]
# Returns
- `SeparableVelocitySpecies(Π, Ω, FBeam(vthb, vdb), FRing(vth⊥, vd⊥))`
...

# Example
```julia
```
"""
function RingBeamSpecies(Π, Ω, vthb, vth⊥=vthb, vdb=0.0, vd⊥=0.0;
    minharmonics=DEFAULT_MIN_HARMONICS)
  @assert vthb > 0.0 && vth⊥ > 0.0
  Fz = FBeam(vthb, vdb)
  F⊥ = FRing(vth⊥, vd⊥)
  return SeparableVelocitySpecies(Π, Ω, Fz, F⊥, minharmonics)
end

plasmafrequency(S::AbstractSpecies) = S.Π
cyclotronfrequency(S::AbstractSpecies) = S.Ω

is_normalised(S::AbstractSeparableVelocitySpecies) =
  is_normalised(S.Fz) &&
  is_normalised(S.F⊥)

is_normalised(S::AbstractCoupledVelocitySpecies) =
    is_normalised(S.F)

is_normalised(S::AbstractFluidSpecies) = true

ColdSpecies(s::AbstractKineticSpecies) = ColdSpecies(s.Π, s.Ω)
ColdSpecies(s::WarmSpecies) = ColdSpecies(s.Π, s.Ω)
WarmSpecies(s::ColdSpecies) = WarmSpecies(s.Π, s.Ω, 0.0)
"""
    WarmSpecies(s::T,γ=5/3)where{T<:AbstractSeparableVelocitySpecies{

Create a warm plasma species from the AbstractSeparableVelocitySpecies

...
# Arguments
- `s::T`:
- `γ=5/3where{T<:AbstractSeparableVelocitySpecies{`:
...

"""
function WarmSpecies(s::T, γ=5/3) where {T<:AbstractSeparableVelocitySpecies{
    <:AbstractFParallelAnalytical,
    <:AbstractFPerpendicularAnalytical}}
  (s.Fz.vth != s.F⊥.vth) && throw("Parallel and perp thermal speeds not equal")
  return WarmSpecies(s.Π, s.Ω, s.Fz.vth, γ)
end

minharmonics(s::AbstractKineticSpecies) = s.minharmonics

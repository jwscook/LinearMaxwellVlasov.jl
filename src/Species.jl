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
    WarmSpecies(Π::Float64,Ω::Float64,thermalspeed::Float64,adiabiaticindex::Number)

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
end
(S::CoupledVelocitySpecies)(vz, v⊥) = S.F(vz, v⊥)
function CoupledVelocitySpecies(Π::Float64, Ω::Float64, vthb::Float64,
    vth⊥::Float64=vthb, vzdrift::Float64=0.0, v⊥drift::Float64=0.0)
  return CoupledVelocitySpecies(Π, Ω,
    FCoupledVelocityNumerical(vthb, vth⊥,
      vzdrift, v⊥drift))
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
  F::TF
  function CoupledRelativisticSpecies(Π::TΠ, Ω::TΩ, m::Tm, F::TF
      ) where {TΠ, TΩ, Tm, TF}
    @warn "CoupledRelativisticSpecies not stress tested"
    return new{TΠ,TΩ,Tm,TF}(Π, Ω, m, F)
  end
end
(S::CoupledRelativisticSpecies)(pz, p⊥) = S.F(pz, p⊥)

function CoupledRelativisticSpecies(Π, Ω, m, pthz::Number, pth⊥=pthz, pzdrift=0)
  return CoupledRelativisticSpecies(Π, Ω, m,
    FRelativisticNumerical(pthz, pth⊥, pzdrift))
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
function MaxwellianSpecies(Π, Ω, vthb, vth⊥=vthb, vdb=0.0)
  @assert vthb > 0.0 && vth⊥ > 0.0
  Fz = FBeam(vthb, vdb)
  F⊥ = FPerpendicularMaxwellian(vth⊥)
  return SeparableVelocitySpecies(Π, Ω, Fz, F⊥)
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
function RingBeamSpecies(Π, Ω, vthb, vth⊥=vthb, vdb=0.0, vd⊥=0.0)
  @assert vthb > 0.0 && vth⊥ > 0.0
  Fz = FBeam(vthb, vdb)
  F⊥ = FRing(vth⊥, vd⊥)
  return SeparableVelocitySpecies(Π, Ω, Fz, F⊥)
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
function WarmSpecies(s::T, γ=5/3) where {T<:AbstractSeparableVelocitySpecies{
    <:AbstractFParallelAnalytical,
    <:AbstractFPerpendicularAnalytical}}
  (s.Fz.vth != s.F⊥.vth) && throw("Parallel and perp thermal speeds not equal")
  return WarmSpecies(s.Π, s.Ω, s.Fz.vth, γ)
end

using StaticArrays


struct FRelativisticNumerical{T<:Function, U<:Number,
    } <: AbstractFRelativisticMomentum
  F::T
  normalisation::Tuple{U,U}
end

(f::FRelativisticNumerical)(pb, p⊥) = f.F(pb, p⊥)

function FRelativisticNumerical(mass, pthb::Real, pth⊥::Real=pthb, pbdrift::Real=zero(pthb))
  return FRelativisticNumerical(MaxwellJuttner(mass, pthb, pth⊥, pbdrift),
                               (abs(pbdrift) + pthb, pth⊥))
end

"""
  MaxwellJuttner
Return a drifting Maxwellian function.

# Members
mass::Real - rest mass [kg]
pthz::Real - thermal momentum parallel to magnetic field [kg m/s]
pth⊥::Real - thermal momentum perpendicular to magnetic field [kg m/s]
pzdrift::Real=0.0 - drift parallel to the magnetic field [kg m/s]
lognormalisation::Real - the log of the normalisation constant
"""
struct MaxwellJuttner{T<:Real} <: Function
  mass::T
  pthz::T
  pth⊥::T
  pzdrift::T
  mc::T
  μz::T
  μ⊥::T
  γd::T
  γdm1::T
  βdmc::T
  βd_mc::T
  inv2mc²::T
  lognormalisation::T
  function MaxwellJuttner(mass::T, pthz::T, pth⊥::T=pthz, pzdrift::T=zero(T)
      ) where {T<:Real}
    @assert pthz > 0
    @assert pth⊥ > 0
    mc = mass * c₀
    γd = sqrt(1 + (pzdrift / mc)^2)
    γdm1 = (pzdrift/mc)^2 / (γd + 1)
    βd = pzdrift / (γd * mc)
    μz = 2 * (mc / pthz)^2
    μ⊥ = 2 * (mc / pth⊥)^2
    inv2mc² = 1 / 2 / mc^2

    dummy = new{T}(mass, pthz, pth⊥, pzdrift, mc, μz, μ⊥, γd, γdm1, βd * mc, βd / mc, inv2mc², zero(T))
    lower, upper = [-1.0, 0.0] .+ 1e3eps(), [1.0, 1.0] .- 1e3eps()
    integrand = TransformFromInfinity(
      p -> 2π * p[2] * dummy(p), (pthz, pth⊥))
    normalisation = HCubature.hcubature(integrand, lower, upper,
                                        rtol=1e4eps(), atol=0)[1]
    @assert isfinite(normalisation) && !iszero(normalisation) normalisation
    lognormalisation = log(normalisation)
    return new{T}(mass, pthz, pth⊥, pzdrift, mc, μz, μ⊥, γd, γdm1, βd * mc, βd / mc, inv2mc², lognormalisation)
  end
end

function (f::MaxwellJuttner)(p)
  pz, p⊥ = p[1], p[2] # kg m / s
  s = (pz^2 + p⊥^2) / f.mc^2
  γ = sqrt(1 + s)
  γm1 = s / (γ + 1)
  p̂z = f.γd * (pz - f.βdmc * γ) # boosted by drift
  #Γm1 = f.γd * (γ - f.βd * pz / f.mc) - 1
  #Γm1 = f.γd * γ - f.γd* f.βd * pz / f.mc - 1
  #Γm1 = f.γd * (γ - 1) + f.γd - f.γd* f.βd * pz / f.mc - 1
  #Γm1 = f.γd * γm1 - f.γd * (f.βd * pz / f.mc - 1) - 1
  Γm1 = f.γd * (γm1 - pz * f.βd_mc) + f.γdm1
  #ϵ = Γm1 * f.μ⊥ - p̂z^2 * f.inv2mc² * (f.μz - f.μ⊥)
  ϵ = (p̂z^2 * f.inv2mc² + Γm1) * f.μ⊥ - p̂z^2 * f.inv2mc² * f.μz
  #ϵ_classical = (((pz - f.pzdrift) / f.pthz)^2 + (p⊥ / f.pth⊥)^2)
  return exp(-ϵ - f.lognormalisation)
end
(f::MaxwellJuttner)(pz, p⊥) = f((pz, p⊥))


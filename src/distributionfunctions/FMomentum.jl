using StaticArrays


struct FRelativisticNumerical{T<:Function, U<:Number,
    } <: AbstractFRelativisticMomentum
  F::T
  normalisation::Tuple{U,U}
end

(f::FRelativisticNumerical)(pb, p⊥) = f.F(pb, p⊥)

function FRelativisticNumerical(pthb::Real, pth⊥::Real=pthb, pbdrift::Real=0.0)
  return FRelativisticNumerical(RelativisticMaxwellian(pthb, pth⊥, pbdrift),
                               (abs(pbdrift) + pthb, pth⊥))
end

"""
  RelativisticMaxwellian
Return a drifting Maxwellian function.

# Members
pthz::Real - thermal momentum parallel to magnetic field [kg m/s]
pth⊥::Real - thermal momentum perpendicular to magnetic field [kg m/s]
pzdrift::Real=0.0 - drift parallel to the magnetic field [kg m/s]
lognormalisation::Real - the log of the normalisation constant
"""
struct RelativisticMaxwellian{T<:Real,U<:Real,V<:Real,W<:Real} <: Function
  pthz::T
  pth⊥::U
  pzdrift::V
  lognormalisation::W
  function RelativisticMaxwellian(pthz::T, pth⊥::U=pthz, pzdrift::V=0.0) where
    {T<:Real, U<:Real, V<:Real}
    dummy = new{T,U,V,Float64}(pthz, pth⊥, pzdrift, 0.0)
    lower, upper = [-1.0, 0.0] .+ 1000eps(), [1.0, 1.0] .- 1000eps()
    integrand = TransformFromInfinity(
      p -> 2π * p[2] * dummy(p), [pthz, pth⊥])
    normalisation = HCubature.hcubature(integrand, lower, upper,
      rtol=1000eps(), atol=0.0)[1]
    @assert isfinite(normalisation) && !iszero(normalisation)
    lognormalisation = log(normalisation)
    return new{T,U,V,typeof(lognormalisation)}(pthz, pth⊥, pzdrift,
                                               lognormalisation)
  end
end

function (f::RelativisticMaxwellian)(p)
  ϵ̄ = ((p[1] - f.pzdrift) / f.pthz)^2 + (p[2] / f.pth⊥)^2
  return exp(-ϵ̄ - f.lognormalisation)
end
(f::RelativisticMaxwellian)(pz, p⊥) = f((pz, p⊥))

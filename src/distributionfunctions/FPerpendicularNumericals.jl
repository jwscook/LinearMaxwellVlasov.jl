"""
The probability density at perpendicular velocity, v, of a shifted maxwellian
from thermal velocity vth, and drift velocity vd.
"""
struct ShiftedMaxwellianPerpendicular{T<:Number, U<:Number, V<:Number
    } <: Function
  vth::T
  vd::U
  invvthsquared::T
  twoinvvthsquared::T
  lognormalisation::V
end
function ShiftedMaxwellianPerpendicular(vth::T, vd::U
    ) where {T<:Number, U<:Number}
  @assert vth > 0 "vth = $vth"
  @assert vd ≥ 0 "vd = $vd"
  invvthsquared =  1 / vth^2
  lognormalisation = log(π * (sqrt(π) * vth * vd *
    (1 - erf(-vd / vth)) + exp(-vd^2 / vth^2) * vth^2))
  return ShiftedMaxwellianPerpendicular(vth, vd, invvthsquared, 2invvthsquared,
    lognormalisation)
end
function (smp::ShiftedMaxwellianPerpendicular)(v::Number, ∂F∂v::Bool=false)
  value = exp(-(v - smp.vd)^2 * smp.invvthsquared - smp.lognormalisation)
  return ∂F∂v ? value * (smp.vd - v) * smp.twoinvvthsquared : value
end

"""
Type for an arbitrary perpendicular distribution function that holds, a function
for the distribution function, it's derivative and the limits of integral
"""
struct FPerpendicularNumerical{T<:Function, U<:Function
    } <: AbstractFPerpendicularNumerical
  F::T
  dFdv::U
  lower::Float64
  upper::Float64
  _uniqueid::UInt64
  function FPerpendicularNumerical(F::T, dFdv::U, lower::Float64,
      upper::Float64) where {T<:Function, U<:Function}
    @assert (0 <= lower && lower < upper) "$lower $upper"
    _uniqueid = hash((F, dFdv, lower, upper), hash(:FPerpendicularNumerical))
    output = new{T, U}(F, dFdv, lower, upper, _uniqueid)
    @assert is_normalised(output) "$(integrate(output)), $lower, $upper"
    return output
  end
end
function FPerpendicularNumerical(vth::T, vd::T=0.0) where {T<:Number}
  @assert vth > 0 && vd ≥ 0
  smp = ShiftedMaxwellianPerpendicular(vth, vd)
  F = v -> smp(v, false)
  dFdv = v -> smp(v, true)
  lower = max(vd - default_integral_range * vth, zero(vd - vth))
  upper = vd + default_integral_range * vth
  return FPerpendicularNumerical(F, dFdv, lower, upper)
end
function FPerpendicularNumerical(f::T, lower::Float64, upper::Float64
    ) where {T<:Function}
  @assert 0 ≤ lower < upper
  integrand(v) = 2π * v * f(v)
  _, n = normalise(integrand, lower, upper)
  n == 1 && return FPerpendicularNumerical(f, derivative(f), lower, upper)
  invn = 1 / n
  fn(v) = f(v) * invn
  gn = derivative(fn)
  return FPerpendicularNumerical(fn, gn, lower, upper)
end

function (f::FPerpendicularNumerical)(∂F∂v::Bool=false)
  return v -> ∂F∂v ? f.dFdv(v) : f.F(v)
end
@inline function (f::FPerpendicularNumerical)(v::Number, ∂F∂v::Bool=false)
  return ifelse(∂F∂v, f.dFdv(v), f.F(v))
end

lower(f::FPerpendicularNumerical) = f.lower
upper(f::FPerpendicularNumerical) = f.upper

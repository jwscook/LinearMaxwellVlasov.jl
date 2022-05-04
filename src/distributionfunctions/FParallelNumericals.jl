"""
The probability density at parallel velocity, v, of a shifted maxwellian
from thermal velocity vth, and drift velocity vd.
"""
struct ShiftedMaxwellianParallel{T<:Number, U<:Number} <: Function
  vth::T
  vd::U
  invvthsquared::T
  twoinvvthsquared::T
  logvthsqrt2pi::T
end
function ShiftedMaxwellianParallel(vth::T, vd::U
    ) where {T<:Number, U<:Number}
  @assert vth > 0 "vth = $vth"
  invvthsquared = 1 / vth^2
  twoinvvthsquared = 2 / vth^2
  logvthsqrt2pi = log(vth * sqrt(π))
  return ShiftedMaxwellianParallel(vth, vd, invvthsquared, twoinvvthsquared,
    logvthsqrt2pi)
end
function (sm::ShiftedMaxwellianParallel)(v::T, ∂F∂v::Bool=false
    ) where {T<:Number}
  value = exp(-(v - sm.vd)^2 * sm.invvthsquared - sm.logvthsqrt2pi)
  return ∂F∂v ? value * (sm.vd - v) * sm.twoinvvthsquared : value
end

"""
Type for an arbitrary parallel distribution function that holds, a function
for the distribution function, it's derivative and the limits of integral
"""
struct FParallelNumerical{T<:Function, U<:Function} <: AbstractFParallelNumerical
  F::T
  dFdv::U
  lower::Float64
  upper::Float64
  _uniqueid::UInt64
  function FParallelNumerical(F::T, dFdv::U, lower::Float64, upper::Float64
      ) where {T<:Function, U<:Function}
    _uniqueid = hash((F, dFdv, lower, upper), hash(:FParallelNumerical))
    output = new{T, U}(F, dFdv, lower, upper, _uniqueid)
    @assert is_normalised(output)
    return output
  end
end
function FParallelNumerical(f::T, lower::U, upper::U
    ) where {T<:Function, U<:Real}
  @assert lower < upper
  fn, n = normalise(f, lower, upper)
  return FParallelNumerical(fn, derivative(fn), lower, upper)
end
function FParallelNumerical(vth::Number, vd::Number=0.0)
  @assert vth > 0
  lower = vd - default_integral_range * vth
  upper = vd + default_integral_range * vth
  smp = ShiftedMaxwellianParallel(vth, vd)
  F = v -> smp(v, false)
  dFdv = v -> smp(v, true)
  return FParallelNumerical(F, dFdv, lower, upper)
end

(f::FParallelNumerical)(∂F∂v::Bool=false) = v -> ∂F∂v ? f.dFdv(v) : f.F(v)
function (f::FParallelNumerical)(v::Number, ∂F∂v::Bool=false)
  return ∂F∂v ? f.dFdv(v) : f.F(v)
end

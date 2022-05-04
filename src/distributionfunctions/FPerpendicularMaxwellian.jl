struct FPerpendicularMaxwellian{T<:Number, U<:Function
    } <: AbstractFPerpendicularAnalytical
  vth::T # thermal velocity
  F::U # the distribution function
  _uniqueid::UInt64
  function FPerpendicularMaxwellian(vth::T, F::U
      ) where {T<:Number, U<:Function}
    @assert vth > 0.0
    _uniqueid = hash((vth, F), hash(:FPerpendicularMaxwellian))
    return new{T, U}(vth, F, _uniqueid)
  end
end
function FPerpendicularMaxwellian(vth::T) where {T<:Number}
  @assert vth > 0.0
  smp = ShiftedMaxwellianPerpendicular(vth, 0.0)
  return FPerpendicularMaxwellian(vth, smp)
end

(f::FPerpendicularMaxwellian)(∂F∂v::Bool=false) where {T<:Number} = v -> f.F(v, ∂F∂v)
(f::FPerpendicularMaxwellian)(v::T, ∂F∂v::Bool=false) where {T<:Number} = f.F(v, ∂F∂v)

is_normalised(f::FPerpendicularMaxwellian) = true

lower(f::FPerpendicularMaxwellian) = zero(f.vth)
upper(f::FPerpendicularMaxwellian) = default_integral_range * f.vth

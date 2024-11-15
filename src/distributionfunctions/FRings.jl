struct FRing{T<:Number, U<:Number, V<:Function
    } <: AbstractFPerpendicularNumerical
  vth::T # thermal velocity
  vd::U # drift velocity
  F::V # the distribution function
  _uniqueid::UInt64
  function FRing(vth::T, vd::U, F::V) where {T<:Number, U<:Number, V<:Function}
    @assert vth > 0.0
    _uniqueid = hash((vth, vd, F), hash(:FRing))
    return new{T, U, V}(vth, vd, F, _uniqueid)
  end
end
"""
    FRing(vth::T,vd::U=0.0)where{T<:Number,U<:Number}

Create a suitably normalised ring distribution in perpendicular velocity
proportional to `exp(-(v⊥-vd)^2/vthz^2)`

...
# Arguments
- `vth::T`: thermal velocity [m/s]
- `vd::U=0.0`: bulk perpendicular velocity [m/s]
...

"""
function FRing(vth::T, vd::U=0.0) where {T<:Number, U<:Number}
  @assert vth > 0.0
  return FRing(vth, vd, ShiftedMaxwellianPerpendicular(vth, vd))
end

(f::FRing)(v::T, ∂F∂v::Bool=false) where {T<:Number} = f.F(v, ∂F∂v)
(f::FRing)(∂F∂v::Bool=false) = v -> f.F(v, ∂F∂v)

is_normalised(f::FRing) = true
lower(f::FRing) = max(0.0, f.vd - default_integral_range * f.vth)
upper(f::FRing) = f.vd + default_integral_range * f.vth

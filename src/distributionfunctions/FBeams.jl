struct FBeam{T<:Number, U<:Number, V<:Function
    } <: AbstractFParallelAnalytical
  vth::T # thermal velocity
  vd::U # drift velocity
  F::V # the distribution function
  _uniqueid::UInt64
  function FBeam(vth::T, vd::U, F::V) where {T<:Number, U<:Number, V<:Function}
    @assert vth > 0.0
    _uniqueid = hash((vth, vd, F), hash(:FBeam))
    return new{T, U, V}(vth, vd, F, _uniqueid)
  end
end
function FBeam(vth::T, vd::U=0.0) where {T<:Number, U<:Number}
  @assert vth > 0.0
  smp = ShiftedMaxwellianParallel(vth, vd)
  return FBeam(vth, vd, smp)
end

(f::FBeam)(v::T, ∂F∂v::Bool=false) where {T<:Number} = f.F(v, ∂F∂v)

is_normalised(f::FBeam) = true

lower(f::FBeam) = f.vd - f.vth * default_integral_range
upper(f::FBeam) = f.vd + f.vth * default_integral_range

"""
The Dirac-delta function parallel distribution functions
"""
struct FParallelDiracDelta{T<:Number} <: AbstractFParallelNumerical
  v_drift::T
  _uniqueid::UInt64
  function FParallelDiracDelta(v_drift::T) where {T<:Number}
    _uniqueid = hash((v_drift), hash(:FParallelDiracDelta))
    return new{T}(v_drift, _uniqueid)
  end
end

"""
The Dirac-delta function parallel distribution functions
"""
struct FPerpendicularDiracDelta{T<:Number} <: AbstractFPerpendicularNumerical
  v_drift::T
  inv2πv_drift::T
  _uniqueid::UInt64
  function FPerpendicularDiracDelta(v_drift::T) where {T<:Number}
    @assert v_drift >= 0
    inv2πv_drift = 1 / (2π * v_drift)
    _uniqueid = hash((v_drift), hash(:FPerpendicularDiracDelta))
    return new{T}(v_drift, inv2πv_drift, _uniqueid)
  end
end

"""
The integrals over the distribution functions and associated integral kernel
We avoid doing integrals involving the derivative of the Dirac delta function
by doing integral by parts and knowing that f-> 0 and the integral limits
"""
function integrate(f::FParallelDiracDelta, kernel::T, ∂F∂v::Bool,
    _::Tolerance=Tolerance()) where {T<:Function}
  return ∂F∂v ? -derivative(kernel, f.v_drift) : kernel(f.v_drift)
end

"""
Calculate the parallel integral with principal part and residue
"""
function integrate(f::FParallelDiracDelta, numerator::T, pole::Pole,
    ∂F∂v::Bool, _::Tolerance=Tolerance()
    ) where {T<:Function}
  integrand(v) = numerator(v) / (v - float(pole))
  output = ∂F∂v ? -derivative(integrand, f.v_drift) : integrand(f.v_drift)
  if real(pole) == f.v_drift
    f_residue(v) = ∂F∂v ? -derivative(numerator, v) : numerator(v)
    polefix = wavedirectionalityhandler(pole)
    output += polefix(residue(f_residue, polefix(float(pole))))
  end
  return output
end

"""
The integrals over the distribution functions and associated integral kernel
We avoid doing integrals involving the derivative of the Dirac delta function
by doing integral by parts and knowing that f-> 0 and the integral limits
"""
function integrate(f::FPerpendicularDiracDelta, kernel::T,
    ∂F∂v::Bool, _::Tolerance=Tolerance()) where {T<:Function}
  op(x) = ∂F∂v ? -derivative(kernel, x) : kernel(x)
  output = op(f.v_drift) * f.inv2πv_drift
  @assert !any(isnan, output) "$output, $∂F∂v, $(f.v_drift)"
  return output
end

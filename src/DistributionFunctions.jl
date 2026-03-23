using QuadGK, SpecialFunctions, LinearAlgebra

abstract type AbstractDistributionFunction end
abstract type AbstractFParallel <: AbstractDistributionFunction end
abstract type AbstractFPerpendicular <: AbstractDistributionFunction end

abstract type AbstractFParallelNumerical <: AbstractFParallel end
abstract type AbstractFParallelAnalytical <: AbstractFParallel end
abstract type AbstractFPerpendicularNumerical <: AbstractFPerpendicular end
abstract type AbstractFPerpendicularAnalytical <: AbstractFPerpendicular end

abstract type AbstractCoupledVelocity <: AbstractDistributionFunction end

abstract type AbstractFMomentum <: AbstractDistributionFunction end
abstract type AbstractFCoupledMomentum <: AbstractFMomentum end
abstract type AbstractFRelativisticMomentum <: AbstractFCoupledMomentum end

uniqueid(f::AbstractDistributionFunction) = f._uniqueid

quadnorm(x) = maximum(norm.(x))
function is_normalised(f::AbstractDistributionFunction)
  return isapprox(integrate(f), 1, atol=0, rtol=10000 * eps())
end

"""
Tool to normalize a function f between two integral limits a and b
"""
function normalise(f::T, a::Float64, b::Float64) where {T<:Function}
  n, _ = QuadGK.quadgk(f, a, b, rtol=eps(), atol=0, order=DEFAULT_QUADORDER,
                       norm=quadnorm)
  @assert n > 0 "n = $n"
  n == one(n) && return f, one(n)
  invn = 1 / n
  f_normalised = v -> f(v) * invn
  return f_normalised, n
end

include("./distributionfunctions/FDiracDeltas.jl")
include("./distributionfunctions/FParallelNumericals.jl")
include("./distributionfunctions/FPerpendicularMaxwellian.jl")
include("./distributionfunctions/FPerpendicularNumericals.jl")
include("./distributionfunctions/FBeams.jl")
include("./distributionfunctions/FRings.jl")
include("./distributionfunctions/FMomentum.jl")
include("./distributionfunctions/FCoupledVelocity.jl")

lower(f::AbstractDistributionFunction) = f.lower
upper(f::AbstractDistributionFunction) = f.upper

"""
Integration of Abstract Arbitrary FParallel
"""
function integrate(f::AbstractFParallelNumerical, ∂F∂v::Bool=false,
    tol::Tolerance=Tolerance())
  integrand = f(∂F∂v) # TODO
  return QuadGK.quadgk(integrand, f.lower, f.upper, rtol=tol.rel,
    atol=tol.abs, order=DEFAULT_QUADORDER_PARA, norm=quadnorm)[1]
end
function integrate(f::AbstractFParallelNumerical, numerator_kernel::T,
    ∂F∂v::Bool, tol::Tolerance=Tolerance()) where {T<:Function}
  # No pole on this integral, therefore no residue
  fv = f(∂F∂v) # TODO
  integrand(v) = fv(v) * numerator_kernel(v)
  if f.upper == -f.lower
    return first(QuadGK.quadgk(v->integrand(v) + integrand(-v), 0.0, f.upper,
      rtol=tol.rel, atol=tol.abs, order=DEFAULT_QUADORDER_PARA, norm=quadnorm))
  else
    return first(QuadGK.quadgk(integrand, f.lower, f.upper,
      rtol=tol.rel, atol=tol.abs, order=DEFAULT_QUADORDER_PARA, norm=quadnorm))
  end
end

"""
Integrate over the parallel distribution function multiplied by various kernels.
If the imaginary part of the pole is zero, then do a trick that folds over the
integral from the left of the pole to right.
We need to take into account whether the real part of the pole is negative,
otherwise positive slopes at negative velocities give rise to instability and
not damping.
"""
function integrate(f::AbstractFParallelNumerical, numerator_kernel::T,
    pole::Pole, ∂F∂v::Bool, tol::Tolerance=Tolerance()) where {T<:Function}
  fv = f(∂F∂v) # TODO
  numerator(v) = fv(v) * numerator_kernel(v)
  integrand(v) = numerator(v) / (v - pole)

  limits = [f.lower, f.upper] .+ im * pole.deformation
  Δ = (f.upper - f.lower)
  output = QuadGK.quadgk(integrand, limits[1] - Δ, limits[2] + Δ,
    rtol=tol.rel, atol=tol.abs, order=DEFAULT_QUADORDER_PARA, norm=quadnorm)[1]

  return output + residue(numerator, pole)
end

"""
    integrate(f::AbstractFPerpendicular,kernel::T,∂F∂v::Bool,tol::Tolerance=Tolerance()

Integration of Abstract Arbitrary FPerpendicular with respect to kernel where
the distribution function may be differentiated. If quadrature is involved it
will adaptively calculate the integral until the tolerance is met.

...
# Arguments
- `f::AbstractFPerpendicular`:
- `kernel::T`: the kernel e.g. v⊥^2 * BesselJ(n, v⊥^2 * β^2)^2 etc.
- `∂F∂v::Bool`:  differentiate with respect to velocity, or not.
- `tol::Tolerance=Tolerance`:
...

# Example
```julia
```
"""
function integrate(f::AbstractFPerpendicular, kernel::T,
    ∂F∂v::Bool, tol::Tolerance=Tolerance()
    ) where {T<:Function}
  fv = f(∂F∂v) # TODO
  f_integrand(v) = fv(v) .* kernel(v)
  return first(QuadGK.quadgk(f_integrand, lower(f), upper(f),
    rtol=tol.rel, atol=tol.abs, order=DEFAULT_QUADORDER_PERP, norm=quadnorm))
end
function integrate(f::AbstractFPerpendicularNumerical, ∂F∂v::Bool=false,
    tol::Tolerance=Tolerance())
  return integrate(f, x->2π * x, ∂F∂v, tol)
end

using QuadGK, SpecialFunctions, LinearAlgebra

const default_integral_range = 8.0
const quadorder = 32
const quadorder_para = quadorder
const quadorder_perp = quadorder
const quadrature = QuadGK.quadgk

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
  n, _ = quadrature(f, a, b, rtol=eps(), atol=0, order=quadorder, norm=quadnorm)
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
  return quadrature(integrand, f.lower, f.upper, rtol=tol.rel,
    atol=tol.abs, order=quadorder_para, norm=quadnorm)[1]
end
function integrate(f::AbstractFParallelNumerical, numerator_kernel::T,
    ∂F∂v::Bool, tol::Tolerance=Tolerance()) where {T<:Function}
  # No pole on this integral, therefore no residue
  fv = f(∂F∂v) # TODO
  integrand(v) = fv(v) * numerator_kernel(v)
  return first(quadrature(integrand, f.lower, f.upper,
    rtol=tol.rel, atol=tol.abs, order=quadorder_para, norm=quadnorm))
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
    pole::Pole, ∂F∂v::Bool,
    tol::Tolerance=Tolerance()) where {T<:Function}
  fv = f(∂F∂v) # TODO
  numerator(v) = fv(v) * numerator_kernel(v) # / (v - pole) is implied
  integrand = foldnumeratoraboutpole(numerator, pole)

  Δ_2 = (f.upper - f.lower) / 2

  limits = [f.lower, f.upper]
  overlap = (real(pole) - Δ_2  < f.upper) && (real(pole) + Δ_2 > f.lower)
  if overlap
    limits[1] = min(limits[1], real(pole) - Δ_2)
    limits[2] = max(limits[2], real(pole) + Δ_2)
  end

  limits = unique(sort(limitsfolder(limits, real(pole))))
  @assert 2 <= length(limits) <= 3
  principal = first(quadrature(integrand, limits[1], limits[2],
    rtol=tol.rel, atol=tol.abs, order=quadorder_para, norm=quadnorm))

  if length(limits) == 3 # can't use limits... due to weird type instability
    principal += first(quadrature(integrand, limits[2], limits[3],
      rtol=tol.rel, atol=tol.abs, order=quadorder_para, norm=quadnorm))
  end
  # run out of digits? real(pole) ± Δ_2 == real(pole)
  if !overlap
    limits = limitsfolder(real(pole) .+ [- Δ_2, Δ_2], real(pole))
    limits = unique(sort(limits))
    @assert 2 <= length(limits) <= 3
    principal += first(quadrature(integrand, limits[1], limits[2],
      rtol=tol.rel, atol=tol.abs, order=quadorder_para, norm=quadnorm))
    if length(limits) == 3
      principal += first(quadrature(integrand, limits[2], limits[3],
        rtol=tol.rel, atol=tol.abs, order=quadorder_para, norm=quadnorm))
    end
  end
  polefix = wavedirectionalityhandler(pole)
  residueatpole = polefix(residue(numerator, polefix(pole)))

  output = Complex(principal) + residueatpole
  return output
end

"""
Integration of Abstract Arbitrary FPerpendicular
"""
function integrate(f::AbstractFPerpendicular, kernel::T,
    ∂F∂v::Bool, tol::Tolerance=Tolerance()
    ) where {T<:Function}
  fv = f(∂F∂v) # TODO
  f_integrand(v) = fv(v) .* kernel(v)
  return first(quadrature(f_integrand, lower(f), upper(f),
    rtol=tol.rel, atol=tol.abs, order=quadorder_perp, norm=quadnorm))
end
function integrate(f::AbstractFPerpendicularNumerical, ∂F∂v::Bool=false,
    tol::Tolerance=Tolerance())
  return integrate(f, x->2π * x, ∂F∂v, tol)
end

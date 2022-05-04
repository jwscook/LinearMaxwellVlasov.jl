using CommonSubexpressions, MuladdMacro

struct ParallelKernelNumerator <: Function
  power::UInt64
end
(pk::ParallelKernelNumerator)(v::Number) = -v^pk.power

"""
Compile the kernels of the parallel integral and fetch it from the
DistributionFunctions module. If k parallel is zero then there is no
Landau damping, so this case is made separate, as is easier to deal with.
"""
function parallel(Fz::AbstractFParallelNumerical, ω::T,
    kz::U, n::Integer, Ω::Real, power::Unsigned, ∂F∂v::Bool,
    tol::Tolerance=quadtol()) where {T<:Number, U<:Number}
  nΩ = n * Ω
  pole = Pole((ω - nΩ) / kz, kz)
  V = complex(promote_type(typeof(nΩ - ω), typeof(kz))) # to get type stability
  pkn = ParallelKernelNumerator(power)
  return if iszero(kz)
    integrate(Fz, pkn, ∂F∂v, tol) / V(nΩ - ω)
  else
    integrate(Fz, pkn, pole, ∂F∂v, tol) / V(kz)
  end
end

struct MaxwellianIntegralsParallel{T, U, V}
  vth::T
  vd::U
  ∫⁰¹²::NTuple{3, V}
  vth⁻²2::T
end

function MaxwellianIntegralsParallel(Fz::FBeam, ω, kz, nΩ)
  return MaxwellianIntegralsParallel(Fz.vth, Fz.vd, ω, kz, nΩ)
end

"""
The parallel integral of the Beam only requires an integral
over a drifting Maxwellian subject to the relevant kernels. All this is
calculated here.
"""
function MaxwellianIntegralsParallel(vth, vd, ω, kz, nΩ)
  T = promote_type(typeof.((vth, vd, ω, kz, nΩ))...)
  if iszero(kz) # no need to do anything difficult!
    ∫⁰ = T(1 / (ω - nΩ))
    ∫¹ = T(vd / (ω - nΩ))
    ∫² = T((vth^2 / 2 + vd^2) / (ω - nΩ))
  else
    σ⁻¹ = 1 / (kz * vth)
    @muladd z = (ω - kz * vd - nΩ) * σ⁻¹
    z̄ = wavedirectionalityhandler(z, kz)
    Z0 = plasma_dispersion_function(z̄, UInt64(0))
    Z1 = plasma_dispersion_function(z̄, UInt64(1), Z0)
    Z2 = plasma_dispersion_function(z̄, UInt64(2), Z1)
    @cse @muladd begin
      a = Z0
      b = Z0 * vd + Z1 * vth
      c = Z0 * vd^2 + Z1 * 2vth * vd + Z2 * vth^2
    end
    ∫⁰ = - wavedirectionalityhandler(a, kz) * σ⁻¹
    ∫¹ = - wavedirectionalityhandler(b, kz) * σ⁻¹
    ∫² = - wavedirectionalityhandler(c, kz) * σ⁻¹
  end
  return MaxwellianIntegralsParallel(vth, vd, (∫⁰, ∫¹, ∫²), 2 / vth^2)
end
(mib::MaxwellianIntegralsParallel)(power::Integer) = mib.∫⁰¹²[power + 1]

function (mib::MaxwellianIntegralsParallel)(power::Unsigned, ∂F∂v::Bool)
  output = mib(power)
  if ∂F∂v # can't use muladd macro here for some reason, it causes a bug.
    output *= mib.vd
    output -= mib(power + 1)
    output *= mib.vth⁻²2
  end
  return output
end

"""
Parallel integrals for the Beam - used for testing purposes
"""
function parallel(Fz::FBeam, ω::T, kz::U, n::Integer,
    Ω::Real, power::Unsigned, ∂F∂v::Bool, _::Tolerance=quadtol()
    ) where {T<:Number, U<:Number}
  return MaxwellianIntegralsParallel(Fz, ω, kz, n * Ω)(power, ∂F∂v)
end

"""
Interface
The parallel integral of the distribution function and the relevant kernel
"""
function parallel(S::AbstractKineticSpecies, C::Configuration, n::Integer,
    power::Unsigned, ∂F∂v::Bool)
  return parallel(S.Fz, C.frequency, para(C.wavenumber), n, S.Ω, power, ∂F∂v,
    C.options.quadrature_tol)
end

function parallel(S::T, C::Configuration, n::Int
    ) where {T⊥, Tz<:FBeam, T<:AbstractSeparableVelocitySpecies{Tz, T⊥}}
  ω = C.frequency
  kz = para(C.wavenumber)
  mib = MaxwellianIntegralsParallel(S.Fz.vth, S.Fz.vd, ω, kz, n * S.Ω)
  return paralleltuple(mib)
end
function parallel(S::AbstractKineticSpecies, C::Configuration, n::Int)
  return paralleltuple((power, ∂F∂v)->parallel(S, C, n, power, ∂F∂v))
end
function paralleltuple(f::F) where {F}
    return (f(UInt64(0), false),
            f(UInt64(1), false),
            f(UInt64(2), false),
            f(UInt64(0), true),
            f(UInt64(1), true)) # Don't change this order
end

"""
The function parallel gives the same answers for a given set of inputs.
Calculate the key given the inputs and the (identity) operator for storing
and retrieving values
"""
function (::Type{CacheKeyAndOp{ParallelCache}})(species::AbstractKineticSpecies,
    config::Configuration, n::Int)
  ω = config.frequency
  kz = parallel(config.wavenumber)
  hashsig = (kz, ω - n * species.Ω, uniqueid(config.options))
  key = foldr(hash, hashsig; init=UInt64(2^31 + Int32(typeof(hashsig).hash)))
  return key, CacheOp{ParallelCacheOp}(false)
end

struct ParallelCacheOp <: AbstractCacheOpType end

"""Identity operator for storing and retreiving values from CacheDict"""
(cp::CacheOp{ParallelCacheOp})(x) = x

"""
The non-integer argument to the BesselJ
"""
@inline function Jμν(μν::Pair{<:Integer, <:Integer}, z::T) where {T<:Number}
  isequal(μν[1], μν[2]) && return besselj(μν[1], z)^2
  return besselj(μν[1], z) * besselj(μν[2], z)
end

"""
The kernel of the perpendicular integral
"""
struct PerpendicularKernel{T<:Number, U<:Integer, V<:Integer} <: Function
  k⊥_Ω::T
  μν::Pair{U, U}
  power::V
end
function (pk::PerpendicularKernel)(v::Number)
  return 2π * v^Signed(pk.power) * Jμν(pk.μν, v * pk.k⊥_Ω)
end

"""
Calculate the definite integral kernels of the perpenedicular distribution
function, or derivative thereof, multiplied by the kernel functions
"""
function perpendicular(F⊥::AbstractFPerpendicular, Ω::Real, μν::Pair{T, T},
    k⊥::Number, power::Unsigned, ∂F∂v::Bool, tol::Tolerance=Tolerance()
    ) where {T<:Integer}
  kernel = PerpendicularKernel(k⊥ / Ω, μν, power)
  output = integrate(F⊥, kernel, ∂F∂v, tol)
  @assert !any(isnan, output) "output = $output"
  return output
end

struct MaxwellianIntegralsPerpendicular{T, U, V, W, X}
  vth²₂::T # vth^2 / 2
  k⊥_Ω::U
  n::V
  λ::W
  ∫n₋::X
  ∫n₀::X
  ∫n₊::X
end

"""
Perpendicular integrals for a Maxwellian
"""
function MaxwellianIntegralsPerpendicular(vth::Number, k⊥_Ω::Number, n::Integer)
  vth²₂ = vth^2 / 2
  λ = vth²₂ * k⊥_Ω^2
  ∫n₋ = besselix(n - 1, λ)
  ∫n₊ = n == 0 ? ∫n₋ : besselix(n + 1, λ)
  ∫n₀ = n == 0 ? besselix(n, λ) : λ / 2n * (∫n₋ - ∫n₊) # accurate & fast
  return MaxwellianIntegralsPerpendicular(vth²₂, k⊥_Ω, n, λ, ∫n₋, ∫n₀, ∫n₊)
end

∫Jn²F⊥2πv⊥(mi⊥::MaxwellianIntegralsPerpendicular) = mi⊥.∫n₀
∫Jn²∂F⊥2π(mi⊥::MaxwellianIntegralsPerpendicular) = - mi⊥.∫n₀ / mi⊥.vth²₂
function ∫Jn∂JnF⊥2πv⊥²(mi⊥::MaxwellianIntegralsPerpendicular)
  return - ∫Jn∂Jn∂F⊥2πv⊥(mi⊥) * mi⊥.vth²₂
end
function ∫Jn∂Jn∂F⊥2πv⊥(mi⊥::MaxwellianIntegralsPerpendicular)
  return - mi⊥.k⊥_Ω * ((mi⊥.∫n₋ + mi⊥.∫n₊) / 2 - mi⊥.∫n₀)
end
function ∫∂Jn²F⊥2πv⊥³(mi⊥::MaxwellianIntegralsPerpendicular)
  return - mi⊥.vth²₂ * ∫∂Jn²∂F⊥2πv⊥²(mi⊥)
end
function ∫∂Jn²∂F⊥2πv⊥²(mi⊥::MaxwellianIntegralsPerpendicular{T,U,V}
    ) where {T,U,V}
  return mi⊥.λ * (mi⊥.∫n₊ - 2mi⊥.∫n₀ + mi⊥.∫n₋) +
    mi⊥.n * (mi⊥.∫n₊ - mi⊥.∫n₋) / 2
end

function perpendicular(species::T, config::Configuration, n::Integer,
    ) where {Tz, T⊥<:FPerpendicularMaxwellian,
             T<:AbstractSeparableVelocitySpecies{<:Tz, <:T⊥}}
  k⊥ = perpendicular(config.wavenumber)
  return perpendicular(species.F⊥.vth, k⊥ / species.Ω, n)
end

function perpendicular(vth::Real, k⊥_Ω::Number, n::Integer)
  mi⊥ = MaxwellianIntegralsPerpendicular(vth, k⊥_Ω, n)
  nΩ_k⊥ = n / k⊥_Ω

  ⊥2T11_2⊥2T1_1⊥2T_1_1 = 4 * ∫∂Jn²∂F⊥2πv⊥²(mi⊥)
  V⊥ = typeof(⊥2T11_2⊥2T1_1⊥2T_1_1)

  ⊥3F11_2⊥3F1_1⊥3F_1_1 = 4 * ∫∂Jn²F⊥2πv⊥³(mi⊥)

  ⊥2T_1_1_⊥2T11 = 4nΩ_k⊥ * ∫Jn∂Jn∂F⊥2πv⊥(mi⊥)
  iszero(k⊥_Ω) && (⊥2T_1_1_⊥2T11 = V⊥(- n * 2 * isone(abs(n))))

  ⊥3F_1_1_⊥3F11 = 4nΩ_k⊥ * ∫Jn∂JnF⊥2πv⊥²(mi⊥)
  iszero(k⊥_Ω) && (⊥3F_1_1_⊥3F11 = V⊥(n * 2mi⊥.vth²₂ * isone(abs(n))))

  ⊥1T0_1_⊥1T10 = 2 * ∫Jn∂Jn∂F⊥2πv⊥(mi⊥)

  ⊥2F0_1_⊥2F10 = 2 * ∫Jn∂JnF⊥2πv⊥²(mi⊥)

  ⊥1T0_1⊥1T10 = 2nΩ_k⊥ * ∫Jn²∂F⊥2π(mi⊥)
  iszero(k⊥_Ω) && (⊥1T0_1⊥1T10 = V⊥(0))

  ⊥2F0_1⊥2F10 = 2nΩ_k⊥ * ∫Jn²F⊥2πv⊥(mi⊥)
  iszero(k⊥_Ω) && (⊥2F0_1⊥2F10 = V⊥(0))

  ⊥2T112⊥2T1_1⊥2T_1_1 = 4nΩ_k⊥^2 * ∫Jn²∂F⊥2π(mi⊥)
  iszero(k⊥_Ω) && (⊥2T112⊥2T1_1⊥2T_1_1 = V⊥(-2 * isone(abs(n))))

  ⊥3F112⊥3F1_1⊥3F_1_1 = 4nΩ_k⊥^2 * ∫Jn²F⊥2πv⊥(mi⊥)
  iszero(k⊥_Ω) && (⊥3F112⊥3F1_1⊥3F_1_1 = V⊥(2mi⊥.vth²₂ * isone(abs(n))))

  ⊥1F00 = ∫Jn²F⊥2πv⊥(mi⊥)

  return (⊥2T11_2⊥2T1_1⊥2T_1_1, ⊥3F11_2⊥3F1_1⊥3F_1_1, ⊥2T_1_1_⊥2T11,
    ⊥3F_1_1_⊥3F11, ⊥1T0_1_⊥1T10, ⊥2F0_1_⊥2F10, ⊥1T0_1⊥1T10, ⊥2F0_1⊥2F10,
    ⊥2T112⊥2T1_1⊥2T_1_1, ⊥3F112⊥3F1_1⊥3F_1_1, ⊥1F00)
end

"""
Interface
The perpendicular integral of the distribution function and the relevant kernel
"""
function perpendicular(species::AbstractKineticSpecies,
    config::Configuration, μ::Integer, ν::Integer, power::Unsigned, ∂F∂v::Bool)
  output = perpendicular(species.F⊥, species.Ω, Pair(μ, ν),
    perp(config.wavenumber), power, ∂F∂v, config.options.quadrature_tol)
  @assert !any(isnan, output) "output = $output"
  return output
end

function perpendicular(S::AbstractKineticSpecies, C::Configuration, n::Integer)
  ⊥3F11 = perpendicular(S, C, n + 1, n + 1, UInt64(3), false)
  ⊥3F_1_1 = perpendicular(S, C, n - 1, n - 1, UInt64(3), false)
  ⊥3F1_1 = perpendicular(S, C, n + 1, n - 1, UInt64(3), false)
  ⊥2T11 = perpendicular(S, C, n + 1, n + 1, UInt64(2), true)
  ⊥2T_1_1 = perpendicular(S, C, n - 1, n - 1, UInt64(2), true)
  ⊥2T1_1 = perpendicular(S, C, n + 1, n - 1, UInt64(2), true)
  if iszero(S.Ω * n)
    ⊥1F00 = perpendicular(S, C, n, n, UInt64(1), false)
    ⊥1T10 = perpendicular(S, C, n + 1, n, UInt64(1), true)
    ⊥2F10 = perpendicular(S, C, n + 1, n, UInt64(2), false)
    # when n == 0, besselj(1, x) == -besselj(-1, x) &
    # if Ω == 0, besselj(n, v⊥ k⊥ / Ω) -> besselj(±n, Inf) == ±0
    ⊥2F0_1 = -⊥2F10
    ⊥1T0_1 = -⊥1T10
  else # free to divide by nΩ because it's not zero
    k⊥_2nΩ = perpendicular(C.wavenumber) / (2 * n * S.Ω)
    ⊥2F10 = (⊥3F1_1 + ⊥3F11) * k⊥_2nΩ
    ⊥2F0_1 = (⊥3F1_1 + ⊥3F_1_1) * k⊥_2nΩ
    ⊥1F00 = (⊥2F10 + ⊥2F0_1) * k⊥_2nΩ
    ⊥1T10 = (⊥2T1_1 + ⊥2T11) * k⊥_2nΩ
    ⊥1T0_1 = (⊥2T1_1 + ⊥2T_1_1) * k⊥_2nΩ
  end
  # don't change this order; see test/integrals/Memoisation.jl
  return (⊥3F11, ⊥3F_1_1, ⊥3F1_1, ⊥2T11, ⊥2T_1_1, ⊥2T1_1, ⊥2F10, ⊥2F0_1,
    ⊥1F00, ⊥1T10, ⊥1T0_1)
end

# all the below is needed for memoisation

abstract type PerpendicularCacheOp <: AbstractCacheOpType end
struct GenericPerpendicularCacheOp <: PerpendicularCacheOp end
struct MaxwellianPerpendicularCacheOp <: PerpendicularCacheOp end

function CacheOp{T}(x::Number, n::Signed) where {T<:PerpendicularCacheOp}
  a, b, c = n < 0, real(x) < 0, real(x) * imag(x) < 0
  return CacheOp{T}(a, b, c)
end

"""
Regarding integrals with besselj(n, x)* besselj(n±1, x), e.g.:
where m >= 0
# besselj(-m,x) * besselj(-m-1,x) == -besselj(m,x) * besselj(m+1,x)
# besselj(-m+1,x) * besselj(-m-1,x) == +besselj(m-1,x) * besselj(m+1,x)
"""
function (cp::CacheOp{GenericPerpendicularCacheOp})(x::T)::T where {U,T<:NTuple{11,U}}
  y = if cp.a
    (⊥3F11, ⊥3F_1_1, ⊥3F1_1, ⊥2T11, ⊥2T_1_1, ⊥2T1_1,  ⊥2F10,  ⊥2F0_1, ⊥1F00,
       ⊥1T10,  ⊥1T0_1) = x
    (⊥3F_1_1, ⊥3F11, ⊥3F1_1, ⊥2T_1_1, ⊥2T11, ⊥2T1_1, -⊥2F0_1, -⊥2F10, ⊥1F00,
      -⊥1T0_1, -⊥1T10)
    # don't change this order; see test/integrals/Memoisation.jl
  else
    x
  end
  z = if cp.b
    (a, b, c, d, e, f,  g,  h, i,  j,  k) = y
    (a, b, c, d, e, f, -g, -h, i, -j, -k)
    # don't change this order; see test/integrals/Memoisation.jl
  else
    y
  end
  return cp.c ? conj.(z) : z
end


# the tests prove that these two functions work
makepositive(x::Dual{<:Real}) = Dual(abs2(realpart(x)), abs2(dualpart(x)))
function makepositive(x::Dual{T}) where {T<:Complex}
  r, i = reim(realpart(x))
  return Dual(T(abs2(r), abs2(i)), dualpart(x))
end
makepositive(x::Real) = abs2(x)
makepositive(x::T) where {T<:Complex} = T(abs2(real(x)), abs2(imag(x)))

function (::Type{CacheKeyAndOp{PerpendicularCache}})(
    species::AbstractKineticSpecies, config::Configuration, n::Signed)
  k⊥_Ω = perpendicular(config.wavenumber) / species.Ω
  hashsig = (abs(n), makepositive(k⊥_Ω), uniqueid(config.options))
  key = foldr(hash, hashsig; init=UInt64(2^31 + Int32(typeof(hashsig).hash)))
  return key, CacheOp{GenericPerpendicularCacheOp}(k⊥_Ω, n)
end

function (cp::CacheOp{MaxwellianPerpendicularCacheOp})(x::T)::T where {U,T<:NTuple{11,U}}
  y = if cp.a
    (a, b,  c,  d, e, f,  g,  h, i, j, k) = x
    (a, b, -c, -d, e, f, -g, -h, i, j, k)
    # don't change this order; see test/integrals/Memoisation.jl
  else
    x
  end
  z = if cp.b
    (a, b, c, d,  e,  f,  g,  h, i, j, k) = y
    (a, b, c, d, -e, -f, -g, -h, i, j, k)
    # don't change this order; see test/integrals/Memoisation.jl
  else
    y
  end
  return cp.c ? conj.(z) : z
end

function (::Type{CacheKeyAndOp{PerpendicularCache}})(species::T,
    config::Configuration, n::Signed) where {Tz, T⊥<:FPerpendicularMaxwellian,
             T<:AbstractSeparableVelocitySpecies{<:Tz, <:T⊥}}
  k⊥_Ω = perpendicular(config.wavenumber) / species.Ω
  hashsig = (abs(n), makepositive(k⊥_Ω), species.F⊥.vth)
  key = foldr(hash, hashsig; init=UInt64(2^31 + Int32(typeof(hashsig).hash)))
  return key, CacheOp{MaxwellianPerpendicularCacheOp}(k⊥_Ω, n)
end

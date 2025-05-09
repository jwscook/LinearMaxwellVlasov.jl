using DualNumbers, LinearAlgebra, MuladdMacro, SpecialFunctions, StaticArrays
using CommonSubexpressions, HypergeometricFunctions

derivative(f::T) where {T} = x -> DualNumbers.dualpart(f(Dual(x, 1)))
derivative(f::T, x::Number) where {T} = DualNumbers.dualpart(f(Dual(x, 1)))

import Base.^
^(x::DualNumbers.Dual, p::Complex) = exp((log(abs2(x))/2 + im * angle(x)) * p)

import SpecialFunctions.besselj
function besselj(n::Integer, x::DualNumbers.Dual)
  r, d = realpart(x), dualpart(x)
  return Dual(besselj(n, r), d * (besselj(n - 1, r) - besselj(n + 1, r)) / 2)
end

function isapproxinteger(z::Complex, tol=eps())
  rz, iz = reim(z)
  a = round(Int, rz)
  isapprox(a, rz, rtol=tol, atol=tol) || return false
  isapprox(0, iz, rtol=tol, atol=tol) || return false
  return true
end
function besselj(a::T, z) where {T<:Complex}
  #exp(a * log(z/2)) is faster and more accurate than (z/2)^a
  if !(T <: Dual) && (isinteger(a) || isapproxinteger(a, eps()))
    return promote_type(T, typeof(z))(besselj(round(Int, real(a)), z))
  else
    factor(z::Dual, a) = (z / 2)^a / gamma(a + 1) # Can't do Complex(Dual)
    factor(z, a) = exp(a * log(Complex(z) / 2) - loggamma(a + 1))
    return factor(z, a) * HypergeometricFunctions.pFq(
      (@SArray T[]), (@SArray [a + 1]), -z^2 / 4)
  end
end


import SpecialFunctions.besselix
function besselix(n::Integer, x::DualNumbers.Dual{T}) where T
  r, d = realpart(x), dualpart(x)
  bix = besselix(n, r)
  didr = (besseli(n - 1, r) + besseli(n + 1, r)) / 2
  return Dual(bix,
    d * (didr * exp(-abs(real(r))) - sign(real(r)) * bix))
end



#```math
#\mathcal{Z_0}\left(z\right)=\frac{1}{\sqrt{\pi}}\int_{-\infty}^{\infty}\frac{\exp\left(-x^{2}\right)}{x-z}dx+\sigma i\sqrt{\pi}\exp\left(-z^{2}\right)=i\sqrt{\pi}\exp\left(-z^{2}\right)\mathrm{erfc}\left(-iz\right)
#```
"""
Return the value of the plasma dispersion function
This implementation includes the residue, which is easy to verify
because Z(0) = im sqrt(π).

 - x is the argument to the plasma disperion function
 - power is the moment of the integral

[1] S.D. Baalrud, Phys. Plasmas 20, 012118 (2013) and put ν = -Inf
[2] M. Sampoorna et al., Generalized Voigt functions and their derivatives,
  Journal of Quantitative Spectroscopy & Radiative Transfer (2006),
  doi:10.1016/j.jqsrt.2006.08.011
"""
function plasma_dispersion_function(x::T, power::Unsigned=UInt64(0), Z₋₁=missing
    ) where {T<:Number}
  R = float(real(T))
  @inline function _const(i::Signed)::R
    0 <= i <= 1 && return R(i) # make sure that 0 & 1 are handled quickly
    return iseven(i) ? R(0) : R(prod(1:2:(i - 2)) * R(2)^div(1 - i, 2))
  end
  if ismissing(Z₋₁)
    Z = im * sqrt(R(π)) * erfcx(-im * x)# handle e.g. big float π
    for i ∈ 1:power
      @muladd Z = x * Z + _const(Signed(i))
    end
    return Z
  else
    return @muladd x * Z₋₁ + _const(Signed(power))
  end
end

function plasma_dispersion_function(x::T, power::Int, Z₋₁=missing
    ) where {T<:Number}
  @assert power >= 0 "power, $power, must be >= 0"
  return plasma_dispersion_function(x, Unsigned(power), Z₋₁)
end

"""
Takes a function that when integrated between -Inf and +Inf returns value x,
and returns a new function that returns x when integrated between real(pole)
and +Inf.
"""
function foldnumeratoraboutpole(f::T, pole::Real) where {T}
  folded(v::Number) = (f(v + pole) - f(-v + pole)) / v
  folded(v) = (f(v[1] + pole, v[2]) - f(-v[1] + pole, v[2])) / v[1]
  return folded
end
function foldnumeratoraboutpole(f::T, pole::Number) where {T}
  r, i = reim(pole)
  function folded(v)
    a, b, c = f(r + v), f(r - v), 1 / (v - Complex(0, i))
    return (a - b) * real(c) + (a + b) * Complex(0, imag(c))
  end
  return folded
end
"""
Transform the limits of an integrand
quadrature(foldnumeratoraboutpole(integrand, pole), limitsfolder(limits, pole)...)
"""
function limitsfolder(ls::AbstractVector{T}, pole) where {T}
  U = promote_type(T, typeof(pole))
  return vcat(abs.(U.(ls).- real(pole)), isreal(pole) ? eps(U) : zero(U))
end

"""
Transform a function from domain [-∞, ∞]ⁿ down to [-1, 1]ⁿ
"""
struct TransformFromInfinity{T, U}
  f::T
  scale::U
end
function (tfi::TransformFromInfinity)(x)
  @assert all(-1 .<= x .<= 1)
  n = prod(tfi.scale .* (1 .+ x.^2) ./ (1 .- x.^2).^2)
  return n .* tfi.f(tfi.scale .* x ./ (1 .- x.^2))
end
function coordinate(tfi::TransformFromInfinity, x::Number)
  return (-tfi.scale + sqrt(tfi.scale^2 + 4 * x^2)) / (2 * x)
end

"""
Transform a function to transform an integrand from domain [-∞, 0]×[∞, ∞]
down to [-1, -π/2]×[1, π/2].

Example:
```
julia> f(x) = exp(-(x[1]^2 + x[2]^2)/2) * cos(x[2])^2 * sin(x[1])^2
julia> hcubature(f, [-12.0, 0.0], [12.0, 12.0])
(0.7710130943379178, 1.1482318484139944e-8)
julia> hcubature(UnitSemicircleIntegrandTransform(f, 2.0), [0, -π/2], [1, π/2])
(0.7710130940919337, 1.1464196819507393e-8)
```
"""
struct UnitSemicircleIntegrandTransform{T, U<:Number}
  f::T
  scale::U
end
function (usit::UnitSemicircleIntegrandTransform)(x)
  @cse @muladd begin
    y = 1 / (1 - x[1]^2)
    n = usit.scale * (1 + x[1]^2) * y^2
    r = usit.scale * x[1] * y
    sinθ, cosθ = sincos(x[2])
  end
  return n * r * usit.f((r * sinθ, r * cosθ))
end
function coordinates(usit::UnitSemicircleIntegrandTransform, x)
  r = usit.scale * x[1] / (1 - x[1]^2)
  θ = x[2]
  return (r, θ)
end

"""
Concertina and rescale a function in the sections between its roots
such that all roots lie at -1 or +1. Integration over -1..1 of the resulting
function gives the same answer as integration over original -Inf..Inf domain
"""
function transformaboutroots(f::T, root::Real) where {T}
  function output(x)
    return f(root + x / (x + 1) - 1 / 2) * (1 - x / (x + 1)) / (x + 1) +
           f(root + x / (1 - x) + 1 / 2) * (1 + x / (1 - x)) / (1 - x)
  end
  return output
end

function transformaboutroots(f::T, a::Real, b::Real) where {T}
  @assert -Inf < a < b < Inf
  function output(x)
    return f(a + x / (x + 1) - 1 / 2) * (1 - x / (x + 1)) / (x + 1) +
      f(b + x / (1 - x) + 1 / 2) * (1 + x / (1 - x)) / (1 - x) +
      f(a + (x + 1) / 2 * (b - a)) * (b - a) / 2
  end
  return output
end

parallelperpfrompolar(rθ) = rθ[1] .* sincos(rθ[2])
function polarfromparallelperp(vz⊥)
  return (norm(vz⊥), atan(vz⊥[1], vz⊥[2]))
end

function transformtopolar(f::T) where {T}
  fpolar(rθ::AbstractVector{T}) where {T<:Number} = f(parallelperpfrompolar(rθ))
  return fpolar
end

struct ConcertinaSinpi{F}
  numerator::F
  az::Tuple{Int, Int} # the left and right hand integral bounds in normalised units
  function ConcertinaSinpi(num, az)
    @assert az[1] < az[2]
    return new{typeof(num)}(num, az)
  end
end

"""
    (c::ConcertinaSinpi)(t)

Evaluate and sum the ConcertinaSinpi function at location `(t, v⊥)`
where t ∈ (0, 1). t is mapped to (0, 1) + az[1]:az[2], so that
it takes into account the sign of the denominator, sinpi(t + ax[1]:azp[2]),
before finally dividing by |sinpi(t)|.

...
# Arguments
- `c::ConcertinaSinpi(tv⊥)`:
...
"""
function (c::ConcertinaSinpi)(tv⊥)
  t, v⊥ = tv⊥
  @assert 0 < t < 1
  n = c.numerator(t, v⊥)
  breaks = c.az[1]:c.az[2]
  output = zero(n)
  sgn = sign(sinpi((breaks[1] + breaks[2]) / 2))
  for i in 2:length(breaks)
    term = c.numerator(t + breaks[i-1], v⊥)
    output += sgn * term
    sgn *= -1
  end
  return output / sinpi(t)
end

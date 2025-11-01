using DualNumbers, LinearAlgebra, MuladdMacro, SpecialFunctions, StaticArrays
using CommonSubexpressions, HypergeometricFunctions

derivative(f::T, x::Number) where {T} = DualNumbers.dualpart(f(Dual(x, 1)))

import Base.^
^(x::DualNumbers.Dual, p::Complex) = exp((log(abs2(x))/2 + im * angle(x)) * p)

function isapproxinteger(z::Complex, tol=eps())
  rz, iz = reim(z)
  a = round(Int, rz)
  isapprox(a, rz, rtol=tol, atol=tol) || return false
  isapprox(0, iz, rtol=tol, atol=tol) || return false
  return true
end

import SpecialFunctions.besselix
function besselix(n::Integer, x::DualNumbers.Dual{T}) where T
  r, d = realpart(x), dualpart(x)
  bix = besselix(n, r)
  didr = (besseli(n - 1, r) + besseli(n + 1, r)) / 2
  return Dual(bix,
    d * (didr * exp(-abs(real(r))) - sign(real(r)) * bix))
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

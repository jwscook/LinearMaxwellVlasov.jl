using DualNumbers, LinearAlgebra, MuladdMacro, SpecialFunctions, StaticArrays
using CommonSubexpressions

derivative(f::T) where {T} = x -> DualNumbers.dualpart(f(Dual(x, 1)))
derivative(f::T, x::Number) where {T} = DualNumbers.dualpart(f(Dual(x, 1)))

import SpecialFunctions.besselj
function besselj(n::Integer, x::DualNumbers.Dual)
  r, d = realpart(x), dualpart(x)
  return Dual(besselj(n, r), d * (besselj(n - 1, r) - besselj(n + 1, r)) / 2)
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
function foldnumeratoraboutpole(f::T, pole::Real) where {T<:Function}
  folded(v) = (f(v + pole) - f(-v + pole)) / v
  return folded
end
function foldnumeratoraboutpole(f::T, pole::Number) where {T<:Function}
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
struct TransformFromInfinity{T<:Function, U}
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
struct UnitSemicircleIntegrandTransform{T<:Function, U<:Number}
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

function parallelperpfrompolar(rθ)
  return (rθ[1] * sin(rθ[2]), rθ[1] * cos(rθ[2]))
end
function polarfromparallelperp(vz⊥)
  return (norm(vz⊥), atan(vz⊥[1], vz⊥[2]))
end

function transformtopolar(f::T) where {T<:Function}
  fpolar(rθ::AbstractVector{T}) where {T<:Number} = f(parallelperpfrompolar(rθ))
  return fpolar
end

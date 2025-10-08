struct Pole{T<:Number, U<:Real} <: Number
  pole::T
  causalsign::Int
  deformation::U
end
function Pole(pole::Number, causalsign::Int, cauchydeformationangle=DEFAULT_CAUCHY_DEFORMATION_ANGLE)
  return Pole(pole, causalsign, imagcontourdeformation(pole, causalsign, θ=cauchydeformationangle))
end

function Pole(ω::Number, kz::Number, n::Integer, Ω::Number,
    causalsign::Int, deformation=imagcontourdeformation((ω / Ω - n) * Ω / kz, causalsign))
  @assert real(ω) >= 0
  return Pole((ω / Ω - n) * Ω / kz, causalsign, deformation)
end
Pole(ω, K::Wavenumber, n, Ω) = Pole(ω, parallel(K), n, Ω, K.causalsign)
Pole(ω, K::Wavenumber, n, Ω, deformation) = Pole(ω, parallel(K), n, Ω, K.causalsign, deformation)

for op ∈ (:abs, :conj, :real, :imag, :reim, :isreal, :float, :isfinite, :angle)
  @eval Base.$op(f::Pole) = $op(f.pole)
end
Base.:-(p::Pole) = - p.pole
Base.:-(x::Number, p::Pole) = x - p.pole
Base.:-(p::Pole, x::Number) = p.pole - x
Base.:-(d::Dual, p::Pole) = d - p.pole
Base.:+(x::Number, p::Pole) = x + p.pole
Base.:+(p::Pole, x::Number) = p.pole + x
Base.:*(x::Number, p::Pole) = x * p.pole
Base.:*(p::Pole, x::Number) = p.pole * x
Base.:*(a::Pole, b::Pole) = a.pole * b.pole
Base.:^(p::Pole, n::Integer) = p.pole^n
Base.:/(p::Pole, x::Number) = p.pole / x
pole(p::Pole) = p.pole
import DualNumbers: Dual
Dual(p::Pole, x) = Dual(p.pole, x)

"""
    residue(principalpart,pole::Number)

...
# Arguments
- `principalpart`:
- `pole::Number`:
- `deformation::Real`:
...

"""
function residue(numerator::F, pole::Pole) where {F<:Function}
  principalpart = numerator(pole)
  σ = residuesigma(pole)
  iszero(σ) && return zero(principalpart) # defend against overflow
  return im * (σ * π * principalpart)
end

"""
  residuesigma(pole::Number) = imag(pole) < 0 ? 2 : imag(pole) == 0 ? 1 : 0

Calculate the "sigma" factor of the residue
"""
function residuesigma(pole::Number, causalsign::Real)
  @assert abs(causalsign) == 1
  return if causalsign == 1
    imag(pole) < 0 ? 2 : imag(pole) == 0 ? 1 : 0
  else
    -residuesigma(conj(pole), -causalsign)
  end
end
residuesigma(pole::Pole) = residuesigma(pole.pole - im * pole.deformation, pole.causalsign)

function imagcontourdeformation(pole, causalsign, θ=DEFAULT_CAUCHY_DEFORMATION_ANGLE)
  r, i = reim(pole)
  deformation = if (!iszero(i) && abs(angle(pole)) >= θ)
    zero(r)
  elseif iszero(pole)
    -one(r)
  elseif isfinite(r) && isfinite(i)
    -(tan(θ) * abs(r) + i)
  elseif !isfinite(r) && isfinite(i)
    -(tan(θ) + abs(i))
  else
    -one(r)
  end
  @assert deformation <= 0  "pole, θ, deformation = $pole, $θ, $deformation"
  @assert !iszero(i + deformation) "pole, θ, deformation = $pole, $θ, $deformation"
  return deformation * causalsign
end



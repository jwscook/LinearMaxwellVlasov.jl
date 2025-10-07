struct Pole{T<:Number, U<:Real} <: Number
  pole::T
  causalsign::Int
  deformation::U
end
function Pole(pole::Number, causalsign::Int, cauchydeformationangle=DEFAULT_CAUCHY_DEFORMATION_ANGLE)
  return Pole(pole, causalsign, imagcontourdeformation(pole, causalsign, θ=cauchydeformationangle))
end

function Pole(ω::Number, kz::Number, n::Integer, Ω::Number,
    causalsign::Int, deformation=imagcontourdeformation((ω - n * Ω) / kz, causalsign))
  @assert real(ω) >= 0
  return Pole((ω - n * Ω) / kz, causalsign, deformation)
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


"""
    discretefouriertransform(f::T,n::Int,N=512)where{T<:Function}

calculate the DFT of function f at Fourier mode n with N points in the interval
[0,2π]

...
# Arguments
- `f::T`: DFT of this function
- `n::Int`: DFT of this Fourier mode
- `N=512`: evaluate f at this many points
...
"""
function discretefouriertransform(f::T, n::Int, N=512) where {T<:Function}
  kernel(i) = f(2π * i / N) * cis(i * n / N * 2π)
  output = mapreduce(kernel, +, 0:(N-1)) / N
  output *= iszero(n) ? 1 : 2
  return output
end

"""
Expressing f(x) = ∑ᵢ aᵢ (x - p)ⁱ find a₋₁
"""
function residuepart(f::T, pole::Number, radius::Real=abs(pole) * sqrt(eps()),
    N::Int=64) where {T<:Function}
  kernel(θ) = f(radius * Complex(cos(θ), -sin(θ)) + pole)
  return discretefouriertransform(kernel, -1, N) / 2 * radius
end


"""
    residuepartadaptive(f::T,pole::Number,radius::Real=(isreal(pole)?abs(pole):imag(pole))*sqrt(eps()),

Expressing f(x) = ∑ᵢ aᵢ (x - p)ⁱ find a₋₁ by doing a Laurent transform via discrete Fourier
transform by evaluating `f` on walk around the pole. The number of evaluations grows with
each loop until the tolerance has been met.

...
# Arguments
- `f::T`: function to calculate residue of...
- `pole::Number`: ... at this pole
- `radius::Real=(isreal(pole) ? abs(pole) : imag(pole)) * sqrt(eps())`: radius of the residue
- `N::Int=64`: the number of points to evaluate the residue initially
- `tol::Tolerance=Tolerance()`: Tolerance object
- `maxevals::Integer=typemax(Int)`: Maximum number of sampling points
...

"""
function residuepartadaptive(f::T, pole::Number, radius::Real=abs(pole) * sqrt(eps()),
    N::Int=64, tol::Tolerance=Tolerance(), maxevals::Integer=typemax(Int)) where {T}
  @assert ispow2(N) "N must be a power of 2 but it is $N"
  bitreverser(a, b) = ((bitreverse(i) + 2.0^63) * 2.0^-64 for i in a:b)
  inner(θ) = f(radius * Complex(cos(θ), -sin(θ)) + pole)
  outer(x) = inner(2π * x) * cispi(-2x)
  value = mapreduce(outer, +, bitreverser(0, N-1)) / N * radius
  delta = mapreduce(outer, +, bitreverser(N, 2N-1)) / 2N * radius
  @assert !any(isnan, value) "Initial value in residuepartadaptive must not contain NaNs"
  @assert !any(isnan, delta) "Initial delta in residuepartadaptive must not contain NaNs"
  while !isapprox(value, value / 2 + delta, rtol=tol.rel, atol=tol.abs, nans=true)
    value = value / 2 + delta
    N *= 2
    # stop early and error out rather than give an underaccurate result
    @assert N <= maxevals "Number of samples for adaptive residue exceeds limit"
    delta = mapreduce(outer, +, bitreverser(N, 2N-1)) / 2N * radius
  end
  @assert !any(isnan, value) "Final value in residuepartadaptive must not contain NaNs"
  @assert !any(isnan, delta) "Final delta in residuepartadaptive must not contain NaNs"
  return all(isfinite, delta) ? value / 2 + delta : value
end


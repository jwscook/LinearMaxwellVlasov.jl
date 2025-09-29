struct Pole{T<:Number, U<:Real} <: Number
  pole::T
  multipliersign::Int
  deformation::U
end
function Pole(pole::Number, multipliersign::Int)
  return Pole(pole, multipliersign, imagcontourdeformation(pole))
end

function Pole(ω::Number, kparallel::Number, n::Integer, Ω::Number,
    multipliersign::Int, deformation=imagcontourdeformation((ω - n * Ω) / kz))
  @assert real(ω) >= 0
  @assert kparallel >= 0
  return Pole((ω - n * Ω) / kparallel, multipliersign, deformation)
end
Pole(ω, K::Wavenumber, n, Ω, deformation=nothing) = Pole(ω, parallel(K), n, Ω, K.multipliersign, deformation)

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
Base.:^(p::Pole, n::Integer) = p.pole^n
pole(p::Pole) = p.pole
import DualNumbers: Dual
Dual(p::Pole, x) = Dual(p.pole, x)

"""
    residue(principalpart,pole::Number)

...
# Arguments
- `principalpart`:
- `pole::Number`:
...

"""
function residue(principalpart, pole)
  σ = residuesigma(pole)
  output = im * (σ * π * principalpart)
  iszero(σ) && return zero(output) # defend against overflow
  return output
end

function residue(numerator, pole::Number, deformation::Real)
  principalpart = numerator(pole)
  σ = residuesigma(pole - im * deformation)
  iszero(σ) && return zero(principalpart) # defend against overflow
  return im * (σ * π * principalpart)
end
function residue(numerator::F, p::Pole) where {F<:Function}
  return residue(numerator, p.pole, p.deformation)
end

"""
  residuesigma(pole::Number) = imag(pole) < 0 ? 2 : imag(pole) == 0 ? 1 : 0

Calculate the "sigma" factor of the residue
"""
function residuesigma(pole::Number)
  return imag(pole) < 0 ? 2 : imag(pole) == 0 ? 1 : 0
end
residuesigma(pole::Pole) = residuesigma(pole.pole - im * pole.deformation)

function imagcontourdeformation(x, δ=1.0e-2)
  r, i = reim(x)
  θ = abs(angle(Complex(abs(r), abs(i))))
  deformation = if θ >= δ
    zero(r)
  else
    #-abs(i - δ * (iszero(r) ? one(r) : abs(r)))
    -(atan(δ) * abs(r) + abs(i))
  end
  @assert deformation <= 0 # deformation is always negative or zero
  @assert !iszero(i + deformation) "x, δ = $x, $δ"
  return deformation
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


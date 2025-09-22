struct Pole{T<:Number, U<:Real, V<:Real} <: Number
  pole::T
  realkparallel::U
  deformation::V
  function Pole(pole::T, kparallel::U, deformation::V=0.0) where {T<:Number, U<:Real, V<:Real}
    return new{T, U, V}(pole, kparallel, deformation)
  end
end
function Pole(pole::T, kparallel::Complex) where {T<:Number}
  return Pole(pole, real(kparallel))
end

function Pole(ω::Number, kparallel::Number, n::Integer, Ω::Number)
  T = promote_type(typeof(ω), typeof(kparallel), typeof(n), typeof(Ω))
  op(x) = iszero(imag(ω)) ? real(x) : identity(x)
  return Pole(T((op(ω) - n * Ω) / kparallel), kparallel)
end
Pole(ω, K::Wavenumber, n, Ω) = Pole(ω, parallel(K), n, Ω)

for op ∈ (:abs, :conj, :real, :imag, :reim, :isreal, :float, :isfinite, :angle)
  @eval Base.$op(f::Pole) = $op(f.pole)
end
Base.:-(p::Pole) = - p.pole
Base.:-(x::Number, p::Pole) = x - p.pole
Base.:-(p::Pole, x::Number) = p.pole - x
Base.:+(x::Number, p::Pole) = x + p.pole
Base.:+(p::Pole, x::Number) = p.pole + x
Base.:*(x::Number, p::Pole) = x * p.pole
Base.:*(p::Pole, x::Number) = p.pole * x
Base.:^(p::Pole, n::Integer) = p.pole^n
pole(p::Pole) = p.pole
import DualNumbers: Dual
Dual(p::Pole, x) = Dual(p.pole, x)

"""
    wavedirectionalityhandler(x::Number,kz::Number)

Take into account the size of kz when performing parallel integrals

Most codes, assume that kz is real, but this is not always the case.
Usually, they take the absolute value of kz, but this does not respect
complex kz in convetive instability calculations, for example

...
# Arguments
- `x::Number`:
- `kz::Number`:
...
"""
function wavedirectionalityhandler(x::Number, kz::Number)
  # this way works with DualNumbers
  return real(x) + im * (real(kz) > 0 ? imag(x) : -imag(x))
end

wavedirectionalityhandler(x::Number, pole::Pole) = wavedirectionalityhandler(x, pole.realkparallel)
wavedirectionalityhandler(pole::Pole) = x->wavedirectionalityhandler(x, pole)

"""
    residue(numerator::T,pole::Number)where{T<:Function}

Calculate the residue of a number function at a pole

...
# Arguments
- `numerator::T`:
- `pole::Number`:
...
"""
#function residue(numerator::T, pole::Pole) where {T<:Function}
#  return residue(numerator(pole), pole.pole)
#end
function residue(numerator::T, pole::Number) where {T<:Function}
  return residue(numerator(pole), pole)
end

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

function conditionalconj(x, cond::Bool)
  # this way works with DualNumbers
  return real(x) + im * (cond ? imag(x) : -imag(x))
end

function residue(numerator, pole::Number, realkparallel::Real, deformation::Real)
  cond = (realkparallel > 0)
  principalpart = numerator(pole)
  σ = residuesigma(pole - im * deformation)
  iszero(σ) && return zero(principalpart) # defend against overflow
  return im * (σ * π * principalpart) * (realkparallel > 0 ? 1 : -1)
end
function residue(numerator::F, p::Pole) where {F<:Function}
  return residue(numerator, p.pole, p.realkparallel, p.deformation)
end

"""
  residuesigma(pole::Number) = imag(pole) < 0 ? 2 : imag(pole) == 0 ? 1 : 0

Calculate the "sigma" factor of the residue
"""
function residuesigma(pole::Number)
  return imag(pole) < 0 ? 2 : imag(pole) == 0 ? 1 : 0
#  return realkparallel > 0 ? σ : -σ
end
residuesigma(pole::Pole) = residuesigma(pole.pole - im * pole.deformation)

function imagcontourdeformation(x, δ=1.0e-7)
  r, i = reim(x)
  θ = abs(angle(Complex(abs(r), abs(i))))
  deformation = if θ >= δ
    zero(r)
  else
    δ * (iszero(r) ? one(r) : abs(r)) - i
  end
  @assert δ >= 0 # deformatino is always positive
  @assert !iszero(i + deformation) "x, δ = $x, $δ"
  resultangle = abs(angle(r + im * (i + deformation)))
  @assert resultangle > δ || resultangle ≈ δ "resultangle, x, δ = $resultangle, $x, $δ"
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


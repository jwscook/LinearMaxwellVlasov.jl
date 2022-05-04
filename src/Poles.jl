struct Pole{T<:Number, U<:Real} <: Number
  pole::T
  realkparallel::U
end
function Pole(pole::T, kparallel::Complex) where {T<:Number}
  return Pole(pole, real(kparallel))
end

function Pole(ω, K::Wavenumber, n::Int, Ω::Real)
  kparallel = parallel(K)
  return Pole((ω - n * Ω) / kparallel, kparallel)
end

for op ∈ (:abs, :conj, :real, :imag, :reim, :isreal, :float)
  @eval Base.$op(f::Pole) = $op(f.pole)
end
Base.:-(x::Number, p::Pole) = x - p.pole
Base.:-(p::Pole, x::Number) = p.pole - x
Base.:+(x::Number, p::Pole) = x + p.pole

(pole::Pole)(x::Real) = x
(pole::Pole)(x) = wavedirectionalityhandler(x, pole.realkparallel)
wavedirectionalityhandler(pole::Pole) = x->pole(x)
function wavedirectionalityhandler(x::Number, kz::Number)
  return Complex(real(x), (real(kz) < 0 ? -imag(x) : imag(x)))
end

function residue(numerator::T, pole::Number) where {T<:Function}
  return residue(numerator(pole), pole)
end

function residue(principalpart, pole::Number)
  σ = residuesigma(pole)
  output = im * (σ * π * principalpart)
  iszero(σ) && return zero(output) # defend against overflow
  return output
end
residuesigma(pole::Number) = imag(pole) < 0 ? 2 : imag(pole) == 0 ? 1 : 0

function discretefouriertransform(f::T, n::Int, N=512) where {T<:Function}
  kernel(i) = f(2π * i / N) * exp(im * i * n / N * 2π)
  output = mapreduce(kernel, +, 0:(N-1)) / N
  output *= iszero(n) ? 1 : 2
  return output
end

"""
Expressing f(x) = ∑ᵢ aᵢ (x - p)ⁱ find a₋₁
"""
function principalpart(f::T, pole::Number, radius::Real=abs(pole) * sqrt(eps()),
    N::Int=64) where {T<:Function}
  kernel(θ) = f(radius * Complex(cos(θ), -sin(θ)) + pole)
  return discretefouriertransform(kernel, -1, N) / 2 * radius
end


"""
Expressing f(x) = ∑ᵢ aᵢ (x - p)ⁱ find a₋₁
"""
function principalpartadaptive(f::T, pole::Number,
    radius::Real=(isreal(pole) ? abs(pole) : imag(pole)) * sqrt(eps()),
    N::Int=64, tol::Tolerance=Tolerance(); Nmax=2^20) where {T<:Function}
  @assert ispow2(N) "N must be a power of 2 but it is $N"
  bitreverser(a, b) = ((bitreverse(i) for i in a:b) .+ 2.0^63) / 2.0^64
  inner(θ) = f(radius * Complex(cos(θ), -sin(θ)) + pole)
  outer(x) = inner(2π * x) * exp(- im * 2π * x)
  value = mapreduce(outer, +, bitreverser(0, N-1)) / N * radius
  delta = mapreduce(outer, +, bitreverser(N, 2N-1)) / 2N * radius
  while !isapprox(value, value / 2 + delta, rtol=tol.rel, atol=tol.abs, nans=true)
    value = value / 2 + delta
    N >= 2^52 && break # far far far far too many iterations
    N *= 2
    delta = mapreduce(outer, +, bitreverser(N, 2N-1)) / 2N * radius
  end
  return all(isfinite, delta) ? value / 2 + delta : value
end

using Dates
println("$(now()) $(@__FILE__)")

using Random, Test, InteractiveUtils, QuadGK

using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov
using PlasmaDispersionFunctions


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



@testset "Poles" begin

@testset "CauchyResidues" begin
  for σ ∈ (1, -1)
    @testset "principalpart, wavenumber sign is $σ" begin
      for i ∈ 1:10
        a = 10.0^rand(-5:5) * (rand() - 0.5) + im * 10.0^rand(-5:5) * (rand() - 0.5)
        z = 10.0^rand(-5:5) * (rand() - 0.5) + im * 10.0^rand(-5:5) * (rand() - 0.5)
        pole = LMV.Pole(z, σ)
        f(x) = a / (x - pole)
        RP1 = residuepart(f, pole)
        RP2 = residuepartadaptive(f, pole, 8)
        @test RP1 ≈ a atol=0 rtol=sqrt(eps())
        @test RP2 ≈ a atol=0 rtol=sqrt(eps())
      end
    end
  end
end
@testset "isreal isfinite" begin
  for r in (1.0, Inf), s in (-1, 1)
    @test !isreal(LMV.Pole(r + 1im, s))
    @test isreal(LMV.Pole(r + 0im, s))
  end
  @test !isfinite(LMV.Pole(Inf + 0im, 1))
  @test !isfinite(LMV.Pole(0 + Inf*im, 1))
  @test !isfinite(LMV.Pole(Inf + Inf*im, 1))
end

@testset "integral contour deformation logic" begin
  numerator(x) = exp(-x*x) / sqrt(π)
  foobles(x, z) = numerator(x) ./ (x - z)
  for r in (1.0, ), i in (0.1,-0.1, 0.0), deformation in (0.0, 0.01, -0.01, -0.2, 0.2)
    d = deformation
    z = r + im * i
    iszero(d - i) && continue # otherwise test logic is broken ...
    pole = LMV.Pole(z, 1, d)

    expected = plasma_dispersion_function(z)

    σ = LMV.residuesigma(pole)
    rs = im * π * σ * numerator(z)
    ab = QuadGK.quadgk(x->foobles(x, z), -12 + im * d, 12 + im * d)[1] # ... here
    manualresult = ab + rs

    deformedpole = LMV.Pole(z, 1, d)
    result = QuadGK.quadgk(x->foobles(x, z), -12 + im * d, 12 + im * d)[1] + LMV.residue(numerator, deformedpole)

    @testset "r=$r, i=$i, d=$d" begin
      @test expected ≈ manualresult
      @test expected ≈ result
    end
  end
end
@testset "imagcontourdeformation" begin
 @assert LMV.imagcontourdeformation(Inf + NaN*im, 1, 1.0) <= 0
end
@testset "integral contour deformation" begin
  numerator(x) = exp(-x*x) / sqrt(π)
  foobles(x, z) = numerator(x) ./ (x - z)
  for r in (1.0,), i in (0.1, -0.1, 0.0, 1e-10, -1e-10)
    z = r + im * i

    expected = plasma_dispersion_function(z)

    d = LMV.imagcontourdeformation(z, r >= 0 ? 1 : -1, 1.0)
    deformedpole = LMV.Pole(z, 1, d)
    result = QuadGK.quadgk(x->foobles(x, z), -12 + im * d, 12 + im * d)[1] + LMV.residue(numerator, deformedpole)

    @testset "r=$r, i=$i, d=$d" begin
      @test expected ≈ result
    end
  end
end

@testset "discrete fourier transform" begin
  for n ∈ -5:5
    c = rand()
    f(x) = c * sin(n * x)
    fn = discretefouriertransform(f, n)
    @test c * im * (n != 0) ≈ fn atol=sqrt(eps()) rtol=sqrt(eps())
    g(x) = c * cos(n * x)
    gn = discretefouriertransform(g, n)
    @test c ≈ gn atol=sqrt(eps()) rtol=sqrt(eps())
  end
end

end

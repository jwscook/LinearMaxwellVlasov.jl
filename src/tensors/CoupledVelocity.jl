using CommonSubexpressions, DualNumbers, HCubature, LinearAlgebra, QuadGK
using StaticArrays, SpecialFunctions

function coupledvelocity(S::AbstractCoupledVelocitySpecies,
    C::Configuration, n::Int)

  ω, Ω, nΩ = C.frequency, S.Ω, S.Ω * n
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)
  @assert !iszero(k⊥) "Perpendicular wavenumber must not be zero"

  pole = Pole(C.frequency, C.wavenumber, n, S.Ω)
  polefix = wavedirectionalityhandler(pole)

  function denominatorwithpole(vz)
    output = (pole - vz) * kz
    iszero(output) && (output += Inf)
    return output
  end
  denominator(vz) = iszero(kz) ? (ω - nΩ) : denominatorwithpole(vz)
  invdenominator(vz) = 1 / denominator(vz)

  function numerator(vz⊥)
    @assert length(vz⊥) == 2
    vz, v⊥ = vz⊥

    ξ⊥ = v⊥ * k⊥ / Ω
    dfdvz = DualNumbers.dualpart(S(Dual(vz, 1), v⊥))
    dfdv⊥ = DualNumbers.dualpart(S(vz, Dual(v⊥, 1)))

    Jn₋ = besselj(n - 1, ξ⊥)
    Jn₊ = iszero(n) ? -Jn₋ : besselj(n + 1, ξ⊥)
    Jnd = (Jn₋ - Jn₊) / 2
    Jn = iszero(n) ? besselj(n, ξ⊥) : ξ⊥ / 2n * (Jn₋ + Jn₊)
    nΩJn_k⊥ = iszero(k⊥) ? isone(abs(n)) * typeof(Jn)(v⊥ / 2) : nΩ * Jn / k⊥
    @cse @muladd begin
      U = 2π * (v⊥ * kz * dfdvz + (ω - vz * kz) * dfdv⊥)
      Uv⊥ = U * v⊥
      Wv⊥ = 2π * (v⊥ * (ω - nΩ) * dfdvz + nΩ * vz * dfdv⊥)

      M11 = nΩJn_k⊥^2 * U
      M12 = im * nΩJn_k⊥ * Jnd * Uv⊥
      M13 = nΩJn_k⊥ * Jn * Wv⊥
      M21 = -M12 # Onsager
      M22 = v⊥ * Jnd^2 * Uv⊥
      M23 = -im * v⊥ * Jn * Jnd * Wv⊥
      M31 = M13 # Onsager # M31 = nΩJn_k⊥ / v⊥ * vz * Jn * Uv⊥
      M32 = -M23 # Onsager # M32 = im * vz * Jn * Jnd * Uv⊥
      M33 = vz * Jn^2 * Wv⊥
    end

    output = @SArray [M11 M12 M13; M21 M22 M23; M31 M32 M33]
    @assert !any(isnan, output) "output = $output, vz⊥=$vz⊥, $dfdvz, $dfdv⊥"

    return output
  end

  integrand(vz⊥) = numerator(vz⊥) * invdenominator(vz⊥[1])

  function integral2D()
    ∫dvrdθ(vrθ) = vrθ[1] * integrand(parallelperpfrompolar(vrθ))
    return first(HCubature.hcubature(∫dvrdθ,
      (S.F.lower, -π / 2), (S.F.upper, π / 2), initdiv=64,
      rtol=C.options.quadrature_tol.rel, atol=C.options.quadrature_tol.abs))
  end

  function principalzerokz(v⊥)
    @assert iszero(kz)
    ω == nΩ && throw(DomainError("Singularity detected, kz=0, ω == nΩ"))
    ∫dvz(x) = integrand((x, v⊥))
    output = first(QuadGK.quadgk(∫dvz, -S.F.upper, S.F.upper, order=32,
      atol=C.options.quadrature_tol.abs,
      rtol=C.options.quadrature_tol.rel / 10))
    @assert !any(isnan, output) "v⊥ = $v⊥, output = $output"
    return output
  end

  function principal(v⊥)
    @assert !iszero(kz)
    ∫dvz_kz(x) = - numerator((x, v⊥)) / kz
    ∫dvz_kz_folded = foldnumeratoraboutpole(∫dvz_kz, float(pole))
    output = first(QuadGK.quadgk(∫dvz_kz_folded, S.F.lower, S.F.upper, order=32,
        atol=C.options.quadrature_tol.abs,
        rtol=max(eps(), C.options.quadrature_tol.rel / 10)))
    @assert !any(isnan, output) "v⊥ = $v⊥, output = $output"
    return output
  end

  function coupledresidue(v⊥)
    # this started life in relativistic version - can it be simplified?
    ∫dvz(x) = integrand((x, v⊥))
    ppradius = (iszero(imag(pole)) ? abs(pole) : abs(imag(pole))) * sqrt(eps())
    pp = principalpartadaptive(∫dvz, pole, ppradius, 64,
      C.options.summation_tol, Nmax=2048)
    output = polefix.(residue(pp, polefix(pole)))
    output = sign(real(kz)) .* real(output) .+ im .* imag(output)
    @assert !any(isnan, output) "v⊥ = $v⊥, pp = $pp, pole = $pole"
    return output
  end

  function integralsnested1D(∫dv⊥::T) where T
    return first(QuadGK.quadgk(∫dv⊥, S.F.lower, S.F.upper, order=32,
      atol=C.options.quadrature_tol.abs,
      rtol=C.options.quadrature_tol.rel))
  end

  result = if isreal(pole) && iszero(kz)
    integralsnested1D(principalzerokz)
  elseif isreal(pole)# && !iszero(kz)
    integralsnested1D(principal) + integralsnested1D(coupledresidue)
  else
    integral2D() + integralsnested1D(coupledresidue)
  end
  return result
end

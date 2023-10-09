using CommonSubexpressions, DualNumbers, HCubature, QuadGK, SpecialFunctions
using StaticArrays

function relativisticmomentum(S::CoupledRelativisticSpecies,
    C::Configuration, n::Int)
  ω, nΩ = C.frequency, S.Ω * n
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)
  mc² = S.m * c₀^2
  @assert !iszero(k⊥) "Perpendicular wavenumber must not be zero"
  polesarereal = all(iszero, imag.((ω, kz, k⊥)))

  fγ(pz, p⊥) = sqrt(1 + (pz^2 + p⊥^2) / mc²)

  function denominator(pz⊥)
    pz, p⊥ = pz⊥
    output = fγ(pz, p⊥) * ω - kz * pz / S.m - nΩ
    iszero(output) && (output += convert(eltype(output), Inf))
    return output
  end

  function momentumpoles(p⊥)
    a = 1 - (kz * c₀ / ω)^2
    b = - 2 * nΩ * kz * mc² / ω^2
    c = (p⊥ + mc² * S.m - mc² * S.m * (nΩ / ω)^2)
    pzroot1 = (-b + sqrt(b^2 - 4 * a * c)) / (2a)
    pzroot2 = (-b - sqrt(b^2 - 4 * a * c)) / (2a)
    return (Pole(pzroot1, kz), Pole(pzroot2, kz))
  end

  function numerator(pz⊥)
    pz, p⊥ = pz⊥

    # Following Brambilla's book
    γξ⊥ = p⊥ * k⊥ / S.m / S.Ω

    Jn₋ = besselj(n - 1, γξ⊥)
    Jn₊ = iszero(n) ? -Jn₋ : besselj(n + 1, γξ⊥)
    Jn = iszero(n) ? besselj(n, γξ⊥) : γξ⊥ / 2n * (Jn₋ + Jn₊)
    Jnd = (Jn₋ - Jn₊) / 2

    nJn_γξ⊥ = iszero(γξ⊥) ? typeof(γξ⊥)(isone(abs(n)) / 2) : n * Jn / γξ⊥

    dfdpz = DualNumbers.dualpart(S(Dual(pz, 1), p⊥))
    dfdp⊥ = DualNumbers.dualpart(S(pz, Dual(p⊥, 1)))

    γ = fγ(pz, p⊥)

    @cse @muladd begin
      θF = p⊥ * dfdpz - pz * dfdp⊥

      O⊥p⊥ = 2π * p⊥ * (ω * dfdp⊥ + kz / S.m / γ * θF)
      Ob1p⊥ = 2π * p⊥ * (p⊥ * ω * dfdpz - nΩ / γ * θF)
      Ob2p⊥ = 2π * (p⊥ * pz * ω * dfdpz - nΩ / γ * pz * θF)

      m11 = nJn_γξ⊥^2 * p⊥ * O⊥p⊥
      m12 = im * nJn_γξ⊥ * Jnd * p⊥ * O⊥p⊥
      m13 = nJn_γξ⊥ * Jn * Ob1p⊥
      m21 = -m12 # Onsager
      m22 = Jnd^2 * p⊥ * O⊥p⊥
      m23 = -im * Jn * Jnd * Ob1p⊥
      # m31 = nJn_γξ⊥ * Jn * pz * O⊥p⊥
      m31 = m13 # Onsager
      # m32 = im * Jn * Jnd * pz * O⊥p⊥
      m32 = -m23 # Onsager
      m33 = Jn^2 * Ob2p⊥
    end

    return @SArray [m11 m12 m13; m21 m22 m23; m31 m32 m33]
  end
  integrand(pz⊥) = numerator(pz⊥) ./ denominator(pz⊥)

  bound = 1 - 1000 * eps()

  function integral2D()
    return first(HCubature.hcubature(
      UnitSemicircleIntegrandTransform(integrand, norm(S.F.normalisation)),
      (0, -π/2), (1, π/2), initdiv=16,
      rtol=C.options.quadrature_tol.rel, atol=C.options.quadrature_tol.abs))
  end

  outertol = C.options.quadrature_tol.rel
  innertol = outertol / 10 # inner loop has higher accuracy than outer

  function principal(p⊥) # TODO remove probable type instability
    integrandpz = x -> integrand((x, p⊥))
    pzroots = momentumpoles(p⊥)
    objective = if all(isreal, pzroots)
      transformaboutroots(integrandpz, (real.(pzroots))...)
    else
      transformfrominfinity(integrandpz, S.F.normalisation[1])
    end
    principal = first(QuadGK.quadgk(objective, -bound, bound, rtol=innertol))
    @assert !any(isnan, principal)# "principal = $principal"
    return principal
  end
  function relativisticresidue(p⊥)
    integrandpz(x) = integrand((x, p⊥))
    p⊥roots = momentumpoles(p⊥)
    function localresidue(pole)
      polefix = wavedirectionalityhandler(pole)
      rpradius = (iszero(imag(pole)) ? abs(pole) : abs(imag(pole))) * sqrt(eps())
      rp = residuepartadaptive(integrandpz, pole, rpradius, 64,
        C.options.summation_tol, Nmax=2048)
      output = polefix.(residue(rp, polefix(pole)))
      output = Complex.((real(kz) >= 0 ? 1 : -1) * real(output), imag(output))
      return output
    end
    output = mapreduce(localresidue, +, p⊥roots)
    @assert !any(isnan, output)# "output = $output"
    return output
  end
  function integralsnested1D(∫dpz::T, nrm=1) where {T<:Function}
    p⊥normalisation = S.F.normalisation[2]
    transformfunctor = TransformFromInfinity(∫dpz, p⊥normalisation)
    return first(QuadGK.quadgk(
      transformfunctor,
      coordinate(transformfunctor, p⊥normalisation * 1e-8),
      coordinate(transformfunctor, p⊥normalisation * 1e8),
      atol=max(C.options.quadrature_tol.abs,
               outertol * nrm / 2),
      rtol=outertol))
  end

  result = polesarereal ? integralsnested1D(principal) : integral2D()
  result += integralsnested1D(relativisticresidue, norm(result))
  return result
end

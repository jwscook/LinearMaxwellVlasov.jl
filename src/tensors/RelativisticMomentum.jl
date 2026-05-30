using GeneralBesselj

abstract type AbstractRelativisticIntegrand end

struct NewbergerRelativistic{S,T,U,V}
  species::S
  ω::T
  kz::U
  k⊥::V
  count::Ref{Int}
end
NewbergerRelativistic(s, ω, kz, k⊥) = NewbergerRelativistic(s, ω, kz, k⊥, Ref(0))

fγ(nr::NewbergerRelativistic, pz⊥) = sqrt(1 + sum(x->(x / (nr.species.m * c₀))^2, pz⊥))
function pseudoharmonic(nr::NewbergerRelativistic, pz⊥)
  #return (fγ(nr, pz⊥) - pz⊥[1] * nr.kz / (nr.species.m * nr.ω)) * nr.ω / nr.species.Ω
  γ = fγ(nr, pz⊥)
  t = pz⊥[1] * nr.kz / (nr.species.m * nr.ω)
  # Kahan two-sum for γ - t
  s = γ - t
  e = (γ - s) - t   # rounding error correction
  return (s + e) * nr.ω / nr.species.Ω
end


denominator(nr::NewbergerRelativistic, pz⊥) = sinpi(pseudoharmonic(nr, pz⊥))
function (nr::NewbergerRelativistic)(pz⊥)
  return numerator(nr, pz⊥) ./ denominator(nr, pz⊥)
end
function numerator(nr::NewbergerRelativistic, pz⊥)
  nr.count[] += 1
  pz, p⊥ = pz⊥
  ω = nr.ω
  Ω = nr.species.Ω
  kz = nr.kz
  k⊥ = nr.k⊥
  @assert !iszero(k⊥)
  m = nr.species.m
  nz = kz * c₀ / ω
  n⊥ = k⊥ * c₀ / ω

  γ = fγ(nr, pz⊥)
  a = pseudoharmonic(nr, pz⊥)
  sinπa = sinpi(a)
  γξ⊥ = p⊥ * k⊥ / m / Ω # Swanson's b variable

  dfdpz = DualNumbers.dualpart(nr.species(Dual(pz, 1), p⊥))
  dfdp⊥ = DualNumbers.dualpart(nr.species(pz, Dual(p⊥, 1)))

  if iszero(dfdpz) && iszero(dfdp⊥)
    T = promote_type((eltype.((γξ⊥, a, dfdpz, dfdp⊥)))...)
    return @MArray zeros(T, 3, 3)
  end

  Jadual, J_adual = besselj_v(MVector(a, -a), Dual(γξ⊥, 1))
  Ja, Jad = DualNumbers.realpart(Jadual), DualNumbers.dualpart(Jadual)
  J_a, J_ad = DualNumbers.realpart(J_adual), DualNumbers.dualpart(J_adual)
  @assert !isnan(Ja)
  @assert !isnan(J_a)
  @assert !isnan(Jad)
  @assert !isnan(J_ad)

  θF = p⊥ * dfdpz - pz * dfdp⊥

  Qxx = p⊥ / γξ⊥^2 * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * a * Ja * J_a - sinπa) * a
  Qxy = im * p⊥ / γξ⊥ * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * a * Ja * J_ad + a * sinπa / γξ⊥)
  Qxz = (p⊥ * dfdpz) * (π * a * Ja * J_a) / γξ⊥
  Qxz -= (Ω / (γ * ω) * θF) * (π * a * Ja * J_ad - sinπa) * a / γξ⊥
  Qyx = -Qxy
  Qyy = p⊥ * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * Jad * J_ad + a * sinπa / γξ⊥^2)
  Qyz = -im * p⊥ * dfdpz * (π * Ja * J_ad)
  Qyz += im * (Ω / (γ * ω) * θF) * (π * a * Ja * J_ad + a * sinπa / γξ⊥)
  Qzx = pz / γξ⊥ * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * a * Ja * J_a)
  Qzy = im * pz * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * Ja * J_ad)
  Qzz = pz * dfdpz * (π * Ja * J_a)
  Qzz -= pz / p⊥ * Ω / (γ * ω) * θF * (π * a * Ja * J_a)

  Qij = @MArray [Qxx Qxy Qxz; Qyx Qyy Qyz; Qzx Qzy Qzz]
  if !all(!isnan, Qij)
    @show pz, p⊥, a, γξ⊥, dfdpz, dfdp⊥, γ
    @show Ja, J_a, Jad, J_ad
  end
  @assert all(!isnan, Qij) (Qij, a, Ja, J_a, Jad, J_ad, dfdpz, dfdp⊥)
  common = 2π * p⊥ * ω / Ω
  Qij .*= common
  return Qij
end

function momentumpole(nr::NewbergerRelativistic, p⊥, n, deformation)
  kz = nr.kz
  ω = nr.ω
  Ω = nr.species.Ω
  m = nr.species.m
  a = 1 - (kz * c₀ / ω)^2
  b = - 2 * n * Ω * kz * c₀ / ω^2
  c = (p⊥  / (m * c₀))^2 + (1 - (n * Ω / ω)^2)

  pzroot1, pzroot2 = stablequadraticroots(a, b, c) .* (m * c₀)

  v1 = pseudoharmonic(nr, (pzroot1, p⊥))
  v2 = pseudoharmonic(nr, (pzroot2, p⊥))
  pzroot = abs2(real(v1) - n) < abs2(real(v2) - n) ? pzroot1 : pzroot2

  causalsign = real(kz) >= 0 ? 1 : -1
  return Pole(pzroot, causalsign, deformation)
end

function laurentnumerator(nr::NewbergerRelativistic, pz⊥, n)
  γ = fγ(nr, pz⊥)
  mc₀² = nr.species.m * c₀^2
  Ω∂a∂pz = nr.species.Ω * DualNumbers.dualpart(pseudoharmonic(nr, (Dual(pz⊥[1], 1), pz⊥[2])))
  return (-1)^n * nr.species.Ω * numerator(nr, pz⊥) / π / Ω∂a∂pz
end

function relativisticmomentum(S::CoupledRelativisticSpecies, C::Configuration)
  ω, Ω = C.frequency, S.Ω
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)

  @assert !iszero(k⊥) "Perpendicular wavenumber must not be zero"
  polesarereal = all(iszero, imag.((ω, kz, k⊥)))

  integrand = NewbergerRelativistic(S, ω, kz, k⊥)

  bound = 1 - sqrt(eps())

  cubaatol = C.options.cubature_tol.abs
  cubartol = C.options.cubature_tol.rel
  pchar = norm(S.F.normalisation)
  deformation = imagcontourdeformation(ω / kz, real(kz) >= 0 ? 1 : -1,
    pchar * DEFAULT_INTEGRAL_RANGE / S.m, C.options.cauchydeformationangle) * S.m

  function integral2D()
    t1 = @elapsed output, errorestimate = HCubature.hcubature(x->integrand((x[1] + im * deformation, x[2])),
      (-20pchar, 0), (20pchar, 20pchar), initdiv=16,
      rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals, norm=cubanorm)
    if C.options.erroruponcubaturenonconvergence
      msg = "error / val = $(errorestimate / norm(output))"
      msg *= ", count = $(integrand.count[]), time=$t1 seconds"
      @assert (integrand.count[] < C.options.cubature_maxevals) ||
        errorestimate < max(cubartol * norm(output), cubaatol) msg
    end
    return output
  end

  outertol = C.options.quadrature_tol.rel
  innertol = outertol / 10 # inner loop has higher accuracy than outer

  function relativisticresidue(p⊥, pv)
    causalconj = real(kz) >= 0 ? 1 : -1
    function alllocalresidues(n)
      pole = momentumpole(integrand, p⊥, n, deformation)
      @assert pole.deformation == deformation
      output1 = residue(x->laurentnumerator(integrand, (x, p⊥), n), pole)
      @assert !any(isnan, output1)
      return output1
    end
    output = converge(alllocalresidues, minharmonics(S), C.options.summation_tol)
    @assert !any(isnan, output)# "output = $output"
    return output
  end
  function integralsnested1D(∫dpz::T, pv) where {T<:Function}
    p⊥normalisation = S.F.normalisation[2]
    transformfunctor = TransformFromInfinity(x->∫dpz(x, pv), p⊥normalisation)
    return first(QuadGK.quadgk(
      transformfunctor,
      coordinate(transformfunctor, p⊥normalisation * 1e-16),
      coordinate(transformfunctor, p⊥normalisation * 1e3),
      atol=max(C.options.quadrature_tol.abs, outertol * norm(pv) / 2),
      rtol=outertol, norm=quadnorm))
  end

  result = integral2D()
  if !iszero(kz)
    result += integralsnested1D(relativisticresidue, result)
  end
  return result
end


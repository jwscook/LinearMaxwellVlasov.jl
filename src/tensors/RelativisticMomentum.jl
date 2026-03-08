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

function fγ(nr::NewbergerRelativistic, pz⊥)
  y = pz⊥ ./ (nr.species.m * c₀)
  n = norm(y)
  return sqrt(1 + sum(x->x^2, y ./ n) * n^2)
end
function fa(nr::NewbergerRelativistic, pz⊥)
  return (fγ(nr, pz⊥) - pz⊥[1] * nr.kz / (nr.species.m * nr.ω)) * nr.ω / nr.species.Ω
end

denominator(nr::NewbergerRelativistic, pz⊥) = sinpi(fa(nr, pz⊥))
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
  a = fa(nr, pz⊥)
  sinπa = sinpi(a)
  γξ⊥ = p⊥ * k⊥ / m / Ω

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

"""
Same as `numerator` but for integer `a` and with floating point cancellation problems mitigated.
"""
function numeratorintegera(nr::NewbergerRelativistic, pz⊥)
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

  a = round(Int, fa(nr, pz⊥))
  γξ⊥ = p⊥ * k⊥ / m / Ω

  dfdpz = DualNumbers.dualpart(nr.species(Dual(pz, 1), p⊥))
  dfdp⊥ = DualNumbers.dualpart(nr.species(pz, Dual(p⊥, 1)))

  if iszero(dfdpz) && iszero(dfdp⊥)
    T = promote_type((eltype.((γξ⊥, a, dfdpz, dfdp⊥)))...)
    return @MArray zeros(T, 3, 3)
  end

  Jadual= besselj(a, Dual(γξ⊥, 1))
  Ja, Jad = DualNumbers.realpart(Jadual), DualNumbers.dualpart(Jadual)
  @assert !isnan(Ja)
  @assert !isnan(Jad)
  J_a, J_ad = (-1)^a .* (Ja, Jad)

  θF = (p⊥ * dfdpz - pz * dfdp⊥)

  Qxx = (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * a * Ja * J_a) * a * p⊥ / γξ⊥^2
  Qxy = im * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * a * Ja * J_ad) * p⊥ / γξ⊥
  #Qxz = (p⊥ * dfdpz - Ω / (γ * ω) * θF * (-1)^a * a) * π * a * Ja / γξ⊥ * J_a
  Qxz = ((1 - Ω / (γ * ω) * (-1)^a * a) * dfdpz + Ω / (γ * ω) * (-1)^a * a * pz / p⊥ * dfdp⊥) * π * a * Ja / k⊥ * m * Ω * J_a
  Qyx = -Qxy
  Qyy = p⊥ * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * Jad * J_ad)
  Qyz = im * (Ω / (γ * ω) * θF * a - p⊥ * dfdpz) * π * Ja * J_ad
  Qzx = pz * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * a * Ja * J_a) / p⊥ / k⊥ * m * Ω
  Qzy = im * pz * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π * Ja * J_ad)
  Qzz = pz * dfdpz * (π * Ja * J_a) - pz / p⊥ * Ω / (γ * ω) * θF * (π * a * Ja * J_a)

  Qij = @MArray [Qxx Qxy Qxz; Qyx Qyy Qyz; Qzx Qzy Qzz]
  if any(isnan, Qij)
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
  b = - 2 * n * Ω * kz * m * c₀^2 / ω^2
  c = p⊥^2 + m^2 * c₀^2 * (1 - (n * Ω / ω)^2)

  nrm = maximum(abs, (2a, b, c))
  b /= nrm
  c /= nrm
  a /= nrm

  pzroot1, pzroot2 = if n == 0
    @assert iszero(b)
    sqrt(- c / a) .* (-1, 1)
  else
    absb = abs(b)
    pzroot1 = (-b/absb - sqrt(b^2 / absb^2 - 4 * a * c / absb^2)) / (2a) * absb
    pzroot2 = (-b/absb + sqrt(b^2 / absb^2 - 4 * a * c / absb^2)) / (2a) * absb
    (pzroot1, pzroot2)
  end

  ν1 = fa(nr, (pzroot1, p⊥))
  ν2 = fa(nr, (pzroot2, p⊥))
  causalsign = real(kz) >= 0 ? 1 : -1
  if isapproxinteger(ν1, 100eps())
    @assert isapproxinteger(ν1, 100eps()) ν1
    @assert !isapproxinteger(ν2, 100eps()) ν2
    return Pole(pzroot1, causalsign, deformation)
  else
    @assert !isapproxinteger(ν1, 100eps()) ν1
    @assert isapproxinteger(ν2, 100eps()) ν2
    return Pole(pzroot2, causalsign, deformation)
  end
end

function laurentnumerator(nr::NewbergerRelativistic, pz⊥, n)
  mc = (nr.species.m * c₀)
  factor = pz⊥[1] * nr.ω / sqrt(mc + sum(x->x^2, pz⊥)) - nr.kz / nr.species.m
  return (-1)^n * nr.species.Ω * numeratorintegera(nr, pz⊥) / π / factor
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
    pchar * 1000 / S.m, C.options.cauchydeformationangle) * S.m

  function integral2D()
    integrand.count[] = 0
    #output, errorestimate = HCubature.hcubature(
    #  UnitSemicircleIntegrandTransform(
    #    x->integrand((x[1] + im * deformation, x[2])), pchar/100),
    #  (0, -π/2), (1, π/2), initdiv=2,
    #  rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    t1 = @elapsed output, errorestimate = HCubature.hcubature(x->integrand((x[1] + im * deformation, x[2])),
      (-20pchar, 0), (20pchar, 20pchar), initdiv=16,
      rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    #output, errorestimate = HCubature.hcubature(x->integrand((pchar * x[1] + im * deformation, pchar * x[2])),
    #  (-Inf, 0), (Inf, Inf), initdiv=16,
    #  rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    #output /= pchar
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
      coordinate(transformfunctor, p⊥normalisation * 1e2),
      atol=max(C.options.quadrature_tol.abs, outertol * norm(pv) / 2),
      rtol=outertol))
  end

  result = integral2D()
  if !iszero(kz)
    result += integralsnested1D(relativisticresidue, result)
  end
  return result
end


using GeneralBesselj

abstract type AbstractRelativisticStruct end
const ARS = AbstractRelativisticStruct 

struct NewbergerRelativistic{S,T,U,V} <: ARS
  species::S
  ω::T
  kz::U
  k⊥::V
  count::Ref{Int}
end
NewbergerRelativistic(s, ω, kz, k⊥) = NewbergerRelativistic(s, ω, kz, k⊥, Ref(0))

fγ(ars::ARS, pz⊥) = sqrt(1 + sum(x->x^2, pz⊥) / (ars.species.m * c₀)^2)
function fa(ars::ARS, pz⊥)
  return (fγ(ars, pz⊥) * ars.ω - ars.kz * pz⊥[1] / ars.species.m) / ars.species.Ω
end

function (nr::NewbergerRelativistic)(pz⊥)
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
  sinπa = sin(π * a)
  π_sinπa = π / sinπa
  πa_sinπa = π_sinπa * a
  γξ⊥ = p⊥ * k⊥ / m / Ω

  dfdpz = DualNumbers.dualpart(nr.species(Dual(pz, 1), p⊥))
  dfdp⊥ = DualNumbers.dualpart(nr.species(pz, Dual(p⊥, 1)))

  if iszero(dfdpz) && iszero(dfdp⊥)
    T = promote_type((eltype.((γξ⊥, a, dfdpz, dfdp⊥)))...)
    return @MArray zeros(T, 3, 3)
  end
  @assert !isinteger(a) (a, πa_sinπa, pz)

  Ja, J_a, Jad, J_ad = if real(a) > 0
    Ja, J_a, Ja_1, J_a1 = besselj_v([a, -a, a - 1, -a + 1], γξ⊥)
    (Ja, J_a, Ja_1 - Ja * a / γξ⊥, -J_a * a / γξ⊥ - J_a1)
  else
    Ja, J_a, Ja1, J_a_1 = besselj_v([a, -a, a + 1, -a - 1], γξ⊥)
    (Ja, J_a, Ja * a / γξ⊥ - Ja1, J_a_1 + J_a * a / γξ⊥)
  end

  @assert isfinite(Ja) Ja
  @assert isfinite(J_a) J_a
  @assert isfinite(Jad) Jad
  @assert isfinite(J_ad) J_ad

  θF = p⊥ * dfdpz - pz * dfdp⊥

  Qxx = p⊥ / γξ⊥^2 * (dfdp⊥ + kz / (m * γ * ω) * θF) * (πa_sinπa * Ja * J_a - 1) * a
  Qxy = im * p⊥ / γξ⊥ * (dfdp⊥ + kz / (m * γ * ω) * θF) * (πa_sinπa * Ja * J_ad + a / γξ⊥)
  Qxz = 1 / γξ⊥ * (p⊥ * dfdpz) * (πa_sinπa * Ja * J_a)
  Qxz -= 1 / γξ⊥ * (Ω / (γ * ω) * θF) * (πa_sinπa * Ja * J_ad - 1) * a
  Qyx = -Qxy
  Qyy = p⊥ * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π_sinπa * Jad * J_ad + a / γξ⊥^2)
  Qyz = -im * p⊥ * dfdpz * (π_sinπa * Ja * J_ad)
  Qyz += im * (Ω / (γ * ω) * θF) * (πa_sinπa * Ja * J_ad + a / γξ⊥)
  Qzx = pz / γξ⊥ * (dfdp⊥ + kz / (m * γ * ω) * θF) * (πa_sinπa * Ja * J_a)
  Qzy = im * pz * (dfdp⊥ + kz / (m * γ * ω) * θF) * (π_sinπa * Ja * J_ad)
  Qzz = pz * dfdpz * (π_sinπa * Ja * J_a)
  Qzz -= pz / p⊥ * Ω / (γ * ω) * θF * (πa_sinπa * Ja * J_a)

  Qij = @MArray [Qxx Qxy Qxz; Qyx Qyy Qyz; Qzx Qzy Qzz]
  if !all(isfinite, Qij)
    @show pz, p⊥, a, γξ⊥, dfdpz, dfdp⊥, γ
    @show Ja, J_a, Jad, J_ad
    @show Ja_1 = besselj(a - 1, γξ⊥)
    @show J_a1 = besselj(-a + 1, γξ⊥)
    @show Ja1 = besselj(a + 1, γξ⊥)
    @show J_a_1 = besselj(-a - 1, γξ⊥)
  end
  @assert all(isfinite, Qij) (Qij, a, Ja, J_a, Jad, J_ad, dfdpz, dfdp⊥, π_sinπa, πa_sinπa)
  common = 2π * p⊥ * ω / Ω
  Qij .*= common
  return Qij
end

"""
Relativestic dielectric tensor as a function of cyclotron harmonic
This is only used for the principal part.
"""
struct RelativisticHarmonic{S, T, U, V} <: ARS 
  species::S
  ω::T
  kz::U
  k⊥::V
  n::Int
end

function denominator(rh::RelativisticHarmonic, pz⊥)
  pz, p⊥ = pz⊥
  γ = fγ(rh, pz⊥)
  output = γ * rh.ω - rh.kz * pz / rh.species.m - rh.n * rh.species.Ω
  iszero(output) && (output += convert(eltype(output), Inf))
  return output
end

function numerator(rh::RelativisticHarmonic, pz⊥)
  pz, p⊥ = pz⊥

  kz = rh.kz
  k⊥ = rh.k⊥
  ω = rh.ω
  m = rh.species.m
  Ω = rh.species.Ω
  n = rh.n
  nΩ = n * Ω

  # Following Brambilla's book
  γξ⊥ = p⊥ * k⊥ / m / Ω

  Jn₋ = besselj(n - 1, γξ⊥)
  Jn₊ = iszero(n) ? -Jn₋ : besselj(n + 1, γξ⊥)
  Jn = iszero(n) ? besselj(n, γξ⊥) : γξ⊥ / 2n * (Jn₋ + Jn₊)
  Jnd = (Jn₋ - Jn₊) / 2

  nJn_γξ⊥ = iszero(γξ⊥) ? typeof(γξ⊥)(isone(abs(n)) / 2) : n * Jn / γξ⊥

  dfdpz = DualNumbers.dualpart(rh.species(Dual(pz, 1), p⊥))
  dfdp⊥ = DualNumbers.dualpart(rh.species(pz, Dual(p⊥, 1)))

  γ = fγ(rh, pz⊥)

  @cse @muladd begin
    θF = p⊥ * dfdpz - pz * dfdp⊥

    O⊥p⊥ = 2π * p⊥ * (ω * dfdp⊥ + kz / m / γ * θF)
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
(rh::RelativisticHarmonic)(pz⊥) = numerator(rh, pz⊥) ./ denominator(rh, pz⊥)

function momentumpoles(rh::RelativisticHarmonic, p⊥)
  return momentumpoles(rh, p⊥, rh.n)
end
function momentumpoles(ars::AbstractRelativisticStruct, p⊥, n)
  kz = ars.kz
  ω = ars.ω
  Ω = ars.species.Ω
  m = ars.species.m
  a = 1 - (kz * c₀ / ω)^2
  b = - 2 * n * Ω * kz * m * c₀^2 / ω^2
  c = p⊥^2 + m^2 * c₀^2 * (1 - (n * Ω / ω)^2)

  pzroot1 = (-b + sqrt(b^2 - 4 * a * c)) / (2a)
  pzroot2 = (-b - sqrt(b^2 - 4 * a * c)) / (2a)

  pz1 = Pole(pzroot1, kz)
  output = Vector{typeof(pz1)}()
  if isapproxinteger(fa(ars, (pzroot1, p⊥)), 10000eps())
    push!(output, Pole(pzroot1, kz))
  end
  if isapproxinteger(fa(ars, (pzroot2, p⊥)), 10000eps())
    push!(output, Pole(pzroot2, kz))
  end
  return output
end


function relativisticmomentum(S::CoupledRelativisticSpecies, C::Configuration)
  ω, Ω, m = C.frequency, S.Ω, S.m
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)
  @assert !iszero(k⊥) "Perpendicular wavenumber must not be zero"
  polesarereal = all(iszero, imag.((ω, kz, k⊥)))

  integrand = NewbergerRelativistic(S, ω, kz, k⊥)

  bound = 1 - 1e8 * eps()

  cubaatol = C.options.cubature_tol.abs
  cubartol = C.options.cubature_tol.rel

  function integral2D()
    integrand.count[] = 0
    output, errorestimate = HCubature.hcubature(
      UnitSemicircleIntegrandTransform(integrand, norm(S.F.normalisation)),
      (0, -π/2), (1, π/2), initdiv=2,
      rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    if C.options.erroruponcubaturenonconformance
      @assert (integrand.count[] < C.options.cubature_maxevals) ||
        errorestimate < max(cubartol * norm(output), cubaatol)
    end
    return output
  end

  outertol = C.options.quadrature_tol.rel
  innertol = outertol / 10 # inner loop has higher accuracy than outer

  function principal(p⊥) # TODO remove probable type instability
    function allprincipals(n) # TODO remove probable type instability
      rh = RelativisticHarmonic(S, ω, kz, k⊥, n)
      pchar = S.m * real(ω / kz)
      pzroots = momentumpoles(rh, p⊥, n)
      @assert length(pzroots) == 1
      function integrandpz(x)
        x *= pchar
        return rh((x, p⊥)) .* pchar
      end
      objective = transformaboutroots(integrandpz, real(pole(pzroots[1])/pchar))

      polefix = wavedirectionalityhandler(pzroots[1])
      principal = polefix.(first(QuadGK.quadgk(objective, -bound, bound, order=7,
        atol=C.options.quadrature_tol.abs, rtol=C.options.quadrature_tol.rel)))
      principal = sign(real(kz)) .* real(principal) .+ im .* imag(principal)

      @assert !any(isnan, principal)# "principal = $principal"
      return principal
    end
    return converge(allprincipals, C.options.summation_tol)
  end
  function relativisticresidue(p⊥)
    integrandn = NewbergerRelativistic(S, ω, kz, k⊥)
    function alllocalresidues(n)
      integrandpz(x) = integrandn((x, p⊥))
      p⊥roots = momentumpoles(integrandn, p⊥, n)
      function localresidue(pole)
        polefix = wavedirectionalityhandler(pole)
        rpradius = abs(pole) * sqrt(eps())
        rp = residuepartadaptive(integrandpz, pole, rpradius, 64,
          C.options.summation_tol, C.options.residue_maxevals)
        output1 = polefix.(residue(rp, polefix(pole)))
        output1 = sign(real(kz)) .* real(output1) .+ im .* imag(output1)
        return output1
      end
      return mapreduce(localresidue, +, p⊥roots)
    end
    output = converge(alllocalresidues, C.options.summation_tol)
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


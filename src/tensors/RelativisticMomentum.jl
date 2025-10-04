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

function numerator(nr::NewbergerRelativistic, pz⊥)
  a = fa(nr, pz⊥)
  sinπa = sin(π * a)
  return nr(pz⊥) * sinπa
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

  Jadual, J_adual = besselj_v(MVector(a, -a), Dual(γξ⊥, 1))
  Ja, Jad = DualNumbers.realpart(Jadual), DualNumbers.dualpart(Jadual)
  J_a, J_ad = DualNumbers.realpart(J_adual), DualNumbers.dualpart(J_adual)
  @assert isfinite(Ja)
  @assert isfinite(J_a)
  @assert isfinite(Jad)
  @assert isfinite(J_ad)

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

#"""
#Relativestic dielectric tensor as a function of cyclotron harmonic
#This is only used for the principal part.
#"""
#struct RelativisticHarmonic{S, T, U, V} <: ARS 
#  species::S
#  ω::T
#  kz::U
#  k⊥::V
#  n::Int
#end
#
#function denominator(rh::RelativisticHarmonic, pz⊥)
#  pz, p⊥ = pz⊥
#  γ = fγ(rh, pz⊥)
#  output = γ * rh.ω - rh.kz * pz / rh.species.m - rh.n * rh.species.Ω
#  iszero(output) && (output += convert(eltype(output), Inf))
#  return output
#end
#
#function numerator(rh::RelativisticHarmonic, pz⊥)
#  pz, p⊥ = pz⊥
#
#  kz = rh.kz
#  k⊥ = rh.k⊥
#  ω = rh.ω
#  m = rh.species.m
#  Ω = rh.species.Ω
#  n = rh.n
#  nΩ = n * Ω
#
#  # Following Brambilla's book
#  γξ⊥ = p⊥ * k⊥ / m / Ω
#
#  Jn₋ = besselj(n - 1, γξ⊥)
#  Jn₊ = iszero(n) ? -Jn₋ : besselj(n + 1, γξ⊥)
#  Jn = iszero(n) ? besselj(n, γξ⊥) : γξ⊥ / 2n * (Jn₋ + Jn₊)
#  Jnd = (Jn₋ - Jn₊) / 2
#
#  nJn_γξ⊥ = iszero(γξ⊥) ? typeof(γξ⊥)(isone(abs(n)) / 2) : n * Jn / γξ⊥
#
#  dfdpz = DualNumbers.dualpart(rh.species(Dual(pz, 1), p⊥))
#  dfdp⊥ = DualNumbers.dualpart(rh.species(pz, Dual(p⊥, 1)))
#
#  γ = fγ(rh, pz⊥)
#
#  @cse @muladd begin
#    θF = p⊥ * dfdpz - pz * dfdp⊥
#
#    O⊥p⊥ = 2π * p⊥ * (ω * dfdp⊥ + kz / m / γ * θF)
#    Ob1p⊥ = 2π * p⊥ * (p⊥ * ω * dfdpz - nΩ / γ * θF)
#    Ob2p⊥ = 2π * (p⊥ * pz * ω * dfdpz - nΩ / γ * pz * θF)
#
#    m11 = nJn_γξ⊥^2 * p⊥ * O⊥p⊥
#    m12 = im * nJn_γξ⊥ * Jnd * p⊥ * O⊥p⊥
#    m13 = nJn_γξ⊥ * Jn * Ob1p⊥
#    m21 = -m12 # Onsager
#    m22 = Jnd^2 * p⊥ * O⊥p⊥
#    m23 = -im * Jn * Jnd * Ob1p⊥
#    # m31 = nJn_γξ⊥ * Jn * pz * O⊥p⊥
#    m31 = m13 # Onsager
#    # m32 = im * Jn * Jnd * pz * O⊥p⊥
#    m32 = -m23 # Onsager
#    m33 = Jn^2 * Ob2p⊥
#  end
#
#  return @SArray [m11 m12 m13; m21 m22 m23; m31 m32 m33]
#end
#(rh::RelativisticHarmonic)(pz⊥) = numerator(rh, pz⊥) ./ denominator(rh, pz⊥)

function momentumpoles(ars::AbstractRelativisticStruct, p⊥, n, deformation)
  kz = ars.kz
  ω = ars.ω
  Ω = ars.species.Ω
  m = ars.species.m
  a = 1 - (kz * c₀ / ω)^2
  b = - 2 * n * Ω * kz * m * c₀^2 / ω^2
  c = p⊥^2 + m^2 * c₀^2 * (1 - (n * Ω / ω)^2)

  pzroot1 = (-b + sqrt(b^2 - 4 * a * c)) / (2a)
  pzroot2 = (-b - sqrt(b^2 - 4 * a * c)) / (2a)

  causalsign = real(kz) >= 0 ? 1 : -1
  pz1 = Pole(pzroot1, causalsign, deformation)
  output = Vector{typeof(pz1)}()
  if isapproxinteger(fa(ars, (pzroot1, p⊥)), 10000eps())
    push!(output, Pole(pzroot1, causalsign, deformation))
  end
  if isapproxinteger(fa(ars, (pzroot2, p⊥)), 10000eps())
    push!(output, Pole(pzroot2, causalsign, deformation))
  end
  return output
end


function relativisticmomentum(S::CoupledRelativisticSpecies, C::Configuration)
  ω, Ω = C.frequency, S.Ω
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)

  @assert !iszero(k⊥) "Perpendicular wavenumber must not be zero"
  polesarereal = all(iszero, imag.((ω, kz, k⊥)))

  integrand = NewbergerRelativistic(S, ω, kz, k⊥)

  bound = 1 - 1e8 * eps()

  cubaatol = C.options.cubature_tol.abs
  cubartol = C.options.cubature_tol.rel
  deformation = imagcontourdeformation(S.m * ω / kz, real(kz) >= 0 ? 1 : -1)

  function integral2D()
    integrand.count[] = 0
    output, errorestimate = HCubature.hcubature(
      UnitSemicircleIntegrandTransform(
        x->integrand((x[1] + im * deformation, x[2])),
        norm(S.F.normalisation)),
      (0, -π/2), (1, π/2), initdiv=16,
      rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    if C.options.erroruponcubaturenonconformance
      @assert (integrand.count[] < C.options.cubature_maxevals) ||
        errorestimate < max(cubartol * norm(output), cubaatol)
    end
    return output
  end

  outertol = C.options.quadrature_tol.rel
  innertol = outertol / 10 # inner loop has higher accuracy than outer

  function relativisticresidue(p⊥)
    integrandn = NewbergerRelativistic(S, ω, kz, k⊥)
    function alllocalresidues(n)
      integrandpz(x) = integrandn((x, p⊥))
      p⊥roots = momentumpoles(integrandn, p⊥, n, deformation)
      function localresidue(pole)
        @assert pole.deformation == deformation
        laurentnumerator(x) = -(-1)^n * Ω * numerator(integrandn, (x, p⊥)) / kz / π
        output1 = residue(laurentnumerator, pole)
        @assert !any(isnan, output1)
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

  result = integral2D()
#  result += integralsnested1D(relativisticresidue, norm(result))
  return result
end


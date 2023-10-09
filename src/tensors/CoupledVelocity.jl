using CommonSubexpressions, DualNumbers, HCubature, LinearAlgebra, QuadGK
using StaticArrays, SpecialFunctions

abstract type AbstractCoupledIntegrand end

(igrand::AbstractCoupledIntegrand)(vzv⊥) = igrand(vzv⊥[1], vzv⊥[2])

struct HarmonicSum{S,T,U,V,W} <: AbstractCoupledIntegrand
  species::S
  ω::T
  Ω::U
  kz::V
  k⊥::W
  n::Int
end

function numerator(harmonicsum::HarmonicSum, vz⊥)
  @assert length(vz⊥) == 2
  vz, v⊥ = vz⊥

  S = harmonicsum.species
  ω = harmonicsum.ω
  Ω = harmonicsum.Ω
  kz = harmonicsum.kz
  k⊥ = harmonicsum.k⊥
  n = harmonicsum.n

  nΩ = n * Ω

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

  output = @MArray [M11 M12 M13; M21 M22 M23; M31 M32 M33]
  @assert !any(isnan, output)# "output = $output, vz⊥=$vz⊥, $dfdvz, $dfdv⊥"

  return output
end

#=
# deprecated

function denominator(h::HarmonicSum, vz)
  output = h.ω - h.n * h.Ω - vz * h.kz
  iszero(output) && (output += Inf)
  return output
end

(h::HarmonicSum)(vz⊥) = numerator(h, vz⊥) / denominator(h, vz⊥[1])

function coupledvelocity(S::AbstractCoupledVelocitySpecies,
    C::Configuration, n::Int)
  #@warn "Deprecation: calling coupledvelocity(species, config) instead will
  #be significantly faster for most cases"

  ω, Ω, nΩ = C.frequency, S.Ω, S.Ω * n
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)
  @assert !iszero(k⊥) "Perpendicular wavenumber must not be zero"

  pole = Pole(C.frequency, C.wavenumber, n, S.Ω)
  polefix = wavedirectionalityhandler(pole)

  integrand = HarmonicSum(S, ω, Ω, kz, k⊥, n)

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
    @assert !any(isnan, output)# "v⊥ = $v⊥, output = $output"
    return output
  end

  function principal(v⊥)
    @assert !iszero(kz)
    ∫dvz_kz(x) = - numerator(integrand, (x, v⊥)) / kz
    ∫dvz_kz_folded = foldnumeratoraboutpole(∫dvz_kz, float(pole))
    output = first(QuadGK.quadgk(∫dvz_kz_folded, S.F.lower, S.F.upper, order=32,
        atol=C.options.quadrature_tol.abs,
        rtol=max(eps(), C.options.quadrature_tol.rel / 10)))
    @assert !any(isnan, output)# "v⊥ = $v⊥, output = $output"
    return output
  end

  function coupledresidue(v⊥)
    # this started life in relativistic version - can it be simplified?
    ∫dvz(x) = integrand((x, v⊥))
    rpradius = (iszero(imag(pole)) ? abs(pole) : abs(imag(pole))) * sqrt(eps())
    rp = residuepartadaptive(∫dvz, pole, rpradius, 64,
      C.options.summation_tol, Nmax=2048)
    output = polefix.(residue(rp, polefix(pole)))
    output = sign(real(kz)) .* real(output) .+ im .* imag(output)
    @assert !any(isnan, output)# "v⊥ = $v⊥, rp = $rp, pole = $pole"
    return output
  end

  function integralsnested1D(∫dv⊥::T, nrm=1) where T
    return first(QuadGK.quadgk(∫dv⊥, S.F.lower, S.F.upper, order=32,
      atol=max(C.options.quadrature_tol.abs,
               C.options.quadrature_tol.rel * nrm / 2),
      rtol=C.options.quadrature_tol.rel))
  end

  result = if isreal(pole) && iszero(kz)
    integralsnested1D(principalzerokz)
  elseif isreal(pole)# && !iszero(kz)
    pp = integralsnested1D(principal)
    pp .+ integralsnested1D(coupledresidue, norm(pp))
  elseif !iszero(kz) # && !isreal(pole)
    i2d = integral2D()
    i2d .+ integralsnested1D(coupledresidue, norm(i2d))
  else # iszero(kz) && !isreal(pole)
    integral2D()
  end
  return result
end
=#


struct Newberger{S,T,U,V,W} <: AbstractCoupledIntegrand
  species::S
  ω::T
  Ω::U
  kz::V
  k⊥::W
end
function (newberger::Newberger)(vz, v⊥)
  S = newberger.species
  ω = newberger.ω
  Ω = newberger.Ω
  kz = newberger.kz
  k⊥ = newberger.k⊥
  dfdvz = DualNumbers.dualpart(S(Dual(vz, 1), v⊥))
  dfdv⊥ = DualNumbers.dualpart(S(vz, Dual(v⊥, 1)))

  T = promote_type(typeof.((dfdvz, dfdv⊥, ω, Ω, kz, k⊥))...)
  (iszero(dfdvz) && iszero(dfdv⊥)) && return @MArray zeros(T, 3, 3)

  a = (ω - kz * vz) / Ω
  π_sinπa = π / sinpi(a)
  @assert isfinite(π_sinπa)
  z = k⊥ * v⊥ / Ω

  Ja = besselj(a, z)
  J_a = besselj(-a, z)
  Jad = (besselj(a - 1, z) - besselj(a + 1, z)) / 2
  Jad = isfinite(Jad) ? Jad : besselj(a - 1, z) - besselj(a, z) * a / z
  Jad = isfinite(Jad) ? Jad : besselj(a, z) * a / z - besselj(a + 1, z)
  Jad = isfinite(Jad) ? Jad : DualNumbers.dualpart(besselj(a, Dual(z, 1)))
  J_ad = (besselj(-a - 1, z) - besselj(-a + 1, z)) / 2
  J_ad = isfinite(J_ad) ? J_ad : besselj(-a - 1, z) + besselj(-a, z) * a / z
  J_ad = isfinite(J_ad) ? J_ad : -besselj(-a, z) * a / z - besselj(-a + 1, z)
  J_ad = isfinite(J_ad) ? J_ad : DualNumbers.dualpart(besselj(-a, Dual(z, 1)))

  isfinite(Ja + J_a + Jad + J_ad) || @show a, z, Ja, J_a, Jad, J_ad

  @assert isfinite(Ja)
  @assert isfinite(J_a)
  @assert isfinite(Jad)
  @assert isfinite(J_ad)
  Q_a = π_sinπa * J_a * Ja # Eq 33
  Qd_a = π_sinπa * (J_ad * Ja + J_a * Jad) # Eq 33
  Xzz = 2π * Ω * vz * (v⊥ * dfdvz - vz * dfdv⊥) / Ω # Part of Eq 34 (x'ed by ω/Ω)
  U = (kz * v⊥ * dfdvz + (ω - kz * vz) * dfdv⊥) / Ω # Eq 4 (multiplied by ω/Ω)
  T11 = a / (k⊥ / Ω)^2 * (a * Q_a - 1)
  T12 = im / 2z * a * Qd_a * v⊥^2
  T13 = (a * Q_a - 1) / (k⊥ / Ω) * vz
  T22 = (π_sinπa * J_ad * Jad * v⊥^2 + a / (k⊥ / Ω)^2)
  T23 = - vz * im / 2 * Qd_a * v⊥
  T33 = Q_a * vz^2
  T21 = -T12
  T31 = T13
  T32 = -T23
  T = @MArray [T11 T12 T13; T21 T22 T23; T31 T32 T33]
  @assert all(isfinite, T) "$T"
  Xij = (2π * U) .* T # Eq 34, part
  Xij[3, 3] += Xzz
  return Xij # Eq 34 (U is multiplied by ω)
end

function coupledvelocity(S::AbstractCoupledVelocitySpecies, C::Configuration)
  ω, Ω = C.frequency, S.Ω
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)
  @assert !iszero(k⊥) "Perpendicular wavenumber must not be zero"

  integrand = Newberger(S, ω, Ω, kz, k⊥)

  lower = max(S.F.lower, eps())

  function integral2D()
    return first(HCubature.hcubature(integrand,
      (-S.F.upper, lower), (S.F.upper, S.F.upper), initdiv=16,
      rtol=C.options.quadrature_tol.rel, atol=C.options.quadrature_tol.abs))
  end

  function principal()
    @assert !iszero(kz)
    @assert iszero(imag(kz)) || iszero(imag(ω))
    function allprincipals(n)
      harmonicsum = HarmonicSum(S, ω, Ω, kz, k⊥, n)
      pole = Pole(C.frequency, C.wavenumber, n, S.Ω)
      @assert iszero(imag(pole))
      principalintegrand(vz, v⊥) = - numerator(harmonicsum, (vz, v⊥)) / kz
      folded = foldnumeratoraboutpole(principalintegrand, real(float(pole)))
      output = first(HCubature.hcubature(folded,
        (0.0, S.F.lower), (S.F.upper, S.F.upper), initdiv=16,
        rtol=C.options.quadrature_tol.rel, atol=C.options.quadrature_tol.abs))
      @assert !any(isnan, output)
      return output
    end
    return converge(allprincipals, C.options.summation_tol)
  end

  function coupledresidue(v⊥)
    # this started life in relativistic version - can it be simplified?
    function allresidues(n)
      pole = Pole(C.frequency, C.wavenumber, n, S.Ω)
      polefix = wavedirectionalityhandler(pole)
      residuesigma(polefix(pole)) == 0 && return 0.0
      rpradius = (iszero(imag(pole)) ? abs(pole) : abs(imag(pole))) * sqrt(eps())
      output = residuepartadaptive(vz->integrand((vz, v⊥)),
        pole, rpradius, 64, C.options.quadrature_tol, Nmax=2^12)
      output = polefix.(residue(output, polefix(pole)))
      output = sign(real(kz)) .* real(output) .+ im .* imag(output)
      @assert !any(isnan, output)# "v⊥ = $v⊥, pp = $pp, pole = $pole"
      return output
    end
    return converge(allresidues, C.options.summation_tol)
  end

  function perpendicularintegral(∫dv⊥::T, nrm=1) where T
    return first(QuadGK.quadgk(∫dv⊥, S.F.lower, S.F.upper, order=7,
      atol=max(C.options.quadrature_tol.abs,
               C.options.quadrature_tol.rel * nrm / 2),
      rtol=C.options.quadrature_tol.rel))
  end

  # if logic here is a confusing!
  result = if iszero(kz) || (!isreal(ω) || !isreal(kz))
    i2d = integral2D()
    if !iszero(kz) && !iszero(imag(ω))
      i2d .+= perpendicularintegral(coupledresidue, norm(i2d))
    end
    i2d
  else
    @assert isreal(ω) && isreal(kz)
    @assert !iszero(kz) # obviously
    @warn "Coupled species calculations with zero imaginary pole coupled species are slow and inaccurate"
    pp = principal()
    pp .+ perpendicularintegral(coupledresidue, norm(pp))
  end
  return result
end

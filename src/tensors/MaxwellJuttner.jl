using QuadGK, SeriesAccelerators, StaticArrays

struct MaxwellJuttnerIntegrand{T1, T2, T3<:Real, T4, T5, T6}
  m::T1
  ω::T2
  Ω::T3
  kz::T4
  k⊥::T5
  μ::T6
  invbesselkx2μ::T6
end

const COEFFS_K2 = (
                   1.0,
                   1.875,
                   0.8203125,
                   -0.3076171875,
                   0.317230224609375,
                   -0.5154991149902344,
                   1.127654314041138,
                   -3.080912679433823,
                   10.06110546877608,
                   -38.14835823577596,
                   164.5147948917838,
                   -794.5316798750923,
                   4245.778664332524
                  )

const COEFFS_K3 = (
                   1.0,
                   4.375,
                   7.3828125,
                   3.3837890625,
                   -1.374664306640625,
                   1.546497344970703,
                   -2.73858904838562,
                   6.504148989915848,
                   -19.20756498584524,
                   67.49324918637285,
                   -274.1913248196397,
                   1261.903256272205,
                   -6480.399013981222
                  )

"""
    besselk_ratio(ν, sqrtR, μ)

Return `K_ν(√R) / K₂(μ)` computed via exponentially-scaled Bessel functions so
that neither numerator nor denominator underflows for large arguments.

Mathematically equivalent to `besselk(ν, sqrtR) / besselk(2, μ)`, but evaluable
when `μ` is large enough that `besselk(2, μ)` itself underflows to zero.
"""
@inline function besselk_ratio(ν::Integer, sqrtR::Number, rmi)
  @assert ν == 2 || ν == 3
  besselkxνsqrtR = if abs(sqrtR) > 156
    evalpoly(1 / sqrtR, ν == 2 ? COEFFS_K2 : COEFFS_K3) * sqrt(pi/2sqrtR)
  else
    besselkx(ν, sqrtR)
  end
  output = exp(rmi.μ - sqrtR) * besselkxνsqrtR * rmi.invbesselkx2μ
  @assert isfinite(output) (sqrtR, exp(rmi.μ - sqrtR), besselkxνsqrtR, rmi.invbesselkx2μ)
  return output
end

# See Swanson section 4.7.1, Eqs 4.243-245
function (rmi::MaxwellJuttnerIntegrand)(ξ)
  m, ω, Ω, kz, k⊥ = rmi.m, rmi.ω, rmi.Ω, rmi.kz, rmi.k⊥
  sinξ, cosξ = sincos(ξ)
  T1 = @SArray [cosξ sinξ 0; -sinξ cosξ 0; 0 0 1] # transpose of Swanson Eq 4.244
  @cse begin
    Qxx = k⊥^2 * sinξ^2
    Qxy = k⊥^2 * sinξ * (1 - cosξ) # note negative of Swanson Eq 4.245 entry
    Qxz = k⊥ * kz * ξ * sinξ
    Qyy = -k⊥^2 * (1 - cosξ)^2
    Qyz = -k⊥ * kz * ξ * (1 - cosξ) # note negative of Swanson Eq 4.245 entry
    Qzz = kz^2 * ξ^2
  end
  T2 = (c₀ / Ω)^2 * @SArray [Qxx Qxy Qxz; -Qxy Qyy Qyz; Qxz -Qyz Qzz]
  R = ((rmi.μ * Ω - im * ξ * ω)^2 + 2 * (k⊥ * c₀)^2 * (1 - cosξ) + (kz * c₀ * ξ)^2) / Ω^2
  @assert !iszero(R)
  sqrtR = sqrt(R)
  output = (besselk_ratio(2, sqrtR, rmi) * T1 - besselk_ratio(3, sqrtR, rmi) / sqrtR * T2) / R
  @assert all(isfinite, output) (ξ, R, sqrtR)
  return output
end

function maxwelljuttner(species::MaxwellJuttnerSpecies, config::Configuration)
  ω = config.frequency
  kz = para(config.wavenumber)
  k⊥ = perp(config.wavenumber)
  μ = species.m * c₀^2 / q₀ / species.TeV
  if imag(ω) < 0 && iszero(kz)
    throw(ArgumentError("Negative imaginary frequency and zero kz not supported by this method"))
  end

  besselkx2μ = besselkx(2, μ)
  rmi = MaxwellJuttnerIntegrand(species.m, ω, species.Ω, kz, k⊥, μ, 1 / besselkx2μ)

  σ = sign(real(species.Ω)) # Swanson's derivation uses ξ = |Ω|t
  atol = config.options.quadrature_tol.abs
  chunkrtol = config.options.quadrature_tol.rel / 2
  scale = Ref(0.0)
  function period(n)
    atol = max(config.options.quadrature_tol.abs, chunkrtol * scale[])
    o = QuadGK.quadgk(rmi, σ * 2n * π, σ * (n + 1) * 2π;
                      atol=atol, rtol=chunkrtol, order=25, norm=quadnorm)[1]
    n == 0 && (scale[] = quadnorm(o))
    return o
  end

  igrand, isconverged = SeriesAccelerators.shanks(period, 4, 10_000;
                                                  atol=config.options.quadrature_tol.abs,
                                                  rtol=config.options.quadrature_tol.rel)
  isconverged || @warn "Maxwell-Juttner time integral did not converge" maxlog=1
  @assert all(isfinite, igrand)
  factor = im * ω / species.Ω * μ^2
  @assert isfinite(factor)
  output = igrand * factor
  return output
end



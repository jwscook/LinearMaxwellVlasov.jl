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
    @assert !any(isnan, output)# "output = $output, vz⊥=$vz⊥, $dfdvz, $dfdv⊥"

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
    @assert !any(isnan, output)# "v⊥ = $v⊥, output = $output"
    return output
  end

  function principal(v⊥)
    @assert !iszero(kz)
    ∫dvz_kz(x) = - numerator((x, v⊥)) / kz
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
    ppradius = (iszero(imag(pole)) ? abs(pole) : abs(imag(pole))) * sqrt(eps())
    pp = principalpartadaptive(∫dvz, pole, ppradius, 64,
      C.options.summation_tol, Nmax=2048)
    output = polefix.(residue(pp, polefix(pole)))
    output = sign(real(kz)) .* real(output) .+ im .* imag(output)
    @assert !any(isnan, output)# "v⊥ = $v⊥, pp = $pp, pole = $pole"
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


function coupledvelocity(S::AbstractCoupledVelocitySpecies, C::Configuration)
  @show "Newberger!"

  ω, Ω = C.frequency, S.Ω
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)
  @assert !iszero(k⊥) "Perpendicular wavenumber must not be zero"
  ψ = angle(C.wavenumber)

  integrand(vz⊥) = integrand(vz⊥[1], vz⊥[2])
  function integrand(vz, v⊥)
    a = (ω - kz * vz) / Ω
    b = k⊥ * v⊥ / Ω
    c = kz * vz / ω
    Ja = besselj(a, b)
    J_a = besselj(-a, b)
    Jad = (besselj(a - 1, b) - besselj(a + 1, b)) / 2
    J_ad = (besselj(-a - 1, b) - besselj(-a + 1, b)) / 2
    π_sinπa = π / sinpi(a)
    dfdvz = DualNumbers.dualpart(S(Dual(vz, 1), v⊥))
    dfdv⊥ = DualNumbers.dualpart(S(vz, Dual(v⊥, 1)))
    F⊥ = v⊥ * (dfdv⊥ * (1 - c) + c * dfdvz)
    ω_Ω = ω / Ω
    A0 = ω_Ω * π_sinπa * Jad * J_ad * F⊥
    A1 = ω_Ω * π_sinπa / b^2 * Ja * J_a * F⊥
    A2 = ω_Ω * π_sinπa * a / b^2 * Ja * J_ad * F⊥
    A⊥ = ω_Ω * a / b^2 * F⊥
    A3 = ω_Ω * π_sinπa * Ja * J_a
    A4 = ω_Ω * π_sinπa * a * Ja * J_a * F⊥
    A5 = ω_Ω * π_sinπa * Ja * J_ad * F⊥
    A6 = A4
    A7 = A5
    sinψ, cosψ = sincos(ψ)
    sin²ψ = sinψ^2
    cos²ψ = 1 - sin²ψ
    sin2ψ, cos2ψ = sincos(2ψ)
    Kxx = A0 * sin²ψ + A1 * cos²ψ - A⊥ * cos2ψ
    # Kxx = A1 + sin²ψ * A0
    Kyy = A0 * cos²ψ + A1 * sin²ψ + A⊥ * cos2ψ
    Kxy = im * A2 + (A1 - A0) * sin2ψ / 2 + A⊥ * (im - sin2ψ)
    Kxz = cosψ * A4 + sinψ * A5
    Kyz = -im * A2 + (A1 - A0) * sin2ψ / 2 - A⊥ * (im + sin2ψ)
    Kzz = A3
    return @MArray [Kxx Kxy Kxz; -Kxy Kyy Kyz; -Kxz -Kyz Kzz]
  end

  pole = Pole(C.frequency, C.wavenumber, 0, S.Ω)
  polefix = wavedirectionalityhandler(pole)
  function integral2D()
    ∫dvrdθ(vrθ) = vrθ[1] * integrand(parallelperpfrompolar(vrθ))
    return first(HCubature.hcubature(∫dvrdθ,
      (S.F.lower, -π / 2), (S.F.upper, π / 2), initdiv=64,
      rtol=C.options.quadrature_tol.rel, atol=C.options.quadrature_tol.abs))
  end

  function principalzerokz(v⊥)
    @assert iszero(kz)
    ∫dvz(x) = integrand((x, v⊥))
    output = first(QuadGK.quadgk(∫dvz, -S.F.upper, S.F.upper, order=32,
      atol=C.options.quadrature_tol.abs,
      rtol=C.options.quadrature_tol.rel / 10))
    @assert !any(isnan, output)# "v⊥ = $v⊥, output = $output"
    return output
  end

  function principal(v⊥)
    @assert !iszero(kz)
    ∫folded = foldnumeratoraboutpole(x->integrand((x, v⊥)), float(pole))
    output = first(QuadGK.quadgk(∫dvz_kz_folded, S.F.lower, S.F.upper, order=32,
        atol=C.options.quadrature_tol.abs,
        rtol=max(eps(), C.options.quadrature_tol.rel / 10)))
    @assert !any(isnan, output)# "v⊥ = $v⊥, output = $output"
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
    @assert !any(isnan, output)# "v⊥ = $v⊥, pp = $pp, pole = $pole"
    return output
  end

  function integralsnested1D(∫dv⊥::T, nrm=1) where T
    return first(QuadGK.quadgk(∫dv⊥, S.F.lower, S.F.upper, order=32,
      atol=max(C.options.quadrature_tol.abs,
               C.options.quadrature_tol.rel * nrm / 2),
      rtol=C.options.quadrature_tol.rel))
  end

  @show isreal(pole), iszero(kz)
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

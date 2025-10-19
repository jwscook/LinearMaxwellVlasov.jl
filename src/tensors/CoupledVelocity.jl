using CommonSubexpressions, DualNumbers, HCubature, LinearAlgebra, QuadGK
using StaticArrays, SpecialFunctions
using GeneralBesselj

(igrand::AbstractCoupledIntegrand)(vzv⊥) = igrand(vzv⊥[1], vzv⊥[2])

struct NewbergerClassical{S,T,W<:Wavenumber} <: AbstractCoupledIntegrand
  species::S
  ω::T
  k::W
  count::Ref{Int}
end
NewbergerClassical(s, ω, k::Wavenumber) = NewbergerClassical(s, ω, k, Ref(0))

(nc::NewbergerClassical)(vz, v⊥) = numerator(nc, vz, v⊥) / denominator(nc, vz, v⊥)

function denominator(nc::NewbergerClassical, vz, v⊥)
  a = pseudoharmonic(nc, vz)
  sinπa = sinpi(a)
  @assert !isnan(sinπa) "a = $a, vz = $vz"
  return sinπa
end

function pseudoharmonic(nc::NewbergerClassical, vz)
  ω = nc.ω
  Ω = nc.species.Ω
  kz = para(nc.k)
  a = (ω / Ω - kz * vz / Ω)
  return a
end

function numerator(nc::NewbergerClassical, vz, v⊥)
  nc.count[] += 1

  a = pseudoharmonic(nc, vz)
  sinπa = sinpi(a)

  kz = para(nc.k)

  S = nc.species
  Ω = S.Ω
  k⊥ = perp(nc.k)

  z = k⊥ * v⊥ / Ω

  f = S(real(vz), v⊥)
  T = promote_type(typeof.((f, a, z))...)
  # if the real part of f is zero, then there are no particles at (vz, v⊥)
  iszero(f) && return @MArray zeros(T, 3, 3)

  dfdvz = DualNumbers.dualpart(S(Dual(vz, 1), v⊥))
  @assert !isnan(dfdvz)
  dfdv⊥ = DualNumbers.dualpart(S(vz, Dual(v⊥, 1)))
  @assert !isnan(dfdv⊥)
  (iszero(dfdvz) && iszero(dfdv⊥)) && return @MArray zeros(T, 3, 3)

  Jadual, J_adual = besselj_v(MVector(a, -a), Dual(z, 1); maxiters=2^20)
  Ja, Jad = DualNumbers.realpart(Jadual), DualNumbers.dualpart(Jadual)
  J_a, J_ad = DualNumbers.realpart(J_adual), DualNumbers.dualpart(J_adual)
  @assert !isnan(Ja)
  @assert !isnan(J_a)
  @assert !isnan(Jad)
  @assert !isnan(J_ad)

  #@cse begin
    Q_a = π * J_a * Ja # Eq 33
    Qd_a = π * (J_ad * Ja + J_a * Jad) # Eq 33
    Xzz = 2π * Ω * vz * (v⊥ * dfdvz - vz * dfdv⊥) / Ω # Part of Eq 34 (x'ed by ω/Ω)
    U = (v⊥ * kz / Ω * dfdvz + a * dfdv⊥) # Eq 4 (multiplied by ω/Ω)
    T11 = a * (Ω / k⊥)^2 * (a * Q_a - sinπa)
    T12 = im / 2z * a * Qd_a * v⊥^2
    T13 = (a * Q_a - sinπa) * (Ω / k⊥) * vz 
    T22 = (π * J_ad * Jad * v⊥^2 + sinπa * a * (Ω / k⊥)^2)
    T23 = - vz * im / 2 * Qd_a * v⊥
    T33 = Q_a * vz^2
    T21 = -T12
    T31 = T13
    T32 = -T23
  #end
  Tij = @MArray [T11 T12 T13; T21 T22 T23; T31 T32 T33]
  Xij = (2π * U) .* Tij # Eq 34, part
  Xij[3, 3] += Xzz * sinπa
  @assert !any(isnan, Xij) (Xij, vz, v⊥)
  return Xij # Eq 34 (U is multiplied by ω)
end

function coupledvelocity(S::AbstractCoupledVelocitySpecies, C::Configuration)
  ω, Ω = C.frequency, S.Ω
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)

  cubaatol = C.options.cubature_tol.abs
  cubartol = C.options.cubature_tol.rel
  nc = NewbergerClassical(S, ω, C.wavenumber)

  deformation = imagcontourdeformation(ω / kz,
                                       real(kz) >= 0 ? 1 : -1,
                                       C.options.cauchydeformationangle)

  function robustintegral2D()
    nc.count[] = 0

    t1 = @elapsed output, integral2Derrorestimate = if S.F.lower == 0
      HCubature.hcubature(vz⊥ -> nc((vz⊥[1] + im * deformation, vz⊥[2])),
        (-S.F.upper, 0.0), (S.F.upper, S.F.upper), initdiv=32,
        rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    else
      @assert S.F.lower > 0
      ∫dvrdθ(vrθ) = vrθ[1] * nc(parallelperpfrompolar(vrθ) .+ (im * deformation, 0))
      HCubature.hcubature(∫dvrdθ,
        (S.F.lower, -π / 2), (S.F.upper, π / 2), initdiv=32,
        rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    end

    if C.options.erroruponcubaturenonconformance
      msg = "error / val = $(integral2Derrorestimate / norm(output))"
      msg *= ", count = $(nc.count[]), time=$t1 seconds"
      @assert ((nc.count[] < C.options.cubature_maxevals) ||
        integral2Derrorestimate < max(cubartol * norm(output), cubaatol)) msg
    end
    return output, deformation
  end

  function residueperharmonic(n, firstpart::T0)::T0 where T0
    pole = Pole(C.frequency, C.wavenumber, n, Ω, deformation)
    @assert pole.deformation == deformation
    function inner(v⊥)
      laurentnumerator(x) = -(-1)^n * Ω * numerator(nc, x, v⊥) / kz / π
      output = residue(laurentnumerator, pole)
      @assert !any(isnan, output)
      return output
    end
    lv⊥ = sqrt(max(S.F.lower^2 - real(pole)^2, 0.0))
    uv⊥ = sqrt(max(S.F.upper^2 - real(pole)^2, 0.0))
    lv⊥ >= uv⊥ && return zero(T0)
    return first(QuadGK.quadgk(inner, lv⊥, uv⊥, order=DEFAULT_QUADORDER_PERP,
      atol=max(cubaatol, cubartol * norm(firstpart)), rtol=cubartol))
  end

  function robustresidue(firstpart)
    res = converge(n->residueperharmonic(n, firstpart), minharmonics(S), C.options.cubature_tol)
  end

  t1 = @elapsed result, deformation = robustintegral2D()
  t2 = @elapsed if !iszero(kz)
    result += robustresidue(result)
  end

  return result
end

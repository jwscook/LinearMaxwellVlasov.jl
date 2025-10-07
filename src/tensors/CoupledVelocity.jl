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
  @assert !isnan(sinπa) (a, vz)
  return sinπa
end

function pseudoharmonic(nc::NewbergerClassical, vz)
  ω = nc.ω
  Ω = nc.species.Ω
  kz = para(nc.k)
  a = (ω - kz * vz) / Ω
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

  dfdvz = DualNumbers.dualpart(S(Dual(vz, 1), v⊥))
  @assert !isnan(dfdvz) (dfdvz, vz, v⊥)
  dfdv⊥ = DualNumbers.dualpart(S(vz, Dual(v⊥, 1)))
  @assert !isnan(dfdv⊥) (dfdv⊥, vz, v⊥)

  z = k⊥ * v⊥ / Ω

  T = promote_type(typeof.((dfdvz, dfdv⊥, a, z))...)
  (iszero(dfdvz) && iszero(dfdv⊥)) && return @MArray zeros(T, 3, 3)

  Jadual, J_adual = besselj_v(MVector(a, -a), Dual(z, 1))
  Ja, Jad = DualNumbers.realpart(Jadual), DualNumbers.dualpart(Jadual)
  J_a, J_ad = DualNumbers.realpart(J_adual), DualNumbers.dualpart(J_adual)
  @assert !isnan(Ja)
  @assert !isnan(J_a)
  @assert !isnan(Jad)
  @assert !isnan(J_ad)

  @cse begin
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
  end
  Tij = @MArray [T11 T12 T13; T21 T22 T23; T31 T32 T33]
  @assert all(!isnan, Tij) Tij
  Xij = (2π * U) .* Tij # Eq 34, part
  Xij[3, 3] += Xzz * sinπa
  return Xij # Eq 34 (U is multiplied by ω)
end

function numeratorintegera(nc::NewbergerClassical, vz, v⊥)
  nc.count[] += 1

  a = round(Int, pseudoharmonic(nc, vz))

  kz = para(nc.k)

  S = nc.species
  Ω = S.Ω
  k⊥ = perp(nc.k)

  dfdvz = DualNumbers.dualpart(S(Dual(vz, 1), v⊥))
  @assert !isnan(dfdvz) (dfdvz, vz, v⊥)
  dfdv⊥ = DualNumbers.dualpart(S(vz, Dual(v⊥, 1)))
  @assert !isnan(dfdv⊥) (dfdv⊥, vz, v⊥)

  z = k⊥ * v⊥ / Ω

  T = promote_type(typeof.((dfdvz, dfdv⊥, a, z))...)
  (iszero(dfdvz) && iszero(dfdv⊥)) && return @MArray zeros(T, 3, 3)

  Jadual= besselj(a, Dual(z, 1))
  Ja, Jad = DualNumbers.realpart(Jadual), DualNumbers.dualpart(Jadual)
  @assert !isnan(Ja)
  @assert !isnan(Jad)
  J_a, J_ad = (-1)^a .* (Ja, Jad)

  @cse begin
    Q_a = π * J_a * Ja # Eq 33
    Qd_a = π * (J_ad * Ja + J_a * Jad) # Eq 33
    Xzz = 2π * Ω * vz * (v⊥ * dfdvz - vz * dfdv⊥) / Ω # Part of Eq 34 (x'ed by ω/Ω)
    U = (v⊥ * kz / Ω * dfdvz + a * dfdv⊥) # Eq 4 (multiplied by ω/Ω)
    T11 = a * (Ω / k⊥)^2 * (a * Q_a)
    T12 = im / 2z * a * Qd_a * v⊥^2
    T13 = (a * Q_a) * (Ω / k⊥) * vz 
    T22 = (π * J_ad * Jad * v⊥^2)
    T23 = - vz * im / 2 * Qd_a * v⊥
    T33 = Q_a * vz^2
    T21 = -T12
    T31 = T13
    T32 = -T23
  end
  Tij = @MArray [T11 T12 T13; T21 T22 T23; T31 T32 T33]
  @assert all(!isnan, Tij) Tij
  Xij = (2π * U) .* Tij # Eq 34, part
  return Xij # Eq 34 (U is multiplied by ω)
end

function coupledvelocity(S::AbstractCoupledVelocitySpecies, C::Configuration)
  ω, Ω = C.frequency, S.Ω
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)

  cubaatol = C.options.cubature_tol.abs
  cubartol = C.options.cubature_tol.rel
  nc = NewbergerClassical(S, ω, C.wavenumber)

  trickiestn = round(Int, real(ω / Ω))
  deformation = imagcontourdeformation((ω - trickiestn * Ω) / kz, real(kz) >= 0 ? 1 : -1,
                                       C.options.cauchydeformationangle)

  function robustintegral2D()
    nc.count[] = 0

    t1 = @elapsed output, integral2Derrorestimate = if S.F.lower == 0
      HCubature.hcubature(vz⊥ -> nc((vz⊥[1] * (1 + 0im) + im * deformation, vz⊥[2])),
        (-S.F.upper, 0.0), (S.F.upper, S.F.upper), initdiv=16,
        rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    else
      @assert S.F.lower > 0
      ∫dvrdθ(vrθ) = vrθ[1] * nc(parallelperpfrompolar(vrθ) .+ (im * deformation, zero(vrθ[2])))
      HCubature.hcubature(∫dvrdθ,
        (S.F.lower, -π / 2), (S.F.upper, π / 2), initdiv=16,
        rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    end

    if C.options.erroruponcubaturenonconformance
      msg = "error / val = $(integral2Derrorestimate / norm(output))"
      msg *= ", count = $(nc.count[]), time=$t1"
      @assert ((nc.count[] < C.options.cubature_maxevals) ||
        integral2Derrorestimate < max(cubartol * norm(output), cubaatol)) msg
    end
    return output
  end

  function trickintegrand2D(yv⊥, n, part)
    y, v⊥ = yv⊥
    vz = (ω - (y + n) * Ω) / kz * (1 + 0im) + im * deformation
    num = numerator(nc, vz, v⊥)
    den = denominator(nc, vz, v⊥)
    return (num - part) / den
  end
  function trickintegrandpart(v⊥, n)
    vzn = (ω - n * Ω) / kz * (1 + 0im) + im * deformation
    return numeratorintegera(nc, vzn, v⊥)
  end
  trickintegrand1D(v⊥, n)= -trickintegrandpart(v⊥, n) * im * (-1)^n

  function robustintegral2Dnew()
    nc.count[] = 0

    # x = (ω - kv) / Ω
    hmax = ceil(Int, max(abs(real((ω - kz * S.F.upper) / Ω)),
                         abs(real((ω + kz * S.F.upper) / Ω))))

#    t1 = @elapsed output, integral2Derrorestimate = HCubature.hcubature(
#      vz⊥ -> sum(trickintegrand2D(vz⊥, n) for n in -hmax:hmax),
#      (-0.5, 0.0), (0.5, S.F.upper), initdiv=16,
#      rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
#
#    if C.options.erroruponcubaturenonconformance
#      msg = "error / val = $(integral2Derrorestimate / norm(output))"
#      msg *= ", count = $(nc.count[]), time=$t1"
#      @assert ((nc.count[] < C.options.cubature_maxevals) ||
#        integral2Derrorestimate < max(cubartol * norm(output), cubaatol)) msg
#    end
#
#    output += QuadGK.quadgk(v⊥->sum(trickintegrand1D(v⊥, n) for n in -hmax:hmax),
#      0.0, S.F.upper, order=DEFAULT_QUADORDER, atol=cubaatol, rtol=cubartol)[1]

    t1 = @elapsed output = converge(n->
      HCubature.hcubature(vz⊥->trickintegrand2D(vz⊥, n, trickintegrandpart(vz⊥[2], n)),
      (-0.5, 0.0), (0.5, S.F.upper), initdiv=1,
      rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)[1] +
      QuadGK.quadgk(v⊥->trickintegrand1D(v⊥, n),
      S.F.lower, S.F.upper, order=7, atol=cubaatol, rtol=cubartol)[1],
      hmax, C.options.summation_tol)

    return output * Ω / kz
  end


  function perpendicularintegral(∫dv⊥::T, nrm=1) where T
    return first(QuadGK.quadgk(∫dv⊥, S.F.lower, S.F.upper, order=7,
      atol=max(C.options.quadrature_tol.abs, C.options.quadrature_tol.rel * nrm / 2),
      rtol=C.options.quadrature_tol.rel))
  end

  function robustcoupledresidue(v⊥, ::Type{T0}, deformation)::T0 where T0
    function allresidues(n)
      pole = Pole(C.frequency, C.wavenumber, n, Ω, deformation)
      @assert pole.deformation == deformation
      laurentnumerator(x) = -(-1)^n * Ω * numerator(nc, x, v⊥) / kz / π
      output = residue(laurentnumerator, pole)
      @assert !any(isnan, output)
      return output
    end
    return converge(allresidues, minharmonics(S), C.options.summation_tol)
  end

  t1 = @elapsed robustintegral = robustintegral2D()
  t2 = @elapsed robustintegralnew = robustintegral2Dnew()
  @show norm((robustintegralnew .- robustintegral) ./ robustintegral), t2 / t1
  result = robustintegral
  t2 = @elapsed if !iszero(kz)
    res = perpendicularintegral(
      v⊥->robustcoupledresidue(v⊥, typeof(robustintegral), deformation),
      norm(robustintegral))
    result += res
  end

  return result
end

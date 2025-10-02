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
  # Don't multiply by the sign changer vz *= nc.k.multipliersign
  a = pseudoharmonic(nc, vz)
  sinπa = sinpi(a)
  @assert isfinite(sinπa) "a = $a, vz = $vz"
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

  # now change the sign if required
#  ms = nc.k.multipliersign
  kz = para(nc.k)
  @assert kz >= 0
  
  a = pseudoharmonic(nc, vz) # vz sign unchanged
  sinπa = sinpi(a) # vz sign unchanged

  S = nc.species
  Ω = S.Ω
  k⊥ = perp(nc.k)

  dfdvz = DualNumbers.dualpart(S(Dual(vz, 1), v⊥))
  @assert isfinite(dfdvz)
  dfdv⊥ = DualNumbers.dualpart(S(vz, Dual(v⊥, 1)))
  @assert isfinite(dfdv⊥)

  z = k⊥ * v⊥ / Ω

  T = promote_type(typeof.((dfdvz, dfdv⊥, a, z))...)
  (iszero(dfdvz) && iszero(dfdv⊥)) && return @MArray zeros(T, 3, 3)

  Jadual, J_adual = besselj_v(MVector(a, -a), Dual(z, 1))
  Ja, Jad = DualNumbers.realpart(Jadual), DualNumbers.dualpart(Jadual)
  J_a, J_ad = DualNumbers.realpart(J_adual), DualNumbers.dualpart(J_adual)
  @assert isfinite(Ja)
  @assert isfinite(J_a)
  @assert isfinite(Jad)
  @assert isfinite(J_ad)

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
  @assert all(isfinite, Tij) Tij
  Xij = (2π * U) .* Tij # Eq 34, part
  Xij[3, 3] += Xzz * sinπa
  return Xij # Eq 34 (U is multiplied by ω)
end

function coupledvelocity(S::AbstractCoupledVelocitySpecies, C::Configuration)
  ω, Ω = C.frequency, S.Ω
  @assert !iszero(Ω)
  kz, k⊥ = para(C.wavenumber), perp(C.wavenumber)

  cubaatol = C.options.cubature_tol.abs
  cubartol = C.options.cubature_tol.rel
  nc = NewbergerClassical(S, ω, C.wavenumber)

  deformation = imagcontourdeformation(ω / kz)

  function robustintegral2D()
    nc.count[] = 0
    ms = C.wavenumber.multipliersign

    t1 = @elapsed output, integral2Derrorestimate = if S.F.lower == 0
      HCubature.hcubature(vz⊥ -> nc((vz⊥[1]#= * ms=# + im * deformation, vz⊥[2])),
        (-4S.F.upper, 0.0), (4S.F.upper, 4S.F.upper), initdiv=16,
        rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    else
      @assert S.F.lower > 0
      ∫dvrdθ(vrθ) = vrθ[1] * nc(parallelperpfrompolar(vrθ)#= .* (ms, 1)=# + (im * deformation, zero(vrθ[2])))
      HCubature.hcubature(∫dvrdθ,
        (S.F.lower/4, -π / 2), (4S.F.upper, π / 2), initdiv=16,
        rtol=cubartol, atol=cubaatol, maxevals=C.options.cubature_maxevals)
    end

    if C.options.erroruponcubaturenonconformance
      msg = "error / val = $(integral2Derrorestimate / norm(output))"
      msg *= ", count = $(nc.count[]), time=$t1"
      @assert ((nc.count[] < C.options.cubature_maxevals) ||
        integral2Derrorestimate < max(cubartol * norm(output), cubaatol)) msg
    end
    return output, deformation
  end

  function perpendicularintegral(∫dv⊥::T, nrm=1) where T
    return first(QuadGK.quadgk(∫dv⊥, S.F.lower/4, 4S.F.upper, order=7,
      atol=max(C.options.quadrature_tol.abs, C.options.quadrature_tol.rel * nrm / 2),
      rtol=C.options.quadrature_tol.rel))
  end

  function robustcoupledresidue(v⊥, ::Type{T0}, deformation)::T0 where T0
    ms = C.wavenumber.multipliersign
    function allresidues(n)
      pole = Pole(C.frequency, C.wavenumber, n, S.Ω, deformation)
      @assert pole.deformation == deformation
      laurentnumerator(x) = -(-1)^n * Ω * numerator(nc, x, v⊥) / kz / π
      output = residue(laurentnumerator, pole)
      @assert !any(isnan, output)
      return output
    end
    return converge(allresidues, C.options.summation_tol)
  end

  t1 = @elapsed robustintegral, deformation = robustintegral2D()
  result = robustintegral
  t2 = @elapsed if !iszero(kz)
    res = perpendicularintegral(
      v⊥->robustcoupledresidue(v⊥, typeof(robustintegral), deformation),
      norm(robustintegral))
    result += res
  end

  return result
end

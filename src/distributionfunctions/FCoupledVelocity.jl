using StaticArrays, HCubature

struct FCoupledVelocityNumerical{T<:Function, U<:Number
    } <: AbstractCoupledVelocity
  F::T
  normalisation::Tuple{U,U}
  lower::Float64 # minimum speed for integration bounds
  upper::Float64 # maximum speed for integration bounds
  function FCoupledVelocityNumerical(F::T, normalisation::Tuple{U,U},
      lower=0.0, upper=4default_integral_range * norm(normalisation);
      autonormalise::Bool=false) where {T, U}
    output = new{T,U}(F, normalisation, lower, upper)
    if autonormalise
      n = integrate(output)
      newf(vz, v⊥) = F(vz, v⊥) / n
      return new{typeof(newf),U}(newf, normalisation, lower, upper)
    else
      is_normalised(output) || throw(ArgumentError("F is not normalised"))
      return output
    end
  end
end

(f::FCoupledVelocityNumerical)(v) = f.F(v[1], v[2])
(f::FCoupledVelocityNumerical)(vz, v⊥) = f.F(vz, v⊥)

lower(f::FCoupledVelocityNumerical) = f.lower
upper(f::FCoupledVelocityNumerical) = f.upper

struct ShiftedMaxwellianCoupled{T,U,V,W,X} <: Function
  Fz::ShiftedMaxwellianParallel{T,U}
  F⊥::ShiftedMaxwellianPerpendicular{V,W,X}
end
(f::ShiftedMaxwellianCoupled)(v) = f.Fz(v[1]) * f.F⊥(v[2])
(f::ShiftedMaxwellianCoupled)(vz, v⊥) = f.Fz(vz) * f.F⊥(v⊥)

function FCoupledVelocityNumerical(vthz::Real, vth⊥::Real=vthz,
    vdz::Number=0.0, vd⊥::Number=0.0)
  Fz = ShiftedMaxwellianParallel(vthz, vdz)
  F⊥ = ShiftedMaxwellianPerpendicular(vth⊥, vd⊥)
  F = ShiftedMaxwellianCoupled(Fz, F⊥)
  normalisation = (vthz + abs(vdz), vth⊥ + abs(vd⊥))
  return FCoupledVelocityNumerical(F, normalisation)
end

function is_normalised(F::AbstractCoupledVelocity)
  return isapprox(integrate(F), 1, rtol=1e5eps(), atol=0)
end

function integrate(F::AbstractCoupledVelocity)
  objective(v) = v[2] * F(v)
  ∫dvrdθ(vrθ) = vrθ[1] * objective(parallelperpfrompolar(vrθ))
  return 2π * first(HCubature.hcubature(∫dvrdθ,
    [lower(F), -π / 2], [upper(F), π / 2], initdiv=50, maxevals=1_000_000,
    rtol=1e5eps(), atol=0.0))
end

function FShell(vth::Real, vshell::Real)
  @assert vshell >= 0 "Shell speed, vshell, must be >= 0"
  @assert vth >= 0 "Thermal spread of the shell, vth, must be > 0"
  inv2vth² = 1 / 2 / vth^2
  funnormalised = v -> exp(-(sqrt(sum(x->x^2, v)) - vshell)^2 * inv2vth²)
  funnormalisedpolar = transformtopolar(funnormalised)
  function fpolar_for_normalisation(vrθ)
    _, v⊥ = parallelperpfrompolar(vrθ)
    return 2π * v⊥ * vrθ[1] * funnormalisedpolar(vrθ)
  end
  lower = max(0.0, vshell - default_integral_range * vth)
  upper = vshell + default_integral_range * vth
  normalisation = HCubature.hcubature(fpolar_for_normalisation,
    [lower, -π/2], [upper, π/2], rtol=1e5eps(), atol=0.0, initdiv=10)[1]
  @assert isfinite(normalisation) && normalisation > 0
  F(vz, v⊥) = funnormalised(@SArray [vz, v⊥]) / normalisation
  characteristics = (vth + vshell, vth + vshell)
  return FCoupledVelocityNumerical(F, characteristics, lower, upper)
end

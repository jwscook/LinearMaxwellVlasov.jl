using StaticArrays, HCubature

"""

FCoupledVelocityNumerical

A disribution function where vz and v⊥ are coupled, i.e. non-separable.

...
# Arguments
- `F::T`: the distrubtion function
- `normalisation::Tuple{U,U}`: the speeds used for normalisation in parallel and perp [m/s]
- `lower::Float64`: minimum speed for integration bounds [m/s]
- `upper::Float64`: maximum speed for integration bounds [m/s]
...

# Example
```julia
vth = 1e4
vshell = 1e5
fshell = FShell(vth, vshell)
lower = vshell - 12 * vth # where the shell is zero
upper = vshell + 12 * vth # where the shell is zero
FCoupledVelocityNumerical(fshell, (vshell, vshell), lower, upper)
```
"""
struct FCoupledVelocityNumerical{T<:Function, U<:Number
    } <: AbstractCoupledVelocity
  F::T
  normalisation::Tuple{U,U}
  lower::Float64 # minimum speed for integration bounds
  upper::Float64 # maximum speed for integration bounds
  function FCoupledVelocityNumerical(F::T, normalisation::Tuple{U,U},
      lower=0.0, upper=default_integral_range * norm(normalisation);
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
# the following line gets missed by the coverage tool, but it is needed!
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


"""
    FShell(vth::Real,vshell::Real)

The shell distribution function, as though f is only non-zero on or
close to the surface of a sphere.

...
# Arguments
- `vth::Real`: the thermal velocity of the shell [m/s]
- `vshell::Real`: the speed of the shell [m/s]
...

# Example
```julia
```
"""
function FShell(vth::Real, vshell::Real)
  @assert vshell >= 0 "Shell speed, vshell, must be >= 0"
  @assert vth >= 0 "Thermal spread of the shell, vth, must be > 0"
  inv2vth² = 1 / 2 / vth^2
  funnormalised = v -> exp(-(sqrt(sum(x->x^2, v)) - vshell)^2 * inv2vth²)

  lower = max(0.0, vshell - default_integral_range * vth)
  upper = vshell + default_integral_range * vth

  funnormalisedpolar = transformtopolar(funnormalised)
  function fpolar_for_normalisation(vrθ)
    _, v⊥ = parallelperpfrompolar(vrθ)
    return 2π * v⊥ * vrθ[1] * funnormalisedpolar(vrθ)
  end
  normalisation = HCubature.hcubature(fpolar_for_normalisation,
    [lower, -π/2], [upper, π/2], rtol=1e5eps(), atol=0.0, initdiv=10)[1]

  @assert isfinite(normalisation) && normalisation > 0
  F(vz, v⊥) = funnormalised(@SArray [vz, v⊥]) / normalisation
  characteristics = (vth + vshell, vth + vshell)
  return FCoupledVelocityNumerical(F, characteristics, lower, upper)
end

"""
    FSlowingDown(vbeam::Real,vcrit::Real,vcutoffwidth::Real)

The slowing down distribution

...
# Arguments
- `vbeam::Real`: the beam speed [m/s]
- `vcrit::Real`: the critival velocity [m/s]
- `vcutoffwidth::Real`: the width of the error function used to smooth
the distribution function at vbeam [m/s]
...

# Example
```julia
```
"""
function FSlowingDown(vbeam::Real, vcrit::Real, vcutoffwidth::Real)
  @assert 0 < vbeam
  @assert 0 < vcrit < vbeam
  @assert 0 < vcutoffwidth
  funnormalised(vz⊥) = funnormalised(vz⊥[1], vz⊥[2])
  function funnormalised(vz, v⊥)
    u = sqrt(vz^2 + v⊥^2)
    return v⊥ / (u^3 + vcrit^3) * (1 + erf((vbeam - u) / vcutoffwidth)) / 2
  end

  upper = vbeam + default_integral_range * vcutoffwidth

  funnormalisedpolar = transformtopolar(funnormalised)
  function fpolar_for_normalisation(vrθ)
    _, v⊥ = parallelperpfrompolar(vrθ)
    return 2π * v⊥ * vrθ[1] * funnormalisedpolar(vrθ)
  end
  normalisation = HCubature.hcubature(fpolar_for_normalisation,
    [0.0, -π/2], [upper, π/2], rtol=1e5eps(), atol=0.0, initdiv=10)[1]

  @assert isfinite(normalisation) && normalisation > 0
  F(vz, v⊥) = funnormalised(@SArray [vz, v⊥]) / normalisation
  characteristics = (vbeam, vbeam)
  return FCoupledVelocityNumerical(F, characteristics, 0.0, upper)
end

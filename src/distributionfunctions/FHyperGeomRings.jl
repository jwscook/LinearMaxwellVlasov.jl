using HypergeometricFunctions
# f = (v/u)^p * exp(-v^2/t^2) / n
# dfdv = (p*(v/u)^(p-1)/u - 2 v / t^2* v^p) / n
# dfdv = (p / v - 2 v/t^2) * (v/u)^p / n
# dfdv = 0 = (p / v - 2 v/t^2)
# p = 2 vd^2 / t^2
# doesn't matter what u is

#using SymPy
#t = SymPy.symbols("t", real=true, positive=true)
#u = SymPy.symbols("u", real=true, positive=true)
#p = SymPy.symbols("p", real=true, positive=true)
#v = SymPy.symbols("v", real=true, positive=true)
# f = (v/u)^p * exp(-v^2/t^2)# / n
#I = SymPy.integrate(2 * SymPy.PI * v * f, (v, 0, Inf))
#   2  p  -p  ⎛p    ⎞
#π⋅t ⋅t ⋅u  ⋅Γ⎜─ + 1⎟
#             ⎝2    ⎠

struct FHyperGeomRing{T<:Number, U<:Number, V<:Number
    } <: AbstractFPerpendicularNumerical
  vth::T # thermal velocity
  vd::U # drift velocity
  α::T
  p::V
  vchar::V
  invnormconst::V
  _uniqueid::UInt64
  function FHyperGeomRing(vth::T, vd::U) where {T<:Number, U<:Number}
    _uniqueid = hash((vth, vd), hash(:FHyperGeomRing))
    @assert vth > 0.0
    vchar = sqrt(vd^2 + vth^2)
    @show t = sqrt(2) * vth # factor of √2 requied to make this similar the FRing
    @show α = 1 / t^2
    p = 2 * vd^2 / t^2
    #invnormconst = 1 / exp(log(π * t^2) + p * log(t/vchar) + lgamma((p + 1)/2))
    invnormconst = 1 / (π * t^2 * (t/vchar)^p * gamma((p + 2)/2))
    @assert isfinite(invnormconst)
    @assert invnormconst > 0
    V = promote_type(T, U)
    return new{T, U, V}(vth, vd, α, p, vchar, invnormconst, _uniqueid)
  end
end

#a*v^p * exp(y) = exp(log(a*v^p * exp(y))) = exp(log(a*v^p) + log(exp(y)))
#exp(log(a*v^p) + log(exp(y))) = exp(log(a*v^p) + y)) =
#exp(log(a*v^p) + y)) = exp(log(a) + p * log(v) + y)) =
function (f::FHyperGeomRing)(v::T, ∂F∂v::Bool=false) where {T<:Number}
  output = exp(-f.α * v^2 + f.p * log(v / f.vchar) + log(f.invnormconst))
  if ∂F∂v
    output *= (f.p / v - 2 * f.α * v)
  end
  return output
end
(f::FHyperGeomRing)(∂F∂v::Bool=false) = v -> f(v, ∂F∂v)

is_normalised(f::FHyperGeomRing) = true
lower(f::FHyperGeomRing) = max(0.0, f.vd - default_integral_range * f.vth)
upper(f::FHyperGeomRing) = f.vd + default_integral_range * f.vth

# Abramovtiz & Stegun Eq. 6.633.5
function AbramSteg66335(λ, μ::Unsigned, ν::Unsigned, α, β)
  @assert real(ν + λ + μ) > 0
  @assert real(α) > 0
  halfsum = (ν+λ+μ)/2
  output = 1 / 2.0^(ν+μ+1) * α^(-halfsum) * β^(ν + μ)
  output *= exp(lgamma(halfsum) - lgamma(μ + 1) - lgamma(ν + 1))
  output *= HypergeometricFunctions.pFq(
    ((ν + μ + 1)/2, (ν + μ + 2)/2, halfsum),
    (μ + 1, ν + 1, μ + ν + 1), -β^2 / α)
  return output
end

function integrate(f::FHyperGeomRing, kernel::T,
    ∂F∂v::Bool, tol::Tolerance=Tolerance()
    ) where {T<:Function}
  β = kernel.k⊥_Ω # units of 1 / velocity, use a normalisation
  μ, ν = kernel.μν
  sgnμ, μ = (μ < 0) ? (isodd(μ) ? -1 : 1, -μ) : (1, μ)
  sgnν, ν = (ν < 0) ? (isodd(ν) ? -1 : 1, -ν) : (1, ν)
  sgn = sgnμ * sgnν
  λ = kernel.power + f.p + 1
  α = f.α
  @show α, β, λ
  #@show α = f.α / β^2 # noramlised α (α has units of 1/velocity^2)
  #@show β /= β # noramlised β
  output = if ∂F∂v
    # f = (v/u)^p * exp(-v^2/t^2)
    # df/dv = (p / v - 2 * v / t^2) * (v/u)^p * exp(-v^2/t^2)
    i1 = AbramSteg66335(λ - 1, Unsigned(μ), Unsigned(ν), α, β) * f.p
    i2 = AbramSteg66335(λ + 1, Unsigned(μ), Unsigned(ν), α, β) * 2f.α # α = 1 / t^2
    (i1 - i2)
  else
    AbramSteg66335(λ, Unsigned(μ), Unsigned(ν), α, β)
  end
  @show output, f.invnormconst, 1 / f.vchar^f.p, sgn
  @show output * f.invnormconst * 1 / f.vchar^f.p * sgn
  return 2π * output * f.invnormconst / f.vchar^f.p * sgn
  #return output * exp(log(f.invnormconst) - f.p * log(f.vchar)) * sgn
end


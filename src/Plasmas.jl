
abstract type AbstractPlasma end

struct NeutralPlasma{T<:AbstractVector{<:AbstractSpecies}} <: AbstractPlasma
  species::T
  function NeutralPlasma(plasma::T) where {T<:AbstractVector{<:AbstractSpecies}}
    isneutral(plasma) || error(ArgumentError("Plasma is not neutral!"))
    return new{T}(plasma)
  end
end

struct Plasma{T<:AbstractVector{<:AbstractSpecies}} <: AbstractPlasma
  species::T
end

import Base.iterate, Base.length
Base.iterate(p::AbstractPlasma) = iterate(p.species)
Base.iterate(p::AbstractPlasma, counter) = iterate(p.species, counter)
Base.length(p::AbstractPlasma) = length(p.species)

function isneutral(vs::AbstractVector{T}, atol=100eps()
    ) where {T<:AbstractSpecies}
  inner(s) = s.Π^2 / s.Ω
  charge_proxy = mapreduce(inner, +, vs) / mapreduce(x->abs(inner(x)), +, vs)
  return isapprox(0, charge_proxy, atol=atol)
end

isneutral(p::AbstractPlasma) = isneutral(p.species)

alfvenspeed(p::AbstractPlasma) = alfvenspeed((s.Π / s.Ω for s ∈ p))

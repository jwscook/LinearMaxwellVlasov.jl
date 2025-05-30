
abstract type AbstractPlasma end

# container may be a vector or a tuple:w
const PlasmaContainerType = Union{AbstractVector{<:AbstractSpecies},
  NTuple{N, AbstractSpecies} where {N}}

"""
NeutralPlasma checks that the plasma is neutral
"""
struct NeutralPlasma{T<:PlasmaContainerType} <: AbstractPlasma
  species::T
  function NeutralPlasma(plasma::T) where {T}
    isneutral(plasma) || error(ArgumentError("Plasma is not neutral!"))
    return new{T}(plasma)
  end
end

"""
Plasma doesn't check if the plasma is neutral
"""
struct Plasma{T<:PlasmaContainerType} <: AbstractPlasma
  species::T
end

import Base.iterate, Base.length
Base.iterate(p::AbstractPlasma) = iterate(p.species)
Base.iterate(p::AbstractPlasma, counter) = iterate(p.species, counter)
Base.length(p::AbstractPlasma) = length(p.species)

"""
    isneutral(vs::AbstractVector{T},atol=100eps()

Calculate if a plasma is neutral

...
# Arguments
- `vs::AbstractVector{T}`: A vector of species
- `atol=100eps()`: tolerance for deciding upon neutrality
...

"""
function isneutral(vs::AbstractVector{T}, atol=100eps()
    ) where {T<:AbstractSpecies}
  inner(s) = s.Π^2 / s.Ω
  charge_proxy = mapreduce(inner, +, vs) / mapreduce(x->abs(inner(x)), +, vs)
  return isapprox(0, charge_proxy, atol=atol)
end

isneutral(p::AbstractPlasma) = isneutral(p.species)

alfvenspeed(p::AbstractPlasma) = alfvenspeed((s.Π / s.Ω for s ∈ p))

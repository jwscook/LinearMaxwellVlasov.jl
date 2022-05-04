mutable struct Configuration{T,U,V,W}
  frequency::T
  wavenumber::Wavenumber{U,V}
  options::Options{W}
end

function Configuration(args...)
  frequency = filter(x->typeof(x) <: Number, args)
  wavenumber = filter(x->typeof(x) <: Wavenumber, args)
  options = filter(x->typeof(x) <: Options, args)
  ω = isempty(frequency) ? ComplexF64(0) : first(frequency)
  k = isempty(wavenumber) ? Wavenumber() : first(wavenumber)
  o = isempty(options) ? Options() : first(options)
  return Configuration(ω, k, o)
end

function Base.hash(c::Configuration)
  return mapreduce(x->hash(getproperty(c, x)), hash, propertynames(c))
end
Base.hash(c::Configuration, h::UInt64) = hash(hash(c), h)

import Base.==
function ==(a::Configuration, b::Configuration)
  return mapreduce(x->getproperty(a, x) == getproperty(b, x), &,
    propertynames(a))
end


# Everything below, except the includes, is for memoising the parallel
# and perpendicular integrals

abstract type AbstractCacheType end
struct ParallelCache <: AbstractCacheType end
struct PerpendicularCache <: AbstractCacheType end

abstract type AbstractCacheOpType end

struct CacheOp{T<:AbstractCacheOpType} <: Function
  a::Bool
  b::Bool
  c::Bool
  CacheOp{T}(a::Bool, b::Bool=false, c::Bool=false) where {T} = new{T}(a, b, c)
end

struct CacheKeyAndOp{T<:AbstractCacheType} <: Function end

include("./integrals/Parallel.jl")
include("./integrals/Perpendicular.jl")


struct CacheDict{C,V} <: AbstractDict{UInt64, V}
  data::Dict{UInt64,V}
  function CacheDict{C}(keyval::Pair{UInt64,V}) where {C<:AbstractCacheType, V}
    return new{C,V}(Dict{UInt64,V}(keyval))
  end
end

mutable struct Cache
  parallel::Dict{UInt64,CacheDict{ParallelCache}}
  perpendicular::Dict{UInt64,CacheDict{PerpendicularCache}}
end

Cache() = Cache(Dict{UInt64,CacheDict{ParallelCache}}(),
                Dict{UInt64,CacheDict{PerpendicularCache}}())

typehash(::Type{ParallelCache}, s) = uniqueid(s.Fz)::UInt64
typehash(::Type{PerpendicularCache}, s) = uniqueid(s.F⊥)::UInt64

"""
    Base.get!(f::F,data::Dict{UInt64,CacheDict{C}},species,config)where{F<:Function,C<:AbstractCacheType}

Return the correct cache dictionary for the given species and configuration.

This creates a CacheKeyOp with the arguments given so that it can create a
CacheDict with the correct type, which then populates the CacheDict if it's empty.
This process acts as a function barrier to get better type inference.

...
# Arguments
- `f::F`: The function to cache
- `data::Dict{UInt64`: the dictionary of caches
- `species`: The species argument to f
- `config`: The configuration argument to f
...
"""
function Base.get!(f::F, data::Dict{UInt64,CacheDict{C}}, species, config
    ) where {F<:Function, C<:AbstractCacheType}
  key = typehash(C, species)
  args = (species, config, 0)
  subkey, cacheop = CacheKeyAndOp{C}(args...)
  return get!(()->CacheDict{C}(subkey=>cacheop(f(args...))), data, key)
end

"""
    Base.get!(f::F,cache::CacheDict{C,V},i...)::Vwhere{F<:Function,C<:AbstractCacheType,V}

description

...
# Arguments
- `f::F`:
- `cache::CacheDict`:
- `i...`:
...

"""
function Base.get!(f::F, cache::CacheDict{C,V}, i...
    )::V where {F<:Function, C<:AbstractCacheType, V}
  key, cacheop = CacheKeyAndOp{C}(i...)
  function g()::V return cacheop(f(i...)) end # to get better inference
  return cacheop(get!(g, cache.data, key))
end

"""
    parallel_integral(species::S,config::Configuration,parallelcaches=Dict{UInt64,CacheDict{ParallelCache}}()) where {S}

Memoise the inputs to the integrals; export the memoised function, and
the dictionary of inputs -> outputs.
This uses a parallel cache which holds the logic for the integrals and
ensure that the values are only calculated once.
It only matters what n-m is, not the values of n and m individually
So change (n, m)->(n-m, 0) to do fewer calculations

...
# Arguments
- `species::S`:  the species to do parallel integrals over
- `config::Configuration`:  the configuration holding frequency and wavenumber, and options.
- `parallelcaches=Dict{UInt64,CacheDict{ParallelCache}}`:  The dict to store parallel integrals
...

# Example
```julia
```
"""
function parallel_integral(species::S, config::Configuration,
    parallelcaches=Dict{UInt64,CacheDict{ParallelCache}}()) where {S}
  if config.options.memoiseparallel
    return constructparallel_memoised(species, config, parallelcaches)
  end
  return parallel
end

@inline function constructparallel_memoised(species, config::Configuration,
    parallelcaches=Dict{UInt64,CacheDict{ParallelCache}}())
  cache = get!(parallel, parallelcaches, species, config)
  parallel_memoised(i...) = get!(parallel, cache, i...)
  return parallel_memoised
end

"""
Memoise the inputs to the integrals; export the memoised function, and
the dictionary of inputs -> outputs
besselj harmonic numbers n and l can be any order
Use this to do only half the calculatons (n, l)->(min(n, l), max(n, l))
Also know that besselj(-n, x) = (-1)^n * besselj(n, x)
"""
function perpendicular_integral(species, config::Configuration,
    perpendicularcaches=Dict{UInt64,CacheDict{PerpendicularCache}}())
  if config.options.memoiseperpendicular
    return constructperpendicular_memoised(species, config, perpendicularcaches)
  end
  return perpendicular
end

@inline function constructperpendicular_memoised(species, config::Configuration,
    perpendicularcaches=Dict{UInt64,CacheDict{PerpendicularCache}}())
  cache = get!(perpendicular, perpendicularcaches, species, config)
  perpendicular_memoised(i...) = get!(perpendicular, cache, i...)
  return perpendicular_memoised
end


# Everything below, except the includes, is for memoising the parallel
# and perpendicular integrals

abstract type AbstractCacheType end
struct ParallelCache <: AbstractCacheType end
struct PerpendicularCache <: AbstractCacheType end
struct GeneralCache <: AbstractCacheType end

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
include("./integrals/GeneralCache.jl")


struct CacheDict{C,V} <: AbstractDict{UInt64, V}
  data::Dict{UInt64,V}
  function CacheDict{C}(keyval::Pair{UInt64,V}) where {C<:AbstractCacheType, V}
    return new{C,V}(Dict{UInt64,V}(keyval))
  end
end

mutable struct Cache
  parallel::Dict{UInt64,CacheDict{ParallelCache}}
  perpendicular::Dict{UInt64,CacheDict{PerpendicularCache}}
  general::Dict{UInt64, CacheDict{GeneralCache}}
end

Cache() = Cache(Dict{UInt64,CacheDict{ParallelCache}}(),
                Dict{UInt64,CacheDict{PerpendicularCache}}(),
                Dict{UInt64,CacheDict{GeneralCache}}())

typehash(::Type{ParallelCache}, s) = uniqueid(s.Fz)::UInt64
typehash(::Type{PerpendicularCache}, s) = uniqueid(s.F⊥)::UInt64
# General contributions don't depend on the plasma frequency; it multiplies
# the result afterwards. So hash every field except Π, and seed with the type
# name rather than the type, because the type carries Π's type parameter.
function typehash(::Type{GeneralCache}, s)
  h = hash(nameof(typeof(s)))
  for name in fieldnames(typeof(s))
    name === :Π && continue
    h = hash(getfield(s, name), h)
  end
  return h::UInt64
end

Base.empty!(c::Cache) = (map(empty!, (c.parallel, c.perpendicular, c.general)); return c)

# The arguments that the cached function is called with when first populating
# a CacheDict. Parallel and perpendicular integrals are per cyclotron harmonic,
# whereas general caches hold the whole tensor for a (species, config) pair.
cacheinitargs(::Type{<:AbstractCacheType}, species, config) = (species, config, 0)
cacheinitargs(::Type{GeneralCache}, species, config) = (species, config)

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
  args = cacheinitargs(C, species, config)
  subkey, cacheop = CacheKeyAndOp{C}(args...)
  cache = get!(()->CacheDict{C}(subkey=>cacheop(f(args...))), data, key)
  V = Base.promote_op(f, typeof.(args)...)
  return isconcretetype(V) ? cache::CacheDict{C,V} : cache
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

"""
    memoise(f, species, config::Configuration, generalcaches=Dict{UInt64,CacheDict{GeneralCache}}())

Return `f(species, config)`, memoised in `generalcaches` if
`config.options.memoisegeneral` is true.

This is for contributions that are calculated in one go, rather than summed
over cyclotron harmonics e.g. coupled velocity, relativistic momentum and
Maxwell-Juttner species. The outer dictionary is keyed on the species and the
inner `CacheDict` is keyed on the frequency, wavenumber and options.
"""
function memoise(f::F, species, config::Configuration,
    generalcaches=Dict{UInt64,CacheDict{GeneralCache}}()) where {F<:Function}
  config.options.memoisegeneral || return f(species, config)
  cache = get!(f, generalcaches, species, config) # cache for this species (except Π) + config
  return get!(f, cache, species, config) # memoised function f with the correct data dict
end

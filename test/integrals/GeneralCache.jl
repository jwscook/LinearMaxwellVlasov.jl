using Dates
println("$(now()) $(@__FILE__)")

using Test, Random
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

Random.seed!(0)

@testset "GeneralCache" begin
  memopts(; kwargs...) = Options(; memoisegeneral=true, kwargs...)
  # number of configs cached for a given species
  nentries(cache, s) = length(cache.general[LMV.typehash(LMV.GeneralCache, s)].data)

  @testset "memoise logic with a counting function" begin
    ncalls = Ref(0)
    # cheap stand-in for coupledvelocity etc. that depends on every key input
    # but, like them, not on the plasma frequency
    f(s, c) = (ncalls[] += 1;
      c.frequency * s.Ω .* [1 parallel(c.wavenumber); perpendicular(c.wavenumber) 1])

    species = CoupledVelocitySpecies(rand(4)...)
    config = Configuration(ComplexF64(1.0, 0.1), Wavenumber(kz=0.5, k⊥=0.25), memopts())
    cache = Cache()

    a = LMV.memoise(f, species, config, cache.general)
    @test a == f(species, config)
    ncalls[] = 0 # reset the counter
    for _ in 1:3
      @test LMV.memoise(f, species, config, cache.general) == a
    end
    @test ncalls[] == 0 # hits only
    @test length(cache.general) == 1

    # changing each part of the key forces a recalculation and a new answer
    for newconfig in (
        Configuration(ComplexF64(2.0, 0.1), config.wavenumber, memopts()),
        Configuration(config.frequency, Wavenumber(kz=-0.5, k⊥=0.5), memopts()),
        Configuration(config.frequency, Wavenumber(kz=0.5, k⊥=0.5), memopts()),
        Configuration(config.frequency, Wavenumber(kz=0.5, k⊥=0.75), memopts()),
        Configuration(config.frequency, config.wavenumber, memopts(cubature_rtol=1e-3)))
      ncalls[] = 0
      b = LMV.memoise(f, species, newconfig, cache.general)
      @test ncalls[] == 1
      @test b == f(species, newconfig) # adds a call directory
      @test LMV.memoise(f, species, newconfig, cache.general) == b # adds no call
      @test ncalls[] == 2 # still one
    end
    @test length(cache.general) == 1 # still one species
    @test nentries(cache, species) == 6 # six distinct configs

    # a species that differs only in plasma frequency shares the cache entry.
    # Reuse the same F because FCoupledVelocityNumerical's normalised closure
    # captures a Core.Box, so separately constructed F's never hash equal.
    samebutΠ = CoupledVelocitySpecies(2 * species.Π, species.Ω, species.F,
                                      species.minharmonics)
    ncalls[] = 0
    @test LMV.memoise(f, samebutΠ, config, cache.general) == a
    @test ncalls[] == 0
    @test length(cache.general) == 1 # just one "species"

    # a different species gets its own inner cache and doesn't see the old value
    other = CoupledVelocitySpecies(rand(4)...)
    LMV.memoise(f, other, config, cache.general)
    @test length(cache.general) == 2 # another species

    # memoisegeneral=false bypasses the cache entirely
    nomemo = Configuration(config.frequency, config.wavenumber, Options())
    ncalls[] = 0
    LMV.memoise(f, species, nomemo, cache.general)
    LMV.memoise(f, species, nomemo, cache.general)
    @test ncalls[] == 2

    empty!(cache)
    @test isempty(cache.general)
  end
end

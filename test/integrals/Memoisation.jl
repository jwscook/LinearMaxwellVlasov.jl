using Dates
println("$(now()) $(@__FILE__)")

using Test, Random, SpecialFunctions
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

import LinearMaxwellVlasov.CacheOp
import LinearMaxwellVlasov.GenericPerpendicularCacheOp
import LinearMaxwellVlasov.MaxwellianPerpendicularCacheOp

Random.seed!(0)

@testset "Testing perp memoisation logic" begin
  @testset "Maxwellian" begin
    kr, ki = rand(2)
    vth = rand()
    for n in -2:2
      a = perpendicular(vth, kr + im*ki, abs(n))
      for x in (kr+im*ki, kr-im*ki, -kr+im*ki, -kr-im*ki)
        cp = CacheOp{MaxwellianPerpendicularCacheOp}(
           n < 0,
           real(x)<0,
           imag(x)*real(x)<0)
        b = perpendicular(vth, x, n)
        if !all(b .≈ cp(a))
          @show n, sign(real(x)), sign(imag(x))
          @show round.(b ./ cp(a))
        end
        @test all(b .≈ cp(a))
        @test all(b .== cp(cp(b)))
      end
      a = perpendicular(vth, kr, abs(n))
      for x in (-kr, kr)
        cp = CacheOp{MaxwellianPerpendicularCacheOp}(
           n < 0,
           real(x)<0,
           imag(x)*real(x)<0)
        b = perpendicular(vth, x, n)
        if !all(b .≈ cp(a))
          @show n, sign(real(x)), sign(imag(x))
          @show round.(b ./ cp(a))
        end
        @test all(b .≈ cp(a))
        @test all(b .== cp(cp(b)))
      end
    end
  end
  @testset "Generic" begin
    species = RingBeamSpecies(rand(6)...)
    kr, ki = rand(2)
    for n in -2:2
      x = kr + im*ki
      config = Configuration(Wavenumber(kz=1.0, k⊥=x))
      a = perpendicular(species, config, abs(n))
      for x in (kr+im*ki, -kr+im*ki, kr-im*ki, -kr-im*ki)
        cp = CacheOp{GenericPerpendicularCacheOp}(
          n<0, real(x)<0, imag(x)*real(x)<0)
        config = Configuration(Wavenumber(kz=1.0, k⊥=x))
        b = perpendicular(species, config, n)
        @test all(b .≈ cp(a))
        @test all(b .== cp(cp(b)))
      end
      config = Configuration(Wavenumber(kz=1.0, k⊥=kr))
      a = perpendicular(species, config, abs(n))
      for x in (-kr, kr)
        config = Configuration(Wavenumber(kz=1.0, k⊥=x))
        cp = CacheOp{GenericPerpendicularCacheOp}(n<0, x < 0)
        b = perpendicular(species, config, n)
        @test all(b .== cp(a))
        @test all(b .== cp(cp(b)))
      end
    end
  end
end
@testset "Testing perp memoisation" begin
  for species ∈ (MaxwellianSpecies(rand(5)...), RingBeamSpecies(rand(6)...))
    for kr ∈ (rand(), -rand()), ki ∈ (rand(), -rand())
      for n ∈ -2:2
        K = Wavenumber(kz=1.0, k⊥=(rand(Bool) ? kr+im*ki : kr))
        config = Configuration(K)
        perpmem = LMV.generate_perpendicular_memoised(species, config)
        a = perpendicular(species, config, n)
        for _ ∈ 1:2
          b = perpmem(species, config, n)
          c = perpmem(species, config, n)
          @test all(b .≈ a)
          @test all(b .== c)
        end
      end
    end
  end
end
@testset "@inferred" begin
  for species ∈ (MaxwellianSpecies(rand(5)...), RingBeamSpecies(rand(6)...))
    n = 0
    K = Wavenumber(kz=1.0, k⊥=rand())
    config = Configuration(K)
    perpmem = LMV.generate_perpendicular_memoised(species, config)
    try
      @inferred perpendicular(species, config, n)
      @test true
    catch
      @test false
    end
    try
      @inferred perpmem(species, config, n)
      @test true
    catch
      @test false
    end
  end
end


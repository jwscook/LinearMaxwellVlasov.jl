using Dates
println("$(now()) $(@__FILE__)")

using DualNumbers, Test, Random, SpecialFunctions
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

import LinearMaxwellVlasov.CacheOp
import LinearMaxwellVlasov.GenericPerpendicularCacheOp
import LinearMaxwellVlasov.MaxwellianPerpendicularCacheOp

Random.seed!(0)

@testset "Memoisation" begin
  fs = (identity, y->Dual(y, 1), y->Dual(y, -1))
  @testset "Testing perp memoisation logic" begin
    @testset "Maxwellian" begin
      kr, ki = rand(2)
      vth = rand()
      for n in -2:2
        for f in fs
          a = perpendicular(vth, f(kr + im*ki), abs(n))
          for x in f.((kr+im*ki, kr-im*ki, -kr+im*ki, -kr-im*ki))
            cp = CacheOp{MaxwellianPerpendicularCacheOp}(x, n)
            b = perpendicular(vth, x, n)
            if !all(isapprox.(b, cp(a), rtol=sqrt(eps())))
              @show n, sign(real(x)), sign(imag(x))
              @show round.(b ./ cp(a))
            end
            @test all(isapprox.(b, cp(a), rtol=sqrt(eps())))
            @test all(isapprox.(b, cp(cp(b)), rtol=sqrt(eps())))
          end
          a = perpendicular(vth, f(kr), abs(n))
          for x in f.((-kr, kr))
            cp = CacheOp{MaxwellianPerpendicularCacheOp}(x, n)
            b = perpendicular(vth, x, n)
            if !all(isapprox.(b, cp(a), rtol=sqrt(eps())))
              @show n, sign(real(x)), sign(imag(x))
              @show round.(b ./ cp(a))
            end
            @test all(isapprox.(b, cp(a), rtol=sqrt(eps())))
            @test all(isapprox.(b, cp(cp(b)), rtol=sqrt(eps())))
          end
        end
      end
    end
    @testset "Generic" begin
      species = RingBeamSpecies(rand(6)...)
      kr, ki = rand(2) # the positive parts
      for n in -2:2
        for f in fs
          x = f(kr + im*ki)
          config = Configuration(Wavenumber(kz=1.0, k⊥=x))
          a = perpendicular(species, config, abs(n))
          for x in f.((kr+im*ki, -kr+im*ki, kr-im*ki, -kr-im*ki)) # permutations
            cp = CacheOp{GenericPerpendicularCacheOp}(x, n)
            config = Configuration(Wavenumber(kz=1.0, k⊥=x))
            b = perpendicular(species, config, n)
            @test all(isapprox.(b, cp(a), rtol=sqrt(eps())))
            @test all(b .== cp(cp(b)))
          end
          config = Configuration(Wavenumber(kz=1.0, k⊥=f(kr))) # real positive
          a = perpendicular(species, config, abs(n))
          for x in f.((-kr, kr)) # real only permutations
            config = Configuration(Wavenumber(kz=1.0, k⊥=x))
            cp = CacheOp{GenericPerpendicularCacheOp}(x, n)
            b = perpendicular(species, config, n)
            @test all(b .== cp(a))
            @test all(b .== cp(cp(b)))
          end
        end
      end
    end
  end
  @testset "Testing perp memoisation" begin
    for species ∈ (MaxwellianSpecies(rand(5)...), RingBeamSpecies(rand(6)...))
      for f in fs
        for kr ∈ (rand(), -rand()), ki ∈ (rand(), -rand())
          for B in (true, false)
            K = Wavenumber(kz=1.0, k⊥=f((B ? kr+im*ki : kr)))
            config = Configuration(K)
            for n ∈ -2:2
              a = perpendicular(species, config, n)
              perpmem = LMV.constructperpendicular_memoised(species, config)
              for _ ∈ 1:2
                b = perpmem(species, config, n)
                c = perpmem(species, config, n)
                @test all(isapprox.(b, a, rtol=sqrt(eps())))
                @test all(b .== c)
              end
            end
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
      perpmem = LMV.constructperpendicular_memoised(species, config)
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
end


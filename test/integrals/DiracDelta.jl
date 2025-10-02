using Dates
println("$(now()) $(@__FILE__)")

using InteractiveUtils, Random, Test, UnitTestDesign
using LinearMaxwellVlasov

const LMV = LinearMaxwellVlasov

Random.seed!(0)
@testset "DiracDelta vs Analytical vs Numerical in spiky limit" begin
reltol = 1.0e-5
smallness = 1e-5

myrand(T::Type{Float64}) = 2 * rand(T) - 1
function myrand(T::Type{ComplexF64})
  xy = 2 .* rand(2) .- 1
  sort!(xy, by=abs, rev=true)
  return xy[1] .+ im * xy[2]
end

@testset "Parallel" begin
  function paralleltest(ω, K, n, Ω, pow, dFdv, vd, vth)
    Bb = FBeam(vth, vd)
    Db = FParallelDiracDelta(vd)
    resonance = (real(ω) - n * Ω) / real(K.parallel) == vd
    @testset "kz=$(K.parallel * K.multipliersign), ω=$ω,pow=$(pow), dFdv=$(dFdv), resonance=$resonance" begin
    try
      resultBb = LMV.parallel(Bb, ω, K, n, Ω, pow, dFdv)
      isnan(resultBb) && return
      resultDb = @inferred LMV.parallel(Db, ω, K, n*Ω, pow, dFdv)
      if isfinite(resultBb)
        @test resultBb ≈ resultDb rtol=reltol atol=0
      end
    catch e
      @warn e
      @test false
    end
    end
  end

  ωTs = (Float64, ComplexF64)
  Ks = (0.0, -1.0, 1.0)
  ns = -1:1
  pows = Unsigned.(0:2)
  for (ωT, kz, n, pow, dFdv) ∈ all_pairs(ωTs, Ks, ns, pows, (true, false))
    vth = 2 * rand() * smallness
    vd = myrand(Float64)
    ω = myrand(ωT)
    Ω = myrand(Float64)
    K = LMV.Wavenumber(parallel=kz, perpendicular=NaN)
    pow == 2  && dFdv && continue
    paralleltest(ω, K, n, Ω, pow, dFdv, vd, vth)
  end
  Ω = rand()
  vth = 2 * smallness
  for kz in (-2.0, 2.0), ω in (1.5Ω + im, 1.5Ω - im)
    K = LMV.Wavenumber(parallel=kz, perpendicular=NaN)
    vd = real(ω) / kz
    # this is what this test is about, so error if wrong:
    p = (ω - 0 * Ω) / kz 
    @assert real(p) == vd "p, vd = $p, $vd"
    paralleltest(ω, K, 0, Ω, UInt64(0), false, vd, vth)
    paralleltest(ω, K, 0, Ω, UInt64(0), true, vd, vth)
    paralleltest(ω, K, 0, Ω, UInt64(1), false, vd, vth)
    paralleltest(ω, K, 0, Ω, UInt64(1), true, vd, vth)
    paralleltest(ω, K, 0, Ω, UInt64(2), false, vd, vth)
  end
end

@testset "Perpendicular" begin
  Ks = ComplexF64.((-1.0, 0.0, 1.0))
  ns = -1:1
  ls = -1:1
  pows = Unsigned.(0:2)
  for (K, n, l, pow, dFdv) ∈ all_pairs(Ks, ns, ls, pows, (true, false))
    vth⊥, vd⊥ = 2 .* rand(2)
    Ω = myrand(Float64)
    R⊥ = FRing(abs(vd⊥) * smallness, vd⊥)
    D⊥ = FPerpendicularDiracDelta(vd⊥)
    nl = Pair(n, l)
    try
      resultR⊥ = @inferred LMV.perpendicular(R⊥, Ω, nl, K, pow, dFdv)
      resultD⊥ = @inferred LMV.perpendicular(D⊥, Ω, nl, K, pow, dFdv)
      @test resultR⊥ ≈ resultD⊥ rtol=reltol atol=0
    catch e
      @test false
    end
  end
end


end

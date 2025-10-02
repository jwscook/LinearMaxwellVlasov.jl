using Dates
println("$(now()) $(@__FILE__)")

using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

import LinearMaxwellVlasov: parallel
function parallel(fz::LMV.FBeam, ω, kz, n::Integer, Ω::Number, p::Integer, diffbool::Bool, tol=LMV.Tolerance())
    return parallel(fz, ω, kz, n * Ω, p, diffbool, tol)
end
function parallel(fz::LMV.FBeam, ω, kz::Number, nΩ::Number, p::Integer, diffbool::Bool, tol=LMV.Tolerance())
  k = LMV.Wavenumber(kz, NaN)
  return parallel(fz, ω, k, nΩ, p, diffbool, tol)
end

function parallel(fz::LMV.FBeam, ω, k::LMV.Wavenumber, nΩ::Number, p::Unsigned, diffbool::Bool, tol=LMV.Tolerance())
  outputtuple = LMV.parallel(fz, ω, k, nΩ)
  if !((Unsigned(p), diffbool) in LMV.PARALLEL_TUPLE_ORDER)
      @show Unsigned(p), diffbool
      @show LMV.PARALLEL_TUPLE_ORDER
  end
  @assert (Unsigned(p), diffbool) in LMV.PARALLEL_TUPLE_ORDER
  ind = findfirst(i == (Unsigned(p), diffbool) for i in LMV.PARALLEL_TUPLE_ORDER)
  return outputtuple[ind]
end

function LinearMaxwellVlasov.parallel(fz::LMV.AbstractFParallelNumerical, ω, kz, n, Ω, p::Integer, diffbool::Bool, tol=LMV.Tolerance())
  return parallel(fz, ω, kz, n * Ω, p, diffbool, tol)
end
function LinearMaxwellVlasov.parallel(fz::LMV.AbstractFParallelNumerical, ω, kz, nΩ, p::Integer, diffbool::Bool, tol=LMV.Tolerance())
  k = kz isa Number ? LMV.Wavenumber(kz, NaN) : kz
  return LMV.parallel(fz, ω, k, nΩ, p, diffbool, tol)
end

using Test
@testset "Integrals tests" begin
  include("integrals/DiracDelta.jl")
  include("integrals/Inferred.jl")
  include("integrals/Integrals.jl")
  include("integrals/Parallel.jl")
  include("integrals/ParallelBasic.jl")
  include("integrals/ParallelSignOfKpara.jl")
  include("integrals/Parallel_unittest1.jl")
  include("integrals/Parallel_unittest2.jl")
  include("integrals/Parallel_unittest3.jl")
  include("integrals/Parallel_unittest4.jl")
  include("integrals/Perpendicular.jl")
  include("integrals/PerpendicularKeyAndFactor.jl")
  include("integrals/Perpendicular_maxwellian.jl")
  include("integrals/Memoisation.jl")
end

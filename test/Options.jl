using Dates
println("$(now()) $(@__FILE__)")

using LinearMaxwellVlasov
using Test


@testset "Options" begin
  @testset "memoiseparallel" begin
    o = Options(memoiseparallel=false)
    @test o.memoiseparallel == false
  end

  @testset "memoiseperpendicular" begin
    o = Options(memoiseperpendicular=false)
    @test o.memoiseperpendicular == false
  end

  @testset "tols" begin
    tol = Tolerance(rtol=0.0, atol=1.0)
    o = Options(tols=tol)
    @test tol == o.quadrature_tol
    @test tol == o.summation_tol
  end

  @testset "cubature_maxevals" begin
    o = Options(cuba_evals=101)
    @test o.cubature_maxevals == 101
    o = Options(maxevals=102)
    @test o.cubature_maxevals == 102
    o = Options(cubature_maxevals=103)
    @test o.cubature_maxevals == 103
  end

end

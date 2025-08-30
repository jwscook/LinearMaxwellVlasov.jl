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
    @test tol == o.cubature_tol
    @test tol == o.summation_tol
  end

  @testset "cubature_maxevals" begin
    o = Options(cubature_maxevals=101)
    @test o.cubature_maxevals == 101
    o = Options(cuba_evals=102)
    @test o.cubature_maxevals == 102
  end

  @testset "residue_maxevals" begin
    o = Options(residue_maxevals=101)
    @test o.residue_maxevals == 101
    o = Options(res_evals=102)
    @test o.residue_maxevals == 102
  end

end

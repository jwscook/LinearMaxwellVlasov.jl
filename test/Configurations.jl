using Dates
println("$(now()) $(@__FILE__)")

using LinearMaxwellVlasov, Test


@testset "hash is unique when changing wavenumber" begin
  C = Configuration(Wavenumber(1.0, 0.1))
  h1 = hash(C)
  C.wavenumber = Wavenumber(2.0, 0.2)
  h2 = hash(C)
  @test h1 != h2
end

@testset "hash is unique when changing frequencies" begin
  C = Configuration(ComplexF64(1.0, 0.1))
  h1 = hash(C)
  C.frequency = ComplexF64(2.0, 0.2)
  h2 = hash(C)
  @test h1 != h2
end

@testset "hash is unique when changing options" begin
  C = Configuration(ComplexF64(1.0, 0.1))
  h1 = hash(C)
  O1 = deepcopy(C.options)
  O2 = Options(memoiseparallel=false)
  @assert O1 != O2
  C.options = O2
  h2 = hash(C)
  @test h1 != h2
end

@testset "hash takes second arg" begin
  C = Configuration(Wavenumber(1.0, 0.1), Options())
  h1 = hash(C)
  @test h1 != hash(C, h1)
end

@testset "different ctor" begin
  try
    C = Configuration(ComplexF64(1.0, 0.1), Options())
    @test true
  catch
    @test false
  end
end

@testset "deepcopy ==" begin
  C = Configuration(Wavenumber(1.0, 0.1))
  @test deepcopy(C) == C
end


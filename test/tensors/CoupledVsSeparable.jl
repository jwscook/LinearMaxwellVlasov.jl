using Dates
println("$(now()) $(@__FILE__)")

using Test, Random
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

Random.seed!(0)

@testset "Separable vs Coupled velocity tensors" begin
  mₑ = LMV.mₑ
  mi = 1836*mₑ
#  realratios  = Float64[]
#  imagratios  = Float64[]
  for m ∈ (60mₑ,)#, _ ∈ 1:2
    B0 = rand() * 4
    n0 = rand() * 1e20
    Ω = cyclotronfrequency(B0, m, 1)
    Π = plasmafrequency(n0, m, 1)
    Va = sqrt(B0^2 / LMV.μ₀ / n0 / mi)
    ϵV = rand() * 1.0e4
    vth = thermalspeed(ϵV, m)
    λD = vth / Π

    argsM = 2 * rand(3) * vth
    coupledMaxwellian = CoupledVelocitySpecies(Π, Ω, argsM...)
    separableMaxwellian = MaxwellianSpecies(Π, Ω, argsM...)
    argsR = 2 * rand(4) * vth
    coupledRingBeam = CoupledVelocitySpecies(Π, Ω, argsR...)
    separableRingBeam = RingBeamSpecies(Π, Ω, argsR...)

    for (coupled, separable) ∈ (
                                (coupledMaxwellian, separableMaxwellian),
                                #(coupledRingBeam, separableRingBeam),
                               )
      k = abs(Ω / Va / 2)
      ωr = real(abs(vth * abs(k))) # real ωr must be > 0
      rtol=1e-4
      atol=eps()
      σs = (0, -1, 1)#-0.1, -0.01, -1e-3,) #1e-3, 1e-2, 1e-1, 1.0)#, 0
      kzs = (2k, k, k/2, 0, -k/2, -k, -2k)#-k, -k/10, -k/100)
      k⊥s = (k, k/2, k, 2k)
      for σ ∈ σs, kz in kzs, k⊥ in k⊥s
        F = ComplexF64(ωr, σ * ωr / 100)
        K = Wavenumber(kz=kz, k⊥=k⊥)
        iszero(K) && continue
        config = Configuration(F, K)
        config.options = Options(quadrature_rtol=1.0e-9, summation_rtol=1e-8)
        #outputS = LMV.contribution(separable, config)
        #outputC = LMV.contribution(coupled, config)

        #@test separable(0.0, 0.0) ≈ coupled(0.0, 0.0)

        #@testset "real, $σ, $kz" begin
        #  @test real(outputC[1,1])≈real(outputS[1,1]) rtol=rtol atol=atol
        #  @test real(outputC[1,2])≈real(outputS[1,2]) rtol=rtol atol=atol
        #  @test real(outputC[1,3])≈real(outputS[1,3]) rtol=rtol atol=atol
        #  @test real(outputC[2,1])≈real(outputS[2,1]) rtol=rtol atol=atol
        #  @test real(outputC[2,2])≈real(outputS[2,2]) rtol=rtol atol=atol
        #  @test real(outputC[2,3])≈real(outputS[2,3]) rtol=rtol atol=atol
        #  @test real(outputC[3,1])≈real(outputS[3,1]) rtol=rtol atol=atol
        #  @test real(outputC[3,2])≈real(outputS[3,2]) rtol=rtol atol=atol
        #  @test real(outputC[3,3])≈real(outputS[3,3]) rtol=rtol atol=atol
        #end
        #@testset "imag, $σ, $kz" begin
        #  @test imag(outputC[1,1])≈imag(outputS[1,1]) rtol=rtol atol=atol
        #  @test imag(outputC[1,2])≈imag(outputS[1,2]) rtol=rtol atol=atol
        #  @test imag(outputC[1,3])≈imag(outputS[1,3]) rtol=rtol atol=atol
        #  @test imag(outputC[2,1])≈imag(outputS[2,1]) rtol=rtol atol=atol
        #  @test imag(outputC[2,2])≈imag(outputS[2,2]) rtol=rtol atol=atol
        #  @test imag(outputC[2,3])≈imag(outputS[2,3]) rtol=rtol atol=atol
        #  @test imag(outputC[3,1])≈imag(outputS[3,1]) rtol=rtol atol=atol
        #  @test imag(outputC[3,2])≈imag(outputS[3,2]) rtol=rtol atol=atol
        #  @test imag(outputC[3,3])≈imag(outputS[3,3]) rtol=rtol atol=atol
        #end
        #(iszero(imag(F)) && iszero(kz)) && continue

        config.options = Options(quadrature_rtol=1.0e-15, summation_rtol=4eps())

        t0 = @elapsed summed = LMV.contribution(separable, config)
        #try
          config.options = Options(quadrature_rtol=1.0e-5, quadrature_atol=1e-64,
                                   summation_rtol=1e-12)
          t1 = @elapsed newberger = LMV.coupledvelocity(coupled, config)
          if !isapprox(summed, newberger, rtol=1e-4)
            ratio = newberger ./ summed
            for i in eachindex(ratio)
              @show ratio[i], abs(ratio[i])
            end
            @show σ, k, t1 / t0
          end
          @testset "newberger, $σ, $kz, $k⊥" begin
            @test isapprox(summed, newberger, rtol=1e-4)
          end
          #t0 = @elapsed LMV.converge(n->LMV.coupledvelocity(coupled, config, n))
          #@show t1 / t0
        #catch err
        #  @show err
        #  @testset "newberger, $σ, $kz, fail" begin
        #    @test false
        #  end
        #end
      end
    end
  end
end

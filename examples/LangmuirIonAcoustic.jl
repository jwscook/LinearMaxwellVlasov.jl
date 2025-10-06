using Dates
println("Starting Langmuir Ion Acoustic waves at ", now())
# addprocs(4-nprocs())
# @everywhere begin
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

using Plots, WindingNelderMead, LinearAlgebra
function run()
  mₑ = LinearMaxwellVlasov.mₑ
  q₀ = LinearMaxwellVlasov.q₀
  mi = 183.6*mₑ # reduced mass ratio
  n0 = 1.0e19
  B0 = 1.0
  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωi = cyclotronfrequency(B0, mi, 1)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πi = plasmafrequency(n0, mi, 1)
  ϵV = 1.0e3
  vthe = thermalspeed(ϵV, mₑ)
  vthi = thermalspeed(ϵV, mi)
  λD = vthe / Πe
 
  electron = MaxwellianSpecies(Πe, Πe*1.0e-3, vthe, vthe)
  proton = MaxwellianSpecies(Πi, Πi*1.0e-3, vthi, vthi)
 
  S = Plasma([electron, proton])
  N = 128
  ks = Vector{Float64}(range(0.01, stop=5, length=N))
  O = Options(memoiseperpendicular=true)
 
  solutions1 = Vector{Any}()
  solutions2 = Vector{Any}()
  cache = LinearMaxwellVlasov.Cache()
  function f!(x::Vector{Float64}, C::Configuration)
    C.frequency = ComplexF64(x[1], x[2])
    output = electrostatic(S, C, cache)
#    @show output
    return output
  end
  for i = 1:N
    K = Wavenumber(ks[i] * 2π/λD, 0.0)
    @assert K.perpendicular == 0
    C = Configuration(K, O)
    for (icnorm, sols) in ((1.25Πe, solutions1), (Πi, solutions2))
      for _ in 1:10
        ic = [2rand(), randn()/10] .* icnorm
        neldermeadsol = WindingNelderMead.optimise(x->norm(f!(x, C)),
          ic, [1, 1] .* icnorm / 100; stopval=0.1, timelimit=1000,
          maxiters=1000, ftol_rel=0, ftol_abs=0, xtol_rel=0, xtol_abs=1e-15*icnorm)
        simplex, windingnumber, returncode, numiterations = neldermeadsol
        #@show windingnumber, returncode, numiterations, simplex.vertices[1].value
        if (windingnumber == 1 && returncode == :XTOL_REACHED) || returncode == :STOPVAL_REACHED
          c = deepcopy(C)
          push!(sols, c)
        end
      end
    end
  end

  ω1 = [solution.frequency for solution in solutions1]./Πe
  ω2 = [solution.frequency for solution in solutions2]./Πe
  @show length(ω1), length(ω2)

  h = Plots.scatter(ks, real(ω1))
  Plots.scatter!(h, ks, real(ω2))
  #Plots.scatter!(h, ks, imag(ω1))
  #Plots.scatter!(h, ks, imag(ω2))
  
  #Plots.plot!(h, ks, @. sqrt(Πe^2 + 3*vthe^2*(ks*2*pi/λD)^2)/Πe)
  ## Plots.plot!(h, ks, @. sqrt(Πe^2 + (2.99e8)^2*(ks*2*pi/λD)^2)/Πe)
  #Plots.plot!(h, ks, @. (ks*2*pi/λD).*sqrt(LMV.q₀*(ϵV + 3*ϵV)/mi)/Πe)
  ##vphase1 = real(ω1)./ks * Πe/(2π/λD)
  ##vphase2 = real(ω2)./ks * Πe/(2π/λD)
  ##Plots.plot!(h, ks, -pi*Πe^2./(ks*2π/λD).^2 .* electron.Fb.dFdv(vphase1) / Πe)
  ## Plots.plot!(h, ks, -pi*Πi^2./(ks*2π/λD).^2 .* proton.Fb.dFdv(vphase2) / Πe)
  
#  Plots.savefig(h, "LangmuirIonAcoustic.pdf")
  h
end
h = run()
Plots.display(h)
println("Ending at ", now())

# rmprocs(3)











#

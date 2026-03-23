using Dates
println("Starting Cold Plasma waves at ", now())
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

using NelderMead, Plots, Random, LinearAlgebra

const mₑ = LMV.mₑ
const c0 = LMV.c₀
const mi = 18.36*mₑ# reduced mass
const n0 = 1.0e19
const B0 = 1.0
const Va = B0/sqrt(n0*LMV.μ₀*mi)
const Πe = plasmafrequency(n0, mₑ, -1)
const Πi = plasmafrequency(n0, mi, 1)
const Ωe = cyclotronfrequency(B0, mₑ, -1)
const Ωi = cyclotronfrequency(B0, mi, 1)
const electron = ColdSpecies(Πe, Ωe)
const proton = ColdSpecies(Πi, Ωi)

#electron = NumericalSpecies(Πe, Ωe,
#  FParallelDiracDelta(),
#  FPerpendicularDiracDelta())
#proton = NumericalSpecies(Πi, Ωi,
#  FParallelDiracDelta(),
#  FPerpendicularDiracDelta())

const k0 = Ωi / Va
const k1 = Πe / c0
#const  ω0 = Πe
const ω0 = Ωi
const S = Plasma([electron, proton])

solutions1 = Vector{Any}()
solutions2 = Vector{Any}()
function f!(x, C::Configuration)
  C.frequency = ComplexF64(x[1], 0.0)
  return abs(det(tensor(S, C)))
end

N = 1000
ks = Vector{Float64}(range(-2.0, stop=2.0, length=N))
Random.seed!(0)
for i = 1:N
  K = Wavenumber(k=ks[i] * k1, θ=π/4)
  C = Configuration(K, Options())
  for j in 1:20
    ic = rand()*(10.0 - 0.01) + 0.01
    t1 = @elapsed neldermeadsol = NelderMead.optimise(x->norm(f!(x, C)),
        [ic] * Πe, [Πe/100]; stopval=1.0e-8, timelimit=1200, ftol_rel=10eps(),
        maxiters=200)
    minimizer, val, returncode, numiterations = neldermeadsol
    if returncode == :STOPVAL_REACHED
      push!(solutions1, deepcopy(C))
    end
  end
  for j in 1:20
    ic = rand()*(10.0 - 0.01) + 0.01
    t1 = @elapsed neldermeadsol = NelderMead.optimise(x->norm(f!(x, C)),
        [ic] * Ωi, [Ωi/100]; stopval=1.0e-8, timelimit=1200, ftol_rel=10eps(),
        maxiters=200)
    minimizer, val, returncode, numiterations = neldermeadsol
    if returncode == :STOPVAL_REACHED
      push!(solutions1, deepcopy(C))
    end
  end
end
# Profile.print(format=:flat, combine=true, sortedby=:count)
k1s = [perp(solution.wavenumber) for solution in solutions1]./k0
k2s = [perp(solution.wavenumber) for solution in solutions2]./k0
ω1s = [solution.frequency for solution in solutions1]./ω0
ω2s = [solution.frequency for solution in solutions2]./ω0

h = Plots.scatter(k1s, real(ω1s), markercolor=:green, markershape=:circle)
Plots.scatter!(k2s, real(ω2s), markercolor=:green, markershape=:square)
Plots.display(h)

println("Ending at ", now())

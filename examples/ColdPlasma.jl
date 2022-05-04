using Dates
println("Starting Cold Plasma waves at ", now())
using LinearMaxwellVlasov
const LMV = LinearMaxwellVlasov

using Optim, Plots, Random
Plots.plotly()

const mₑ = LMV.mₑ
const c0 = LMV.c0
const mi = 18.36*mₑ
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
const S = [electron, proton]
O = Options()
O.solution_tol.abs = 1.0e-3
O.solution_tol.rel = 1.0e-6

# Profile.clear()
# Profile.init(n = 10^7, delay = 0.01)
# end # @everywhere
solutions1 = Vector{Solutions.Solution}()
solutions2 = Vector{Solutions.Solution}()
function f!(x, C::Configuration)
  C.frequency = ComplexF64(x[1], 0.0)
  return abs(det(tensor(S, C)))
end

N = 1000
ks = Vector{Float64}(range(-2.0, stop=2.0, length=N))
for method ∈ (:LN_COBYLA, :LN_NELDER_MEAD)
duration = @elapsed begin
  Random.seed!(0)
  for i = 1:N
    K = Wavenumber(k=ks[i] * k1, θ=π/4)
    C = Configuration(K, O)
    for j in 1:20
      ic = rand()*(10.0 - 0.01) + 0.01
      push!(solutions1,
            Solvers.solve(f!, C, [ic]*Πe, [0.01]*Πe, [10.0]*Πe, method)...)
    end
  end
  for i = 1:N
    K = Wavenumber(k=ks[i] * k0, θ=π/4)
    C = Configuration(K, O)
    for j in 1:20
      ic = rand()*(10.0-0.01) + 0.01
      push!(solutions2,
            Solvers.solve(f!, C, [ic]*Ωi, [0.01]*Ωi, [10]*Ωi, method)...)
    end
  end
end
# Profile.print(format=:flat, combine=true, sortedby=:count)
values1 = [solution.value for solution in solutions1]
values2 = [solution.value for solution in solutions2]
mask1 = values1 .< 2
mask2 = values2 .< 2
k1s = [perp(solution.wavenumber) for solution in solutions1]./k0
k2s = [perp(solution.wavenumber) for solution in solutions2]./k0
ω1s = [solution.frequency for solution in solutions1]./ω0
ω2s = [solution.frequency for solution in solutions2]./ω0

h = Plots.scatter(k1s[mask1], real(ω1s[mask1]),# c=log10.(values1),
              markercolor=:green, markershape=:circle)
Plots.scatter!(k2s[mask2], real(ω2s[mask2]), #c=log10.(values2))
              markercolor=:green, markershape=:square)
Plots.title!("$method, $duration")
Plots.display(h)

println("Ending at ", now())
end

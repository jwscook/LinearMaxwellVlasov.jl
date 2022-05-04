using Dates
println("Starting Warm Plasma waves at ", now())
# addprocs(4-nprocs())
# @everywhere begin
using LinearMaxwellVlasov
,.. Solvers

using ,..Integrals

using Optim, Plots
Plots.plotly()

mₑ = mₑ
c0 = c0
mi = 18.36*mₑ
n0 = 1.0e19
B0 = 1.0
Va = B0/sqrt(n0*μ₀*mi)
Πe = plasmafrequency(n0, mₑ, -1)
Πi = plasmafrequency(n0, mi, 1)
Ωe = cyclotronfrequency(B0, mₑ, -1)
Ωi = cyclotronfrequency(B0, mi, 1)
ϵV = 1.0e2
vthe = thermalspeed(ϵV, mₑ)
vthi = thermalspeed(ϵV, mi)
electron = WarmSpecies(Πe, Ωe, 5/3*vthe, 0.0)
proton = WarmSpecies(Πi, Ωi, vthi, 0.0)

#electron = NumericalSpecies(Πe, Ωe,
#  FParallelDiracDelta(),
#  FPerpendicularDiracDelta())
#proton = NumericalSpecies(Πi, Ωi,
#  FParallelDiracDelta(),
#  FPerpendicularDiracDelta())

k0 = Ωi / Va
k1 = Πe / c0
# ω0 = Πe
ω0 = Ωi
S = [electron, proton]
O = Options()
#O.quadrature_tol.abs = eps()
#O.quadrature_tol.rel = sqrt(eps())
#O.summation_tol.rel = 1.0e-8
#O.solution_tol.abs = sqrt(eps())
#O.solution_tol.rel = sqrt(eps())

# Profile.clear()
# Profile.init(n = 10^7, delay = 0.01)
# end # @everywhere
solutions1 = Vector{Solutions.Solution}()
solutions2 = Vector{Solutions.Solution}()
function f!(x, C::Configuration)
  C.frequency = ComplexF64(x[1], 0.0)
  return abs(det(tensor(S, C)))
end

N = 100
ks = Vector{Float64}(range(0.01, stop=2.0, length=N))
@timev begin
  for i = 1:N
    K1 = Wavenumber(k=ks[i] * k1, θ=π/4)
    C1 = Configuration(K1, O)
    push!(solutions1, Solvers.solve(f!, C1, [0.1]*Πe, [3.0]*Πe, :NLopt)...)

    K2 = Wavenumber(k=ks[i] * k0, θ=π/4)
    C2 = Configuration(K2, O)
    push!(solutions2, Solvers.solve(f!, C2, [0.01]*Ωi, [1.1]*Ωi, :NLopt)...)
  end
end
# Profile.print(format=:flat, combine=true, sortedby=:count)
k1 = [abs(solution.wavenumber) for solution in solutions1]./k0
k2 = [abs(solution.wavenumber) for solution in solutions2]./k0
ω1 = [solution.frequency for solution in solutions1]./ω0
ω2 = [solution.frequency for solution in solutions2]./ω0
values1 = [solution.value for solution in solutions1]
values2 = [solution.value for solution in solutions2]

# Plots.pygui(true)
# Plots.show()


h = Plots.scatter(k1, real(ω1), markercolor=:red)
Plots.scatter!(k2, real(ω2), markercolor=:blue)
Plots.plot!(k1, k1*k0*c0 / ω0)
Plots.plot!(k2, k2*k0*c0 / ω0)
Plots.plot!(k1, k1*k0*Va / ω0)
Plots.plot!(k2, k2*k0*Va / ω0)
Plots.display(h)

println("Ending at ", now())

# rmprocs(3)











#

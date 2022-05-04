using Dates
println("Starting TemperatureAnisotropy at ", now())

using LinearMaxwellVlasov

using ,..Species, ..Solutions, ..Solvers

using Optim, Plots, Base.Threads, Suppressor, Random, LinearAlgebra
using StatProfilerHTML, Profile
Plots.plotly()
Random.seed!(0) # seed rand
const mₑ = mₑ
const md = 2*1836*mₑ
const n0 = 1.0e20
const B0 = 1.0
const nd = n0
const θ = 0.0 / 180.0 * π
const Va = sqrt(B0^2/μ₀/nd/md)

const Ωe = cyclotronfrequency(B0, mₑ, -1)
const Ωd = cyclotronfrequency(B0, md, 1)
const Πe = plasmafrequency(n0, mₑ, -1)
const Πd = plasmafrequency(nd, md, 1)
const ϵV = 1.0e3
const Tratio = 4
const vthe = thermalspeed(ϵV, mₑ)
const vthdb = thermalspeed(ϵV, md)
const vthd⊥ = thermalspeed(ϵV * Tratio, md)
const electron_cold = ColdSpecies(Πe, Ωe)
const electron_warm = WarmSpecies(Πe, Ωe, vthe)
const electron_maxw = MaxwellianSpecies(Πe, Ωe, vthe, vthe)
const electron_hot = NumericalSpecies(Πe, Ωe,
  FParallelNumerical(vthe),
  FPerpendicularNumerical(vthe))

const deuteron_hot = NumericalSpecies(Πd, Ωd,
  FParallelNumerical(vthdb),
  FPerpendicularNumerical(vthd⊥))
const deuteron_cold = ColdSpecies(Πd, Ωd)
const deuteron_warm = WarmSpecies(Πd, Ωd, vthdb, vthd⊥)
const deuteron_maxw = MaxwellianSpecies(Πd, Ωd, vthdb, vthd⊥)

const Smaxw = [electron_maxw, deuteron_maxw]
const Srb = [electron_maxw, deuteron_maxw]
const Shot = [electron_maxw, deuteron_maxw]
const Scold = [electron_cold, deuteron_cold]
const Swarm = [electron_warm, deuteron_warm]
const Sgyro = [electron_maxw, deuteron_maxw]

f0 = abs(Ωd)
k0 = f0 / abs(Va)

O = Options()

function icsandbounds_ω(ω0, γ0)
  γmax = abs(f0)
  ub = [max(f0, ω0*2), γmax]
  lb = [min(f0*0.5, ω0/2), -γmax]
  ic = [(ub[1] + lb[1])/2, -γmax/10]
  @assert all(lb .< ic .< ub)
  return (ic, lb, ub)
end

function solve_given_ks_ω(ks, f)
  solutions = Vector{Solutions.Solution}()
  for k in ks
    K = Wavenumber(k=k, θ=θ)
    C = Configuration(K, O)
    ω1 = abs(k * Va)
    lb, ub = [0.01 * ω1], [2.0 * ω1]
    solution = Solvers.solve(f, C, lb, ub, :NLopt)
    push!(solutions, solution...)
  end
  return solutions
end

#method = :Newton
#method = :NLopt
method = :LN_NELDERMEAD
#method = :LN_NEWUOA_BOUND
#method = :NelderMead
#method = :GN_DIRECT
function solve_given_ks(ks, f, icsandbounds)
  solutions_out = Vector{Solutions.Solution}()
  sl = Threads.SpinLock()
  #@suppress @threads for k in ks
  for (i, k) in enumerate(ks)
    K = Wavenumber(k=k, θ=θ)
    ωVaS = Solutions.slowmagnetoacousticfrequency(Va, vthdb, K)
    for ω in (ωVaS, Ωd*0.1, Ωd*0.25, Ωd*0.5, Ωd*0.75, Ωd*0.9)
      F = ComplexF64(ω, 0.0)
      C = Configuration(F, K, O)
      ω0, γ0 = reim(F.ω)
      ic, lb, ub = icsandbounds(ω0, γ0)
      t = @elapsed solution = Solvers.solve(f, C, ic, lb, ub, method)
      lock(sl)
      push!(solutions_out, solution...)
      unlock(sl)
    end
  end
  return solutions_out
end

function f2Dω!(x::Vector{Float64}, C::Configuration, S)
  C.frequency = ComplexF64(x[1], x[2])
  t = @elapsed output = det(tensor(S, C)[1:2, 1:2])
  return output
end

f2Dωmaxw!(x::Vector, C::Configuration) = f2Dω!(x, C, Smaxw)
f2Dωwarm!(x::Vector, C::Configuration) = f2Dω!(x, C, Swarm)
f2Dωrb!(x::Vector, C::Configuration) = f2Dω!(x, C, Srb)
f2Dωhot!(x::Vector, C::Configuration) = f2Dω!(x, C, Shot)

ks_positive = range(0.001, stop=2.0, length=2^10) * k0
ks = sort(vcat(-ks_positive, ks_positive))

@time solutions_maxw = solve_given_ks(ks, f2Dωmaxw!, icsandbounds_ω)
@time solutions_warm = solve_given_ks(ks, f2Dωwarm!, icsandbounds_ω)

hs = Vector{Any}()
for sols in (solutions_maxw, solutions_warm)
  ω1 = [sol.frequency for sol in sols]./f0
  kb1 = [para(sol.wavenumber) for sol in sols]./k0
  k1 = [sol.wavenumber.k for sol in sols]./k0
  Ks = [deepcopy(sol.wavenumber) for sol in sols]
  converged1 = [sol.converged for sol in sols]
  notconverged1 = [!sol.converged for sol in sols]
  values1 = [sol.value for sol in sols]

  validω = @. (abs(ω1) < 100) & (imag(ω1) .> -50) & (real(ω1) .> 0.0)
  validω = @. (values1 < 0.1) & validω
  unstable1 = @. (imag(ω1) > 0) & validω
  stable1 = @. (!unstable1) & validω

  imag_scale = 1

  h = Plots.scatter(k1[unstable1], imag(ω1[unstable1])*imag_scale,
             markercolor=:red, markershape=:utriangle,
             xticks = -100:10:100, yticks=-50:1:100)
  Plots.scatter!(k1[unstable1], real(ω1[unstable1]),
             markercolor=:red, markershape=:square)
  Plots.scatter!(k1[stable1], imag(ω1[stable1])*imag_scale,
             markercolor=:blue, markershape=:circle)
  Plots.scatter!(k1[stable1], real(ω1[stable1]),
             markercolor=:cyan, markershape=:diamond)
  push!(hs, h)
end
Plots.display.(hs)
println("Ending at ", now())











#

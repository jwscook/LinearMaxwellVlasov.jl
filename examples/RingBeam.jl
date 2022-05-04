using Dates
println("Starting RingBeam plotter at ", now())

using LinearMaxwellVlasov

using ,..Species, ..Solutions, ..Solvers

using Optim, Plots, Base.Threads, Suppressor, Random, LinearAlgebra
using StatProfilerHTML, Profile
Plots.plotly()
Random.seed!(0) # seed rand
const mₑ = mₑ
const md = 2*1836*mₑ
const mα = 2*md
const n0 = 1.0e19
const B0 = 2.1
const ξ = 1.0e-2
const nd = n0 / (1.0 + 2*ξ)
const nα = ξ*nd
@assert n0 ≈ 2*nα + nd
const θ = 86.0 / 180.0 * π
const Va = sqrt(B0^2/μ₀/nd/md)

const Ωe = cyclotronfrequency(B0, mₑ, -1)
const Ωd = cyclotronfrequency(B0, md, 1)
const Πe = plasmafrequency(n0, mₑ, -1)
const Πd = plasmafrequency(nd, md, 1)
const ϵV = 1.0e3
const vthe = thermalspeed(ϵV, mₑ)
const vthd = thermalspeed(ϵV, md)
const Ωα = cyclotronfrequency(B0, mα, 2)
const Πα = plasmafrequency(nα, mα, 2)
const vα = thermalspeed(3.52e6, mα)
const θpseudopitch = -135.0 * π/180
const vα⊥ = vα * abs(sin(θpseudopitch))
const vαb = vα * cos(θpseudopitch)
const vαth = thermalspeed(10000.0, mα)

const alpha_hot = NumericalSpecies(Πα, Ωα,
  FParallelNumerical(vαth, vαb),
  FPerpendicularNumerical(vαth, vα⊥))
const alpha_ringbeam = RingBeamSpecies(Πα, Ωα, vαth, vαth, vαb, vα⊥)
const alpha_beam = RingBeamSpecies(Πα, Ωα, vαth, vαth, vαb, 0.0)
const alpha_maxw = MaxwellianSpecies(Πα, Ωα, vαth, vαth, vαb)

const alpha_gyro = NumericalSpecies(Πα, Ωα,
  FParallelNumerical(vαth, vαb),
  FPerpendicularNumerical(vαth, vα⊥),
  FGyroangleCosine(x->1/2/pi + 1/20/pi*cos(8*x)))

k = Va / (8 * Ωα)
K = Wavenumber(k=k, θ=θ)

f0 = abs(Ωα)
k0 = f0 / abs(Va)

function fastwavephasespeeds()
  electron_cold = ColdSpecies(Πe, Ωe)
  electron_warm = WarmSpecies(Πe, Ωe, vthe)
  electron_maxw = MaxwellianSpecies(Πe, Ωe, vthe, vthe)

  deuteron_cold = ColdSpecies(Πd, Ωd)
  deuteron_warm = WarmSpecies(Πd, Ωd, vthd)
  deuteron_maxw = MaxwellianSpecies(Πd, Ωd, vthd, vthd)

  Scold = [electron_cold, deuteron_cold]
  Swarm = [electron_warm, deuteron_warm]
  Smaxw = [electron_maxw, deuteron_maxw]

  O = Options()

  function icsandbounds_ω(ω0, γ0)
    γmax = abs(Ωd) * 0.25
    ub = [ω0*2.5, γmax]
    lb = [ω0*0.9, -γmax]
    ic = [ω0, -γmax/10]
    @assert all(lb .< ic .< ub)
    return (ic, lb, ub)
  end

  function solve_given_ks_ω(ks, f)
    solutions = Vector{Solutions.Solution}()
    for k in ks
      K = Wavenumber(k, θ)
      C = Configuration(K, O)
      ωVa = abs(perp(K) * Va)
      for i in (0.25, 0.5, 1.0, 1.5)
        lb, ub = [0.1 * ωVa], [i * ωVa]
        solution = Solvers.solve(f, C, lb, ub, :NLopt)
        for s in solution
          real(s.frequency) ≈ lb[1] && continue
          real(s.frequency) ≈ ub[1] && continue
          push!(solutions, s)
        end
      end
    end
    return solutions
  end

  method = :NelderMead
  function solve_given_ks(ks, f, icsandbounds)
    solutions_out = Vector{Solutions.Solution}()
    sl = Threads.SpinLock()
    for k in ks
      K = Wavenumber(k, θ)
      ωVa = abs(perp(K) * Va)
      for m in [0.5, 1.0]
        F = ComplexF64(max(eps(), m*ωVa), 0.0)
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
    return det(tensor(S, C))
  end
  function f1Dω!(x::Vector{Float64}, C::Configuration, S)
    C.frequency = ComplexF64(x[1], imag(C.frequency))
    return det(tensor(S, C))
  end

  f2Dωmaxw!(x::Vector, C::Configuration) = f2Dω!(x, C, Smaxw)
  f1Dωcold!(x::Vector, C::Configuration) = f1Dω!(x, C, Scold)
  f1Dωwarm!(x::Vector, C::Configuration) = f1Dω!(x, C, Swarm)

  ks = range(0.5, stop=300.0, length=2^5) * k0
  ks = sort(vcat(-ks, ks))
  solutions = solve_given_ks_ω(ks, f1Dωcold!)
  solutions = solve_given_ks(ks, f2Dωmaxw!, icsandbounds_ω)
  hs = Vector{Any}()
  frequencies = []
  wavenumbers = []
  for sols in (solutions, )
    ω1 = [sol.frequency for sol in sols]./f0
    kb1 = [para(sol.wavenumber) for sol in sols]./k0
    converged1 = [sol.converged for sol in sols]
    notconverged1 = [!sol.converged for sol in sols]
    values1 = [sol.value for sol in sols]
    validω = @. (abs(ω1) < 150) & (imag(ω1) .> -50) & (real(ω1) .> 0.0)
    k = [sol.wavenumber for sol in sols][validω]
    ω = [sol.frequency for sol in sols][validω]
    push!(frequencies, ω...)
    push!(wavenumbers, k...)
  end
  return (frequencies, wavenumbers)
end

@time ωϕ, kϕ = fastwavephasespeeds()
kbϕ = real(para.(kϕ))
vbϕ = real(ωϕ ./ kbϕ)
ωranges = []
kbranges = []
colours = []
for i in 1:2, j in 1:2
  ii = j == 1 ? vbϕ .> 0 : vbϕ .< 0
  ranges = extrema(vbϕ[ii])
  indexes = findall(x -> x==ranges[i], vbϕ)
  @assert length(indexes) >= 1
  index = indexes[1]
  push!(ωranges, ωϕ[index])
  push!(kbranges, kbϕ[index])
  push!(colours, i == 1 ? :red : :green)
end

nplot = 2^10
fα0 = zeros(Float64, nplot, nplot)
vbs = range(vαb - 6*vαth, stop=vαb + 6*vαth, length=nplot)
v⊥s = range(vα⊥ - 6*vαth, stop=vα⊥ + 6*vαth, length=nplot)
@time for (i, vb) in enumerate(vbs)
  for (j, v⊥) in enumerate(v⊥s)
    fα0[i, j] = alpha_hot(vb, v⊥)
  end
end

for n in 50:60
  h = Plots.contour(v⊥s/Va, vbs/Va, fα0, nlevels=20,
                  contour_labels=true,
                  xlabel="v perpendicular [Va]",
                  ylabel="v parallel [Va]")
  xab = [minimum(v⊥s[:]), maximum(v⊥s[:])] ./ Va
  for i in eachindex(ωranges)
    vϕi = real((ωranges[i] .- n*Ωα) ./ kbranges[i])
    yab = [vϕi, vϕi] ./ Va
    Plots.plot!(xab, yab, color=colours[i],
               xlabel="v perpendicular [Va]",
               ylabel="v parallel [Va]",
               title="n = $n")
  end
  Plots.display(h)
end
println("Ending at ", now())








#

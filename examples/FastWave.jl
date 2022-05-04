
using LinearMaxwellVlasov

using ,..Species, ..Solutions, ..Solvers

using Optim, Plots, Random, LinearAlgebra

Plots.gr()
Random.seed!(0) # seed rand
function fastwave()
  mₑ = mₑ
  md = 2*1836*mₑ
  mα = 2*md
  n0 = 1.0e19
  B0 = 2.1
  nd = n0
  θ = 85.0 / 180.0 * π
  Va = sqrt(B0^2/μ₀/nd/md)

  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωd = cyclotronfrequency(B0, md, 1)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πd = plasmafrequency(nd, md, 1)
  ϵV = 1.0e3
  vthe = thermalspeed(ϵV, mₑ)
  vthd = thermalspeed(ϵV, md)

  electron_cold = ColdSpecies(Πe, Ωe)
  electron_warm = WarmSpecies(Πe, Ωe, vthe)
  electron_maxw = MaxwellianSpecies(Πe, Ωe, vthe, vthe)

  deuteron_cold = ColdSpecies(Πd, Ωd)
  deuteron_warm = WarmSpecies(Πd, Ωd, vthd)
  deuteron_maxw = MaxwellianSpecies(Πd, Ωd, vthd, vthd)

  Scold = [electron_cold, deuteron_cold]
  Swarm = [electron_warm, deuteron_warm]
  Smaxw = [electron_maxw, deuteron_maxw]

  f0 = abs(Ωd)
  k0 = f0 / abs(Va)

  O = Options(rtols=1.0e-14)

  function icsandbounds_ω(ω0)
    γmax = abs(Ωd) * 0.1
    ub = [ω0*3, γmax]
    lb = [ω0/10, -γmax]
    ic = [ω0, -γmax/20]
    @assert all(lb .< ic .< ub)
    return (ic, lb, ub)
  end

  method = :LN_NELDERMEAD
  function solve_given_ks(ks, f, icsandbounds)
    solutions_out = Vector{Solutions.Solution}()
    for k in ks
      K = Wavenumber(k=k, θ=θ)
      ωVa = fastmagnetoacousticfrequency(Va, vthd, K)
      F = ComplexF64(ωVa, 0.0)
      C = Configuration(F, K, O)
      ic, lb, ub = icsandbounds(ωVa)
      t = @elapsed solution = Solvers.solve(f, C, ic, lb, ub, method)
      ismissing(solution) && continue
      push!(solutions_out, solution)
    end
    return solutions_out
  end

  function f2Dω!(x::Vector{Float64}, C::Configuration, S)
    C.frequency = ComplexF64(x[1], x[2])
    return abs(det(tensor(S, C)))
  end
  function f1Dω!(x::Vector{Float64}, C::Configuration, S)
    C.frequency = ComplexF64(x[1], imag(C.frequency))
    return abs(det(tensor(S, C)))
  end

  f2Dωmaxw!(x::Vector, C::Configuration) = f2Dω!(x, C, Smaxw)
  f1Dωcold!(x::Vector, C::Configuration) = f1Dω!(x, C, Scold)
  f1Dωwarm!(x::Vector, C::Configuration) = f1Dω!(x, C, Swarm)

  ks = range(1e-5, stop=250.0, length=2^9) * k0
  ks = sort(vcat(-ks, ks))
  @time solutions_maxw = solve_given_ks(ks, f2Dωmaxw!, icsandbounds_ω)
  phasespeeds = Vector{Float64}()
  for sols in (solutions_maxw, )

    ω1 = [sol.frequency for sol in sols]./f0
    kb1 = [para(sol.wavenumber) for sol in sols]./k0
    k1 = [perp(sol.wavenumber) for sol in sols]./k0
    converged1 = [sol.converged for sol in sols]
    notconverged1 = [!sol.converged for sol in sols]
    values1 = [sol.value for sol in sols]

    unstable1 = @. (imag(ω1) > 0)
    stable1 = @. (!unstable1)

    push!(phasespeeds, real(ω1 ./ k1)...)

    xlabel = "\$\\mathrm{Wavenumber} \\enskip [\\Omega_{ci} / V_A]\$"
    ylabel = "\$\\mathrm{Frequency} \\enskip [\\Omega_{ci}]\$"
    h = Plots.scatter(k1[unstable1], imag(ω1[unstable1]),
      markersize=3, markerstrokewidth=0, alpha=0.8, framestyle=:box, lims=:round,
      markercolor=:red, markershape=:utriangle, xlabel=xlabel, ylabel=ylabel)
    Plots.scatter!(k1[unstable1], real(ω1[unstable1]),
      markersize=3, markerstrokewidth=0, alpha=0.8, framestyle=:box, lims=:round,
      markercolor=:red, markershape=:square)
    Plots.scatter!(k1[stable1], imag(ω1[stable1]),
      markersize=3, markerstrokewidth=0, alpha=0.8, framestyle=:box, lims=:round,
      markercolor=:blue, markershape=:circle)
    Plots.scatter!(k1[stable1], real(ω1[stable1]),
      markersize=3, markerstrokewidth=0, alpha=0.8, framestyle=:box, lims=:round,
      markercolor=:cyan, markershape=:diamond, xticks=-260:20:260,
      yticks=-5:5:150)
  end
  Plots.plot!(legend=false)
  Plots.savefig("FastWave.pdf")
end

fastwave()

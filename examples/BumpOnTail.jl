using Dates

println("Starting BumpOnTail at ", now())
using LinearMaxwellVlasov
using Random, Profile, StatProfilerHTML
using Optim, Plots, LinearAlgebra, NelderMead
using BenchmarkTools, InteractiveUtils
Plots.gr()
Random.seed!(0)

function run()
  mₑ = LinearMaxwellVlasov.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  Πe = plasmafrequency(n0, mₑ, -1)
  Ωe = 1.0e-8 * Πe
  ϵV = 1.0e2
  vthe = thermalspeed(ϵV, mₑ)
  vd = -6.0 * vthe
  λD = vthe / Πe
  kD = 2π/λD
  ratio = 1.0e-3

  fu(v)  = exp.(-v.^2/vthe^2) + ratio * exp.(-(v-vd).^2/vthe^2)
  f_parallel = FParallelNumerical(fu, -10*vthe + vd, vd + 10*vthe)

  electron_and_beam = SeparableVelocitySpecies(Πe, Ωe, f_parallel,
      FPerpendicularNumerical(vthe))
  Snum = Plasma([electron_and_beam])

  electron_max = MaxwellianSpecies(Πe*sqrt(1.0-ratio), Ωe, vthe, vthe)
  beam_max = MaxwellianSpecies(Πe*sqrt(ratio), Ωe, vthe, vthe, vd)
  Smax = Plasma([electron_max, beam_max])

  electron_rb = RingBeamSpecies(Πe*sqrt(1.0-ratio), Ωe, vthe, vthe)
  beam_rb = RingBeamSpecies(Πe*sqrt(ratio), Ωe, vthe, vthe, vd)
  Srb = Plasma([electron_rb, beam_rb])

  plot_it = false
  if plot_it
    v = Vector(range(-12.0*vthe, stop=12.0*vthe, length=10000))
    Plots.plot(v/(Πe/kD), electron_and_beam.Fb.F(v)/vthe, linestyle=:cyan)
    Plots.plot!(v/(Πe/kD), electron_and_beam.Fb.dFdv(v), linestyle=:black)
  end
  function f!(C::Configuration, x::Vector{T}, S) where {T<:Number}
    C.frequency = ComplexF64(x[1], x[2])
    return abs(tensor(S, C)[3, 3])
  end

  N = 64
  ks = Vector{Float64}(range(1.0e-3, stop=1.0, length=N)) .* 0.08
  ks = sort(vcat(-ks, ks))
  O = Options()
  F = (1.05 + im * 0.1)*Πe
  K = Wavenumber()
  C = Configuration(F, K, O)

  function solve(args)
    solutions = Vector()
    @show N, args[end]
    for i ∈ eachindex(ks)
      objective, config, lb, ub = args
      c = deepcopy(config)
      c.wavenumber = Wavenumber(wavenumber=ks[i] * kD,
        propagationangle=1.0e-8*π)
      neldermeadsol = NelderMead.optimise(x->norm(objective(c, x)),
        (ub .+ lb) / 2, (ub .- lb) ./ 1000; stopval=1.0e-4, timelimit=10)
      val, minimizer, returncode, numiterations = neldermeadsol
      if returncode == :STOPVAL_REACHED
        objective(c, minimizer)
        push!(solutions, c)
      end
    end
    return solutions
  end
  f!max(C, x) = f!(C, x, Smax)
  f!rb(C, x) = f!(C, x, Srb)
  f!num(C, x) = f!(C, x, Snum)
  args = Vector{Any}()
  push!(args, (f!max, C, [0.9, -0.1]*Πe, [1.5, 0.2]*Πe))
  push!(args, (f!rb, C, [0.9, -0.1]*Πe, [1.5, 0.2]*Πe))
  push!(args, (f!num, C, [0.9, -0.1]*Πe, [1.5, 0.2]*Πe))

  @time for (arg, title) in zip(args, ("Maxwellian", "Ring-Beam", "Quadrature"))
    @time solutions = solve(arg)
    ω = [solution.frequency for solution in solutions]
    Ks = [para(solution.wavenumber) for solution in solutions]

    ωanalytical = sqrt.(Πe^2 .+ 3*(Ks*vthe).^2)
    dFdv = f_parallel.(ωanalytical ./ Ks, true)
    γanalytical = π/2 * ωanalytical * Πe.^2 ./ Ks.^2 .* sign.(Ks) .* dFdv

    h = Plots.plot(Ks / kD, ωanalytical ./ Πe,
                  markercolor=:blue, markershape=:circle,
                  xlabel="\$Wavenumber [2\\pi / \\lambda_D]\$",
                  ylabel="\$Frequency [\\Omega_{pe}]\$")
    Plots.plot!(Ks / kD, γanalytical ./ Πe,
               markercolor=:blue, markershape=:square)
    Plots.scatter!(Ks / kD, real.(ω) ./ Πe,
                  markercolor=:green, markershape=:circle)
    Plots.scatter!(Ks / kD, imag.(ω) ./ Πe,
                  markercolor=:green, markershape=:square)
    Plots.plot!(legend=false)
    Plots.savefig("$title.pdf")
  end
end
run()
println("Ending at ", now())

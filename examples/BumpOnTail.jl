using Dates

println("Starting BumpOnTail at ", now())
using LinearMaxwellVlasov
using Random, Plots, LinearAlgebra, WindingNelderMead, NelderMead
Plots.gr()
Random.seed!(0)
using ThreadSafeDicts

const LMV = LinearMaxwellVlasov

function run(DIRECTION)
  mₑ = LinearMaxwellVlasov.mₑ
  mi = 1836*mₑ
  n0 = 1.0e19
  Πe = plasmafrequency(n0, mₑ, -1)
  Ωe = Πe #/ 100
  ϵV = 1.0e2
  vthe = thermalspeed(ϵV, mₑ)
  vd = 6.0 * vthe * DIRECTION
  λD = vthe / Πe
  kD = 2π/λD
  ratio = 1.0e-3

  fu(v) = exp.(-v.^2/vthe^2) + ratio * exp.(-(v-vd).^2/vthe^2)
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

  fring = LMV.FRing(vthe)
  coupled(vz, v⊥) = f_parallel(vz) * fring(v⊥)
  fcoupled = FCoupledVelocityNumerical(coupled, (vd, vthe); autonormalise=true)
  electron_cp = CoupledVelocitySpecies(Πe, Ωe, fcoupled)
  Scp = Plasma([electron_cp])

  plot_it = false
  if plot_it
    v = Vector(range(-12.0*vthe, stop=12.0*vthe, length=10000))
    Plots.plot(v/(Πe/kD), electron_and_beam.Fb.F(v)/vthe, linestyle=:cyan)
    Plots.plot!(v/(Πe/kD), electron_and_beam.Fb.dFdv(v), linestyle=:black)
  end
  function f!(C::Configuration, x::Union{Tuple, AbstractVector}, S)
    C.frequency = ComplexF64(x[1], x[2])
    return tensor(S, C)[3, 3]
  end

  N = 32
  ks = (collect(1/2N:1/N:1-1/2N) .* 0.9 .+ 0.1) .* 0.08
  ks = sort(vcat(-ks, ks))
  freq = (1.05 + im * 0.04)*Πe

  function solve(args)
    solutions = ThreadSafeDict()
    Base.Threads.@threads :dynamic for i ∈ eachindex(ks)
      #@show ks[i]
      objective, lb, ub, opts = args
      wavenumber = Wavenumber(wavenumber=ks[i] * kD, propagationangle=1e-4)
      c = Configuration(freq, wavenumber, opts)
      try
        t1 = @elapsed neldermeadsol = NelderMead.optimise(x->norm(objective(c, x)),
          (ub .+ lb) / 2, (ub .- lb) ./ 1000; stopval=1.0e-3, timelimit=1200, ftol_rel=10eps(),
          maxiters=200)
        minimizer, val, returncode, numiterations = neldermeadsol
        if returncode == :STOPVAL_REACHED
          solutions[ks[i]] = c
        end
        @show returncode, numiterations, val, t1
      catch err
        @warn err
        rethrow(err)
      end
      #neldermeadsol = WindingNelderMead.optimise(x->objective(c, x),
      #  (ub .+ lb) / 2, (ub .- lb) ./ 1000; stopval=1.0e-4, timelimit=10,
      #  maxiters=200, ftol_rel=0, ftol_abs=0,
      #  xtol_rel=norm(ub .- lb) / 1000, xtol_abs=norm(ub .- lb) / 1000)
      #simplex, windingnumber, returncode, numiterations = neldermeadsol
      #if (windingnumber == 1 && returncode == :XTOL_REACHED)
      #  c = deepcopy(config)
      #  minimiser = if windingnumber == 0
      #    WindingNelderMead.position(WindingNelderMead.bestvertex(simplex))
      #  else
      #    WindingNelderMead.centre(simplex)
      #  end
      #  unitobjective!(c, minimiser)
      #  push!(output, c)
      #end
    end
    return solutions
  end
  f!max(C, x) = f!(C, x, Smax)
  f!rb(C, x) = f!(C, x, Srb)
  f!num(C, x) = f!(C, x, Snum)
  f!cp(C, x) = f!(C, x, Scp)
  args = Vector{Any}()

  otheropts = Options()
  coupledopts = Options(cuba_rtol=1e-4, cuba_evals=1_000_000, error_cuba=false, quad_rtol=1e-4, sum_rtol=1e-4)
  push!(args, ((f!cp, [0.9, -0.11]*Πe, [1.5, 0.1]*Πe, coupledopts), "Coupled"))
  push!(args, ((f!max, [0.9, -0.11]*Πe, [1.5, 0.1]*Πe, otheropts), "Maxwellian"))
  push!(args, ((f!rb, [0.9, -0.11]*Πe, [1.5, 0.1]*Πe, otheropts), "Ring-Beam"))
  push!(args, ((f!num, [0.9, -0.11]*Πe, [1.5, 0.1]*Πe, otheropts), "Quadrature"))

  @time for (arg, title) in args
    @time solutions = solve(arg)
    ω = [solution.frequency for solution in values(solutions)]
    Ks = [para(solution.wavenumber) * solution.wavenumber.multipliersign for solution in values(solutions)]
    inds = sortperm(Ks)
    ω = ω[inds]
    Ks = Ks[inds]

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
    Plots.savefig("$(title)_$(DIRECTION).pdf")
  end
end
run(-1)
run(1)
println("Ending at ", now())

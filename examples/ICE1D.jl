using Dates, Distributed
println("Starting ICE at ", now())

propagationangle = try; parse(Float64, ARGS[1]); catch; 89.5; end
pitch = try; parse(Float64, ARGS[2]); catch; 0.5; end
name_extension = try; string(ARGS[1], "_", ARGS[2]); catch;
  string(propagationangle); end

#nprocsadded = Int(Sys.CPU_THREADS / 2)
nprocsadded = 0
addprocs(nprocsadded)
@everywhere using ProgressMeter
@everywhere begin
  using LinearMaxwellVlasov
  
  using ,..Species
  

  using Plots, Base.Threads, Random, LinearAlgebra
  using StatProfilerHTML, Profile, ProgressMeter
  using NelderMead
  Plots.gr()
  #Plots.plotly()
  Plots.default(show=false)
  Random.seed!(1) # seed rand

  mₑ = mₑ
  md = 2*1836*mₑ
  mα = 2*md
  n0 = 1e19
  B0 = 2.1
  ξ = 1.0e-3
  nd = n0 / (1.0 + 2*ξ)
  nα = ξ*nd
  @assert n0 ≈ 2*nα + nd
  θ = Float64(@fetchfrom 1 propagationangle) * π / 180
  Va = sqrt(B0^2/μ₀/nd/md)

  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωd = cyclotronfrequency(B0, md, 1)
  Ωα = cyclotronfrequency(B0, mα, 2)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πd = plasmafrequency(nd, md, 1)
  Πα = plasmafrequency(nα, mα, 2)
  vthe = thermalspeed(1e3, mₑ)
  vthd = thermalspeed(1e3, md)
  vα = thermalspeed(3.52e6, mα)
  pitch = Float64(@fetchfrom 1 pitch)
  vα⊥ = vα * sqrt(1 - pitch^2)
  vαb = vα * pitch
  vαth = thermalspeed(1.0e2, mα)

  electron_cold = ColdSpecies(Πe, Ωe)
  electron_warm = WarmSpecies(Πe, Ωe, vthe)
  electron_maxw = MaxwellianSpecies(Πe, Ωe, vthe, vthe)
  electron_hot = SeparableVelocitySpecies(Πe, Ωe,
    FParallelNumerical(vthe),
    FPerpendicularNumerical(vthe))

  deuteron_hot = SeparableVelocitySpecies(Πd, Ωd,
    FParallelNumerical(vthd),
    FPerpendicularNumerical(vthd))
  deuteron_cold = ColdSpecies(Πd, Ωd)
  deuteron_warm = WarmSpecies(Πd, Ωd, vthd)
  deuteron_maxw = MaxwellianSpecies(Πd, Ωd, vthd, vthd)

  alpha_cold = ColdSpecies(Πα, Ωα)
  alpha_hot = SeparableVelocitySpecies(Πα, Ωα,
    FParallelNumerical(vαth, vαb),
    FPerpendicularNumerical(vαth, vα⊥))
  alpha_ringbeam = SeparableVelocitySpecies(Πα, Ωα,
    FBeam(vαth, vαb),
    FPerpendicularNumerical(vαth, vα⊥))
  alpha_maxw = MaxwellianSpecies(Πα, Ωα, vαth, vαth, vαb)
  alpha_delta = SeparableVelocitySpecies(Πα, Ωα,
    FParallelDiracDelta(vαb),
    FPerpendicularDiracDelta(vα⊥))
  alpha_beamdelta = SeparableVelocitySpecies(Πα, Ωα,
    FBeam(vαth, vαb),
    FPerpendicularDiracDelta(vα⊥))
  alpha_shell = CoupledVelocitySpecies(Πα, Ωα,
    FShell(vαth, vα))

  Smaxw = [electron_maxw, deuteron_maxw, alpha_maxw]
  Srb = [electron_maxw, deuteron_maxw, alpha_ringbeam]
  Smmh = [electron_maxw, deuteron_maxw, alpha_hot]
  Shhh = [electron_hot, deuteron_hot, alpha_hot]
  Sbd = [electron_maxw, deuteron_maxw, alpha_beamdelta]
  Scold = [electron_cold, deuteron_warm]
  Smcr = [electron_maxw, deuteron_cold, alpha_ringbeam]
  Smcs = [electron_maxw, deuteron_cold, alpha_shell]

  f0 = abs(Ωα)
  k0 = f0 / abs(Va)

  O = Options(rtols=1.0e-3)
  cache = Cache()
  method = :LN_NELDERMEAD # :NelderMead #
  function solve_given_ks(ks, f_input)
    function icsandbounds(ω0)
      γmax = abs(Ωα) / 2
      ub = [max(ω0 * 1.1, ω0 + abs(Ωd)), γmax]
      lb = [max(eps(), min(ω0 * 0.6, ω0 - abs(Ωd))), -γmax]
      ic = [(ub[1] + lb[1]) / 2, abs(Ωα) / 4]
      @assert all(lb .< ic .< ub)
      return (ic, lb, ub)
    end
    function insert_cache(cache)
      g(a, b) = f_input(a, b, cache)
      return g, cache
    end
    function boundify(f::T, lb, ub) where T
      function bounded(x)
        output = f(x)
        U = eltype(output)
        return output .+ (all(lb .<= x .<= ub) ? zero(U) : U(Inf))
      end
      return bounded
    end

    solutions = Vector{Any}()
    for k in ks
      cache = Cache() # cache local to loop iteration
      objective, cache = insert_cache(cache) # generate cached objective
      K = Wavenumber(wavenumber=k, propagationangle=θ)
      ω0 = fastmagnetoacousticfrequency(Va, vthd, K)
      config = Configuration(K, O)
      ic, lb, ub = icsandbounds(ω0)
      objectivecold(x, C) = f1Dω!(x, C, Scold)

      neldermeadsol = NelderMead.optimise(
        boundify(x->objectivecold(x, deepcopy(config)),
                 [ω0 * 0.8], [min(ω0 * 2.0, f0 * 150)]),
        [ω0 * 0.9], [ω0 * 1e-4]; stopval=0.5, timelimit=30,
        ftol_rel=0, ftol_abs=0, xtol_rel=eps(), xtol_abs=0)
      val, minimiser, returncode, numiterations = neldermeadsol
      if returncode == :STOPVAL_REACHED
        c = deepcopy(config)
        objectivecold(c, minimiser)
        guessω = real(config.frequency)
        (lb[1] < guessω < ub[1]) || continue
        ic[1] = guessω
      end

      neldermeadsol = NelderMead.optimise(
        boundify(objective, lb, ub), ic, [ω0 * 1e-4, ω0 * 1e-8];
        stopval=1e-4, timelimit=30,
        ftol_rel=0, ftol_abs=0, xtol_rel=eps(), xtol_abs=0)
      val, minimiser, returncode, numiterations = neldermeadsol
      if returncode == :STOPVAL_REACHED
        c = deepcopy(config)
        objectivecold(c, minimiser)
        push!(solutions, c)
      end
    end
    return solutions
  end
  function f1Dω!(x::AbstractVector{T}, config::Configuration,
                 species) where {T<:Number}
    config.frequency = ComplexF64(x[1], 0.0)
    return abs(det(tensor(Plasma(species), config, cache)))
  end
  function f2Dω!(x::Vector{T}, C::Configuration, S, cache
                 ) where {T<:Number}
    C.frequency = ComplexF64(x[1], x[2])
    t = @elapsed output = abs(det(tensor(S, C, cache)))
    return output
  end

  f2Dωmaxw!(x::Vector, C::Configuration, c) = f2Dω!(x, C, Smaxw, c)
  f2Dωrb!(x::Vector, C::Configuration, c) = f2Dω!(x, C, Srb, c)
  f2Dωmmh!(x::Vector, C::Configuration, c) = f2Dω!(x, C, Smmh, c)
  f2Dωmch!(x::Vector, C::Configuration, c) = f2Dω!(x, C, Smch, c)
  f2Dωmcr!(x::Vector, C::Configuration, c) = f2Dω!(x, C, Smcr, c)
  f2Dωhhh!(x::Vector, C::Configuration, c) = f2Dω!(x, C, Shhh, c)
  f2Dωbd!(x::Vector, C::Configuration, c) = f2Dω!(x, C, Sbd, c)

  #ks_positive = range(0.5, stop=20.0, length=2^8) * k0
  #ks = sort(vcat(-ks_positive, ks_positive)) # ks_positive #
#  ks = range(0.5, stop=41.0, length=2^8) * k0
  ks = range(0.5, stop=22.0, length=2^8) * k0
  #@time solutions_maxw = solve_given_ks(ks, f2Dωmaxw!)
  distributedmmh(x) = solve_given_ks(x, f2Dωmmh!)
  distributedmch(x) = solve_given_ks(x, f2Dωmch!)
  distributedhhh(x) = solve_given_ks(x, f2Dωhhh!)
  distributedrb(x) = solve_given_ks(x, f2Dωrb!)
  distributedbd(x) = solve_given_ks(x, f2Dωbd!)
  distributedmcr(x) = solve_given_ks(x, f2Dωmcr!)
  function findsolutions(objective)
    return vcat((@showprogress pmap(objective, ks, batch_size=1))...)
  end
end #@everywhere

function plotsolutions(sols, title)
  @show length(sols)
  ω1 = [sol.frequency for sol in sols]./f0
  kb1 = [para(sol.wavenumber) for sol in sols]./k0
  k1 = [abs(sol.wavenumber) for sol in sols]./k0

  validω = @. (abs(ω1) < 100) & (imag(ω1) .> -50) & (real(ω1) .> 0.0)
  #validω = @. validω & (abs.(real(ω1) ./ k1) .< Va * k0 / f0)
  #validω = @. validω & (abs.(imag(ω1)) .< γmax / f0 * 0.99999)

  unstable1 = @. (imag(ω1) > 1.0e-3 * abs(real(ω1))) & validω
  stable1 = @. (!unstable1) & validω

  imag_scale = 10

  h = Plots.scatter(k1[unstable1], imag(ω1[unstable1])*imag_scale,
             markerstrokecolor=:red,
             markercolor=:red, markershape=:utriangle, markersize=2,
             xticks = -200:10:200, yticks=-50:1:100,
             xlabel="\$ \\mathrm{Wavenumber} \\enskip [\\Omega_{ci} / V_A]\$",
             ylabel="\$ \\mathrm{Frequency} \\enskip [\\Omega_{ci}]\$")
  Plots.scatter!(k1[unstable1], real(ω1[unstable1]), markersize=2,
             markerstrokecolor=:red,
             markercolor=:red, markershape=:square)
  Plots.scatter!(k1[stable1], imag(ω1[stable1])*imag_scale, markersize=2,
             markerstrokecolor=:blue,
             markercolor=:blue, markershape=:circle)
  Plots.scatter!(k1[stable1], real(ω1[stable1]), markersize=2,
             markerstrokecolor=:cyan,
             markercolor=:cyan, markershape=:diamond)
  Plots.annotate!(1, 1, text("\$ \\gamma \times $imag_scale)\$", 7))
  Plots.plot!(legend=false)
  Plots.savefig(title)
end

@time solutions = vcat(map(distributedmcr, ks)...)
#@time solutions = vcat(pmap(distributedhot, ks, batch_size=8)...)
#@time solutions = findsolutions(distributedhhh)
#@time solutions = findsolutions(distributedmmh)
@time solutions = findsolutions(distributedmcr)
rmprocs(nprocsadded)
plotsolutions(solutions, "ICE_$name_extension.pdf")
println("Ending at ", now())

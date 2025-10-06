using Dates
using Plots, Random, ImageFiltering, Statistics
using Dierckx, JLD2

println("Starting at ", now())
# thermal width of ring as a fraction of its speed # Dendy PRL 1993
const filecontents = [i for i in readlines(open(@__FILE__))]
const name_extension = "TwoStream"

using ProgressMeter # for some reason must be up here on its own
using StaticArrays

using LinearMaxwellVlasov, LinearAlgebra, WindingNelderMead

mₑ = LinearMaxwellVlasov.mₑ
n0 = 1e19
Πe = plasmafrequency(n0, mₑ, -1)
vbeam = thermalspeed(1e3, mₑ)

left = SeparableVelocitySpecies(Πe, eps(),
  FParallelDiracDelta(-vbeam), FPerpendicularDiracDelta(eps()))
right = SeparableVelocitySpecies(Πe, eps(),
  FParallelDiracDelta(vbeam), FPerpendicularDiracDelta(eps()))
Sdelta = [left, right]

f0 = abs(Πe) * sqrt(2)
k0 = f0 / abs(vbeam)

options = Options(memoiseparallel=false, memoiseperpendicular=true)

function expectedgrowthrate(k)
  x = k * vbeam / Πe
  return imag(sqrt(Complex(x^2+1-sqrt(4*x^2+1)))) * Πe
end
function solve_given_ks(K, objective!)
  lb = @SArray [f0 * -0.5, f0 * -0.1]
  ub = @SArray [f0 * 3.0, f0 * 0.6]

  function boundify(f::T, lb, ub) where {T}
    isinbounds(x) = all(i->lb[i] <= x[i] <= ub[i], eachindex(x))
    maybeaddinf(x::U, addinf::Bool) where {U} = addinf ? x + U(Inf) : x
    bounded(x) = maybeaddinf(f(x), !isinbounds(x))
    return bounded
  end

  config = Configuration(K, options)

  function unitobjective!(c, x::T) where {T}
    return objective!(c, T([x[i] * (ub[i] - lb[i]) + lb[i] for i in eachindex(x)]))
  end
  unitobjectivex! = x -> unitobjective!(config, x)
  boundedunitobjective! = boundify(unitobjectivex!,
                                   zeros(2), ones(2))
  xtol_abs = f0 .* (@SArray [1e-3, 1e-3]) ./ (ub .- lb)
  γ = expectedgrowthrate(parallel(K))
  output = []
  # we know the answer but don't give it as starting point, because
  # we want it to be a meaningful test
  @elapsed for ic ∈ ([-0.1*f0, γ * 0.9], [f0, 0.0])
    @assert all(i->lb[i] <= ic[i] <= ub[i], eachindex(ic))
    neldermeadsol = WindingNelderMead.optimise(
      boundedunitobjective!, ((ic .- lb) ./ (ub .- lb)),
      1e-2 * (@SArray ones(2)); stopval=1e-15, timelimit=30,
      maxiters=200, ftol_rel=0, ftol_abs=0, xtol_rel=0, xtol_abs=xtol_abs)
    simplex, windingnumber, returncode, numiterations = neldermeadsol
    if (windingnumber == 1 && returncode == :XTOL_REACHED)
      c = deepcopy(config)
      minimiser = if windingnumber == 0
        WindingNelderMead.position(WindingNelderMead.bestvertex(simplex))
      else
        WindingNelderMead.centre(simplex)
      end
      unitobjective!(c, minimiser)
      push!(output, c)
    end
  end
  return output
end

function f2Dω!(config::Configuration, x::AbstractArray,
               species, cache=Cache())
  config.frequency = Complex(x[1], x[2])
  output = tensor(Plasma(species), config, cache)[3,3]
  return output
end

function findsolutions(species)
  ngridpoints = 2^9
  kzs = collect(range(-2, stop=2, length=ngridpoints)) * k0
  push!(kzs, Πe / vbeam)
  push!(kzs, -Πe / vbeam)
  sort!(kzs)
  # change order for better distributed scheduling
  cache = Cache()
  objective! = (C, x) -> f2Dω!(C, x, species, cache)
  solutions = Vector()
  for (ikz, kz) ∈ enumerate(kzs)
    K = Wavenumber(parallel=kz, perpendicular=0.0)
    output = solve_given_ks(K, objective!)
    isnothing(output) && continue
    push!(solutions, output...)
  end
  return solutions
end

Plots.gr()
function plotit(sols, file_extension=name_extension, fontsize=9)
  sols = sort(sols, by=s->imag(s.frequency))
  ωs = [sol.frequency for sol in sols]./f0
  kzs = [para(sol.wavenumber) for sol in sols]./k0
  p = sortperm(imag.(ωs))
  ωs = ωs[p]
  kzs = kzs[p]

  xlabel = "\$\\mathrm{Wavenumber} \\quad [\\Pi_{0} / v_b]\$"
  ylabel = "\$\\mathrm{Frequency} \\quad [\\Pi_{0}]\$"
  Plots.scatter(kzs, real.(ωs), zcolor=imag.(ωs),
    markersize=2, markerstrokewidth=0, alpha=0.8, markershape=:square,
    xlabel=xlabel, ylabel=ylabel, label="\$\\omega\$")
  Plots.scatter!(kzs, imag.(ωs), zcolor=real.(ωs),
    framestyle=:box, lims=:round, markersize=3,
    markerstrokewidth=0, markershape=:circle, c=Plots.cgrad(),
    xlabel=xlabel, ylabel=ylabel, legend=:topleft, label="\$\\gamma\$")
  Plots.scatter!([sqrt(3/8)], [1/sqrt(8)], markershape=:star, markersize=5,
                 label=nothing, color=0)
  Plots.annotate!([(sqrt(3/8), 0.5,
    text("\$k=\\sqrt{3/8}\\Pi_0/v_b\$", fontsize, :black))])
  Plots.annotate!([(sqrt(3/8), 0.65,
    text("\$\\gamma=1/\\sqrt{8}\\Pi_0\$", fontsize, :black))])
  Plots.plot!(framestyle=:box, lims=:round)
  Plots.savefig("TwoStream_DispersionRelation.pdf")
end

@time plasmasols = findsolutions(Sdelta)
@time plotit(plasmasols)
@save "TwoStream.jld" filecontents plasmasols f0 k0

println("Ending at ", now())



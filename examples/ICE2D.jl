using Distributed, Dates
using Plots, Random, ImageFiltering, Statistics
using Dierckx, Contour, JLD2

println("Starting at ", now())
# pitchanglecosine of 0 is all vperp, and (-)1 is (anti-)parallel to B field
const pitchanglecosine = try; parse(Float64, ARGS[1]); catch; -0.64; end
@assert -1 <= pitchanglecosine <= 1
# thermal width of ring as a fraction of its speed # Dendy PRL 1993
const vthermalfractionz = try; parse(Float64, ARGS[2]); catch; 0.01; end
const vthermalfraction⊥ = try; parse(Float64, ARGS[3]); catch; 0.01; end
const name_extension = if length(ARGS) >= 4
  ARGS[4]
else
  "$(pitchanglecosine)_$(vthermalfractionz)_$(vthermalfraction⊥)"
end
const filecontents = [i for i in readlines(open(@__FILE__))]
const nprocsadded = div(Sys.CPU_THREADS, 2)

addprocs(nprocsadded)

@everywhere using ProgressMeter # for some reason must be up here on its own
@everywhere using StaticArrays
@everywhere using FastClosures
@everywhere begin
  using LinearMaxwellVlasov, LinearAlgebra, WindingNelderMead

  mₑ = LinearMaxwellVlasov.mₑ
  md = 2*1836*mₑ
  mα = 2*md

  # Fig 18 Cottrell 1993
  n0 = 1.7e19 # central electron density 3.6e19
  B0 = 2.07 # OR 2.191 for B field at 4m when 2.8 T on axis R0 3.13m
  ξ = 1.5e-4 # nα / ni = 1.5 x 10^-4

  nd = n0 / (1.0 + 2*ξ)
  nα = ξ*nd
  @assert n0 ≈ 2*nα + nd
  Va = sqrt(B0^2/LinearMaxwellVlasov.μ₀/nd/md)

  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωd = cyclotronfrequency(B0, md, 1)
  Ωα = cyclotronfrequency(B0, mα, 2)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πd = plasmafrequency(nd, md, 1)
  Πα = plasmafrequency(nα, mα, 2)
  vthe = thermalspeed(1e3, mₑ)
  vthd = thermalspeed(1e3, md)
  vα = thermalspeed(3.6e6, mα)
  # pitchanglecosine = cos(pitchangle)
  # acos(pitchanglecosine) = pitchangle
  pitchanglecosine = Float64(@fetchfrom 1 pitchanglecosine)
  vα⊥ = vα * sqrt(1 - pitchanglecosine^2) # perp speed
  vαz = vα * pitchanglecosine # parallel speed
  vthermalfractionz = Float64(@fetchfrom 1 vthermalfractionz)
  vthermalfraction⊥ = Float64(@fetchfrom 1 vthermalfraction⊥)
  vαthz = vα * vthermalfractionz
  vαth⊥ = vα * vthermalfraction⊥

  electron_cold = ColdSpecies(Πe, Ωe)
  electron_warm = WarmSpecies(Πe, Ωe, vthe)
  electron_maxw = MaxwellianSpecies(Πe, Ωe, vthe, vthe)

  deuteron_cold = ColdSpecies(Πd, Ωd)
  deuteron_warm = WarmSpecies(Πd, Ωd, vthd)
  deuteron_maxw = MaxwellianSpecies(Πd, Ωd, vthd, vthd)

  alpha_cold = ColdSpecies(Πα, Ωα)
  alpha_maxw = MaxwellianSpecies(Πα, Ωα, vαthz, vαth⊥, vαz)
  alpha_ringbeam = SeparableVelocitySpecies(Πα, Ωα,
    FBeam(vαthz, vαz),
    FRing(vαth⊥, vα⊥))
  alpha_delta = SeparableVelocitySpecies(Πα, Ωα,
    FParallelDiracDelta(vαz),
    FPerpendicularDiracDelta(vα⊥))
  alpha_beamdelta = SeparableVelocitySpecies(Πα, Ωα,
    FBeam(vαthz, vαz),
    FPerpendicularDiracDelta(vα⊥))
  alpha_deltaring = SeparableVelocitySpecies(Πα, Ωα,
    FParallelDiracDelta(vαz),
    FRing(vαth⊥, vα⊥))
  alpha_shell = CoupledVelocitySpecies(Πα, Ωα,
    FShell(sqrt(vαthz^2 + vαth⊥^2), vα))

  Smmr = Plasma([electron_maxw, deuteron_maxw, alpha_ringbeam])
  Smmd = Plasma([electron_maxw, deuteron_maxw, alpha_delta])
  #Smmb = Plasma([electron_maxw, deuteron_maxw, alpha_beamdelta])

  f0 = abs(Ωα)
  k0 = f0 / abs(Va)

  γmax = abs(Ωα) * 0.15
  γmin = -abs(Ωα) * 0.075
  function bounds(ω0)
    lb = @SArray [ω0 * 0.5, γmin]
    ub = @SArray [ω0 * 1.2, γmax]
    return (lb, ub)
  end

  options = Options(memoiseparallel=false, memoiseperpendicular=true)

  function solve_given_ks(K, objective!)
    ω0 = fastzerobetamagnetoacousticfrequency(Va, K, Ωd)

    lb, ub = bounds(ω0)

    function boundify(f::T) where {T}
      @inbounds isinbounds(x) = all(i->0 <= x[i] <= 1, eachindex(x))
      maybeaddinf(x::U, addinf::Bool) where {U} = addinf ? x + U(Inf) : x
      bounded(x) = maybeaddinf(f(x), !isinbounds(x))
      return bounded
    end

    config = Configuration(K, options)

    ics = ((@SArray [ω0*0.8, f0*0.08]),
           (@SArray [ω0*0.9, f0*0.04]),
           (@SArray [ω0*1.0, f0*0.01]))


    function unitobjective!(c, x::T) where {T}
      return objective!(c,
        T([x[i] * (ub[i] - lb[i]) + lb[i] for i in eachindex(x)]))
    end
    unitobjectivex! = x -> unitobjective!(config, x)
    boundedunitobjective! = boundify(unitobjectivex!)
    xtol_abs = f0 .* (@SArray [1e-4, 1e-5]) ./ (ub .- lb)
    @elapsed for ic ∈ ics
      @assert all(i->lb[i] <= ic[i] <= ub[i], eachindex(ic))
      neldermeadsol = WindingNelderMead.optimise(
        boundedunitobjective!, SArray((ic .- lb) ./ (ub .- lb)),
        1.1e-2 * (@SArray ones(2)); stopval=1e-15, timelimit=30,
        maxiters=200, ftol_rel=0, ftol_abs=0, xtol_rel=0, xtol_abs=xtol_abs)
      simplex, windingnumber, returncode, numiterations = neldermeadsol
      if (windingnumber == 1 && returncode == :XTOL_REACHED)# || returncode == :STOPVAL_REACHED
        c = deepcopy(config)
        minimiser = if windingnumber == 0
          WindingNelderMead.position(WindingNelderMead.bestvertex(simplex))
        else
          WindingNelderMead.centre(simplex)
        end
        unitobjective!(c, minimiser)
        return c
      end
    end
    return nothing
  end

  function f2Dω!(config::Configuration, x::AbstractArray, plasma, cache)
    config.frequency = Complex(x[1], x[2])
    return det(tensor(plasma, config, cache))
  end

  function findsolutions(plasma)
    ngridpoints = 2^9
    kzs = range(-2.0, stop=2.0, length=ngridpoints) * k0
    k⊥s = range(0.0, stop=15.0, length=ngridpoints) * k0

    # change order for better distributed scheduling
    k⊥s = shuffle(vcat([k⊥s[i:nprocs():end] for i ∈ 1:nprocs()]...))
    @assert length(k⊥s) == ngridpoints
    solutions = @sync @showprogress @distributed (vcat) for k⊥ ∈ k⊥s
      cache = Cache()
      objective! = @closure (C, x) -> f2Dω!(C, x, plasma, cache)
      innersolutions = Vector()
      for (ikz, kz) ∈ enumerate(kzs)
        K = Wavenumber(parallel=kz, perpendicular=k⊥)
        output = solve_given_ks(K, objective!)
        isnothing(output) && continue
        push!(innersolutions, output)
      end
      innersolutions
    end
    return solutions
  end
end #@everywhere

"Select the largest growth rates if multiple solutions found for a wavenumber"
function selectlargeestgrowthrate(sols)
    imagfreq(s) = imag(s.frequency)
    d = Dict{Any, Vector{eltype(sols)}}()
    for sol in sols
      if haskey(d, sol.wavenumber)
        push!(d[sol.wavenumber], sol)
      else
        d[sol.wavenumber] = [sol]
      end
    end
    output = Vector()
    sizehint!(output, length(sols))
    for (_, ss) in d
      push!(output, ss[findmax(map(imagfreq, ss))[2]])
    end
    return output
end

function selectpropagationrange(sols, lowangle=0, highangle=180)
  function propangle(s)
    kz = para(s.wavenumber)
    k⊥ = perp(s.wavenumber)
    θ = atan.(k⊥, kz) * 180 / pi
  end
  output = Vector()
  for s in sols
    (lowangle <= propangle(s) <= highangle) || continue
    push!(output, s)
  end
  return output
end

Plots.pyplot()
function plotit(sols, file_extension=name_extension, fontsize=9)
  sols = sort(sols, by=s->imag(s.frequency))
  ωs = [sol.frequency for sol in sols]./f0
  kzs = [para(sol.wavenumber) for sol in sols]./k0
  k⊥s = [perp(sol.wavenumber) for sol in sols]./k0
  xk⊥s = sort(unique(k⊥s))
  ykzs = sort(unique(kzs))

  function make2d(z1d)
    z2d = Array{Union{Float64, Missing}, 2}(zeros(Missing, length(ykzs),
                                            length(xk⊥s)))
    for (j, k⊥) in enumerate(xk⊥s), (i, kz) in enumerate(ykzs)
      index = findlast((k⊥ .== k⊥s) .& (kz .== kzs))
      isnothing(index) || (z2d[i, j] = z1d[index])
    end
    return z2d
  end

  ks = [abs(sol.wavenumber) for sol in sols]./k0
  kθs = atan.(k⊥s, kzs)
  extremaangles = collect(extrema(kθs))

  _median(patch) = median(filter(x->!ismissing(x), patch))
  realωssmooth = make2d(real.(ωs))
  try
    realωssmooth = mapwindow(_median, realωssmooth, (5, 5))
    realωssmooth = imfilter(realωssmooth, Kernel.gaussian(3))
  catch
    @warn "Smoothing failed"
  end

  realspline = nothing
  imagspline = nothing
  try
    smoothing = length(ωs) * 1e-4
    realspline = Dierckx.Spline2D(xk⊥s, ykzs, realωssmooth'; kx=4, ky=4, s=smoothing)
    imagspline = Dierckx.Spline2D(k⊥s, kzs, imag.(ωs); kx=4, ky=4, s=smoothing)
  catch err
    @warn "Caught $err. Continuing."
  end

  function plotangles(;writeangles=true)
    for θdeg ∈ vcat(collect.((30:5:80, 81:99, 100:5:150))...)
      θ = θdeg * π / 180
      xs = sin(θ) .* [0, maximum(ks)]
      ys = cos(θ) .* [0, maximum(ks)]
      linestyle = mod(θdeg, 5) == 0 ? :solid : :dash
      Plots.plot!(xs, ys, linecolor=:grey, linewidth=0.5, linestyle=linestyle)
      writeangles || continue
      if atan(maximum(k⊥s), maximum(kzs)) < θ < atan(maximum(k⊥s), minimum(kzs))
        xi, yi = xs[end], ys[end]
        isapprox(yi, maximum(kzs), rtol=0.01, atol=0.01) && (yi += 0.075)
        isapprox(xi, maximum(k⊥s), rtol=0.01, atol=0.01) && (xi += 0.1)
        isapprox(xi, minimum(k⊥s), rtol=0.01, atol=0.01) && (xi -= 0.2)
        Plots.annotate!([(xi, yi, text("\$ $(θdeg)^{\\circ}\$", fontsize, :black))])
      end
    end
  end
  function plotcontours(spline, contourlevels, skipannotation=x->false)
    isnothing(spline) && return nothing
    x, y = sort(unique(k⊥s)), sort(unique(kzs))
    z = evalgrid(spline, x, y)
    for cl ∈ Contour.levels(Contour.contours(x, y, z, contourlevels))
      lvl = try; Int(Contour.level(cl)); catch; Contour.level(cl); end
      for line ∈ Contour.lines(cl)
          xs, ys = Contour.coordinates(line)
          θs = atan.(xs, ys)
          ps = sortperm(θs, rev=true)
          xs, ys, θs = xs[ps], ys[ps], θs[ps]
          mask = minimum(extremaangles) .< θs .< maximum(extremaangles)
          any(mask) || continue
          xs, ys = xs[mask], ys[mask]
          Plots.plot!(xs, ys, color=:grey, linewidth=0.5)
          skipannotation(ys) && continue
          yi, index = findmax(ys)
          xi = xs[index]
          if !(isapprox(xi, minimum(k⊥s), rtol=0.01, atol=0.5) || 
               isapprox(yi, maximum(kzs), rtol=0.01, atol=0.01))
            continue
          end
          isapprox(xi, maximum(k⊥s), rtol=0.1, atol=0.5) && continue
          isapprox(yi, maximum(kzs), rtol=0.1, atol=0.5) && (yi += 0.075)
          isapprox(xi, minimum(k⊥s), rtol=0.1, atol=0.5) && (xi = -0.1)
          Plots.annotate!([(xi, yi, text("\$\\it{$lvl}\$", fontsize-1, :black))])
      end
    end
  end

  msize = 2
  mshape = :square
  function plotter2d(z, xlabel, ylabel, colorgrad,
      climmin=minimum(z[@. !ismissing(z)]), climmax=maximum(z[@. !ismissing(z)]))
    zcolor = make2d(z)
    dx = (xk⊥s[2] - xk⊥s[1]) / (length(xk⊥s) - 1)
    dy = (ykzs[2] - ykzs[1]) / (length(ykzs) - 1)
    h = Plots.heatmap(xk⊥s, ykzs, zcolor, framestyle=:box, c=colorgrad,
      xlims=(minimum(xk⊥s) - dx/2, maximum(xk⊥s) + dx/2),
      ylims=(minimum(ykzs) - dy/2, maximum(ykzs) + dy/2),
      clims=(climmin, climmax), xticks=0:Int(round(maximum(xk⊥s))),
      xlabel=xlabel, ylabel=ylabel)
  end
  xlabel = "\$\\mathrm{Perpendicular\\ Wavenumber} \\ [\\Omega_{i} / V_A]\$"
  ylabel = "\$\\mathrm{Parallel\\ Wavenumber} \\ [\\Omega_{i} / V_A]\$"
  #colorgrad = Plots.cgrad([:cyan, :lightblue, :blue, :darkblue, :black,
  #                        :darkred, :red, :orange, :yellow])
  zs = real.(ωs)
  climmax = maximum(zs)
  plotter2d(zs, xlabel, ylabel, Plots.cgrad(), 0.0, climmax)
  Plots.title!(" ")
  plotangles(writeangles=false)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_real_$file_extension.pdf")

  #ω0s = [fastmagnetoacousticfrequency(Va, vthd, sol.wavenumber) for
  ω0s = [fastzerobetamagnetoacousticfrequency(Va, sol.wavenumber, Ωd) for
    sol in sols] / f0
  zs = real.(ωs) ./ ω0s
  climmin = minimum(zs)
  climmax = maximum(zs)
  plotter2d(zs, xlabel, ylabel, Plots.cgrad(), climmin, climmax)
  Plots.title!(" ")
  plotangles(writeangles=false)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_real_div_guess_$file_extension.pdf")

  zs = iseven.(Int64.(floor.(real.(ωs))))
  climmax = maximum(zs)
  plotter2d(zs, xlabel, ylabel, Plots.cgrad(), 0.0, climmax)
  Plots.title!(" ")
  plotangles(writeangles=false)
  plotcontours(realspline, collect(1:50), y -> y[end] < 0)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_evenfloorreal_real_$file_extension.pdf")

  zs = imag.(ωs)
  climmax = maximum(zs)
  colorgrad = Plots.cgrad([:cyan, :black, :darkred, :red, :orange, :yellow])
  plotter2d(zs, xlabel, ylabel, colorgrad, -climmax / 4, climmax)
  Plots.title!(" ")
  Plots.plot!(legend=false)
  plotcontours(realspline, collect(1:50), y -> y[end] < 0)
  plotangles(writeangles=false)
  Plots.savefig("ICE2D_imag_$file_extension.pdf")

  colorgrad = Plots.cgrad()

  xlabel = "\$\\mathrm{Frequency} \\ [\\Omega_{i}]\$"
  ylabel = "\$\\mathrm{Parallel\\ Wavenumber} \\ [\\Omega_{i} / V_A]\$"

  imaglolim = 1e-5

  mask = shuffle(findall(@. (imag(ωs) > imaglolim) & (real(ωs) <= 12)))
  @warn "Scatter plots rendering with $(length(mask)) points."
  perm = sortperm(imag.(ωs[mask]))
  h0 = Plots.scatter(real.(ωs[mask][perm]), kzs[mask][perm],
     zcolor=imag.(ωs[mask][perm]), framestyle=:box, lims=:round,
    markersize=msize+1, markerstrokewidth=0, markershape=:circle,
    c=colorgrad, xticks=(0:12), yticks=unique(Int.(round.(ykzs))),
    xlabel=xlabel, ylabel=ylabel, legend=:topleft)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_KF12_$file_extension.pdf")


  ylabel = "\$\\mathrm{Growth\\ Rate} \\ [\\Omega_{i}]\$"
#  mask = shuffle(findall(imag.(ωs) .>= 0))
#  h = Plots.scatter(real.(ωs[mask]), imag.(ωs[mask]),
#    zcolor=kθs[mask] .* 180 / π, framestyle=:box, lims=:round,
#    markersize=msize, markerstrokewidth=0,
#    c=colorgrad,
#    xlabel=xlabel, ylabel=ylabel, legend=:topleft)
#  Plots.plot!(legend=false)
#  Plots.savefig("ICE2D_F_$file_extension.pdf")

  mask = shuffle(findall(@. (imag(ωs) > imaglolim) & (real(ωs) <= 12)))
  h1 = Plots.scatter(real.(ωs[mask]), imag.(ωs[mask]),
    zcolor=kθs[mask] .* 180 / π, framestyle=:box, lims=:round,
    markersize=msize+1, markerstrokewidth=0, markershape=:circle,
    c=colorgrad, xticks=(0:12),
    xlabel=xlabel, ylabel=ylabel, legend=:topleft)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_F12_$file_extension.pdf")

  colorgrad1 = Plots.cgrad([:cyan, :red, :blue, :orange, :green,
                            :black, :yellow])
  mask = shuffle(findall(@. (imag(ωs) > imaglolim) & (real(ωs) <= 12)))
  h2 = Plots.scatter(real.(ωs[mask]), imag.(ωs[mask]),
    zcolor=(real.(ωs[mask]) .- vαz/Va .* kzs[mask]), framestyle=:box, lims=:round,
    markersize=msize+1, markerstrokewidth=0, markershape=:circle,
    c=colorgrad1, clims=(0, 13), xticks=(0:12),
    xlabel=xlabel, ylabel=ylabel, legend=:topleft)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_F12_Doppler_$file_extension.pdf")

#  xlabel = "\$\\mathrm{Wavenumber} \\ [\\Omega_{i} / V_A]\$"
#  ylabel = "\$\\mathrm{Growth\\ Rate} \\ [\\Omega_{i}]\$"
#  mask = shuffle(findall(imag.(ωs) .>= 0))
#  h = Plots.scatter(ks[mask], real.(ωs[mask]), zcolor=kθs[mask] * 180 / π,
#    markersize=msize, markerstrokewidth=0, alpha=0.8, framestyle=:box, lims=:round,
#    xlabel=xlabel, ylabel=ylabel)
#  imagscale = 10
#  h3 = Plots.scatter!(ks[mask], imag.(ωs[mask]) * imagscale,
#    zcolor=kθs[mask] * 180 / π, framestyle=:box,
#    markersize=msize, markerstrokewidth=0, alpha=0.8,
#    xlabel=xlabel, ylabel=ylabel)
#  Plots.annotate!(minimum(ks) + 1, 1,
#                  text("\$ \\gamma \\times $imagscale\$", fontsize))
#  Plots.plot!(legend=false)
#  Plots.savefig("ICE2D_DR_$file_extension.pdf")

  xlabel = "\$\\mathrm{Frequency} \\ [\\Omega_{i}]\$"
  ylabel = "\$\\mathrm{Propagation\\ Angle} \\ [^{\\circ}]\$"
#  h = Plots.scatter(real.(ωs), kθs .* 180 / π, zcolor=imag.(ωs), lims=:round,
#    markersize=msize, markerstrokewidth=0, markershape=mshape, framestyle=:box,
#    c=colorgrad,
#    clims=(-maximum(imag.(ωs)), maximum(imag.(ωs))),
#    xlabel=xlabel, ylabel=ylabel)
#  Plots.plot!(legend=false)
#  Plots.savefig("ICE2D_θF_$file_extension.pdf")

  colorgrad = Plots.cgrad([:cyan, :black, :darkred, :red, :orange, :yellow])
  mask = findall(@. (imag(ωs) > imaglolim) & (real(ωs) <= 12))
  h4 = Plots.scatter(real.(ωs[mask]), kθs[mask] .* 180 / π,
    zcolor=imag.(ωs[mask]), lims=:round,
    markersize=msize, markerstrokewidth=0, markershape=mshape, framestyle=:box,
    c=Plots.cgrad([:black, :darkred, :red, :orange, :yellow]),
    clims=(0, maximum(imag.(ωs[mask]))),
    yticks=(0:10:180), xticks=(0:12), xlabel=xlabel, ylabel=ylabel)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_TF12_$file_extension.pdf")

  function relative(p, rx, ry)
    xlims = Plots.xlims(p)
    ylims = Plots.ylims(p)
    return xlims[1] + rx * (xlims[2]-xlims[1]), ylims[1] + ry * (ylims[2] - ylims[1])
   end
  Plots.xlabel!(h1, "")
  Plots.xticks!(h1, 0:-1)
  Plots.annotate!(h1, [(relative(h1, 0.02, 0.95)..., text("(a)", fontsize, :black))])
  Plots.annotate!(h0, [(relative(h0, 0.02, 0.95)..., text("(b)", fontsize, :black))])
  Plots.plot(h1, h0, link=:x, layout=@layout [a; b])
  Plots.savefig("ICE2D_Combo_$file_extension.pdf")
end

if true
#  for (plasma, ext) in ((Smmr, "ringbeam_"),
#                        (Smmd, "deltas_"),
#                        (Smmb, "beamdelta_"))
#    @time plasmasols = findsolutions(Smmr)
#    plasmasols = selectlargeestgrowthrate(plasmasols)
#    jldext = ext * "$name_extension"
#    @time plotit(plasmasols, jldext)
#    @save "solutions2D_$jldext.jld" filecontents plasmasols f0 k0
#  end
  @time plasmasols = findsolutions(Smmr)
  plasmasols = selectlargeestgrowthrate(plasmasols)
  @show length(plasmasols)
  @time plotit(plasmasols)
  @save "solutions2D_$name_extension.jld" filecontents plasmasols f0 k0
  rmprocs(nprocsadded)
else
  rmprocs(nprocsadded)
  @load "solutions2D_$name_extension.jld" filecontents solutions f0 k0
  @time plotit(solutions)
end

println("Ending at ", now())



using Distributed, Dates
using StatProfilerHTML

using FastClosures, Plots, Random, ImageFiltering, Statistics, StaticArrays
using Dierckx, Contour, JLD2, StaticArrays
using OwnTime, Profile, Coverage, PProf
Plots.pyplot()

println("Starting at ", now())
# pitch of 0 is all vperp, and (-)1 is (anti-)parallel to B field
const pitch = try; parse(Float64, ARGS[1]); catch; -0.64; end
@assert -1 <= pitch <= 1
# thermal width of ring as a fraction of its speed
const vthermalfraction = try; parse(Float64, ARGS[2]); catch; 0.01; end
const name_extension = length(ARGS) < 3 ? now() : ARGS[3]
const filecontents = [i for i in readlines(open(@__FILE__))]
const nprocsadded = 0Sys.CPU_THREADS # Int(floor(Sys.CPU_THREADS / 2))
addprocs(nprocsadded)

using LinearMaxwellVlasov, LinearAlgebra, WindingNelderMead
@everywhere using ProgressMeter # for some reason must be up here on its own
@everywhere begin
  using LinearMaxwellVlasov, LinearAlgebra, WindingNelderMead
  const LMV = LinearMaxwellVlasov

  mₑ = LMV.mₑ
  md = 2*1836*mₑ
  mα = 2*md

  # Fig 18 Cottrell 1993
  n0 = 1.7e19 # central electron density 3.6e19
  B0 = 2.07 # OR 2.191 for B field at 4m when 2.8 T on axis R0 3.13m
  ξ = 1.5e-4 # nα / ni = 1.5 x 10^-4

  nd = n0 / (1.0 + 2*ξ)
  nα = ξ*nd
  @assert n0 ≈ 2*nα + nd
  Va = sqrt(B0^2/LMV.μ₀/nd/md)

  Ωe = cyclotronfrequency(B0, mₑ, -1)
  Ωd = cyclotronfrequency(B0, md, 1)
  Ωα = cyclotronfrequency(B0, mα, 2)
  Πe = plasmafrequency(n0, mₑ, -1)
  Πd = plasmafrequency(nd, md, 1)
  Πα = plasmafrequency(nα, mα, 2)
  vthe = thermalspeed(1e3, mₑ)
  vthd = thermalspeed(1e3, md)
  vα = thermalspeed(3.6e6, mα)
  pitch = Float64(@fetchfrom 1 pitch) # pitchangle = acos(pitch)
  vα⊥ = vα * sqrt(1 - pitch^2) # perp speed
  vαb = vα * pitch # parallel speed
  vαth = thermalspeed(5.0e2, mα) # thermal speed
  vthermalfraction = Float64(@fetchfrom 1 vthermalfraction)
  vαth = vα * vthermalfraction # Dendy PRL 1993

  electron_cold = ColdSpecies(Πe, Ωe)
  electron_warm = WarmSpecies(Πe, Ωe, vthe)
  electron_maxw = MaxwellianSpecies(Πe, Ωe, vthe, vthe)

  deuteron_cold = ColdSpecies(Πd, Ωd)
  deuteron_warm = WarmSpecies(Πd, Ωd, vthd)
  deuteron_maxw = MaxwellianSpecies(Πd, Ωd, vthd, vthd)

  alpha_cold = ColdSpecies(Πα, Ωα)
  alpha_maxw = MaxwellianSpecies(Πα, Ωα, vαth, vαth, vαb)
  alpha_ringbeam = SeparableVelocitySpecies(Πα, Ωα,
    FBeam(vαth, vαb),
    FRing(vαth, vα⊥))
  alpha_delta = SeparableVelocitySpecies(Πα, Ωα,
    FParallelDiracDelta(vαb),
    FPerpendicularDiracDelta(vα⊥))
  alpha_beamdelta = SeparableVelocitySpecies(Πα, Ωα,
    FBeam(vαth, vαb),
    FPerpendicularDiracDelta(vα⊥))
  alpha_shell = CoupledVelocitySpecies(Πα, Ωα,
    FShell(vαth, vα))

  Smmr = [electron_maxw, deuteron_maxw, alpha_ringbeam]

  f0 = abs(Ωα)
  k0 = f0 / abs(Va)

  options = Options(memoiseparallel=false, memoiseperpendicular=true)

  function solve_given_ks(K, objective!)
    function bounds(ω0)
      γmax = abs(Ωα) * 0.15
      γmin = -abs(Ωα) * 0.075
      lb = @SArray [ω0 * 0.5, γmin]
      ub = @SArray [ω0 * 1.2, γmax]
      return (lb, ub)
    end
    solutions_out = Vector()
    atol = 1e-4
    ω0 = fastzerobetamagnetoacousticfrequency(Va, K, Ωd)

    lb, ub = bounds(ω0)

    function boundify(f::T) where T
      isinbounds(x) = all(i->0 <= x[i] <= 1, eachindex(x))
      maybeaddinf(x::U, addinf::Bool) where {U} = addinf ? x + U(Inf) : x
      bounded(x) = maybeaddinf(f(x), isinbounds(x))
      return bounded
    end

    config = Configuration(K, options)

    ics = ((@SArray [ω0*0.8, f0*0.08]),
           (@SArray [ω0*0.8, f0*0.04]),
           (@SArray [ω0, f0*0.01]))

    function unitobjective!(c, x::T) where {T}
      return objective!(c,
        T(x[i] * (ub[i] - lb[i]) + lb[i] for i in eachindex(x)))
    end
    boundedunitobjective! = boundify(x->unitobjective!(config, x))
    xtol_abs = f0 .* (@SArray [1e-3, 1e-4]) ./ (ub .- lb)
    @elapsed for ic ∈ ics
      @assert all(i->lb[i] <= ic[i] <= ub[i], eachindex(ic))
      neldermeadsol = WindingNelderMead.optimise(
        boundedunitobjective!, ((ic .- lb) ./ (ub .- lb)),
        (@SArray ones(2)) * 1.1e-2; stopval=atol, timelimit=30,
        maxiters=200, ftol_rel=0, ftol_abs=0, xtol_rel=0, xtol_abs=xtol_abs)
      simplex, windingnumber, returncode, numiterations = neldermeadsol
      if (abs(windingnumber) == 1 && returncode == :XTOL_REACHED) || returncode == :STOPVAL_REACHED
        c = deepcopy(config)
        minimiser = if windingnumber == 0
          WindingNelderMead.position(WindingNelderMead.bestvertex(simplex))
        else
          WindingNelderMead.centre(simplex)
        end
        unitobjective!(c, minimiser)
        return c
        break
      end
    end

    return nothing
  end

  function f2Dω!(config::Configuration, x::AbstractArray,
                 species, cache=Cache())
    config.frequency = ComplexF64(x[1], x[2])
    return det(tensor(Plasma(species), config, cache))
  end

  function findsolutions(species)
    ngridpoints = 2^9
    kzs = range(-2.0, stop=2.0, length=Int(ngridpoints)) * k0
    k⊥s = collect(range(0.0, stop=15.0, length=Int(ngridpoints))) * k0
    # change order for better distributed scheduling
    shuffle!(k⊥s)
    function inner(k⊥, _objective!::T) where {T}
      innersolutions = Vector()
      for (ikz, kz) ∈ enumerate(kzs)
        K = Wavenumber(parallel=kz, perpendicular=k⊥)
        output = solve_given_ks(K, _objective!)
        isnothing(output) && continue
        push!(innersolutions, output)
      end
      return innersolutions
    end

    #solutions = @sync @showprogress @distributed (vcat) for k⊥ ∈ k⊥s
    solutions = Vector()
    for k⊥ ∈ k⊥s
      cache = Cache()#paralleltype=ComplexF64, perpendiculartype=Float64)
      objective! = (C, x) -> f2Dω!(C, x, species, cache)
      #k⊥ == k⊥s[2] && @profilehtml [inner(k⊥, objective!) for _ in 1:10]
      #k⊥ == k⊥s[2] && @profilehtml inner(k⊥, objective!)
      #k⊥ == k⊥s[2] && throw(error("asdgfdfds"))
      #inner(k⊥, objective!)
      innersolutions = Vector()
      for (ikz, kz) ∈ enumerate(kzs)
        K = Wavenumber(parallel=kz, perpendicular=k⊥)
        iszero(K) && continue
        output = solve_given_ks(K, objective!)
        isnothing(output) && continue
        push!(innersolutions, output)
      end
      push!(solutions, innersolutions...)
    end
    return solutions
  end
end #@everywhere

function plotit(sols, fontsize=9)
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

  smoothing = length(ωs) * 1e-4
  _median(patch) = median(filter(x->!ismissing(x), patch))
  realωssmooth = make2d(real.(ωs))
  realωssmooth = mapwindow(_median, realωssmooth, (5, 5))
  realωssmooth = imfilter(realωssmooth, Kernel.gaussian(3))

  realspline = nothing
  imagspline = nothing
  try
    realspline = Dierckx.Spline2D(xk⊥s, ykzs, realωssmooth'; kx=4, ky=4, s=smoothing)
    imagspline = Dierckx.Spline2D(k⊥s, kzs, imag.(ωs); kx=4, ky=4, s=smoothing)
  catch
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
      clims=(climmin, climmax),
      xlabel=xlabel, ylabel=ylabel)
  end
  xlabel = "\$\\mathrm{Perpendicular\\ Wavenumber} \\quad [\\Omega_{i} / V_A]\$"
  ylabel = "\$\\mathrm{Parallel\\ Wavenumber} \\quad [\\Omega_{i} / V_A]\$"
  #colorgrad = Plots.cgrad([:cyan, :lightblue, :blue, :darkblue, :black,
  #                        :darkred, :red, :orange, :yellow])
  zs = real.(ωs)
  climmax = maximum(zs)
  plotter2d(zs, xlabel, ylabel, Plots.cgrad(), 0.0, climmax)
  Plots.title!(" ")
  plotangles(writeangles=false)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_real_$name_extension.pdf")

  #ω0s = [fastmagnetoacousticfrequency(Va, vthd, sol.wavenumber) for
  ω0s = [fastzerobetamagnetoacousticfrequency(Va, sol.wavenumber, Ωd) for
    sol in sols] / f0
  zs = real.(ωs) ./ ω0s
  climmax = maximum(zs)
  plotter2d(zs, xlabel, ylabel, Plots.cgrad(), 0.0, climmax)
  Plots.title!(" ")
  plotangles(writeangles=false)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_guess_div_real_$name_extension.pdf")

  zs = iseven.(Int64.(floor.(real.(ωs))))
  climmax = maximum(zs)
  plotter2d(zs, xlabel, ylabel, Plots.cgrad(), 0.0, climmax)
  Plots.title!(" ")
  plotangles(writeangles=false)
  plotcontours(realspline, collect(1:50), y -> y[end] < 0)
  Plots.plot!(legend=false)
  Plots.savefig("ICE2D_evenfloorreal_real_$name_extension.pdf")

  zs = imag.(ωs)
  climmax = maximum(zs)
  colorgrad = Plots.cgrad([:cyan, :black, :darkred, :red, :orange, :yellow])
  plotter2d(zs, xlabel, ylabel, colorgrad, -climmax / 4, climmax)
  Plots.title!(" ")
  Plots.plot!(legend=false)
  plotcontours(realspline, collect(1:50), y -> y[end] < 0)
  plotangles(writeangles=false)
  Plots.savefig("ICE2D_imag_$name_extension.pdf")

  colorgrad = Plots.cgrad()
  xlabel = "\$\\mathrm{Frequency} \\quad [\\Omega_{i}]\$"
  ylabel = "\$\\mathrm{Growth\\ Rate} \\quad [\\Omega_{i}]\$"
#  mask = shuffle(findall(imag.(ωs) .>= 0))
#  h = Plots.scatter(real.(ωs[mask]), imag.(ωs[mask]),
#    zcolor=kθs[mask] .* 180 / π, framestyle=:box, lims=:round,
#    markersize=msize, markerstrokewidth=0,
#    c=colorgrad,
#    xlabel=xlabel, ylabel=ylabel, legend=:topleft)
#  Plots.plot!(legend=false)
#  Plots.savefig("ICE2D_F_$name_extension.pdf")
  
  if true # lots of data takes too long to plot
    mask = shuffle(findall(@. (imag(ωs) > 0) & (real(ωs) <= 12)))
    h = Plots.scatter(real.(ωs[mask]), imag.(ωs[mask]),
      zcolor=kθs[mask] .* 180 / π, framestyle=:box, lims=:round,
      markersize=msize+1, markerstrokewidth=0,
      c=colorgrad,
      xlabel=xlabel, ylabel=ylabel, legend=:topleft)
    Plots.plot!(legend=false)
    Plots.savefig("ICE2D_F12_$name_extension.pdf")
  end

#  xlabel = "\$\\mathrm{Wavenumber} \\quad [\\Omega_{i} / V_A]\$"
#  ylabel = "\$\\mathrm{Growth\\ Rate} \\quad [\\Omega_{i}]\$"
#  mask = shuffle(findall(imag.(ωs) .>= 0))
#  h = Plots.scatter(ks[mask], real.(ωs[mask]), zcolor=kθs[mask] * 180 / π,
#    markersize=msize, markerstrokewidth=0, alpha=0.8, framestyle=:box, lims=:round,
#    xlabel=xlabel, ylabel=ylabel)
#  imagscale = 10
#  h = Plots.scatter!(ks[mask], imag.(ωs[mask]) * imagscale,
#    zcolor=kθs[mask] * 180 / π, framestyle=:box,
#    markersize=msize, markerstrokewidth=0, alpha=0.8,
#    xlabel=xlabel, ylabel=ylabel)
#  Plots.annotate!(minimum(ks) + 1, 1,
#                  text("\$ \\gamma \\times $imagscale\$", fontsize))
#  Plots.plot!(legend=false)
#  Plots.savefig("ICE2D_DR_$name_extension.pdf")

  xlabel = "\$\\mathrm{Frequency} \\quad [\\Omega_{i}]\$"
  ylabel = "\$\\mathrm{Propagation\\ Angle} \\quad [^{\\circ}]\$"
#  h = Plots.scatter(real.(ωs), kθs .* 180 / π, zcolor=imag.(ωs), lims=:round,
#    markersize=msize, markerstrokewidth=0, markershape=mshape, framestyle=:box,
#    c=colorgrad,
#    clims=(-maximum(imag.(ωs)), maximum(imag.(ωs))),
#    xlabel=xlabel, ylabel=ylabel)
#  Plots.plot!(legend=false)
#  Plots.savefig("ICE2D_θF_$name_extension.pdf")

  if false # lots of data takes too long to plot
    colorgrad = Plots.cgrad([:cyan, :black, :darkred, :red, :orange, :yellow])
    mask = findall(@. ((real(ωs) < 12) & (imag(ωs) > -f0 / 10000)))
    h = Plots.scatter(real.(ωs[mask]), kθs[mask] .* 180 / π,
      zcolor=imag.(ωs[mask]), lims=:round,
      markersize=msize, markerstrokewidth=0, markershape=mshape, framestyle=:box,
      c=colorgrad,
      clims=(-maximum(imag.(ωs[mask])) / 4, maximum(imag.(ωs[mask]))),
      xlabel=xlabel, ylabel=ylabel)
    Plots.plot!(legend=false)
    Plots.savefig("ICE2D_θF12_$name_extension.pdf")
  end
end
# Select the largest growth rates if multiple solutions found for a wavenumber
function selectlargeestgrowthrate!(sols)
    matchwavenumber(a, sols) = sols[findall(s->a.wavenumber == s.wavenumber, sols)]
    imagfreq(s) = imag(complexfrequency(s.frequency))
    sols = filter(s->imagfreq(s) == maximum(imagfreq.(matchwavenumber(s, sols))),
                  sols)
    return nothing
end
function selectlargeestgrowthrate(sols)
    imagfreq(s) = imag(complexfrequency(s.frequency))
    d = Dict{Any, Vector{eltype(sols)}}()
    for sol in sols
      if haskey(d, sol.wavenumber)
        push!(d[sol.wavenumber], sol)
      else
        d[sol.wavenumber] = [sol]
      end
    end
    sols = Vector()
    sizehint!(sols, length(sols))
    for (_, ss) in d
      push!(sols, ss[findmax(map(imagfreq, ss))[2]])
    end
    return sols
end

if true
  #@time solutions = findsolutions(Smmr)
  Profile.clear_malloc_data()
  Profile.init(n=10^8, delay=0.001)
  @profilehtml findsolutions(Smmr)
  #@profile findsolutions(Smmr); #Profile.print(format=:flat, sortedby=:count)
  #pprof(); sleep(1000000)
  #@show owntime()
  rmprocs(nprocsadded)
  #solutions = selectlargeestgrowthrate(solutions)
  #@save "solutions2D_$name_extension.jld" filecontents solutions f0 k0
else
  rmprocs(nprocsadded)
  #@load "solutions2D_$name_extension.jld" filecontents solutions f0 k0
end

#@time plotit(solutions)
println("Ending at ", now())

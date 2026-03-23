using Dates
println("Starting Ion Bernstein at ", now())

using LinearMaxwellVlasov
using Plots, Random, LinearAlgebra, WindingNelderMead

Random.seed!(0) # seed rand

const mₑ = LinearMaxwellVlasov.mₑ
const md = 2*1836*mₑ
const B0 = 1.0
const Ωd = cyclotronfrequency(B0, md, 1)
const Ωe = -Ωd * md / mₑ
const ratio = 10
const Πe = abs(Ωe) * sqrt(ratio)
const n0 = Πe^2 * LinearMaxwellVlasov.ϵ₀ * mₑ / LinearMaxwellVlasov.q₀^2
@assert (Πe / Ωe)^2 ≈ ratio
const Πd = plasmafrequency(n0, md, 1)
const ϵV = 1e3
const vthe = thermalspeed(ϵV, mₑ)
const vthd = thermalspeed(ϵV, md)
const ρd = vthd / Ωd
const electron_cold = ColdSpecies(Πe, Ωe)
const electron_warm = WarmSpecies(Πe, Ωe, vthe)
const electron_maxw = MaxwellianSpecies(Πe, Ωe, vthe, vthe)
const ratio⊥ = 10
const deuteron_maxw = MaxwellianSpecies(Πd, Ωd, vthd, vthd * ratio⊥)
const deuteron_numerical = SeparableVelocitySpecies(Πd, Ωd,
        FParallelNumerical(vthd),
        FPerpendicularNumerical(vthd * ratio⊥))
const deuteron_ringbeam = RingBeamSpecies(Πd, Ωd, vthd, vthd * ratio⊥)

const Smaxw = Plasma((electron_cold, deuteron_maxw))
const Srb = Plasma((electron_cold, deuteron_ringbeam))
const Snum = Plasma((electron_cold, deuteron_numerical))

const k0 = 1 / ρd
const f0 = abs(Ωd)

const options = Options(memoiseperpendicular=true)

function icsandbounds(ω0, γ0)
  ub = [ω0 + Ωd, abs(γ0) + Ωd/10]
  lb = [ω0 - Ωd, -abs(γ0) - Ωd]
  ic = [ω0, iszero(γ0) ? - Ωd/1000 : γ0]
  @assert all(lb .< ic .< ub)
  return (ic, lb, ub)
end

function local_solve(ks, ωs, f)
  solutions_out = Vector{Any}()
  kωs = Vector()
  for kc in ks, ωc in ωs push!(kωs, (kc, ωc)) end
  for (kc, ωc) in kωs
    K = Wavenumber(0.0, kc * k0)
    @assert K.parallel == 0
    C = Configuration(K, options)
    ic, lb, ub = icsandbounds(ωc * Ωd, 0.0)
    bounded(x) = all(lb .<= x .<= ub) ? f(x, C) : Inf + Inf*im 
    t1 = @elapsed neldermeadsol = WindingNelderMead.optimise(bounded,
      ic, [1, 1] .* Ωd * 1e-2; stopval=1e-3, timelimit=1000,
      maxiters=1000, ftol_rel=0, ftol_abs=0, xtol_rel=0, xtol_abs=1e-15*Ωd)
    simplex, windingnumber, returncode, numiterations = neldermeadsol
    @show t1, windingnumber, returncode, numiterations
    if (windingnumber == 1 && returncode == :XTOL_REACHED) || returncode == :STOPVAL_REACHED
      c = deepcopy(C)
      push!(solutions_out, deepcopy(C))
    end
  end
  return solutions_out
end

function f2Dω!(x::Vector{Float64}, C::Configuration, S, cache)
  C.frequency = ComplexF64(x[1], x[2])
  return electrostatic(S, C, cache)
end

const cache = LinearMaxwellVlasov.Cache()
f2Dωmaxw!(x::Vector{Float64}, C::Configuration) = f2Dω!(x, C, Smaxw, cache)
f2Dωrb!(x::Vector{Float64}, C::Configuration) = f2Dω!(x, C, Srb, cache)
f2Dωnum!(x::Vector{Float64}, C::Configuration) = f2Dω!(x, C, Snum, cache)

const ks = reverse(collect(range(0.01, stop=2.0, length=128)))
#ks = sort(vcat(-ks, ks))
const ωs = range(0.0, stop=10, length=20)
solutions_rb = Any[]
solutions_num = Any[]
solutions_maxw = Any[]
#@time solutions_num = local_solve(ks, ωs, f2Dωnum!)
@time solutions_maxw = local_solve(ks, ωs, f2Dωmaxw!)
@time solutions_rb = local_solve(ks, ωs, f2Dωrb!)

for (solutions, str) in ((solutions_maxw, "_maxw"), (solutions_rb, "_rb"), (solutions_num, "_num"))
  @show isempty(solutions)
  isempty(solutions) && continue
  ω1 = [solution.frequency for solution in solutions]./f0
  kb1 = [perp(solution.wavenumber) for solution in solutions]./k0
  k1 = [perp(solution.wavenumber) for solution in solutions]./k0

  validω = @. abs(real(ω1)) < 500
  unstable1 = @. (imag(ω1) > 0) & validω
  stable1 = @. (!unstable1) & validω

  h = Plots.scatter(k1[unstable1], imag(ω1[unstable1]))
  Plots.scatter!(k1[unstable1], real(ω1[unstable1]))
  Plots.scatter!(k1[stable1], imag(ω1[stable1]))
  Plots.scatter!(k1[stable1], real(ω1[stable1]))
  Plots.display(h)
  Plots.savefig(h, "IonBernstein_$str.png")
end

println("Ending at ", now())









#

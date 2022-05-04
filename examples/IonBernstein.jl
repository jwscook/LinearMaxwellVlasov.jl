using Dates
println("Starting Ion Bernstein at ", now())

using LinearMaxwellVlasov
, ..Solvers
using ,..Species, ..Solutions
using Optim, Plots, Base.Threads, Suppressor, Random
Plots.plotly()

Random.seed!(0) # seed rand

const mₑ = mₑ
const md = 2*1836*mₑ
const mα = 2*md
const n0 = 1.0e19
const B0 = 1.0e-3
const nd = n0
const θ = 88.0 * π/180
const Va = sqrt(B0^2/μ₀/nd/md)
const Ωe = cyclotronfrequency(B0, mₑ, -1)
const Ωd = cyclotronfrequency(B0, md, 1)
const Πe = plasmafrequency(n0, mₑ, -1)
const Πd = plasmafrequency(nd, md, 1)
const ϵV = 1.0e3
const vthe = thermalspeed(ϵV, mₑ)
const vthd = thermalspeed(ϵV, md)
const ratio = 10
electron_cold = ColdSpecies(Πe, Ωe)
electron_warm = WarmSpecies(Πe, Ωe, vthe)
electron_maxw = MaxwellianSpecies(Πe, Ωe, vthe, vthe)
deuteron_maxw = MaxwellianSpecies(Πd, Ωd, vthd, vthd * ratio)
deuteron_numerical = NumericalSpecies(Πd, Ωd,
  FParallelNumerical(vthd),
  FPerpendicularNumerical(vthd * ratio))
oneover2pi = 1.0 / (2*pi)
gyroperturbation = Dict{Int, Float64}(0=>oneover2pi, 1=>oneover2pi/2)
deuteron_ringbeam = RingBeamSpecies(Πd, Ωd, vthd, vthd * ratio)
deuteron_gyroringbeam = RingBeamSpecies(Πd, Ωd, vthd, vthd * ratio, 0.0,
  0.0, FGyroangleCosine(gyroperturbation))
deuteron_gyronum = NumericalSpecies(Πd, Ωd, vthd, vthd * ratio, 0.0, 0.0,
  FGyroangleCosine(gyroperturbation))

Smaxw = [electron_maxw, deuteron_maxw]
Srb = [electron_maxw, deuteron_ringbeam]
Snum = [electron_maxw, deuteron_numerical]
Sgyro = [electron_maxw, deuteron_gyronum]

k0 = abs(Ωd / Va)
f0 = abs(Ωd)

options = Options()

function icsandbounds_ω(ω0, γ0)
  ub = [ω0 + Ωd, abs(γ0) + Ωd/10]
  lb = [ω0 - Ωd, -abs(γ0) - Ωd]
  ic = [ω0, iszero(γ0) ? - Ωd/1000 : γ0]
  @assert all(lb .< ic .< ub)
  return (ic, lb, ub)
end

#method = :NLopt
method = :Optim
#method = :WindingNumber
function local_solve(ks, ωs, f, icsandbounds)
  solutions_out = Vector{Solutions.Solution}()
  sl = Threads.SpinLock()
  kωs = Vector()
  for kc in ks, ωc in ωs push!(kωs, (kc, ωc)) end
  #@suppress @threads for (kc, ωc) in kωs
  for (kc, ωc) in kωs
    K = Wavenumber(k=kc * k0, θ=θ)
    C = Configuration(K, options)
    ic, lb, ub = icsandbounds(ωc * Ωd, 0.0)
    t = @elapsed solution = Solvers.solve(f, C, ic, lb, ub, method)
    lock(sl)
    push!(solutions_out, solution...)
    unlock(sl)
  end
  return solutions_out
end
function local_solve(solutions_in, f, icsandbounds)
  solutions_out = Vector{Solutions.Solution}()
  sl = Threads.SpinLock()
  #@suppress @threads for sol in solutions_in
  for sol in solutions_in
    C = Configuration(sol.frequency, sol.wavenumber, options)
    ic, lb, ub = icsandbounds(reim(sol.frequency)...)
    t = @elapsed solution = Solvers.solve(f, C, ic, lb, ub, method)
    lock(sl)
    push!(solutions_out, solution...)
    unlock(sl)
  end
  return solutions_out
end


function f2Dω!(x::Vector{Float64}, C::Configuration, S)
  C.frequency = ComplexF64(x[1], x[2])
  t = @elapsed output = det(tensor(S, C))
  return output
end

f2Dωmaxw!(x::Vector{Float64}, C::Configuration) = f2Dω!(x, C, Smaxw)
f2Dωrb!(x::Vector{Float64}, C::Configuration) = f2Dω!(x, C, Srb)
f2Dωnum!(x::Vector{Float64}, C::Configuration) = f2Dω!(x, C, Snum)
f2Dωgyro!(x::Vector{Float64}, C::Configuration) = f2Dω!(x, C, Sgyro)

ks = reverse(collect(range(0.1, stop=10.0, length=4)))
#ks = sort(vcat(-ks, ks))
ωs = range(1.5, stop=5.5, length=4)
solutions_gyro = Vector{Solutions.Solution}()
solutions_rb = Vector{Solutions.Solution}()
solutions_num = Vector{Solutions.Solution}()
solutions_maxw = Vector{Solutions.Solution}()
#@time solutions_gyro = local_solve(ks, ωs, f2Dωgyro!, icsandbounds_ω)
#@time solutions_num = local_solve(ks, ωs, f2Dωnum!, icsandbounds_ω)
@time solutions_maxw = local_solve(ks, ωs, f2Dωmaxw!, icsandbounds_ω)
#@time solutions_rb = local_solve(ks, ωs, f2Dωrb!, icsandbounds_ω)

#@time solutions_rb = local_solve(solutions_maxw, f2Dωrb!, icsandbounds_ω)
#@time solutions_num = local_solve(solutions_maxw, f2Dωnum!, icsandbounds_ω)
#@time solutions_gyro = local_solve(solutions_maxw, f2Dωgyro!, icsandbounds_ω)
for solutions in (solutions_maxw, solutions_rb, solutions_num, solutions_gyro)
  isempty(solutions) && continue
  ω1 = [solution.frequency for solution in solutions]./f0
  kb1 = [para(solution.wavenumber) for solution in solutions]./k0
  k1 = [para(solution.wavenumber) for solution in solutions]./k0
  converged1 = [solution.converged for solution in solutions]
  notconverged1 = [!solution.converged for solution in solutions]
  values1 = [solution.value for solution in solutions]

  validω = @. abs(real(ω1)) < 500
  unstable1 = @. (imag(ω1) > 0) & validω
  stable1 = @. (!unstable1) & validω

  h = Plots.scatter(k1[unstable1], imag(ω1[unstable1]))#, c=values1)
  Plots.scatter!(k1[unstable1], real(ω1[unstable1]))#, c=values1)
  Plots.scatter!(k1[stable1], imag(ω1[stable1]))#, c=values1)
  Plots.scatter!(k1[stable1], real(ω1[stable1]))#, c=values1)
  Plots.display(h)
end

println("Ending at ", now())











#

using Dates
println("Starting Langmuir Ion Acoustic waves at ", now())
# addprocs(4-nprocs())
# @everywhere begin
using LinearMaxwellVlasov


using ,..Integrals, ..Solvers

using Optim, PyPlot
mₑ = mₑ
mi = 1836*mₑ
n0 = 1.0e19
B0 = 0.0
Ωe = cyclotronfrequency(B0, mₑ, -1)
Ωi = cyclotronfrequency(B0, mi, 1)
Πe = plasmafrequency(n0, mₑ, -1)
Πi = plasmafrequency(n0, mi, 1)
ϵV = 1.0e2
vthe = thermalspeed(ϵV, mₑ)
vthi = thermalspeed(ϵV, mi)
λD = vthe / Πe
electron = NumericalSpecies(Πe, Ωe, FParallelNumerical(vthe), FPerpendicularNumerical(vthe))
proton = NumericalSpecies(Πi, Ωi, FParallelNumerical(vthi), FPerpendicularNumerical(vthi))

#electron = MaxwellianSpecies(Πe, Πe*1.0e-3, vthe, vthe)
#proton = MaxwellianSpecies(Πi, Πi*1.0e-3, vthi, vthi)

S = [electron, proton]
N = 20
ks = Vector{Float64}(range(0.01, stop=0.25, length=N))
O = Options()
#O.quadrature_tol.abs = eps()
#O.quadrature_tol.rel = sqrt(eps())
O.summation_tol.rel = 1.0e-3 #sqrt(eps())
O.summation_tol.abs = 1.0e-3 #sqrt(eps())
O.solution_tol.abs  = 1.0e-3 #sqrt(eps())
O.solution_tol.rel  = 1.0e-3 #sqrt(eps())

parallelisation = :serial
# end # @everywhere
method = :Optim #:NLopt
solutions1 = Vector{Solutions.Solution}()
solutions2 = Vector{Solutions.Solution}()
function f!(x::Vector{Float64}, C::Configuration)
  C.frequency = ComplexF64(x[1], x[2])
  return abs_tensor(S, C)
end
if parallelisation == :serial
  @timev begin
    for i = 1:N
      K = Wavenumber(k=ks[i] * 2π/λD, θ=0.01)
      C = Configuration(K, O)
Profile.clear()
Profile.init(n=10^7, delay=0.05)
      @profile push!(solutions1, Solvers.solve(f!, C, [1.5, -0.1]*Πe, method)...)
      @profile push!(solutions2, Solvers.solve(f!, C, [1.0, -0.05]*Πi, method)...)
Profile.print(format=:flat, combine=true, sortedby=:count)
    end
  end
elseif parallelisation == :threaded
  @timev begin
    Threads.@threads for i = 1:N
    end
  end
elseif parallelisation == :shared
  @assert nprocs() > 1
  ω1 = SharedArray(ω1)
  ω2 = SharedArray(ω2)
  @timev begin
    @parallel for i = 1:N
      #@async
      #@async
    end
  end
end
@show ω1 = [solution.frequency for solution in solutions1]./Πe
@show ω2 = [solution.frequency for solution in solutions2]./Πe
@show values1 = [solution.value for solution in solutions1]
@show values2 = [solution.value for solution in solutions2]


# PyPlot.pygui(true)
# PyPlot.show()

PyPlot.figure()
#PyPlot.plot(ks, real(ω1), "co", ks, real(ω2), "ko")
#PyPlot.plot(ks, imag(ω1), "cs", ks, imag(ω2), "ks")

PyPlot.scatter(ks, real(ω1), c=values1)
PyPlot.scatter(ks, real(ω2), c=values2)
PyPlot.scatter(ks, imag(ω1), c=values1)
PyPlot.scatter(ks, imag(ω2), c=values2)
PyPlot.colorbar()

PyPlot.plot(ks, sqrt.(Πe.^2 + 3*vthe^2*(ks*2*pi/λD).^2)/Πe)
# PyPlot.plot(ks, sqrt(Πe.^2 + (2.99e8)^2*(ks*2*pi/λD).^2)/Πe)
PyPlot.plot(ks, (ks*2*pi/λD).*sqrt(q0*(ϵV + 3*ϵV)/mi)/Πe)
vphase1 = real(ω1)./ks * Πe/(2π/λD)
vphase2 = real(ω2)./ks * Πe/(2π/λD)
#PyPlot.plot(ks, -pi*Πe^2./(ks*2π/λD).^2 .* electron.Fb.dFdv(vphase1) / Πe)
# PyPlot.plot(ks, -pi*Πi^2./(ks*2π/λD).^2 .* proton.Fb.dFdv(vphase2) / Πe)

PyPlot.pygui(true)
PyPlot.show()

println("Ending at ", now())

# rmprocs(3)











#

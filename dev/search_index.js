var documenterSearchIndex = {"docs":
[{"location":"#LinearMaxwellVlasov.jl-Documentation","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.jl Documentation","text":"","category":"section"},{"location":"","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.jl Documentation","text":"CurrentModule = LinearMaxwellVlasov","category":"page"},{"location":"","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.jl Documentation","text":"Modules = [LinearMaxwellVlasov]","category":"page"},{"location":"#Core.Type-Tuple{LinearMaxwellVlasov.AbstractKineticSpecies, Configuration, Int64}","page":"LinearMaxwellVlasov.jl Documentation","title":"Core.Type","text":"The function parallel gives the same answers for a given set of inputs. Calculate the key given the inputs and the (identity) operator for storing and retrieving values\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.CacheOp-Tuple{Any}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.CacheOp","text":"Identity operator for storing and retreiving values from CacheDict\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.ColdSpecies","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.ColdSpecies","text":"ColdSpecies <: AbstractFluidSpecies\n\nCold plasma species. Fields:\n\nΠ Plasma frequency [rad / s]\nΩ Cyclotron Frequency [rad / s]\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.CoupledRelativisticSpecies","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.CoupledRelativisticSpecies","text":"Kinetic plasma species defined by one coupled distribution function in momentum space such that the relativistic dielectric tensor can be calculated. Fields:\n\nΠ Plasma frequency [rad / s]\nΩ Cyclotron Frequency [rad / s]\nmass Species particle mass [kg]\nF :: AbstractFRelativisticMomentum Distribution function in momentum space parallel and perpendicular to the background magnetic field (normalised)\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.CoupledVelocitySpecies","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.CoupledVelocitySpecies","text":"Kinetic plasma species defined by one coupled distribution function in velocity space, parallel and perpendicular to the background magnetic field Fields:\n\nΠ Plasma frequency [rad / s]\nΩ Cyclotron Frequency [rad / s]\nF :: AbstractCoupledVelocity Distribution function in velocity space parallel and perpendicular to the background magnetic field (normalised)\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.FCoupledVelocityNumerical","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.FCoupledVelocityNumerical","text":"FCoupledVelocityNumerical\n\nA disribution function where vz and v⊥ are coupled, i.e. non-separable.\n\n...\n\nArguments\n\nF::T: the distrubtion function\nnormalisation::Tuple{U,U}: the speeds used for normalisation in parallel and perp [m/s]\nlower::Float64: minimum speed for integration bounds [m/s]\nupper::Float64: maximum speed for integration bounds [m/s]\n\n...\n\nExample\n\nvth = 1e4\nvshell = 1e5\nfshell = FShell(vth, vshell)\nlower = vshell - 12 * vth # where the shell is zero\nupper = vshell + 12 * vth # where the shell is zero\nFCoupledVelocityNumerical(fshell, (vshell, vshell), lower, upper)\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.FParallelDiracDelta","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.FParallelDiracDelta","text":"The Dirac-delta function parallel distribution functions\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.FParallelNumerical","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.FParallelNumerical","text":"Type for an arbitrary parallel distribution function that holds, a function for the distribution function, it's derivative and the limits of integral\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.FPerpendicularDiracDelta","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.FPerpendicularDiracDelta","text":"The Dirac-delta function parallel distribution functions\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.FPerpendicularNumerical","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.FPerpendicularNumerical","text":"Type for an arbitrary perpendicular distribution function that holds, a function for the distribution function, it's derivative and the limits of integral\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.MaxwellianIntegralsParallel-NTuple{5, Any}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.MaxwellianIntegralsParallel","text":"The parallel integral of the Beam only requires an integral over a drifting Maxwellian subject to the relevant kernels. All this is calculated here.\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.MaxwellianIntegralsPerpendicular-Tuple{Number, Number, Integer}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.MaxwellianIntegralsPerpendicular","text":"Perpendicular integrals for a Maxwellian\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.PerpendicularKernel","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.PerpendicularKernel","text":"The kernel of the perpendicular integral\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.RelativisticMaxwellian","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.RelativisticMaxwellian","text":"RelativisticMaxwellian Return a drifting Maxwellian function.\n\nMembers\n\npthz::Real - thermal momentum parallel to magnetic field [kg m/s] pth⊥::Real - thermal momentum perpendicular to magnetic field [kg m/s] pzdrift::Real=0.0 - drift parallel to the magnetic field [kg m/s] lognormalisation::Real - the log of the normalisation constant\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.SeparableVelocitySpecies","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.SeparableVelocitySpecies","text":"Kinetic plasma species with separable distribution functions parallel and perpendicular to the magnetic field. Fields:\n\nΠ Plasma frequency [rad / s]\nΩ Cyclotron Frequency [rad / s]\nFz :: AbstractFParallel Distribution function parallel to magnetic field (normalised)\nF⊥ :: AbstractFPerpendicular Distribution function perpendicular to magnetic field (normalised)\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.ShiftedMaxwellianParallel","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.ShiftedMaxwellianParallel","text":"The probability density at parallel velocity, v, of a shifted maxwellian from thermal velocity vth, and drift velocity vd.\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.ShiftedMaxwellianPerpendicular","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.ShiftedMaxwellianPerpendicular","text":"The probability density at perpendicular velocity, v, of a shifted maxwellian from thermal velocity vth, and drift velocity vd.\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.TransformFromInfinity","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.TransformFromInfinity","text":"Transform a function from domain [-∞, ∞]ⁿ down to [-1, 1]ⁿ\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.UnitSemicircleIntegrandTransform","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.UnitSemicircleIntegrandTransform","text":"Transform a function to transform an integrand from domain [-∞, 0]×[∞, ∞] down to [-1, -π/2]×[1, π/2].\n\nExample:\n\njulia> f(x) = exp(-(x[1]^2 + x[2]^2)/2) * cos(x[2])^2 * sin(x[1])^2\njulia> hcubature(f, [-12.0, 0.0], [12.0, 12.0])\n(0.7710130943379178, 1.1482318484139944e-8)\njulia> hcubature(UnitSemicircleIntegrandTransform(f, 2.0), [0, -π/2], [1, π/2])\n(0.7710130940919337, 1.1464196819507393e-8)\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.WarmSpecies","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.WarmSpecies","text":"Warm plasma species, with speeds of sound parallel and perpendicular to the magnetic field. Fields:\n\nΠ Plasma frequency [rad / s]\nΩ Cyclotron Frequency [rad / s]\nsoundspeed Sound speed [m/s]\n\n\n\n\n\n","category":"type"},{"location":"#LinearMaxwellVlasov.WarmSpecies-NTuple{4, Any}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.WarmSpecies","text":"WarmSpecies(Π::Float64,Ω::Float64,thermalspeed::Float64,adiabiaticindex::Number)\n\nWarm plasma species - accept thermalspeed and ratio of specific heats to get sound speed\n\n...\n\nArguments\n\nΠ: Plasma frequency [rad / s]\nΩ: Cyclotron Frequency [rad / s]\nthermalspeed: Thermal speed of Maxwellian distribution\nadiabiaticindex: Equation of state Gruneisen gamma (ratio of specific heats)\n\n...\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.Wavenumber","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.Wavenumber","text":"Wavenumber decomposed into parallel and perpendicular components Construct with parallel and wavenumber components, or by keyword arguement pairs\n\nparallel & perpendicular\nwavenumber & propagationangle\n\n\n\n\n\n","category":"type"},{"location":"#Base.Filesystem.cp-Union{Tuple{T}, Tuple{U}} where {U, T<:NTuple{11, U}}","page":"LinearMaxwellVlasov.jl Documentation","title":"Base.Filesystem.cp","text":"Regarding integrals with besselj(n, x)* besselj(n±1, x), e.g.: where m >= 0\n\nbesselj(-m,x) * besselj(-m-1,x) == -besselj(m,x) * besselj(m+1,x)\n\nbesselj(-m+1,x) * besselj(-m-1,x) == +besselj(m-1,x) * besselj(m+1,x)\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.FShell-Tuple{Real, Real}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.FShell","text":"FShell(vth::Real,vshell::Real)\n\nThe shell distribution function, as though f is only non-zero on or close to the surface of a sphere.\n\n...\n\nArguments\n\nvth::Real: the thermal velocity of the shell [m/s]\nvshell::Real: the speed of the shell [m/s]\n\n...\n\nExample\n\n\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.FSlowingDown-Tuple{Real, Real, Real}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.FSlowingDown","text":"FSlowingDown(vbeam::Real,vcrit::Real,vcutoffwidth::Real)\n\nThe slowing down distribution\n\n...\n\nArguments\n\nvbeam::Real: the beam speed [m/s]\nvcrit::Real: the critival velocity [m/s]\nvcutoffwidth::Real: the width of the error function used to smooth\n\nthe distribution function at vbeam [m/s] ...\n\nExample\n\n\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.Jμν-Union{Tuple{T}, Tuple{Pair{var\"#s38\", var\"#s37\"} where {var\"#s38\"<:Integer, var\"#s37\"<:Integer}, T}} where T<:Number","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.Jμν","text":"The non-integer argument to the BesselJ\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.MaxwellianSpecies","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.MaxwellianSpecies","text":"MaxwellianSpecies(Π,Ω,vthb,vth⊥=vthb,vdb=0.0)\n\nKinetic Maxwellian Plasma species that can optionally have a drift along the magnetic field\n\n...\n\nArguments\n\nΠ: Plasma frequency [rad / s]\nΩ: Cyclotron Frequency [rad / s]\nvthb: parallel thermal speed [m/s]\nvth⊥=vthb: perpendicular thermal speed [m/s]\nvdb=0.0: parallel beam speed [m/s]\n\nReturns\n\nSeparableVelocitySpecies(Π, Ω, FBeam(vthb, vdb), FPerpendicularMaxwellian(vth⊥))\n\n...\n\nExample\n\n\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.RingBeamSpecies","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.RingBeamSpecies","text":"Create a kinetic plasma species with separable distribution functions parallel f(v_parallel) and perpendicular f(v_perpendicular) to the magnetic field, which are defined as a drifting beam and a ring respectively. ...\n\nArguments\n\nΠ: Plasma frequency [rad / s]\nΩ: Cyclotron Frequency [rad / s]\nvthb: parallel thermal speed [m/s]\nvth⊥=vthb: perpendicular thermal speed [m/s]\nvdb=0.0: parallel beam speed [m/s]\nvd⊥=0.0: perpendicular ring speed [m/s]\n\nReturns\n\nSeparableVelocitySpecies(Π, Ω, FBeam(vthb, vdb), FRing(vth⊥, vd⊥))\n\n...\n\nExample\n\n\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.alfvenspeed-Tuple{Any}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.alfvenspeed","text":"Calculate the Alfven speed given Π_Ωs, which is an iterable container of the ratio of the plasma to the cyclotron frequency of all the species of the plasma\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.conductivity","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.conductivity","text":"The conducivity tensor for a given species\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.contribution","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.contribution","text":"Calculate the converged unitless susceptibility tensor contribution for a maxwellian distribution function, having summed over bessel indices n.\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.contribution-2","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.contribution","text":"Calculate the unitless susceptibility tensor for a warm plasma species Swanson 3.63\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.contribution-3","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.contribution","text":"Calculate the unitless susceptibility tensor for a cold plasma species\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.contribution-4","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.contribution","text":"Calculate the converged unitless susceptibility tensor contribution for a maxwellian distribution function, having summed over bessel indices n.\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.contribution-Union{Tuple{T}, Tuple{T⊥}, Tuple{Tz}, Tuple{V}, Tuple{U}, Tuple{T, Configuration, Int64}, Tuple{T, Configuration, Int64, U}, Tuple{T, Configuration, Int64, U, V}} where {U<:Function, V<:Function, Tz, T⊥<:FPerpendicularMaxwellian, T<:(LinearMaxwellVlasov.AbstractSeparableVelocitySpecies{var\"#s5\", var\"#s4\"} where {var\"#s5\"<:Tz, var\"#s4\"<:T⊥})}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.contribution","text":"Calculate the unitless susceptibility tensor for a maxwellian distribution function given the bessel indices n.\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.dielectric","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.dielectric","text":"The dielectric tensor for a given plasma\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.dielectriccontribution","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.dielectriccontribution","text":"The contribution to the dielectric tensor for a given species\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.electrostatictensor","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.electrostatictensor","text":"The electrostatic dielectric tensor for a given plasma, a zero valued determinant of which represents a solution to the linear poisson-vlasov system of equations\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.foldnumeratoraboutpole-Union{Tuple{T}, Tuple{T, Real}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.foldnumeratoraboutpole","text":"Takes a function that when integrated between -Inf and +Inf returns value x, and returns a new function that returns x when integrated between real(pole) and +Inf.\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.integrate","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.integrate","text":"Integration of Abstract Arbitrary FParallel\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.integrate-Union{Tuple{T}, Tuple{FParallelDiracDelta, T, Bool}, Tuple{FParallelDiracDelta, T, Bool, Tolerance}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.integrate","text":"The integrals over the distribution functions and associated integral kernel We avoid doing integrals involving the derivative of the Dirac delta function by doing integral by parts and knowing that f-> 0 and the integral limits\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.integrate-Union{Tuple{T}, Tuple{FParallelDiracDelta, T, LinearMaxwellVlasov.Pole, Bool}, Tuple{FParallelDiracDelta, T, LinearMaxwellVlasov.Pole, Bool, Tolerance}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.integrate","text":"Calculate the parallel integral with principal part and residue\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.integrate-Union{Tuple{T}, Tuple{FPerpendicularDiracDelta, T, Bool}, Tuple{FPerpendicularDiracDelta, T, Bool, Tolerance}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.integrate","text":"The integrals over the distribution functions and associated integral kernel We avoid doing integrals involving the derivative of the Dirac delta function by doing integral by parts and knowing that f-> 0 and the integral limits\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.integrate-Union{Tuple{T}, Tuple{LinearMaxwellVlasov.AbstractFParallelNumerical, T, LinearMaxwellVlasov.Pole, Bool}, Tuple{LinearMaxwellVlasov.AbstractFParallelNumerical, T, LinearMaxwellVlasov.Pole, Bool, Tolerance}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.integrate","text":"Integrate over the parallel distribution function multiplied by various kernels. If the imaginary part of the pole is zero, then do a trick that folds over the integral from the left of the pole to right. We need to take into account whether the real part of the pole is negative, otherwise positive slopes at negative velocities give rise to instability and not damping.\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.integrate-Union{Tuple{T}, Tuple{LinearMaxwellVlasov.AbstractFPerpendicular, T, Bool}, Tuple{LinearMaxwellVlasov.AbstractFPerpendicular, T, Bool, Tolerance}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.integrate","text":"Integration of Abstract Arbitrary FPerpendicular\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.limitsfolder-Union{Tuple{T}, Tuple{AbstractVector{T}, Any}} where T","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.limitsfolder","text":"Transform the limits of an integrand quadrature(foldnumeratoraboutpole(integrand, pole), limitsfolder(limits, pole)...)\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.normalise-Union{Tuple{T}, Tuple{T, Float64, Float64}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.normalise","text":"Tool to normalize a function f between two integral limits a and b\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.parallel-Tuple{LinearMaxwellVlasov.AbstractKineticSpecies, Configuration, Integer, Unsigned, Bool}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.parallel","text":"Interface The parallel integral of the distribution function and the relevant kernel\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.parallel-Union{Tuple{U}, Tuple{T}, Tuple{FBeam, T, U, Integer, Real, Unsigned, Bool}, Tuple{FBeam, T, U, Integer, Real, Unsigned, Bool, Tolerance}} where {T<:Number, U<:Number}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.parallel","text":"Parallel integrals for the Beam - used for testing purposes\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.parallel-Union{Tuple{U}, Tuple{T}, Tuple{LinearMaxwellVlasov.AbstractFParallelNumerical, T, U, Integer, Real, Unsigned, Bool}, Tuple{LinearMaxwellVlasov.AbstractFParallelNumerical, T, U, Integer, Real, Unsigned, Bool, Tolerance}} where {T<:Number, U<:Number}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.parallel","text":"Compile the kernels of the parallel integral and fetch it from the DistributionFunctions module. If k parallel is zero then there is no Landau damping, so this case is made separate, as is easier to deal with.\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.parallel_integral-Union{Tuple{S}, Tuple{S, Configuration}, Tuple{S, Configuration, Any}} where S","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.parallel_integral","text":"Memoise the inputs to the integrals; export the memoised function, and the dictionary of inputs -> outputs It only matters what n-m is, not the values of n and m individually So change (n, m)->(n-m, 0) to do fewer calculations\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.perpendicular-Tuple{LinearMaxwellVlasov.AbstractKineticSpecies, Configuration, Integer, Integer, Unsigned, Bool}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.perpendicular","text":"Interface The perpendicular integral of the distribution function and the relevant kernel\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.perpendicular-Union{Tuple{T}, Tuple{LinearMaxwellVlasov.AbstractFPerpendicular, Real, Pair{T, T}, Number, Unsigned, Bool}, Tuple{LinearMaxwellVlasov.AbstractFPerpendicular, Real, Pair{T, T}, Number, Unsigned, Bool, Tolerance}} where T<:Integer","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.perpendicular","text":"Calculate the definite integral kernels of the perpenedicular distribution function, or derivative thereof, multiplied by the kernel functions\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.perpendicular_integral","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.perpendicular_integral","text":"Memoise the inputs to the integrals; export the memoised function, and the dictionary of inputs -> outputs besselj harmonic numbers n and l can be any order Use this to do only half the calculatons (n, l)->(min(n, l), max(n, l)) Also know that besselj(-n, x) = (-1)^n * besselj(n, x)\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.plasma_dispersion_function-Union{Tuple{T}, Tuple{T, Unsigned}, Tuple{T, Unsigned, Any}} where T<:Number","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.plasma_dispersion_function","text":"Return the value of the plasma dispersion function This implementation includes the residue, which is easy to verify because Z(0) = im sqrt(π).\n\nx is the argument to the plasma disperion function\npower is the moment of the integral\n\n[1] S.D. Baalrud, Phys. Plasmas 20, 012118 (2013) and put ν = -Inf [2] M. Sampoorna et al., Generalized Voigt functions and their derivatives,   Journal of Quantitative Spectroscopy & Radiative Transfer (2006),   doi:10.1016/j.jqsrt.2006.08.011\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.principalpart-Union{Tuple{T}, Tuple{T, Number}, Tuple{T, Number, Real}, Tuple{T, Number, Real, Int64}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.principalpart","text":"Expressing f(x) = ∑ᵢ aᵢ (x - p)ⁱ find a₋₁\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.principalpartadaptive-Union{Tuple{T}, Tuple{T, Number}, Tuple{T, Number, Real}, Tuple{T, Number, Real, Int64}, Tuple{T, Number, Real, Int64, Tolerance}} where T<:Function","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.principalpartadaptive","text":"Expressing f(x) = ∑ᵢ aᵢ (x - p)ⁱ find a₋₁\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.tensor","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.tensor","text":"Calculate the tensor representing the linear Maxwell-Vlasov set of equations. The determinant is zero when the wavenumber and frequency represent a solution to the linear Maxwell-Vlasov system of equations for these species.\n\n\n\n\n\n","category":"function"},{"location":"#LinearMaxwellVlasov.transformaboutroots-Union{Tuple{T}, Tuple{T, Real}} where T","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.transformaboutroots","text":"Concertina and rescale a function in the sections between its roots such that all roots lie at -1 or +1. Integration over -1..1 of the resulting function gives the same answer as integration over original -Inf..Inf domain\n\n\n\n\n\n","category":"method"},{"location":"#LinearMaxwellVlasov.zerobetamagnetoacousticfrequency-Tuple{Number, Wavenumber, Number, Int64}","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.zerobetamagnetoacousticfrequency","text":"Eq. 25 of R. O. Dendy, C. N. Lashmore-Davies, and K. G. McClements, G. A. Cottrell, The excitation of obliquely propagating fast Alfven waves at fusion ion cyclotron harmonics, Phys. Plasmas 1 (6), June 1994\n\n\n\n\n\n","category":"method"},{"location":"#Index","page":"LinearMaxwellVlasov.jl Documentation","title":"Index","text":"","category":"section"},{"location":"","page":"LinearMaxwellVlasov.jl Documentation","title":"LinearMaxwellVlasov.jl Documentation","text":"","category":"page"}]
}

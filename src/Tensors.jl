using CommonSubexpressions, LinearAlgebra, MuladdMacro, StaticArrays

include("tensors/RelativisticMomentum.jl")

function contribution(species::CoupledRelativisticSpecies,
    config::Configuration, cyclotronharmonic::Int)
  return relativisticmomentum(species, config, cyclotronharmonic)
end

include("tensors/CoupledVelocity.jl")

"""
Calculate the unitless susceptibility tensor contribution for a
coupled velocity distribution function, without the need to sum over harmonics.
"""
function contribution(species::AbstractCoupledVelocitySpecies,
    config::Configuration, _::Cache=Cache())
  return coupledvelocity(species, config)
end

function contribution(species::AbstractSeparableVelocitySpecies,
    config::Configuration, cyclotronharmonic::Int,
    parallelintegral::U=parallel, perpendicularintegral::V=perpendicular
    ) where {U<:Function, V<:Function}

  ω = config.frequency
  kz, k⊥ = para(config.wavenumber), perp(config.wavenumber)

  # z /⊥ stands for parallel/perp integral
  # 2nd character q means that v_para^q in the integrand
  # F / T as third character in variable name stands for value / derivative
  (z0F, z1F, z2F, z0T, z1T) = parallelintegral(species, config, cyclotronharmonic)

  # ⊥ stands for perpendicular integral
  # F / T as third character in variable name stands for value / derivative
  # 2nd character q means that v_perp ^ q in the integrand
  # 4th-7th character (0, 1, _1) -> l + (0, 1, -1)
  ∫⊥s = perpendicularintegral(species, config, cyclotronharmonic)
  ⊥3F11, ⊥3F_1_1, ⊥3F1_1, ⊥2T11, ⊥2T_1_1, ⊥2T1_1, ⊥2F10, ⊥2F0_1, ⊥1F00, ⊥1T10, ⊥1T0_1 = ∫⊥s

  @cse @muladd begin
    m11 = (kz*(-(⊥2T11 + 2*⊥2T1_1 + ⊥2T_1_1)*z1F + (⊥3F11 + 2*⊥3F1_1 + ⊥3F_1_1)*z0T) + (⊥2T11 + 2*⊥2T1_1 + ⊥2T_1_1)*ω*z0F)/4
    m21 = im*(kz*((-⊥2T11 + ⊥2T_1_1)*z1F + (⊥3F11 - ⊥3F_1_1)*z0T) + (⊥2T11 - ⊥2T_1_1)*ω*z0F)/4
    m31 = (kz*(-(⊥1T0_1 + ⊥1T10)*z2F + (⊥2F0_1 + ⊥2F10)*z1T) + (⊥1T0_1 + ⊥1T10)*ω*z1F)/2
    m12 = -m21 # Onsager -im*(kz*((-⊥2T11 + ⊥2T_1_1)*z1F + (⊥3F11 - ⊥3F_1_1)*z0T) + (⊥2T11 - ⊥2T_1_1)*ω*z0F)/4
    m22 = (kz*(-(⊥2T11 - 2*⊥2T1_1 + ⊥2T_1_1)*z1F + (⊥3F11 - 2*⊥3F1_1 + ⊥3F_1_1)*z0T) + (⊥2T11 - 2*⊥2T1_1 + ⊥2T_1_1)*ω*z0F)/4
    m32 = -im*(kz*((⊥1T0_1 - ⊥1T10)*z2F + (⊥2F10 - ⊥2F0_1)*z1T) + (⊥1T10 - ⊥1T0_1)*ω*z1F)/2
    m13 = m31 # Onsager (k⊥*((⊥2T11 + 2*⊥2T1_1 + ⊥2T_1_1)*z1F + (-⊥3F11 - 2*⊥3F1_1 - ⊥3F_1_1)*z0T) + 2*(⊥2F0_1 + ⊥2F10)*ω*z0T)/4
    m23 = -m32 # Onsager im*(k⊥*(⊥2T11 - ⊥2T_1_1)*z1F + (k⊥*(⊥3F_1_1 - ⊥3F11) + 2*(⊥2F10 - ⊥2F0_1)*ω)*z0T)/4
    m33 = (k⊥*(⊥1T0_1 + ⊥1T10)*z2F + (-k⊥*(⊥2F0_1 + ⊥2F10) + 2*⊥1F00*ω)*z1T)/2
  end

  return @SArray [m11 m12 m13; m21 m22 m23; m31 m32 m33]
end

"""
Calculate the unitless susceptibility tensor for a maxwellian distribution
function given the bessel indices n.
"""
function contribution(species::T, config::Configuration, cyclotronharmonic::Int,
    parallelintegral::U=parallel, perpendicularintegral::V=perpendicular
    ) where {U<:Function, V<:Function,
             Tz, T⊥<:FPerpendicularMaxwellian,
             T<:AbstractSeparableVelocitySpecies{<:Tz, <:T⊥}}

  ω = config.frequency
  kz, k⊥ = para(config.wavenumber), perp(config.wavenumber)

  # z /⊥ stands for parallel/perp integral
  # 2nd character q means that v_para^q in the integrand
  # F / T as third character in variable name stands for value / derivative
  (z0F, z1F, z2F, z0T, z1T) = parallelintegral(species, config, cyclotronharmonic)

  # ⊥ stands for perpendicular integral
  # F / T as third character in variable name stands for value / derivative
  # 2nd character q means that v_perp ^ q in the integrand
  # 3rd (& 4th) charcter (0, 1, _1) -> l + (0, 1, -1)
  # final (& penultimate) character (0, 1, _1) -> n + (0, 1, -1)
  ∫⊥s = perpendicularintegral(species, config, cyclotronharmonic)
  ⊥2T11_2⊥2T1_1⊥2T_1_1, ⊥3F11_2⊥3F1_1⊥3F_1_1, ⊥2T_1_1_⊥2T11,
    ⊥3F_1_1_⊥3F11, ⊥1T0_1_⊥1T10, ⊥2F0_1_⊥2F10, ⊥1T0_1⊥1T10, ⊥2F0_1⊥2F10,
    ⊥2T112⊥2T1_1⊥2T_1_1, ⊥3F112⊥3F1_1⊥3F_1_1, ⊥1F00 = ∫⊥s

  @cse @muladd begin
    m11 = (kz*(-(⊥2T112⊥2T1_1⊥2T_1_1)*z1F + (⊥3F112⊥3F1_1⊥3F_1_1)*z0T) + (⊥2T112⊥2T1_1⊥2T_1_1)*ω*z0F)/4
    m21 = im*(kz*((⊥2T_1_1_⊥2T11)*z1F - (⊥3F_1_1_⊥3F11)*z0T) - (⊥2T_1_1_⊥2T11)*ω*z0F)/4
    m31 = (kz*(-(⊥1T0_1⊥1T10)*z2F + (⊥2F0_1⊥2F10)*z1T) + (⊥1T0_1⊥1T10)*ω*z1F)/2
    m12 = -m21 # Onsager -im*(kz*((⊥2T_1_1_⊥2T11)*z1F - (⊥3F_1_1_⊥3F11)*z0T) - (⊥2T_1_1_⊥2T11)*ω*z0F)/4
    m22 = (kz*(-(⊥2T11_2⊥2T1_1⊥2T_1_1)*z1F + (⊥3F11_2⊥3F1_1⊥3F_1_1)*z0T) + (⊥2T11_2⊥2T1_1⊥2T_1_1)*ω*z0F)/4
    m32 = -im*(kz*((⊥1T0_1_⊥1T10)*z2F - (⊥2F0_1_⊥2F10)*z1T) - (⊥1T0_1_⊥1T10)*ω*z1F)/2
    m13 = m31 # Onsager (k⊥*((⊥2T112⊥2T1_1⊥2T_1_1)*z1F - (⊥3F112⊥3F1_1⊥3F_1_1)*z0T) + 2*(⊥2F0_1⊥2F10)*ω*z0T)/4
    m23 = -m32 # Onsager im*(-k⊥*(⊥2T_1_1_⊥2T11)*z1F + (k⊥*(⊥3F_1_1_⊥3F11) - 2*(⊥2F0_1_⊥2F10)*ω)*z0T)/4
    m33 = (k⊥*(⊥1T0_1⊥1T10)*z2F + (-k⊥*(⊥2F0_1⊥2F10) + 2*⊥1F00*ω)*z1T)/2
  end

  return @SArray [m11 m12 m13; m21 m22 m23; m31 m32 m33]
end

"""
Calculate the unitless susceptibility tensor for a warm plasma species
Swanson 3.63
"""
function contribution(species::WarmSpecies, config::Configuration,
    _::Cache=Cache())
  ω = config.frequency
  Ω = species.Ω
  kz, k⊥ = para(config.wavenumber), perp(config.wavenumber)
  θ = angle(config.wavenumber)
  k²Cₛ²_ω² = abs2(config.wavenumber) * species.soundspeed^2 / ω^2

  @cse @muladd begin
    denom = (1 - k²Cₛ²_ω²) - (Ω / ω)^2 * (1 - k²Cₛ²_ω² * cos(θ)^2)

    m11 = - (1 - k²Cₛ²_ω² * cos(θ)^2) / denom
    m21 = im * Ω / ω * (1 - k²Cₛ²_ω² * cos(θ)^2) / denom
    m31 = - k²Cₛ²_ω² * cos(θ) * sin(θ) / denom
    m12 = -m21
    m22 = - (1 - k²Cₛ²_ω²) / denom
    m32 = im * Ω / ω * k²Cₛ²_ω² * cos(θ) * sin(θ) / denom
    m13 = m31
    m23 = -m32
    m33 = - (1 - (Ω / ω)^2 - k²Cₛ²_ω² * sin(θ)^2) / denom
  end

  return @SArray [m11 m12 m13; m21 m22 m23; m31 m32 m33]
end

"""
Calculate the unitless susceptibility tensor for a cold plasma species
"""
function contribution(species::ColdSpecies, config::Configuration, _=missing)
  return contribution(WarmSpecies(species), config)
end

"""
Calculate the converged unitless susceptibility tensor contribution for a
maxwellian distribution function, having summed over bessel indices n.
"""
function contribution(species::AbstractSeparableVelocitySpecies,
    config::Configuration, cache::Cache=Cache())

  ∫para = parallel_integral(species, config, cache.parallel)
  ∫perp = perpendicular_integral(species, config, cache.perpendicular)

  f = n -> contribution(species, config, n, ∫para, ∫perp)

  return converge(f, config.options.summation_tol)
end

"""
Calculate the converged unitless susceptibility tensor contribution for a
maxwellian distribution function, having summed over bessel indices n.
"""
function contribution(species::AbstractKineticSpecies, config::Configuration,
    _::Cache=Cache())

  f = n -> contribution(species, config, n)

  return converge(f, config.options.summation_tol)
end

"""
The contribution to the dielectric tensor for a given species
"""
function dielectriccontribution(species::AbstractSpecies, config::Configuration,
    cache::Cache=Cache())
  return contribution(species, config, cache) * (species.Π / config.frequency)^2
end

"""
The conducivity tensor for a given species
"""
function conductivity(species::AbstractSpecies, config::Configuration,
    cache::Cache=Cache())
  ϵᵢⱼₛ = dielectriccontribution(species, config, cache)
  return -im * config.frequency * ϵ₀ * ϵᵢⱼₛ
end

"""
The dielectric tensor for a given plasma
"""
function dielectric(plasma::AbstractPlasma, config::Configuration,
    cache::Cache=Cache())
  ϵᵢⱼₛ = species -> dielectriccontribution(species, config, cache)
  return mapreduce(ϵᵢⱼₛ, +, plasma) + I
end

"""
Calculate the tensor representing the linear Maxwell-Vlasov set of equations.
The determinant is zero when the wavenumber and frequency represent a solution
to the linear Maxwell-Vlasov system of equations for these species.
"""
function tensor(plasma::AbstractPlasma, config::Configuration,
    cache::Cache=Cache())
  ϵᵢⱼ = dielectric(plasma, config, cache)
  return ϵᵢⱼ + kkT_Ik²(config.wavenumber) * (c₀ / config.frequency)^2
end

"""
The electrostatic dielectric tensor for a given plasma, a zero valued
determinant of which represents a solution to the linear
poisson-vlasov system of equations
"""
function electrostatictensor(plasma::AbstractPlasma, config::Configuration,
    cache::Cache=Cache())
  return dielectric(plasma, config, cache)
end


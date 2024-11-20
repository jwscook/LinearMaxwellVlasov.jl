
"""
Calculate the Alfven speed given Π_Ωs, which is an iterable container
of the ratio of the plasma to the cyclotron frequency of all the species
of the plasma
"""
alfvenspeed(Π_Ωs) = c₀ / sqrt(sum(Π_Ωs.^2))

"""
    plasmafrequency(n::Number,mass::Number,Z::Number=1)

...
# Arguments
- `n::Number`: species number density in m^-3
- `mass::Number`: particle mass in kg
- `Z::Number=1`: Charge number
...
"""
function plasmafrequency(n::Number, mass::Number, Z::Number=1)
  return sqrt((Z * q₀)^2 * n / mass / ϵ₀)
end
"""
    cyclotronfrequency(B::Number,mass::Number,Z::Number=1)

...
# Arguments
- `B::Number`: magnetic field in T
- `mass::Number`:  mass in kg
- `Z::Number=1`: charge number
...

"""
function cyclotronfrequency(B::Number, mass::Number, Z::Number=1)
  return Z * q₀ * B / mass
end

"""
  thermalspeed(ϵ_V::Number, mass::Number)

Thermal speed defined as `sqrt(2 * q₀ * ϵ_V / mass)`
'''
# Arguments
-  `ϵ_V::Number`: energy in electron Volts
- `mass::Number`: particle mass in kg
'''
"""
thermalspeed(ϵ_V::Number, mass::Number) = sqrt(2 * q₀ * ϵ_V / mass)

"""
  thermalmomentum(ϵ_V::Number, mass::Number)

Thermal momentum is defined as `sqrt(2 * q₀ * ϵ_V * mass)`

'''
# Arguments
-  `ϵ_V::Number`: energy in electron Volts
- `mass::Number`: particle mass in kg
'''
"""

thermalmomentum(ϵ_V::Number, mass::Number) = thermalspeed(ϵ_V, mass) * mass


"""
Eq. 25 of R. O. Dendy, C. N. Lashmore-Davies, and K. G. McClements,
G. A. Cottrell, The excitation of obliquely propagating fast Alfven waves at
fusion ion cyclotron harmonics, Phys. Plasmas 1 (6), June 1994
"""
function zerobetamagnetoacousticfrequency(Va::Number, K::Wavenumber,
    Ωi::Number, slowfast::Integer)
  kz = parallel(K)
  k⊥ = perpendicular(K)
  k² = kz^2 + k⊥^2
  X = k² + kz^2 + k² * (kz * Va / Ωi)^2
  return sqrt(Va^2 / 2 * (X + slowfast * sqrt(X^2 - 4 * k² * kz^2)))
end
function slowzerobetamagnetoacousticfrequency(Va::Number,
    K::Wavenumber, Ωi::Number)
  return zerobetamagnetoacousticfrequency(Va, K, Ωi, -1)
end
function fastzerobetamagnetoacousticfrequency(Va::Number,
    K::Wavenumber, Ωi::Number)
  return zerobetamagnetoacousticfrequency(Va, K, Ωi, 1)
end

function magnetoacousticfrequency(Va::Number, Cs::Number,
    K::Wavenumber, slowfast::Integer)
  return abs(K) * magnetoacousticspeed(Va, Cs, K, slowfast)
end
function slowmagnetoacousticfrequency(Va::Number, Cs::Number,
    K::Wavenumber)
  return abs(K) * slowmagnetoacousticspeed(Va, Cs, K)
end
function fastmagnetoacousticfrequency(Va::Number, Cs::Number,
    K::Wavenumber)
  return abs(K) * fastmagnetoacousticspeed(Va, Cs, K)
end

function shearfrequency(Va::Number, K::Wavenumber)
  return abs(parallel(K)) * Va # === abs(K) * shearspeed(Va, K)
end

function magnetoacousticspeed(Va::Number, Cs::Number,
    K::Wavenumber, slowfast::Integer)
  @assert abs(slowfast) == 1
  a = 1
  b = - (Va^2 + Cs^2)
  c = (Va * Cs * cos(angle(K)))^2
  return sqrt((-b + slowfast * sqrt(b^2 - 4 * a * c)) / 2 / a)
end
function slowmagnetoacousticspeed(Va::Number, Cs::Number,
    K::Wavenumber)
  return magnetoacousticspeed(Va, Cs, K, -1)
end
function fastmagnetoacousticspeed(Va::Number, Cs::Number,
    K::Wavenumber)
  return magnetoacousticspeed(Va, Cs, K, 1)
end
function shearspeed(Va::Number, K::Wavenumber)
  return Va * abs(cos(angle(K)))
end

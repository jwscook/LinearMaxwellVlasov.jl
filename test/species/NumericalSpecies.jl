

"""
Create a kinetic plasma species
 - `vthb` parallel thermal speed
 - `vth⊥` perpedicular thermal speed
 - `vdb` parallel drift speed
 - `vd⊥` perpedicular ring speed
"""
function NumericalSpecies(Π::Float64, Ω::Float64,
    vthb::Float64, vth⊥::Float64=vthb, vdb::Float64=0.0, vd⊥::Float64=0.0)
  Fz = FParallelNumerical(vthb, vdb)
  F⊥ = FPerpendicularNumerical(vth⊥, vd⊥)
  return SeparableVelocitySpecies(Π, Ω, Fz, F⊥)
end


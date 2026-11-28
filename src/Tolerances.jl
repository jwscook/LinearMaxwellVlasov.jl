struct Tolerance{T}
  abs::T
  rel::T
  _uniqueid::UInt64
  function Tolerance(atol::T, rtol::T) where {T<:Number}
    _uniqueid = hash((atol, rtol), hash(:Tolerance))
    return new{T}(atol, rtol, _uniqueid)
  end
end
Tolerance(;atol::Number=0.0, rtol::Number=eps()^(3//4)) = Tolerance(atol, rtol)
Tolerance(::Type{T}) where {T} = Tolerance(zero(T), eps(real(T))^(3//4))
Tolerance{T}(a::T, b::T) where {T} = Tolerance(a, b)
uniqueid(t::Tolerance) = t._uniqueid

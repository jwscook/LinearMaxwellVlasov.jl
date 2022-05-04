struct Tolerance{T}
  rel::T
  abs::T
  _uniqueid::UInt64
  function Tolerance(rtol::T, atol::T) where {T<:Number}
    _uniqueid = hash((rtol, atol), hash(:Tolerance))
    return new{T}(rtol, atol, _uniqueid)
  end
end
Tolerance(;rtol::Number=sqrt(eps()), atol::Number=0.0) = Tolerance(rtol, atol)
Tolerance(::Type{T}) where {T} = Tolerance(sqrt(eps(real(T))), zero(T))
Tolerance{T}(a::T, b::T) where {T} = Tolerance(a, b)
Tolerance(n::Number) = Tolerance(typeof(n))
uniqueid(t::Tolerance) = t._uniqueid

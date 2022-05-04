
struct Options{T<:Number}
  quadrature_tol::Tolerance{T}
  summation_tol::Tolerance{T}
  memoiseparallel::Bool
  memoiseperpendicular::Bool
  _uniqueid::UInt64
  function Options(quadrature_tol::Tolerance{T},
      summation_tol::Tolerance{T},
      memoiseparallel::Bool, memoiseperpendicular::Bool) where {T<:Number}
    _uniqueid = hash((quadrature_tol, summation_tol,
      memoiseparallel, memoiseperpendicular), hash(:Options))
    return new{T}(quadrature_tol, summation_tol,
      memoiseparallel, memoiseperpendicular, _uniqueid)
  end
end
uniqueid(o::Options) = o._uniqueid


function defaults(::Type{T}=Float64) where {T}
  output = Dict{Symbol, Any}()
  output[:quadrature_rtol] = eps(T)^(3//4)
  output[:summation_rtol] = eps(T)^(3//4)
  output[:quadrature_atol] = zero(T)
  output[:summation_atol] = zero(T)
  output[:memoiseparallel] = true
  output[:memoiseperpendicular] = true
  return output
end

function Options(::Type{T}=Float64; kwargstuple...) where {T}
  kwargs = Dict(kwargstuple)
  args = defaults(T)
  if haskey(kwargs, :tols)
    return Options(T; atols=kwargs[:tols].abs, rtols=kwargs[:tols].rel)
  end
  specialcasekeys = Dict{Symbol, Any}(
    :quad_atol => :quadrature_atol, :quad_rtol => :quadrature_rtol,
    :sum_atol => :summation_atol, :sum_rtol => :summation_rtol,
    :atols => (:quadrature_atol, :summation_atol),
    :rtols => (:quadrature_rtol, :summation_rtol),)

  for (kwargkey, argkeys) ∈ specialcasekeys
    haskey(kwargs, kwargkey) || continue
    (argkeys isa Tuple) || (argkeys = (argkeys,))
    for argkey ∈ argkeys
      args[argkey] = kwargs[kwargkey]
    end
  end

  for key in filter(!in(keys(specialcasekeys)), keys(kwargs))
    @assert haskey(args, key) "key = $key"
    args[key] = kwargs[key]
  end
  quadrature_tol = Tolerance{T}(args[:quadrature_rtol], args[:quadrature_atol])
  summation_tol = Tolerance{T}(args[:summation_rtol], args[:summation_atol])
  return Options(quadrature_tol, summation_tol,
    args[:memoiseparallel], args[:memoiseperpendicular])
end

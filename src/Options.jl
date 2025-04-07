
struct Options{T<:Number}
  quadrature_tol::Tolerance{T}
  summation_tol::Tolerance{T}
  memoiseparallel::Bool
  memoiseperpendicular::Bool
  cubature_maxevals::Int
  _uniqueid::UInt64
  function Options(quadrature_tol::Tolerance{T},
      summation_tol::Tolerance{T},
      memoiseparallel::Bool, memoiseperpendicular::Bool,
      cubature_maxevals::Int) where {T<:Number}
    _uniqueid = hash((quadrature_tol, summation_tol,
      memoiseparallel, memoiseperpendicular, cubature_maxevals), hash(:Options))
    return new{T}(quadrature_tol, summation_tol,
      memoiseparallel, memoiseperpendicular,
      cubature_maxevals, _uniqueid)
  end
end
uniqueid(o::Options) = o._uniqueid


"""
    defaults(::Type{T}=Float64)where{T}

The default settings for the Options object.

...
# Arguments
- `::Type{T}=Float64where{T}`: the number type to express the tolerances in
...
"""
function defaults(::Type{T}=Float64) where {T}
  output = Dict{Symbol, Any}()
  output[:quadrature_rtol] = eps(T)^(3//4)
  output[:summation_rtol] = eps(T)^(3//4)
  output[:quadrature_atol] = zero(T)
  output[:summation_atol] = zero(T)
  output[:memoiseparallel] = true
  output[:memoiseperpendicular] = true
  output[:cubature_maxevals] = typemax(Int)
  return output
end

"""
    Options(::Type{T}=Float64;kwargstuple...)where{T}

Create an Options object a set of kwargs. Best to read the code to see what counts
as a duplicate. Copying that logic here will only create broken docs in the future.

...
# Arguments
- `::Type{T}=Float64`: the number type to express the tolerances in
Optional args:
- `kwargstuple...`: the keyword arguments that are preferred over the defaults
...

"""
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
    :rtols => (:quadrature_rtol, :summation_rtol),
    :maxevals => :cubature_maxevals,
    :cuba_evals => :cubature_maxevals,
   )

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
    args[:memoiseparallel], args[:memoiseperpendicular],
    args[:cubature_maxevals])
end

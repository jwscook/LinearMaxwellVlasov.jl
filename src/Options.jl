
struct Options{T<:Number}
  quadrature_tol::Tolerance{T}
  cubature_tol::Tolerance{T}
  summation_tol::Tolerance{T}
  memoiseparallel::Bool
  memoiseperpendicular::Bool
  cubature_maxevals::Int
  residue_maxevals::Int
  erroruponcubaturenonconvergence::Bool
  cauchydeformationangle::Float64
  _uniqueid::UInt64
  function Options(::Type{T}=Float64;
      tols::Tolerance{T} = Tolerance(T),
      atols = tols.abs,
      rtols = tols.rel,
      quad_atol=atols,
      quad_rtol=rtols,
      quadrature_atol = quad_atol,
      quadrature_rtol = quad_rtol,
      quadrature_tol::Tolerance = Tolerance(quadrature_atol, quadrature_rtol),
      cuba_atol = atols,
      cuba_rtol = rtols,
      cubature_atol = cuba_atol,
      cubature_rtol = cuba_rtol,
      cubature_tol::Tolerance = Tolerance(cubature_atol, cubature_rtol),
      sum_atol = atols,
      sum_rtol = rtols,
      summation_atol = sum_atol,
      summation_rtol = sum_rtol,
      summation_tol::Tolerance = Tolerance(summation_atol, summation_rtol),
      memoiseparallel::Bool = false,
      memoiseperpendicular::Bool = false,
      cuba_evals::Int = typemax(Int),
      cubature_maxevals::Int = cuba_evals,
      res_evals::Int = typemax(Int),
      residue_maxevals::Int = res_evals,
      error_cuba::Bool = true,
      erroruponcubaturenonconvergence::Bool = error_cuba,
      deformationangle::Float64 = DEFAULT_CAUCHY_DEFORMATION_ANGLE,
      cauchydeformationangle::Float64 = deformationangle
    ) where {T<:Real}

    _uniqueid = hash((quadrature_tol, cubature_tol, summation_tol,
      memoiseparallel, memoiseperpendicular, cubature_maxevals, residue_maxevals,
      erroruponcubaturenonconvergence, cauchydeformationangle),
      hash(:Options))

    return new{T}(quadrature_tol, cubature_tol, summation_tol,
      memoiseparallel, memoiseperpendicular, cubature_maxevals, residue_maxevals,
      erroruponcubaturenonconvergence, cauchydeformationangle, _uniqueid)
  end
end
uniqueid(o::Options) = o._uniqueid


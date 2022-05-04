using LinearAlgebra, StaticArrays

"""
Wavenumber decomposed into parallel and perpendicular components
Construct with parallel and wavenumber components, or by keyword arguement pairs
- parallel & perpendicular
- wavenumber & propagationangle
"""
struct Wavenumber{T<:Number, U<:Number}
  parallel::T
  perpendicular::U
end
function Wavenumber(;kwargstuple...)
  kwargs = Dict(kwargstuple)
  haskey(kwargs, :k) && (kwargs[:wavenumber] = kwargs[:k])
  haskey(kwargs, :θ) && (kwargs[:propagationangle] = kwargs[:θ])
  haskey(kwargs, :kz) && (kwargs[:parallel] = kwargs[:kz])
  haskey(kwargs, :k⊥) && (kwargs[:perpendicular] = kwargs[:k⊥])
  if haskey(kwargs, :parallel) && haskey(kwargs, :perpendicular)
    any(k ∈ keys(kwargs) for k ∈ (:wavenumber, :propagationangle)) && err()
    parallel = kwargs[:parallel]
    perpendicular = kwargs[:perpendicular]
    return Wavenumber(parallel, perpendicular)
  elseif haskey(kwargs, :wavenumber) && haskey(kwargs, :propagationangle)
    any(k ∈ keys(kwargs) for k ∈ (:parallel, :perpendicular)) && err()
    wavenumber = kwargs[:wavenumber]
    propagationangle = kwargs[:propagationangle]
    parallel = wavenumber * cospi(propagationangle / π)
    perpendicular = wavenumber * sinpi(propagationangle / π)
    return Wavenumber(parallel, perpendicular)
  elseif !isempty(kwargs)
    throw(ArgumentError("Invalid Wavenumber construction, kwargs=$kwargstuple"))
  end
  return Wavenumber(0.0, 0.0)
end

@inline parallel(K::Wavenumber) = K.parallel
@inline perpendicular(K::Wavenumber) = K.perpendicular

@inline propagationangle(K::Wavenumber) = atan(K.perpendicular, K.parallel)

Base.angle(K::Wavenumber) = propagationangle(K)


@inline para(K::Wavenumber) = parallel(K)
@inline perp(K::Wavenumber) = perpendicular(K)

@inline cartesian_vector(K::Wavenumber) = @SArray [perp(K), 0.0, para(K)]

@inline function curl_curl(K::Wavenumber)
  k = cartesian_vector(K)
  return k*k' - dot(k, k) * I
end

Base.abs(K::Wavenumber) = sqrt(abs2(K))
Base.abs2(K::Wavenumber) = abs2(para(K)) + abs2(perp(K))
Base.iszero(K::Wavenumber) = iszero(para(K)) && iszero(perp(K))
function Base.:-(a::Wavenumber, b::Wavenumber)
  return Wavenumber(parallel=para(a) - para(b), perpendicular=perp(a) - perp(b))
end

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

"""
julia> using Symbolics
julia> @syms kx ky kz;
julia> k = [kx, ky, kz];
julia> curl = im * [0 -k[3] k[2]; k[3] 0 -k[1]; -k[2] k[1] 0];
julia> curl * curl
3×3 Matrix{Any}:
 ky^2 + kz^2  -kx*ky       -kx*kz
 -kx*ky       kx^2 + kz^2  -ky*kz
 -kx*kz       -ky*kz       kx^2 + ky^2

"""
@inline function curlcurl(K::Wavenumber)
  kx = perp(K) # ky = 0
  kz = para(K)
  return @SArray [ kz^2    0           -kx * kz;
                   0       kx^2 + kz^2  0      ;
                  -kx * kz 0            kx^2   ]
end

Base.abs(K::Wavenumber) = sqrt(abs2(K))
Base.abs2(K::Wavenumber) = abs2(para(K)) + abs2(perp(K))
Base.iszero(K::Wavenumber) = iszero(para(K)) && iszero(perp(K))
function Base.:-(a::Wavenumber, b::Wavenumber)
  return Wavenumber(parallel=para(a) - para(b), perpendicular=perp(a) - perp(b))
end

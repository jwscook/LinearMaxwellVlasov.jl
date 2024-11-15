
"""
    fastisapprox(n,x,y,atol,rtol,nans)



...
# Arguments
- `n`: normalisation
- `x`: compare this
- `y`: with that
- `atol`: absolute tolerance
- `rtol`: relative tolerance
- `nans`: whether to treat NaNs as equal
...

# Example
```julia
```
"""
function fastisapprox(n, x, y, atol, rtol, nans)
  (n <= max(abs2(atol), abs2(rtol) * max(x, y))) && return true
  return nans ? (isnan(x) && isnan(y)) : false
end

function fastisapprox(A::Number, B::Number; atol, rtol, nans)
  n = abs2(A - B)
  x, y = abs2(A), abs2(B)
  return fastisapprox(n, x, y, atol, rtol, nans)
end

function fastisapprox(A::AbstractArray{T, N}, B::AbstractArray{T, N};
    atol, rtol, nans) where {T, N}
  n = x = y = zero(real(T))
  @inbounds @simd for i in eachindex(A, B)
    a, b = A[i], B[i]
    n += abs2(a - b)
    x += abs2(a)
    y += abs2(b)
  end
  return fastisapprox(n, x, y, atol, rtol, nans)
end

"""
    converge(f::T,tol::Tolerance=Tolerance()) where {T}

Sum the output of f starting with argument 0 and (1,-1), (2,-2)... etc
until the convergence criterion is met based on the tolerance

...
# Arguments
- `f::T`: input function or functor
- `tol::Tolerance=Tolerance`: Tolerance object containing relative and absolute tolerances
...

# Example
```julia
```
"""
@inline function converge(f::T, tol::Tolerance=Tolerance()) where {T}
  l = u = n₀ = 0
  value = f(n₀) # f(::Int) is summed outwards from 0 to ± Inf until convergence
  delta = f(l -= 1) + f(u += 1)
  while !fastisapprox(value, value + delta, rtol=tol.rel, atol=tol.abs, nans=true)
    value += delta
    delta = f(l -= 1) + f(u += 1)
  end
  return value + delta
end

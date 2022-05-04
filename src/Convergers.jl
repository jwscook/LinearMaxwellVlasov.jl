
@inline function converge(f::T, tol::Tolerance=Tolerance()) where {T<:Function}
  l = u = n₀ = 0
  value = f(n₀) # f(::Int) is summed outwards from 0 to ± Inf until convergence
  delta = f(l -= 1) + f(u += 1)
  while !isapprox(value, value + delta, rtol=tol.rel, atol=tol.abs, nans=true)
    value += delta
    delta = f(l -= 1) + f(u += 1)
  end
  return value + delta
end

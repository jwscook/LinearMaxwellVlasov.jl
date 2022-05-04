using PyCall

sy = pyimport("sympy")

const simplify = sy.simplify
const symbols = sy.symbols
const SymFunction = sy.Function
coeff(s, x) = s.coeff(x)
subs(s, d) = s.subs(d.first, d.second)
cos = sy.cos
sin = sy.sin
diff = sy.diff
expand = sy.expand
expand_trig = sy.expand_trig
integrate = sy.integrate
exp = sy.exp
factor = sy.factor
iszero(m::PyObject) = m == 0

pie = sy.pi

gyrobunch = false # true # false
kriszero = false # r is for perpendicular
kziszero = false # z is for parallel
@assert !(kriszero && kziszero)
@assert gyrobunch "gyrobunching invalidates this model!!! - don't use it!"

pretty(x) = (display(x); println(" "); println(" "); println(" "); println(" "))
∫(a, b) = simplify(integrate(expand(a), b))

qs, ms, ns, wps = symbols("qs, ms, ns, wps", real=true)
cm, a0 = symbols("cm, a0", real=true)
l, m, n = symbols("l, m, n", integer=true)
w, W, kz, kr = symbols("w, W, kz, kr")
vz, vr, a = symbols("vz, vr, a", real=true)
Ex, Ey, Ez = symbols("Ex, Ey, Ez")
fz = SymFunction("fz")(vz)
fr = SymFunction("fr")(vr)
a1, ac, as = symbols("a1, ac, as")

fa = symbols("fa")
if gyrobunch
  fa = cm * cos(m*a + a0)
  fa = expand_trig(fa)
else
  m = subs(m, m => 0) # sub in as integer 0
  a0 = subs(a0, a0 => 0) # sub in as integer 0
  fa = 1//2 / pie
end

z = vr * kr / W
bj = SymFunction("J")

f0 = fz * fr * fa
f1 = SymFunction("f1")(vz, vr, a)

Cx = vr*cos(a)*kz*diff(f0, vz) + (w - kz*vz)*cos(a)*diff(f0, vr) - (w - kz*vz)/vr*sin(a)*diff(f0, a)
Cy = vr*sin(a)*kz*diff(f0, vz) + (w - kz*vz)*sin(a)*diff(f0, vr) + ((w - kz*vz)/vr*cos(a) - kr)*diff(f0, a)
Cz = (w - kr*vr*cos(a))*diff(f0, vz) + vz*kr*cos(a)*diff(f0, vr) - kr*vz/vr*sin(a)*diff(f0, a)

A = im * (w - vr*kr*cos(a) - vz*kz) / W
B = qs/(ms*w*W) * (Cx * Ex + Cy * Ey + Cz * Ez)

function substitutetrigprodtosum(B)
  iszero(m) && return B
  trigsimp(x) = expand_trig(simplify(x))
  @assert 0 == trigsimp(2*cos(m*a)*cos(a) - (cos((m-1)*a)+cos((m+1)*a)))
  @assert 0 == trigsimp(2*cos(m*a)*sin(a) - (sin((m+1)*a)-sin((m-1)*a)))
  @assert 0 == trigsimp(2*sin(m*a)*sin(a) - (cos((m-1)*a)-cos((m+1)*a)))
  @assert 0 == trigsimp(2*sin(m*a)*cos(a) - (sin((m+1)*a)+sin((m-1)*a)))
  cosmacosa = (cos((m - 1)*a) + cos((m + 1)*a)) // 2
  cosmasina = (sin((m + 1)*a) - sin((m - 1)*a)) // 2
  sinmasina = (cos((m - 1)*a) - cos((m + 1)*a)) // 2
  sinmacosa = (sin((m + 1)*a) + sin((m - 1)*a)) // 2
  B = expand(B)
  B = subs(B, cos(m*a)*cos(a) => cosmacosa)
  B = subs(B, cos(m*a)*sin(a) => cosmasina)
  B = subs(B, sin(m*a)*sin(a) => sinmasina)
  B = subs(B, sin(m*a)*cos(a) => sinmacosa)
  return B
end

B = substitutetrigprodtosum(B)

function substitutetrigforbesseljs(Bexp∫Ada)
  s = -1
  exp_ina = exp(s * im*n*a)
  exp_izsina = bj(n, z) * exp_ina

  cosaexp_izsina =       1//2 * (bj(n-1, z) + bj(n+1, z)) * exp_ina
  sinaexp_izsina = -s * im//2 * (bj(n-1, z) - bj(n+1, z)) * exp_ina

  cosmaexp_izsina =       1//2 * (bj(n-m, z) + bj(n+m, z)) * exp_ina
  sinmaexp_izsina = -s * im//2 * (bj(n-m, z) - bj(n+m, z)) * exp_ina

  cosmplus1aexp_izsina  =       1//2 * (bj(n-m-1, z) + bj(n+m+1, z)) * exp_ina
  cosmminus1aexp_izsina =       1//2 * (bj(n-m+1, z) + bj(n+m-1, z)) * exp_ina
  sinmplus1aexp_izsina  = -s * im//2 * (bj(n-m-1, z) - bj(n+m+1, z)) * exp_ina
  sinmminus1aexp_izsina = -s * im//2 * (bj(n-m+1, z) - bj(n+m-1, z)) * exp_ina

  Bexp∫Ada = expand(Bexp∫Ada)
  Bexp∫Ada = subs(Bexp∫Ada, cos(a*m + a)*exp(-im*z*sin(a)) => cosmplus1aexp_izsina)
  Bexp∫Ada = subs(Bexp∫Ada, cos(a*m - a)*exp(-im*z*sin(a)) => cosmminus1aexp_izsina)
  Bexp∫Ada = subs(Bexp∫Ada, sin(a*m + a)*exp(-im*z*sin(a)) => sinmplus1aexp_izsina)
  Bexp∫Ada = subs(Bexp∫Ada, sin(a*m - a)*exp(-im*z*sin(a)) => sinmminus1aexp_izsina)
  if m != 0
    Bexp∫Ada = subs(Bexp∫Ada, cos(a*m)*exp(-im*z*sin(a)) => cosmaexp_izsina)
    Bexp∫Ada = subs(Bexp∫Ada, sin(a*m)*exp(-im*z*sin(a)) => sinmaexp_izsina)
  end
  Bexp∫Ada = subs(Bexp∫Ada, exp(-im*z*sin(a)) => exp_izsina)
  Bexp∫Ada = expand_trig(simplify(Bexp∫Ada))
  return Bexp∫Ada
end

exp∫Ada = exp(∫(A, a))
Bexp∫Ada = B * exp∫Ada
Bexp∫Ada = substitutetrigforbesseljs(Bexp∫Ada)
∫Bexp∫Adada = - im * W * Bexp∫Ada / (w - vz*kz - n*W)

exp_∫Ada = exp(∫(-A, a))
# substitution by dividing first, and multiplying second
exp_∫Adal = exp_∫Ada / exp(im * z * sin(a)) * bj(l, z) * exp(im * l * a)
exp_∫Adal = expand_trig(simplify(exp_∫Adal))

f1 = exp_∫Adal * ∫Bexp∫Adada
f1 = factor(f1)

function shorthand(f1)
  @show f1
  f1copy = deepcopy(f1)
  f1copy = expand_trig(simplify(f1copy))
  dfzdvz = symbols("dfzdvz")
  f1copy = subs(f1copy, diff(fz, vz) => dfzdvz)

  maxcoeffvr, maxcoeffvz = 0, 0

  f1copy *= (w - kz*vz - n*W) # absorbed into the parallel shorthand variables
  f1copy = expand_trig(simplify(f1copy))
  f1copy = expand(f1copy)
  for i ∈ reverse(0:5), d ∈ ("F", "T")
    s = symbols("z" * string(i) * d)
    iszero(coeff(f1copy, vz^i)) || (maxcoeffvz = max(i, maxcoeffvz))
    if d == "F"
      f1copy = subs(f1copy, vz^i * fz => s)
    else
      f1copy = subs(f1copy, vz^i * dfzdvz => s)
    end
  end

  f1copy /= bj(l, z) # going to substitute for this anyway for perpendicular part
  f1copy = expand(f1copy)

  f1copy /= exp(-im*n*a) # get rid of one part of substitution beforehand
  if !gyrobunch
    f1copy /= exp(im*l*a)
    f1copy = simplify(expand_trig(simplify(f1copy)))
    f1copy = expand_trig(f1copy)
    #l = n
    if iszero(coeff(f1copy, cos(a))) && iszero(coeff(f1copy, sin(a)))
      f1copy *= 2*pie*bj(n, z)
    else
      f1copy = subs(f1copy, cos(a) => pie * (bj(n-1, z) + bj(n+1, z)))
      f1copy = subs(f1copy, sin(a) => -im*pie * (bj(n-1, z) - bj(n+1, z)))
    end
    @assert iszero(coeff(f1copy, cos(a)))
    @assert iszero(coeff(f1copy, sin(a)))
    @assert iszero(coeff(f1copy, exp(im*l*a)))
    # z = vr * kr / W
    #f1copy = subs(f1copy, vr * bj(n-1, z) => W / kr * 2*n*bj(n, z) - vr * bj(n+1, z))
  else
    f1copy = expand_trig(simplify(f1copy))
    f1copy = subs(f1copy, exp(im*l*a)*cos(a) => ac)
    f1copy = subs(f1copy, exp(im*l*a)*sin(a) => as)
    f1copy = subs(f1copy, exp(im*l*a) => a1)
  end

  dfrdvr = symbols("dfrdvr")
  f1copy = subs(f1copy, diff(fr, vr) => dfrdvr)
  f1copy = expand(f1copy)
  δs = gyrobunch ? 0 : -1:1
  for i ∈ reverse(0:5), d ∈ ("F", "T"), M ∈ (0, -1, 1),  Δ = -1:1, δ ∈ δs
    mstr = M == -1 ? "_m" : M == 0 ? "0m" : "m"
    nstr = Δ == -1 ? "_1" : Δ == 0 ? "0" : "1"
    lstr = δ == -1 ? "_1" : δ == 0 ? "0" : "1"
    gyrobunch && (lstr = "") # δ is always 0 if gyrobunched
    gyrobunch || (mstr = "") # M is always 0 if not gyrobunched
    s = symbols("r" * string(i) * d  * lstr * mstr * nstr)
    iszero(coeff(f1copy, vr^i)) || (maxcoeffvr = max(i, maxcoeffvr))
    # if gyrobunching then there is a bj(l, z); goes without saying
    bjbj = gyrobunch ? bj(m*M + n + Δ, z) : bj(n + Δ, z) * bj(n + δ, z)
    f1copy = expand(f1copy)
    if d == "F"
      f1copy = subs(f1copy, vr^i * bjbj * fr => s)
    else
      f1copy = subs(f1copy, vr^i * bjbj * dfrdvr => s)
    end
  end

  if !iszero(coeff(f1copy, vr))
  @show f1copy
  end
  @assert iszero(coeff(f1copy, vr))
  @assert iszero(coeff(f1copy, vz))

  f1copy = factor(f1copy)
  return f1copy
end

f1Ex = coeff(expand(f1), Ex)
f1Ey = coeff(expand(f1), Ey)
f1Ez = coeff(expand(f1), Ez)

Evector = [f1Ex, f1Ey, f1Ez]
othervector = vr .* [vr*cos(a), vr*sin(a), vz]

E0 = symbols("E0") # could get rid of this, but keep for clarity
contributions = Matrix{Any}(undef, 3, 3)
for i ∈ 1:3, j ∈ 1:3
  term = expand(ns * qs / E0 * othervector[i] * Evector[j])
  term = subs(term, qs^2 => wps^2 / (ns / (ms * E0))) # substitute for plasma freq
  contributions[i, j] = expand_trig(simplify(im * term / wps^2 * w))
end

function factor_parts(sol)
  iszero(m) && return sol
  solcopy = deepcopy(sol)
  solcopy = expand(solcopy)
  hasa1 = !iszero(coeff(solcopy, a1))
  hasac = !iszero(coeff(solcopy, ac))
  hasas = !iszero(coeff(solcopy, as))
  @assert hasa1 + hasac + hasas == 1
  solcopy /= hasa1 ? a1 : hasac ? ac : as
  solcopy /= cm
  (m0, m1) = coeffs(solcopy, m)
  output = coeff(m0, cos(a0)) * cos(a0)
  output += coeff(m0, sin(a0)) * sin(a0)
  output += coeff(m1, cos(a0)) * m * cos(a0)
  output += coeff(m1, sin(a0)) * m * sin(a0)
  output *= hasa1 ? a1 : hasac ? ac : as
  output *= cm

  @assert simplify(output / sol) == 1
  return output
end

function writefile()
  m == 0 && @assert !gyrobunch
  filename = "derivation_$(m==0 ? 0 : "m")_kz_$(kziszero)_kr_$(kriszero).txt"
  open(filename, "w") do file
    for i ∈ 1:3, j ∈ 1:3
      sh = shorthand(contributions[i, j])
      sh = factor_parts(sh)
      kriszero && (sh = subs(sh, kr => 0*kr))
      kziszero && (sh = subs(sh, kz => 0*kz))
      iszero(coeff(sh, kr)) || (sh = factor(sh, kr))
      iszero(coeff(sh, kz)) || (sh = factor(sh, kz))
      iszero(coeff(sh, w)) || (sh = factor(sh, w))
      write(file, "m$(i)$(j) = ")
      write(file, "$(sh)")
      write(file, "\n")
    end
  end
  txt = read(filename, String)
  replacements = Dict("kr" => "k⊥", "w" => "ω", "W" => "Ω",
   "a1" => "ϕ1", "as" => "ϕs", "ac" => "ϕc", "cos(a0)" => "cos(ϕ₀)",
   "sin(a0)" => "sin(ϕ₀)", "I" => "im", "pie" => "π")
  for r in replacements
    txt = replace(txt, r)
  end
  for i in 0:5
    txt = replace(txt, "z$(i)F" => "b$(i)F")
    txt = replace(txt, "z$(i)T" => "b$(i)T")
    txt = replace(txt, "r$(i)F" => "⊥$(i)F")
    txt = replace(txt, "r$(i)T" => "⊥$(i)T")
  end

  txt = replace(txt, "T0m" => "Tm")
  txt = replace(txt, "F0m" => "Fm")
  txt = replace(txt, "PyObject" => "")
  txt = replace(txt, "1.0*" => "")
  txt = replace(txt, "0.5*" => "1/2*")
  txt = replace(txt, "0.25*" => "1/4*")
  open(filename, "w") do file
    write(file, txt)
  end
end

writefile()


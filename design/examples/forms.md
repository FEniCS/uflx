# Variational forms

These examples use the core names from [language.md](../language.md) without
imports. Shared helper definitions are in the
[reference library](library.md).

## Poisson

For Poisson's equation $-\Delta u=f$ with homogeneous Dirichlet conditions,
the library supplies a triangle topology `T`, its coordinate function `x`
of value shape `(2,)`, and a real scalar continuous Lagrange space `V` on `T`:

```python
D = image(x, RD(2))
u, v, f = function(V, "u"), function(V, "v"), function(V, "f")
u_D, v_D, f_D = push_forward(u, D), push_forward(v, D), push_forward(f, D)
grad_u, grad_v = gradient(u_D, D), gradient(v_D, D)
i = index("i")

a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), D), (v, u))
L = form(integral(f_D * v_D, D), (v,))
```

The functions live on `T`; pushing them forward gives fields on the physical
mesh `D`. The entries `(v, u)` bind the test and trial arguments of `a`, in
that order. `L` binds only `v`, leaving `f` as the source coefficient.
The solver supplies its values, applies the boundary conditions and solves
the assembled system.

## Constants on a domain

A constant can have no support, or belong to a space on one domain. The
first can act wherever it is integrated; the second retains that domain's
support.

Let `F` and `S` be two meshes in `RD(2)`. The library supplies real spaces
`Q` (scalar) and `V` (shape `(2,)`) on `F.topology`, `W` (shape `(2,)`) on `S.topology`,
`R` of scalar constants with no domain, and `R_S` of scalar functions constant
on `S.topology`.

### Mean constraints

Use multipliers for the mean pressure on `F` and mean volume change on `S`:

```python
p, q = function(Q, "p"), function(Q, "q")
d, w = function(W, "d"), function(W, "w")
p_F, q_F = push_forward(p, F), push_forward(q, F)
d_S, w_S = push_forward(d, S), push_forward(w, S)
div_d, div_w = divergence(d_S, S), divergence(w_S, S)

r, s = function(R, "r"), function(R, "s")
k, z = function(R_S, "k"), function(R_S, "z")
k_S, z_S = push_forward(k, S), push_forward(z, S)
means = integral(p_F * s + q_F * r, F) + integral(k_S * div_w + z_S * div_d, S)
mean_form = form(means, ((q, w, s, z), (p, d, r, k)))
```

`r` and `s` have no support; the integrals decide where they act. `k_S` and
`z_S` have support `S`, so `k_S * q_F` raises `ValueError`. A coefficient
that varies by cell marker instead belongs to a library space constant on
each cell.

### Sensitivities

A constant whose sensitivity is needed must be a `Function`, because that is
what a directional derivative accepts as its variable. With the same spaces:

```python
nu, h = function(R, "nu"), function(R, "h")
mu_s = function(R_S, "mu_s")
u = function(V, "u")
i, j = index("i"), index("j")
two = Integer(2)
grad_d, grad_w = gradient(d_S, S), gradient(w_S, S)
eps_d = (grad_d[i, j] + grad_d[j, i]) / two
eps_w = (grad_w[i, j] + grad_w[j, i]) / two

e = integral(two * push_forward(mu_s, S) * index_sum(eps_d * eps_w, (i, j)), S)
de_dmu = form(directional_derivative(e, mu_s, mu_s, Integer(1)), (w, d))

grad_u = gradient(push_forward(u, F), F)
dissipation = integral(nu * index_sum(grad_u[i, j] * grad_u[i, j], (i, j)), F)
g_nu = form(directional_derivative(dissipation, nu, nu, h), (h,))
```

`de_dmu` is the solid form's sensitivity to `mu_s`. `g_nu` is a linear form
on `R`, a function in `dual_space(R)`. Data that needs no differentiation can
instead be a `runtime_array`.

## Helmholtz with a perfectly matched layer

A perfectly matched layer absorbs outgoing waves by stretching the
Helmholtz equation into the complex plane. For
$s_x=1+i\sigma_x/k_0$ and $s_y=1+i\sigma_y/k_0$, the equation is

$$-\nabla\cdot(A\nabla u)-k_0^2s_xs_yu=f,
\qquad A=\operatorname{diag}(s_y/s_x,\;s_x/s_y).$$

The library supplies a triangle topology `T`, its coordinate function `x` of
value shape `(2,)`, a complex scalar Lagrange space `V`, a real scalar data
space `F` on `T`, and a real scalar space of constants `C` with no domain.
The solver sets `sigma_x` and `sigma_y` to zero outside the layer:

```python
D = image(x, RD(2))
k0 = function(C, "k0")
sigma_x, sigma_y = push_forward(function(F, "sigma_x"), D), push_forward(function(F, "sigma_y"), D)
one, i_unit = Integer(1), ImaginaryUnit()
s_x = one + i_unit * sigma_x / k0
s_y = one + i_unit * sigma_y / k0

u, v, f = function(V, "u"), function(V, "v"), function(F, "f")
u_D, v_D, f_D = push_forward(u, D), push_forward(v, D), push_forward(f, D)
grad_u, grad_v = gradient(u_D, D), gradient(v_D, D)
j = index("j")

A = tensor((s_y / s_x, s_x / s_y))
stiffness = index_sum(A[j] * grad_u[j] * conjugate(grad_v[j]), (j,))
mass = k0 * k0 * s_x * s_y * u_D * conjugate(v_D)
a = form(integral(stiffness - mass, D), (conjugate(v), u))
L = form(integral(f_D * conjugate(v_D), D), (conjugate(v),))
```

The entries `(conjugate(v), u)` make `a` antilinear in `v` and linear in `u`;
`L` has the matching test argument. Binding `(v, u)` would declare the wrong
linearity: construction accepts that declaration, but a linearity check
rejects it. Replacing both stretch factors by one gives the ordinary
Helmholtz form, real when `V` is real.

`DomainMap` is real, so the complex coordinate change is performed
analytically and enters through `A` and $s_xs_y$. The adjoint conjugates the
integrated expression and swaps the arguments:

```python
a_adjoint = form(conjugate(integral(stiffness - mass, D)), (conjugate(u), v))
```

A first-order absorbing boundary condition adds
$\int_{\partial D} i k_0 u\bar v$, using restrictions to the boundary.

## Mixed Poisson with Raviart-Thomas elements

Mixed Poisson couples a flux $\sigma$ and a potential $u$, with
$\sigma=\nabla u$ and $\operatorname{div}\sigma=-f$. The library supplies a
triangle topology `T`, its coordinate function `x` of value shape `(2,)`, a
real Raviart–Thomas space `RT` of value shape `(2,)`, and a real discontinuous
Lagrange space `DG` of degree 0, both on `T`:

```python
D = image(x, RD(2))
sigma, tau = function(RT, "sigma"), function(RT, "tau")
u, v, f = function(DG, "u"), function(DG, "v"), function(DG, "f")
sigma_D, tau_D = push_forward(sigma, D), push_forward(tau, D)
u_D, v_D, f_D = push_forward(u, D), push_forward(v, D), push_forward(f, D)
i = index("i")

e = (
    index_sum(sigma_D[i] * tau_D[i], (i,))
    + u_D * divergence(tau_D, D)
    + divergence(sigma_D, D) * v_D
)
a = form(integral(e, D), ((tau, v), (sigma, u)))
L = form(integral(-f_D * v_D, D), ((tau, v),))
```

Each tuple is one argument in a product space. Thus `a` has argument spaces
`((RT, DG), (RT, DG))`; no mixed-space class is needed. `L` depends on only
the `v` member, which joint linearity permits. Replacing `sigma` and `v` by
`Zero` in the integrand isolates the `tau`–`u` block,
`u_D * divergence(tau_D, D)`.

The space supplies the transformation law:

| Space | Law | Cell value |
| ----- | --- | ---------- |
| `RT` | `transformation_law(("contravariant",), 1)` | $J\hat\sigma/\det J$ |
| `DG` | `transformation_law((), 0)` | The reference value. |

`expand_pull_backs` applies the Piola identity to the divergence, giving
$\widehat{\operatorname{div}}\hat\sigma/\det J$ even on curved cells, without
differentiating $J$. For a Nédélec space, the corresponding covariant law of
weight 0 makes `curl` transform by the contravariant law of weight 1.

## Hyperelasticity and the stress

Differentiate a strain energy to obtain both a stress field and variational
forms. For a Saint Venant–Kirchhoff material,
$S=2\operatorname{sym}(\partial\psi/\partial C)$.

The library supplies a triangle topology `T`, its coordinate function `x` of
value shape `(2,)`, a real vector Lagrange space `V` of value shape `(2,)`, a
real scalar space `M` for the Lamé parameters, both on `T`, and a real space
`R` of constant matrices of value shape `(2,2)` with no domain.

```python
D = image(x, RD(2))
one, zero, two = Integer(1), Integer(0), Integer(2)
I = tensor((tensor((one, zero)), tensor((zero, one))))
i, j, k, l, a = index("i"), index("j"), index("k"), index("l"), index("a")

u, v, w = function(V, "u"), function(V, "v"), function(V, "w")
mu, lmbda = push_forward(function(M, "mu"), D), push_forward(function(M, "lambda"), D)

def strain_energy(C):
    E = (C - I) / two
    return lmbda / two * index_sum(E[a, a], (a,)) ** two + mu * index_sum(E[i, j] * E[i, j], (i, j))

F = I + gradient(push_forward(u, D), D)
C = transpose(F) @ F
psi = strain_energy(C)
```

Differentiate with respect to a matrix symbol, evaluate at `C`, and use
symmetric unit directions:

```python
C_ = function(R, "C")
N = component_tensor((I[i, k] * I[j, l] + I[i, l] * I[j, k]) / two, (i, j))
S = two * component_tensor(directional_derivative(strain_energy(C_), C_, C, N), (k, l))
sigma = F @ S @ transpose(F) / determinant(F)
```

`S` evaluates to $\lambda\operatorname{tr}(E)I+2\mu E$. Both `S` and the
Cauchy stress `sigma` are fields on `D`, which the solver can interpolate or
evaluate for output.

The energy's Gateaux derivatives give the residual and tangent:

```python
Pi = integral(psi, D)
residual = gateaux_derivative(Pi, u, u, v)
tangent = gateaux_derivative(residual, u, u, w)
```

Since `V` is real, the first derivative binds `v`; the second binds `w`
after `v`.

## A classic FEniCS layer for the Poisson problem

A convenience layer can recover the familiar
`inner(grad(u), grad(v)) * dx` syntax while producing ordinary core graphs.
The hypothetical `uflx.legacy` module below is an
[implementation sketch](library.md#implementation-of-the-fenics-layer), not an existing
module. It adds no core subclasses.

### Mesh and markers

Use a two-triangle square. The library supplies `x`, a degree-1 coordinate
function on `square` with vertex values `((0,0), (1,0), (0,1), (1,1))`:

```python
K, I = ReferenceSimplex(2), ReferenceSimplex(1)
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
D = image(x, RD(2))
boundary = StaticTopology(
    (I,), (((0, 1), (1, 3), (3, 2), (2, 0)),), ((("inside", square, 0),),)
)
lower = StaticTopology((K,), (((0, 1, 2),),), ((("inside", square, 0),),))
left = StaticTopology((I,), (((2, 0),),), ((("inside", square, 0),),))
```

The layer wraps these in a mesh and marker dictionaries. Its user code sees
the familiar spaces, coefficients and measures:

```python
from uflx.legacy import (
    Coefficient,
    Constant,
    FacetNormal,
    FunctionSpace,
    Measure,
    Mesh,
    TestFunction,
    TrialFunction,
    VectorFunctionSpace,
    div,
    grad,
    inner,
)

mesh = Mesh(D, boundary)
V = FunctionSpace(mesh, "Lagrange", 1)
W = VectorFunctionSpace(mesh, "Lagrange", 1)
u, v = TrialFunction(V), TestFunction(V)
f, q = Coefficient(V, "f"), Coefficient(W, "q")
kappa = Constant("kappa")
n = FacetNormal(mesh)
dx = Measure("dx", domain=mesh, subdomain_data={1: lower})
ds = Measure("ds", domain=mesh, subdomain_data={1: left})

a = kappa * inner(grad(u), grad(v)) * dx
L = inner(f, v) * dx(1) - inner(div(q), v) * dx + inner(inner(q, n), v) * ds(1)
```

`a` is the Poisson stiffness form with conductivity `kappa`. `L` combines a
source on the lower triangle, a volume source from `div(q)` and a Neumann
flux on the left edge.

### What the layer supplies

| Legacy operation | Construction over the core |
| ---------------- | -------------------------- |
| Test and trial functions | Functions labelled by argument number, pushed forward to the mesh. |
| `inner` | Conjugation of the second operand; this sketch handles scalars and vectors. |
| `expr * dx` | An integral that binds the numbered arguments found in `expr`. |
| `dx(1)`, `ds(1)` | Integrals over side domains of marked topologies. |
| `FacetNormal(mesh)` | A placeholder replaced by the normal of each integration facet. |
| `grad`, `div`, `curl` | Core operations using the expression's support. |

For a real space `S = V.space`, `a` has argument spaces `(S, S)` and `L`
has `(S,)`. For complex `S`, conjugating `inner`'s second operand and binding
the test function as `conjugate(v)` gives `(conjugate_space(S), S)` and
`(conjugate_space(S),)`. An unconjugated test expression bound this way fails
a linearity check. A measure with no arguments returns an `Integral`.

The layer infers an unspecified `dx` domain from the integrand's support;
a constant integrand needs an explicit domain. `ds` needs a mesh because its
boundary is library data. The four measure slots are integral type, domain,
subdomain id and subdomain data; calling a measure fills them.

Restriction pulls back the largest subtrees on the mesh, so `grad(u)` becomes
a restricted gradient, not the derivative of a trace. Subtrees containing the
facet-normal placeholder are rebuilt with the facet's normal. This uses the
core's reconstruction contract.

The compiler receives only core graphs, and all core passes apply unchanged.
For Raviart–Thomas and Nédélec fields, `div` and `curl` retain the Piola
identities; on a surface, `grad` is tangential. A domainless `Constant` has
no gradient in this sketch: matching legacy behaviour would require a mesh
and a correctly shaped `Zero`.

Mixed-space helpers, `derivative` and interior-facet restriction syntax are
not implemented here. They would collect tuple arguments, construct Gateaux
derivatives and pull back through named sides respectively. Boundary
conditions, assembly and solving remain outside the language.

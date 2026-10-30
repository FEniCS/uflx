# Error estimation

These examples use the core names from [language.md](../language.md) without
imports. Shared helper definitions are in the
[reference library](library.md).

## A Bank-Weiser error estimator

The Bank–Weiser estimator solves an independent residual problem on each
cell. A discontinuous space lets us write all these local problems as one
global form, without requiring a global assembly.

For the Poisson problem, seek a quadratic correction `e` with zero vertex
values on each cell:

$$
\begin{aligned}
\int_K \nabla e \cdot \nabla v ={} & \int_K (f + \Delta u_h) \, v \\
& - \frac{1}{2} \int_{\partial K \cap \Omega} n_K \cdot (\nabla u_h|_K - \nabla u_h|_{K'}) \, v \\
& + \int_{\partial K \cap \Gamma_N} (g - \partial_n u_h) \, v,
\end{aligned}
$$

Here `u_h` is the Galerkin solution and the estimator is the energy norm of
`e`, as in Bulle et al.

### Setup

Use a two-triangle square, with its diagonal as the interior facet:

```python
K, I = ReferenceSimplex(2), ReferenceSimplex(1)
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
D = image(x, RD(2))
E = side_domain(D, diagonal)
B = side_domain(D, neumann_facets)
```

The library supplies `x`, a degree-1 coordinate function on `square` with
vertex values `((0,0), (1,0), (0,1), (1,1))`, and `neumann_facets`, a topology
with side `"inside"` into `square`. Only Neumann facets belong to it;
Dirichlet facets contribute nothing. All spaces below are real and scalar,
on `square`:

| Space | Functions |
| ----- | --------- |
| `V` | Continuous piecewise linears for the Galerkin solution. |
| `F` | Data fields. |
| `W` | Discontinuous quadratics with zero vertex values on each cell. |

### Local forms

```python
u_h, f, g = function(V, "u_h"), function(F, "f"), function(F, "g")
e, v = function(W, "e"), function(W, "v")
f_D, g_D, v_D = push_forward(f, D), push_forward(g, D), push_forward(v, D)
i = index("i")
half = Integer(1) / Integer(2)

# Local bilinear form.
grad_e = gradient(push_forward(e, D), D)
grad_v = gradient(v_D, D)
a = form(integral(index_sum(grad_e[i] * grad_v[i], (i,)), D), (v, e))
```

The right-hand side combines the cell residual with facet contributions:

```python
# Cell residual; lap_u vanishes for degree 1.
grad_u = gradient(push_forward(u_h, D), D)
lap_u = divergence(grad_u, D)
r_cell = integral((f_D + lap_u) * v_D, D)

# Interior-facet residual, shared by both cells.
plus, minus = side_map(D, E, "+"), side_map(D, E, "-")
n = normal(plus)
jump_n = index_sum((pull_back(grad_u, plus)[i] - pull_back(grad_u, minus)[i]) * n[i], (i,))
r_facets = integral(-half * jump_n * (pull_back(v_D, plus) + pull_back(v_D, minus)), E)

# The Neumann residual on the boundary.
inside = side_map(D, B, "inside")
n_B = normal(inside)
dn_u = index_sum(pull_back(grad_u, inside)[i] * n_B[i], (i,))
r_boundary = integral((pull_back(g_D, inside) - dn_u) * pull_back(v_D, inside), B)

r = form(r_cell + r_facets + r_boundary, (v,))
```

`W` makes the matrix of `a` block diagonal. The vector of `r` holds each
cell's residual, half of each interior jump and its Neumann contribution.
The basis of `W` stays in the library; it may use a quadratic Lagrange basis
and a change of basis, as in Bulle et al.

The sign of `jump_n` is the same from both sides: reversing the normal also
reverses the gradient jump. Restricting `v` from both sides assigns half the
jump to each cell.

### Assembly and the estimator

The [interior-facet walkthrough](coupling.md#from-the-facet-integral-to-a-kernel) explains
the lowering. This example produces four kernels:

| Term | Kernel output |
| ---- | ------------- |
| `a` | Local cell matrix. |
| `r_cell` | Cell residual. |
| `r_facets` | Contributions to both neighbouring cells, evaluating the jump once. |
| `r_boundary` | Contribution to the single boundary cell. |

A global assembler scatters these into a block-diagonal system. A local
assembler instead calls the cell and facet kernels for each cell, keeps its
rows, solves and computes the estimator, with no global array.

After solving for `e`, the cell-wise squared estimator is:

```python
cell = cell_domain(D, 0)
energy_density = index_sum(grad_e[i] * grad_e[i], (i,))
eta_squared = integral(pull_back(energy_density, cell) * density(cell), K)
```

Its cell index remains free, giving one value per cell. Alternatively, a
topology of `(cell, local facet)` pairs can use sides `"self"` and `"other"`.
Restricting `v` from `"self"` alone then makes the facet kernel return only
that cell's rows, with the pairs already grouped by cell.

## Dual-weighted residual error estimation

Weighting a residual by an adjoint solution estimates the error in a goal
functional. For a linear problem and the exact adjoint $z$,
$J(u)-J(u_h)=R(z)=R(z-z_h)$ for any $z_h$ in the Galerkin space.

Here the library supplies a triangle topology `T`, a coordinate function `x`
of value shape `(2,)`, and these real spaces on `T`:

| Space | Purpose |
| ----- | ------- |
| `V` | Scalar Lagrange degree 1 for the primal solution. |
| `W` | Scalar Lagrange degree 2 for an enriched adjoint approximation. |
| `F` | Scalar source and goal data. |
| `B` | Advection fields of value shape `(2,)`. |

For advection–diffusion, build the primal form, goal and adjoint:

```python
D = image(x, RD(2))
u, v, u_h, z_h = function(V, "u"), function(V, "v"), function(V, "u_h"), function(V, "z_h")
z = function(W, "z")
f, psi = push_forward(function(F, "f"), D), push_forward(function(F, "psi"), D)
b = push_forward(function(B, "b"), D)
i = index("i")

def a_expr(w, y):
    w_D, y_D = push_forward(w, D), push_forward(y, D)
    grad_w, grad_y = gradient(w_D, D), gradient(y_D, D)
    return index_sum(grad_w[i] * grad_y[i], (i,)) + index_sum(b[i] * grad_w[i], (i,)) * y_D

v_D = push_forward(v, D)
a = form(integral(a_expr(u, v), D), (v, u))
L = form(integral(f * v_D, D), (v,))
goal = form(integral(psi * v_D, D), (v,))
a_adjoint = form(integral(a_expr(u, v), D), (u, v))
```

The adjoint swaps the arguments of the real bilinear form. With `u_h` the
computed primal solution, `z` the enriched adjoint solution and `z_h` its
interpolant into `V`, the estimator is:

```python
R = L - action(a, 1, u_h)
eta = action(R, 0, z - z_h)
```

`R` has one argument in `V`. An action accepts any expression with the static
type of a function in `V`, so `z-z_h` from the richer space is allowed.
Substitution into `push_forward(v, D)` uses `V`'s transformation law.
Galerkin orthogonality makes the action on `z_h` zero in value; `eta` is a
functional, one number after assembly.

An assembled load can also be represented in the dual space:

```python
l = function(dual_space(V), "l")
R_h = l - action(a, 1, u_h)
```

The library evaluates actions on `l` as pairings with the assembled data.
For complex forms, constructing the adjoint also conjugates the expression.
To obtain cell-wise indicators, integrate `R` by parts into the cell residual
and facet jumps of the [Bank–Weiser example](#a-bank-weiser-error-estimator),
using `z-z_h` as the weight.

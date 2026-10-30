# Worked examples for the UFLx core

These examples accompany [LANGUAGE.md](LANGUAGE.md) and use its core names
without import lists. Each example states its library-supplied inputs locally;
shared helper definitions are in the [reference library](#reference-library).
The code constructs symbolic expressions; numerical results require later
passes or a solver.

## Examples

### Transpose

A transpose is a fixed composition of index operations. For a matrix expression
`A`, the core's `transpose(A)` builds:

```python
k = Index(("uflx.transpose", A, 0))
l = Index(("uflx.transpose", A, 1))
component_tensor(A[k, l], (l, k))
```

The [derived indices](LANGUAGE.md#derived-indices) keep the labels tied to `A`.
Using user indices instead gives a mathematically equal expression with a
different graph. `Indexed` rejects an operand of array rank other than 2,
including a scalar, with `IndexError`.

### Tensor algebra with `EllipsisIndex`

An ellipsis index lets a library define a tensor operation once for every
array rank. With derived ellipsis indices `r` and `s` and an ordinary derived
index `k`, the definitions are:

| Operation | Fixed composition |
| --------- | ----------------- |
| `outer(u, v)` | `component_tensor(u[r] * v[s], (r, s))` |
| `dot(u, v)` | `component_tensor(index_sum(u[r, k] * v[k, s], (k,)), (r, s))` |
| `inner(A, B)` | `index_sum(A[r] * conjugate(B[r]), (r,))` |

These are library operations. Core tensor arithmetic instead builds one index
per axis. For vectors and matrices, `dot(u, v)` and `u @ v` contract the same
axes but are not structurally equal. Nor does an ellipsis composition reduce
when it meets a binder with separate indices, as in `inner(transpose(A), B)`;
that requires a later rewrite.

### An affinely mapped triangle

Map the unit triangle to the triangle with vertices $(1,0)$, $(3,1)$ and
$(2,3)$:

$$\psi(X) = (1 + 2X_0 + X_1,\; X_0 + 3X_1).$$

Using the library's [UnitSimplex](#the-unit-simplex):

```python
T = UnitSimplex(2)
X = coordinate(T)
one, two, three = Integer(1), Integer(2), Integer(3)
psi = tensor((one + two * X[0] + X[1], X[0] + three * X[1]))
K = domain_map(T, psi, RD(2))
```

The geometric quantities evaluate to:

| Expression | Value |
| ---------- | ----- |
| `jacobian(K)` | $\begin{pmatrix}2 & 1 \\ 1 & 3\end{pmatrix}$ |
| `determinant(jacobian(K))` | $5$ |
| `density(K)` | $5$ |

The physical area is therefore $5/2$. Swapping the last two vertices negates
the determinant but leaves the density and area unchanged.

### Circulation along an edge

Reversing an edge changes circulation's sign but leaves a scalar arc-length
integral unchanged. Take the field $F=(0,x)$ on the unit triangle and its
hypotenuse, traversed in both directions:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
one = Integer(1)
e = domain_map(S, tensor((one - t, t)), T)
e_reversed = domain_map(S, tensor((t, one - t)), T)
F = tensor((Integer(0), X[0]))
l = index("l")

def circulation(edge):
    return integral(index_sum(pull_back(F, edge)[l] * jacobian(edge)[l, 0], (l,)), S)
```

| Expression | Value after evaluation |
| ---------- | ---------------------- |
| `circulation(e)` | $1/2$ |
| `circulation(e_reversed)` | $-1/2$ |
| `integral(pull_back(X[0], e) * density(e), S)` | $\sqrt{2}/2$ |
| `integral(pull_back(X[0], e_reversed) * density(e_reversed), S)` | $\sqrt{2}/2$ |

The tangent is $(-1,1)$ on `e` and $(1,-1)$ on `e_reversed`. The circulation
integrands are respectively $1-t$ and $-t$, while both densities are
$\sqrt{2}$.

### Integrating over an edge

Flux uses a directed normal; arc length uses its magnitude. On the hypotenuse
of the unit triangle, the outward normal scaled by the edge length is
$(\psi'_1,-\psi'_0)$:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
e = domain_map(S, tensor((Integer(1) - t, t)), T)
J = jacobian(e)

flux = integral(pull_back(X[0], e) * J[1, 0] - pull_back(X[1], e) * J[0, 0], S)
length = integral(density(e), S)
weighted_length = integral(pull_back(X[0] * X[1], e) * density(e), S)
```

The values are `flux` $=1$, `length` $=\sqrt{2}$ and `weighted_length`
$=\sqrt{2}/6$. The position field $X$ is tangent to the other two edges, so
all its outward flux passes through the hypotenuse. This agrees with the
divergence theorem: $\operatorname{div}X=2$ and the triangle has area $1/2$.

Here `density(e)` is $\sqrt{\det(J^T J)}$ with a $1\times1$ Gram matrix.
For a decomposed domain, `expand_integrals` inserts this density automatically.

### Integrating over the boundary of a triangle

Represent the three directed edges as one family, then sum their circulations:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
zero, one = Integer(0), Integer(1)
F = tensor((zero, X[0]))
i, k = index("i"), index("k")

starts = tensor((tensor((zero, zero)), tensor((one, zero)), tensor((zero, one))))
ends = tensor((tensor((one, zero)), tensor((zero, one)), tensor((zero, zero))))
psi_edges = component_tensor(starts[i, k] + t * (ends[i, k] - starts[i, k]), (k,))
edges = domain_map(S, psi_edges, T)

integrand = index_sum(pull_back(F, edges)[k] * jacobian(edges)[k, 0], (k,))
result = index_sum(integral(integrand, S), (i,))
```

`result` evaluates to $1/2$, the circulation of $(0,x)$ around the triangle,
as Green's theorem predicts. Without the outer `index_sum`, the result has
free index `i`: one circulation per edge. All edges share the parameter domain
`S`; the order of `starts` and `ends` sets their directions.

### Quadrature

An edge-midpoint rule integrates $f(x,y)=xy$ exactly on the unit triangle:

```python
T = UnitSimplex(2)
X = coordinate(T)
f = X[0] * X[1]
zero, half = Integer(0), Integer(1) / Integer(2)
sixth = Integer(1) / Integer(6)
i = index("i")

points = tensor((tensor((half, zero)), tensor((half, half)), tensor((zero, half))))
weights = tensor((sixth, sixth, sixth))
xs = domain_map(RD(0), points[i], T)
q = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

The points form one family with free index `i`, so `f` occurs once in `q`.
The value is $1/24$, mathematically equal to `integral(f, T)`. Replacing the
integral by quadrature changes its graph and is a pass, with the rule supplied
as an argument.

#### From quadrature to a loop

After lowering, evaluating the point integrals and reducing component tensors,
the expression is:

```python
index_sum(weights[i] * (points[i, 0] * points[i, 1]), (i,))
```

`expand_pull_backs` substitutes the points for `X`; the integral over `RD(0)`
returns its integrand. The compiler's `reduce_component_tensors` changes
`points[i][0]` to `points[i, 0]`, and similarly for component 1. A generated
kernel could be:

```python
weights = (1 / 6, 1 / 6, 1 / 6)
points = ((1 / 2, 0.0), (1 / 2, 1 / 2), (0.0, 1 / 2))

def kernel() -> float:
    q = 0.0
    for i in range(3):
        q += weights[i] * (points[i][0] * points[i][1])
    return q
```

The `IndexSum` becomes the loop, and literals become tables. Runtime arrays
of shapes `(q, 2)` and `(q,)` use the same construction, with an
`AbstractExtent` `q` as the kernel's loop bound, as in the next example.

### Integrating over a cut cell

For a cell cut by $\phi<0$, the solver supplies quadrature points and weights
for the retained region. The number of points is known only at assembly.
These examples use the following library inputs, none of which belongs to the
core:

| Input | Meaning |
| ----- | ------- |
| `RuntimeExtent(label)` | An `AbstractExtent` known at run time. |
| `RuntimeArray(label, value_shape)` | An expression with no support holding quadrature data, coefficients or constants. |
| `RuntimeIndex(label, extent)` | An `AbstractIndexValue` of shape `()`, such as a local facet number; read with `gathered_index(value, ())`. |

With their lower-case construction functions, the cut-cell integral is:

```python
T = UnitSimplex(2)
X = coordinate(T)
f = X[0] * X[1]
i = index("i")
q = runtime_extent("q")
points = runtime_array("points", (q, 2))
weights = runtime_array("weights", (q,))

xs = domain_map(RD(0), points[i], T)
a = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

The level set determines the solver's quadrature data; it does not appear in
the expression, since the core has no conditional. Point counts may differ
between cells. Representing all such rules in one topology-wide expression
requires the still-open [ragged indexing](LANGUAGE.md#extents) convention.

### The centroid of a triangle

An integral is an expression, so a centroid can be written as a ratio of
integrals. For the triangle with vertices $(1,0)$, $(3,1)$ and $(2,3)$:

```python
T = UnitSimplex(2)
X = coordinate(T)
one, two, three = Integer(1), Integer(2), Integer(3)
psi = tensor((one + two * X[0] + X[1], X[0] + three * X[1]))
K = domain_map(T, psi, RD(2))
Y = coordinate(RD(2))
k = index("k")

area = integral(density(K), T)
centroid = component_tensor(integral(pull_back(Y[k], K) * density(K), T) / area, (k,))
```

The density is $5$, the area is $5/2$ and the centroid is $(2,4/3)$.
Since the map is affine, evaluating `psi` at the reference centroid
$(1/3,1/3)$ gives the same value through a different graph.

A decomposed domain permits the shorter vector-integral form:

```python
D = decomposed_domain((K,))
centroid_D = integral(coordinate(D), D) / integral(Integer(1), D)
```

`expand_integrals` splits the vector integrand into scalar components.

### An interior facet integral

The interior-facet terms of the symmetric interior penalty method combine
traces from two cells on their shared facet. Here the cells are two triangles
forming a unit square, and the facet is their diagonal. Gradients are taken
on the cells before restriction to the facet.

#### Mesh and fields

Use the library classes [ReferenceSimplex](#the-unit-simplex-as-a-reference-cell)
and [StaticTopology](#a-static-topology):

```python
K, I = ReferenceSimplex(2), ReferenceSimplex(1)
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
```

The library supplies `X_space`, a vector Lagrange space of degree 1 and value
shape `(2,)`, and `V`, a real scalar discontinuous Lagrange space, both on `square`.
The solver gives the coordinate function `x` the vertex values
`((0, 0), (1, 0), (0, 1), (1, 1))`, in vertex-number order.
`runtime_array` supplies an array of solver data with no support (see
[kernel inputs](#integrating-over-a-cut-cell)).

```python
x = function(X_space, "x")
D = image(x, RD(2))
E = side_domain(D, diagonal)
e = cell_index(diagonal, 0)

u, v = function(V, "u"), function(V, "v")
u_D, v_D = push_forward(u, D), push_forward(v, D)
grad_u, grad_v = gradient(u_D, D), gradient(v_D, D)
```

#### Traces and the facet form

The side names `"+"` and `"-"` select triangles `(0, 1, 2)` and `(1, 3, 2)`.
Pulling back from either side gives an expression on `E`, so the traces can
be subtracted at the same physical point:

```python
plus, minus = side_map(D, E, "+"), side_map(D, E, "-")
u_p, u_m = pull_back(u_D, plus), pull_back(u_D, minus)
v_p, v_m = pull_back(v_D, plus), pull_back(v_D, minus)
du_p, du_m = pull_back(grad_u, plus), pull_back(grad_u, minus)
dv_p, dv_m = pull_back(grad_v, plus), pull_back(grad_v, minus)
```

Use the normal out of the `"+"` cell for both averaged normal derivatives.
`alpha[e]` is the penalty on facet `e`, including the usual cell-size factor
$1/h$.

```python
i = index("i")
n = normal(plus)
half = Integer(1) / Integer(2)
jump_u, jump_v = u_p - u_m, v_p - v_m
avg_dn_u = half * index_sum((du_p[i] + du_m[i]) * n[i], (i,))
avg_dn_v = half * index_sum((dv_p[i] + dv_m[i]) * n[i], (i,))
alpha = runtime_array("alpha", (diagonal.num_cells[0],))

f = -avg_dn_u * jump_v - jump_u * avg_dn_v + alpha[e] * jump_u * jump_v
a_facet = form(integral(f, E), (v, u))
```

`a_facet` binds the functions `v` and `u` on `square`. Its integrand has
support `E` and the free index `e` from the penalty array; the traces themselves
have no free index. Swapping the sides negates the normal and both jumps,
leaving every term unchanged.

For this square, the reference-facet inclusions are
$A_+(t) = (1-t, t)$ and $A_-(t) = (0, t)$. Both select the physical point
$(1-t, t)$. The cell Jacobians are $J_+ = I_2$ and
$J_- = \begin{pmatrix}0 & -1 \\ 1 & 1\end{pmatrix}$, so the normal is
$(1,1)/\sqrt{2}$ and the facet density is $\sqrt{2}$.

#### From the facet integral to a kernel

`expand_integrals` expresses the integral on the reference interval:

```python
facet_cell = cell_domain(E, 0)
form(
    index_sum(integral(pull_back(f, facet_cell) * density(facet_cell), I), (e,)),
    (v, u),
)
```

The side tables supply the neighbouring cell numbers `cells_p`, `cells_m`
and the reference coordinates of the facet vertices `positions_p`,
`positions_m`. The latter define the affine inclusions $A_\pm$ through
[affine_map](#side-tables-of-a-static-topology). For example, the `"+"` trace
lowers to:

```python
g_p = gathered_index(cells_p, (e,))
C_p = cell_map(square, 0, g_p)
t = coordinate(I)[0]
A_p = domain_map(I, tensor((Integer(1) - t, t)), K)
pull_back(pull_back(u, C_p), A_p)
```

The coordinate trace follows the same route, so the facet parameterisation is `pull_back(pull_back(x, C_p), A_p)`;
its Jacobian gives the facet density.

At a reference-facet point $t$, let $\hat u_+$ be the reference-cell field,
$J_+$ the cell Jacobian evaluated at $A_+(t)$, and $\hat m_+$ the outward
reference normal supplied by the side tables. The lowering stage gives:

| Expression on `E` | Value at the facet point |
| ----------------- | ------------------------ |
| `u_p` | $\hat u_+(A_+(t))$ |
| `du_p` | $J_+^{-T}\widehat\nabla\hat u_+(A_+(t))$ |
| `n` | $J_+^{-T}\hat m_+ / \lvert J_+^{-T}\hat m_+\rvert$ |
| `alpha[e]` | Unchanged: it has no support and shares the facet index. |

The `"-"` traces follow the same rules on the other cell. Here
$\hat m_+ = (1,1)/\sqrt{2}$. The compiler's `reduce_component_tensors` pass
then reduces the indexed vertex-coordinate expressions to reads from
`positions_p` and `positions_m`. The resulting form has the structure:

```python
form(index_sum(integral(g, I), (e,)), (v, u))
```

`g` is the lowered integrand on `I`, with the facet index `e` free. The
assembler loops over `e`, packing the two cell numbers, coordinate
coefficients, facet-vertex positions and penalty. For a fixed `e`, the kernel
reads the reference bases of `u`, `v` and `x` on both cells at $A_\pm(t)$.
The vertex positions replace legacy local-facet numbers and determine the
facet direction on each side.

#### Variations

- The penalty can use computed cell sizes instead of an array. With
  `B_p = cell_domain(D, 0, g_p)`, `integral(density(B_p), K)` is the measure
  of the `"+"` cell; a solver's per-cell size array `h` is read as `h[g_p]`.
- Restricting values alone needs no mesh:
  `pull_back(u, side_map(square, diagonal, "+"))` is a trace on the topology.
  Pushing it forward to `E` with the identity law gives `u_p`.

```{dropdown} Legacy UFL
This is the `dS` form
`-avg(dot(grad(u), n)) * jump(v) * dS - jump(u) * avg(dot(grad(v), n)) * dS + alpha / h * jump(u) * jump(v) * dS`.
Here the side maps make restriction explicit; jumps and averages are ordinary
arithmetic on the traces.
```

### A Bank-Weiser error estimator

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

#### Setup

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

#### Local forms

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

#### Assembly and the estimator

The [interior-facet walkthrough](#from-the-facet-integral-to-a-kernel) explains
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

### A fluid and a solid

Couple Stokes flow to linear elasticity through a shared interface. The
fluid velocity and solid displacement live on separate meshes; a multiplier
enforces equality of their velocities over a time step.

#### Domains

For a minimal example, split the unit square into a fluid triangle `(0,1,2)`
and a solid triangle `(1,3,2)`. The library supplies a degree-1 coordinate
function `x` on `square`, with vertex values `((0,0), (1,0), (0,1), (1,1))`:

```python
K, I = ReferenceSimplex(2), ReferenceSimplex(1)
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
D = image(x, RD(2))
```

A marked region is a topology with a side into its parent. Its side domain
is the submesh, and its side table is the entity map:

```python
fluid_cells = StaticTopology((K,), (((0, 1, 2),),), ((("inside", square, 0),),))
solid_cells = StaticTopology((K,), (((1, 3, 2),),), ((("inside", square, 0),),))

F = side_domain(D, fluid_cells)
S = side_domain(D, solid_cells)

interface = StaticTopology(
    (I,), (((1, 2),),), ((("fluid", fluid_cells, 0), ("solid", solid_cells, 0)),)
)
G = side_domain(F, interface)

inlet = StaticTopology((I,), (((2, 0),),), ((("fluid", fluid_cells, 0),),))
G_in = side_domain(F, inlet)
```

`interface` names the fluid and solid topologies. Restricting a field from
`D` would require a third side `("parent", square, 0)`. For a solver mesh,
the library supplies these topologies and side tables from markers and
entity maps; the remaining construction is unchanged.

#### Fields and bulk terms

The library supplies these real Lagrange spaces and the
[runtime-array inputs](#integrating-over-a-cut-cell):

| Space | Topology | Value shape |
| ----- | -------- | ----------- |
| `V`, `Q` | `fluid_cells` | `(2,)`, `()` |
| `W` | `solid_cells` | `(2,)` |
| `L` | `interface` | `(2,)` |

```python
u, v = function(V, "u"), function(V, "v")
p, q = function(Q, "p"), function(Q, "q")
d, w, d_old = function(W, "d"), function(W, "w"), function(W, "d_old")
lam, mu = function(L, "lam"), function(L, "mu")

u_F, v_F, p_F, q_F = (push_forward(f, F) for f in (u, v, p, q))
d_S, w_S, d_old_S = (push_forward(f, S) for f in (d, w, d_old))
lam_G, mu_G = push_forward(lam, G), push_forward(mu, G)

nu, mu_s, lam_s = runtime_array("nu", ()), runtime_array("mu_s", ()), runtime_array("lam_s", ())
dt, g = runtime_array("dt", ()), runtime_array("g", (2,))
i, j = index("i"), index("j")
half, two = Integer(1) / Integer(2), Integer(2)
```

The bulk terms are the usual Stokes and linear-elasticity integrals:

```python
grad_u = gradient(u_F, F)
grad_v = gradient(v_F, F)
div_u, div_v = divergence(u_F, F), divergence(v_F, F)
fluid = integral(
    nu * index_sum(grad_u[i, j] * grad_v[i, j], (i, j)) - p_F * div_v - q_F * div_u, F
)

grad_d = gradient(d_S, S)
grad_w = gradient(w_S, S)
eps_d, eps_w = half * (grad_d[i, j] + grad_d[j, i]), half * (grad_w[i, j] + grad_w[j, i])
div_d, div_w = divergence(d_S, S), divergence(w_S, S)
solid = integral(two * mu_s * index_sum(eps_d * eps_w, (i, j)) + lam_s * div_d * div_w, S)
```

#### Interface coupling

Restrict both fields to `G` before combining them. `lam` supplies the
traction, and `mu` tests the constraint $u=(d-d_{old})/dt$:

```python
from_fluid, from_solid = side_map(F, G, "fluid"), side_map(S, G, "solid")
u_G, v_G = pull_back(u_F, from_fluid), pull_back(v_F, from_fluid)
d_G, w_G = pull_back(d_S, from_solid), pull_back(w_S, from_solid)
d_old_G = pull_back(d_old_S, from_solid)

coupling = integral(
    index_sum(lam_G[i] * (w_G[i] / dt - v_G[i]) + mu_G[i] * (d_G[i] / dt - u_G[i]), (i,)), G
)
a = form(fluid + solid + coupling, ((v, q, w, mu), (u, p, d, lam)))
```

The load contains the inlet traction and old displacement:

```python
traction = integral(index_sum(g[i] * pull_back(v_F, side_map(F, G_in, "fluid"))[i], (i,)), G_in)
l = form(traction + integral(index_sum(mu_G[i] * d_old_G[i], (i,)) / dt, G), ((v, q, w, mu),))
```

The fields on `F` and `S` cannot combine directly: their supports differ.
Their traces share support `G`. The integrals have no support and can be
added; tuple arguments bind functions on all three topologies. No field is
extended by zero onto the parent mesh.

The interface is the image of the coordinate trace from `F`.
`normal(from_fluid)` points out of the fluid, and `normal(from_solid)` is
its negative. Lowering follows the [facet example](#from-the-facet-integral-to-a-kernel):
the coupling kernel reads the fluid and solid reference bases on their
respective cells and the multiplier basis on the facet. The two side tables
tell the assembler which cells to pack.

#### Variations

- Two marker values with the same cell type and sides use separate topologies;
  legacy `ds((3, 4))` is a sum of two integrals. Their kernels agree after
  renaming the topology placeholders. Different cell types use separate
  topology groups.
- A material coefficient varying by cell needs no submesh when all fields
  live on `square`: supply a table `nu_cells` and use `nu_cells[c]`, with
  `c = cell_index(square, 0)`, under an integral over `D`.
- For a moving fluid mesh, use `image(F.coordinate_function + m, RD(2))`,
  where the mesh displacement `m` on `fluid_cells` extends the solid motion.
  Functions on that topology push forward to the moved mesh unchanged.

```{dropdown} Legacy UFL
`dx(1)` and `ds(3)` become integrals over side domains of marked topologies.
A submesh's entity map becomes its side table into the parent.
```

### Constants on a domain

A constant can have no support, or belong to a space on one domain. The
first can act wherever it is integrated; the second retains that domain's
support.

Let `F` and `S` be two meshes in `RD(2)`. The library supplies real spaces
`Q` (scalar) and `V` (shape `(2,)`) on `F.topology`, `W` (shape `(2,)`) on `S.topology`,
`R` of scalar constants with no domain, and `R_S` of scalar functions constant
on `S.topology`.

#### Mean constraints

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

#### Sensitivities

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

### A 1D network in a 3D mesh

Couple pressure in a vessel network to pressure in the surrounding tissue.
Each field lives on its own topology; an exchange term compares them on the
network.

#### A conforming network

The library supplies a tetrahedral topology `T3`, its coordinate function
`x3` of value shape `(3,)`, and a topology `T1` of `ReferenceSimplex(1)`
segments. Each segment is a tetrahedral edge, with side `(("tissue", T3, 0),)`.
`V3` and `V1` are real scalar continuous Lagrange spaces on `T3` and `T1`.
The library supplies `runtime_array` for solver data with no support and
`runtime_extent` for sizes known at run time.

```python
D3 = image(x3, RD(3))
D1 = side_domain(D3, T1)
identity = transformation_law((), 0)

p3, q3 = function(V3, "p3"), function(V3, "q3")
p1, q1 = function(V1, "p1"), function(V1, "q1")
beta = runtime_array("beta", ())
m = index("m")

grad_p3 = gradient(push_forward(p3, D3), D3)
grad_q3 = gradient(push_forward(q3, D3), D3)
grad_p1 = gradient(push_forward(p1, D1), D1)
grad_q1 = gradient(push_forward(q1, D1), D1)
tissue = integral(index_sum(grad_p3[m] * grad_q3[m], (m,)), D3)
vessel = integral(index_sum(grad_p1[m] * grad_q1[m], (m,)), D1)

from_tissue = side_map(T3, T1, "tissue")
p3_1, q3_1 = pull_back(p3, from_tissue), pull_back(q3, from_tissue)
exchange = push_forward(beta * (p1 - p3_1) * (q1 - q3_1), D1, identity)
coupling = integral(exchange, D1)
a = form(tissue + vessel + coupling, ((q3, q1), (p3, p1)))
```

The gradient on `D1` is tangential, of value shape `(3,)`, so the vessel
term uses the same construction as the tissue term. `p3` and `p1` cannot
combine on their different supports. Restricting `p3` to `T1` makes the
exchange possible; pushing it forward then gives a field on `D1`.
The three integrals have no support and add into one form.

For each segment, the side tables select a containing tetrahedron and its
reference-edge inclusion. Continuity of `p3` makes the trace independent of
which neighbouring tetrahedron is selected. The shared vertices also let the
partition of `T3` determine that of `T1`.

Lowering follows the [facet walkthrough](#from-the-facet-integral-to-a-kernel).
The coupling kernel reads one segment and one tetrahedron, and its local
matrix couples their degrees of freedom. The assembler scatters through both
spaces' maps. The reference-edge positions have only twelve possibilities
(six edges, two directions), so a compiler can tabulate the tissue basis once
and select it with a runtime index.

For a nonconforming network, subdivide the segments into parts lying in
individual tetrahedra. A topology of these parts has sides into both meshes;
use its side domain for integration and restrict both fields. The inclusion
points can lie anywhere in the reference tetrahedron, so basis evaluation is
at run time, as for a cut cell.

#### Averaging around the network

A general $H^1(\Omega)$ field has no trace on a curve. A model may instead
average the tissue pressure over a ring of radius $R$ around each network
point. The solver supplies two orthonormal directions perpendicular to the
segment, in a real vector space `N` on `T1` of value shape `(3,)`. It also
supplies cosines, sines and weights summing to one as runtime tables; the core
has no `cos`.

```python
L = runtime_extent("L")
l = index("l")
R = runtime_array("R", ())
n1, n2 = function(N, "n1"), function(N, "n2")
cs, w = runtime_array("cs", (L, 2)), runtime_array("w", (L,))
x1 = D1.coordinate_function

ring = domain_map(
    T1, component_tensor(x1[m] + R * (cs[l, 0] * n1[m] + cs[l, 1] * n2[m]), (m,)), D3
)
p3_avg = index_sum(w[l] * pull_back(push_forward(p3, D3), ring), (l,))
q3_avg = index_sum(w[l] * pull_back(push_forward(q3, D3), ring), (l,))

exchange_avg = push_forward(beta * (p1 - p3_avg) * (q1 - q3_avg), D1, identity)
coupling_avg = integral(exchange_avg, D1)
```

`ring` maps `T1` into `D3`, with free index `l`. Summing over `l` gives a
field on `T1`, ready to combine with `p1`. At zero radius it reduces in value
to the conforming trace.

The ring pull-back remains after core lowering. After quadrature on the
segments, a point-location pass identifies the tetrahedron and reference
position of each ring point. The family has segment, quadrature-point and
ring-point indices; a segment kernel therefore reads a list of tissue cells
rather than one cell. Its interface follows from the expression.

### The Laplace-Beltrami operator

The same `gradient` writes Poisson's equation on a surface. Let `T` be a
topology of `ReferenceSimplex(2)` triangles, `x` a degree-1 coordinate
function of value shape `(3,)`, and `V` a real scalar Lagrange space on `T`.
The library supplies the topology, space and coordinate values:

```python
M = image(x, RD(3))

u, v, f = function(V, "u"), function(V, "v"), function(V, "f")
u_M, v_M, f_M = push_forward(u, M), push_forward(v, M), push_forward(f, M)

i = index("i")
grad_u = gradient(u_M, M)
grad_v = gradient(v_M, M)

a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), M), (v, u))
L = form(integral(f_M * v_M, M), (v,))
```

`M` has dimension 2 in `RD(3)`, so `grad_u` and `grad_v` are tangential
vectors with three components. Neither form needs an orientation.

On each reference cell, with $J$ the $3\times2$ coordinate Jacobian,
lowering gives the stiffness integrand

$$\widehat\nabla u^T (J^T J)^{-1}\widehat\nabla v\,
\sqrt{\det(J^T J)}.$$

The inverse metric comes from the pseudo-inverse in the gradient pull-back;
the square root is the surface density.

For second derivatives, set `P = gradient(coordinate(M), M)` and
`H = gradient(grad_u, M)`. `P` projects onto the tangent plane. `H` need not
be symmetric; the covariant Hessian is $P H P$.

```{dropdown} Legacy UFL
This is `inner(grad(u), grad(v)) * dx` on a mesh of topological dimension 2
in three dimensions. The domain dimensions determine the surface gradient.
```

### A deformed domain

Move a mesh by changing its coordinate function. Functions on the topology
retain their degrees of freedom and push forward to the moved mesh.

For a two-triangle square, let the library supply the coordinate function `x`
with vertex values `((0,0), (1,0), (0,1), (1,1))`, a real vector displacement
space `W` of value shape `(2,)`, and a real scalar Lagrange space `V`, all on
`square`:

```python
K = ReferenceSimplex(2)
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
i = index("i")
s = function(W, "s")

D_s = image(x + s, RD(2))

volume = integral(Integer(1), D_s)

u, v = function(V, "u"), function(V, "v")
grad_u, grad_v = gradient(push_forward(u, D_s), D_s), gradient(push_forward(v, D_s), D_s)
a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), D_s), (v, u))
```

For a cell `C = cell_domain(D_s, 0)`, the reference coordinate gradient is
$J+\widehat\nabla s$. Lowering therefore inserts:

| Quantity | Factor on the reference cell |
| -------- | ---------------------------- |
| Volume density | $\lvert\det(J+\widehat\nabla s)\rvert$ |
| Gradient pull-back | $(J+\widehat\nabla s)^{-T}$ |

Core derivatives hold domains fixed. The [shape-derivative extension](#a-shape-derivative)
first lowers the form, exposing `s` in these factors, then differentiates.
An argument as the direction gives a shape gradient.

A surface uses the same construction. For a triangle topology `T`, coordinate
function `x_T` of value shape `(3,)` and displacement space `W_T` on `T` of
the same value shape:

```python
w = function(W_T, "w")
M_w = image(x_T + w, RD(3))

area = integral(Integer(1), M_w)
```

The cell Jacobian is $3\times2$ and the density is $\sqrt{\det(J^T J)}$.

```{dropdown} Isogeometric analysis and shells
For extraction-based isogeometric analysis, store the NURBS map `F` as a
function on the topology of knot spans and use `image(F, RD(3))`. Multiple
patches use the combined topology; lowering inserts the geometric factors.

For a shell midsurface `image(phi_0, RD(3))`, the cell Jacobian is the
reference gradient of `phi_0`. Its Gram matrix is the first fundamental form,
and its density is the midsurface area element.
```

### A shape derivative

Shape differentiation must expose the geometry before differentiating.
Let `T` be a triangle topology, `X` a real vector Lagrange space of degree 1
and value shape `(2,)`, and `V` a real scalar Lagrange space, both on `T`.
The solver supplies coordinate values for `x`.
`shape_derivative` is the construction function of the
[extension](LANGUAGE.md#shape-derivatives):

```python
x, h, v = function(X, "x"), function(X, "h"), function(X, "v")
D = image(x, RD(2))
grad_u = gradient(push_forward(function(V, "u"), D), D)
i = index("i")
energy = integral(index_sum(grad_u[i] * grad_u[i], (i,)), D)

d_energy = shape_derivative(energy, x, h)
g_energy = form(shape_derivative(energy, x, v), (v,))
```

`d_energy` is the derivative in the fixed direction `h`; `g_energy` binds an
argument and is the shape gradient. The extension runs in this order:

1. Lower the form so `x` appears in Jacobians and densities.
2. Replace `ShapeDerivative(g, x, h)` by `directional_derivative(g, x, x, h)`.
3. Apply derivatives and lower again.

Differentiating `energy` directly before step 1 gives `Zero`: `x` occurs only
inside the domain, which substitution does not enter. The extension needs
only the new node and this lowering rule; core differentiation does the rest.

### The flux through a surface

Flux through a closed surface needs an outward orientation; area does not.
Let `T` be a triangle topology with reference cell `K = ReferenceSimplex(2)`,
and `x` its coordinate function of value shape `(3,)`. The library provides
`T.orientation`, with signs $\pm1$ telling whether each cell's vertex order
agrees with the outward orientation.

```python
K = ReferenceSimplex(2)
M = image(x, RD(3))
cells = cell_domain(M, 0)
J = jacobian(cells)
r, l = index("r"), index("l")

def cross(a, b):
    return tensor((
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ))

nu = cross(component_tensor(J[r, 0], (r,)), component_tensor(J[r, 1], (r,)))

o = pull_back(T.orientation, cell_map(T, 0))

flux = index_sum(
    integral(o * index_sum(pull_back(coordinate(M), cells)[l] * nu[l], (l,)), K),
    (cell_index(T, 0),),
)

area = integral(Integer(1), M)
```

`nu` is the directed cofactor normal: its magnitude already includes the
surface density. Thus `flux` integrates on the reference cells without
another density factor. The signs `o` make every contribution outward, and
the position field's flux is three times the enclosed volume.

Consistently outward vertex ordering makes `o` unnecessary. Arbitrary ordering
without the signs instead lets contributions cancel incorrectly. `area` uses
$\sqrt{\det(J^T J)}$ and is independent of ordering. On a nonorientable
surface such as a Möbius strip, orientation is `None`: area remains defined,
but a global outward flux does not.

### Helmholtz with a perfectly matched layer

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

### Mixed Poisson with Raviart-Thomas elements

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

### Dual-weighted residual error estimation

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

### Hyperelasticity and the stress

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

### Symmetric tensors

Tensor symmetry can live in the data representation. A 2D stiffness tensor
with minor and major symmetries has 16 entries but only six stored values.

Let `T` be a triangle topology, `x` its coordinate function of value shape
`(2,)`, `V` a real vector space of the same value shape on `T`, and `K` a
real space of constants of value shape `(6,)`. The library supplies an
`AbstractIndexValue` table `packed`, of shape `(2,2,2,2)` and extent 6,
selecting the stored value for each tensor entry:

```python
D = image(x, RD(2))
u, v = function(V, "u"), function(V, "v")
i, j, k, l = index("i"), index("j"), index("k"), index("l")
two = Integer(2)
c = function(K, "c")
stiffness = component_tensor(c[gathered_index(packed, (i, j, k, l))], (i, j, k, l))

grad_u, grad_v = gradient(push_forward(u, D), D), gradient(push_forward(v, D), D)
eps_u = (grad_u + transpose(grad_u)) / two
eps_v = (grad_v + transpose(grad_v)) / two
elastic = form(
    integral(index_sum(stiffness[i, j, k, l] * eps_u[k, l] * eps_v[i, j], (i, j, k, l)), D),
    (v, u),
)
```

The kernel receives six values for `c`; gathering reconstructs the full
stiffness tensor. The core needs no symmetry attribute.

```{dropdown} Symmetric directions and Voigt notation
For an energy $C_{01}^2$, unconstrained differentiation gives $2C_{01}$ at
entry `(0,1)` and zero at `(1,0)`. Differentiating in symmetric unit directions
gives $C_{01}$ at both. The results agree in every symmetric direction, but
only the second is symmetric.

Voigt and Mandel representations scale off-diagonal strain entries by 2 and
$\sqrt{2}$ respectively. A library can read either through a table and scale
factors. Differentiation with respect to stored Voigt values differs from
differentiation with respect to the full tensor by a factor of 2 off the diagonal.
```

### A classic FEniCS layer for the Poisson problem

A convenience layer can recover the familiar
`inner(grad(u), grad(v)) * dx` syntax while producing ordinary core graphs.
The hypothetical `uflx.legacy` module below is an
[implementation sketch](#implementation-of-the-fenics-layer), not an existing
module. It adds no core subclasses.

#### Mesh and markers

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

#### What the layer supplies

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

## Reference library

These supporting classes are outside the core. The simplex cell type,
static topology and reference-simplex definitions carry the original draft's
LLM-generated-content designation; the centroid example was also LLM-proposed.

### The unit simplex

Define the unit simplex in $R^n$,

$\{p \in R^n : p_i \ge 0 \text{ for every } i, \ \sum_i p_i \le 1\}$,

as a concrete `AbstractCoordinateDomain`.

```python
class UnitSimplex(AbstractCoordinateDomain):
    identifier = "library.UnitSimplex"

    def __init__(self, n: int):
        if type(n) is not int:
            raise TypeError("Dimension must be an int.")
        if n < 0:
            raise ValueError("Dimension must be non-negative.")
        self._n = n

    @property
    def dimension(self) -> int:
        return self._n
```

`UnitSimplex(1)` is the unit interval $[0, 1]$. `UnitSimplex(2)` is the unit
triangle with vertices $(0, 0)$, $(1, 0)$ and $(0, 1)$.

### The simplex as a cell type

`Simplex(n)` is the simplex of dimension $n \ge 0$ as a cell type, with
$n + 1$ vertices. `Simplex(0)` is a point.

```python
class Simplex(AbstractCellType):
    identifier = "library.Simplex"

    def __init__(self, n: int):
        if type(n) is not int:
            raise TypeError("Dimension must be an int.")
        if n < 0:
            raise ValueError("Dimension must be non-negative.")
        self._n = n

    @property
    def dimension(self) -> int:
        return self._n

    @property
    def num_vertices(self) -> int:
        return self._n + 1
```

### A static topology

`StaticTopology(cell_types, cells, sides)` stores fixed connectivity, mainly
for testing UFLx ideas, including Alfeld splits. It has no coordinates and
does not check its sides; vertex positions come from a mesh's coordinate
function.

| Convention | Meaning |
| ---------- | ------- |
| `cells[j][c][i]` | Topology vertex number of local vertex `i` in cell `c` of type `j`. |
| Vertex pairing | Each side uses the same vertex numbers as its target topology. |
| Inclusions | Affine maps through the paired reference vertices; sides cover sub-entities, not pieces of cells. |
| Two cells in the same group | The lower-numbered cell is on the first side. |
| Two cells in different groups | Each side's group selects its cell. |
| Sub-entities below facet dimension | Select the lowest-numbered containing cell. |

The library derives the containing cells, inclusions and facet normals as
described in [Side tables](#side-tables-of-a-static-topology).

```python
class StaticTopology(AbstractTopology):
    identifier = "library.StaticTopology"

    def __init__(
        self,
        cell_types: tuple[AbstractCellType, ...],
        cells: tuple[tuple[tuple[int, ...], ...], ...],
        sides: tuple[tuple[tuple[str, AbstractTopology, int], ...], ...],
    ):
        if len(cells) != len(cell_types):
            raise ValueError("Need one tuple of cells for each cell type.")
        self._cell_types = cell_types
        self._cells = cells
        self._sides = sides

    @property
    def cell_types(self) -> tuple[AbstractCellType, ...]:
        return self._cell_types

    @property
    def num_cells(self) -> tuple[int, ...]:
        return tuple(len(c) for c in self._cells)

    @property
    def sides(self) -> tuple[tuple[tuple[str, AbstractTopology, int], ...], ...]:
        return self._sides
```

With `K = Simplex(2)` and `I = Simplex(1)`:

```python
# Two triangles forming a square.
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())

# Three triangles meeting at vertex 3 (an Alfeld split).
alfeld = StaticTopology((K,), (((0, 1, 3), (1, 2, 3), (2, 0, 3)),), ())

# Four boundary facets, each seen from one triangle.
boundary = StaticTopology(
    (I,), (((0, 1), (1, 3), (3, 2), (2, 0)),), ((("inside", square, 0),),)
)

# One diagonal facet, seen from two triangles.
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
```

Facets are matched by vertex sets, using the facets of `K`. The square's
diagonal and the Alfeld split's three facets containing vertex 3 each belong
to two cells; their remaining facets belong to one. The four topologies have
`num_cells` `(2,)`, `(3,)`, `(4,)` and `(1,)`, respectively.

#### Side tables of a static topology

The library builds [side tables](LANGUAGE.md#side-tables) from reference
vertices, which the core does not read. For a $k$-dimensional reference cell
$F$ with vertices $V_v$, `affine_map(F, points)` maps the first $k+1$ vertices
to rows $p_v$ of `points`:

$$
X \mapsto p_0 + \sum_{d < k} (p_{d+1} - p_0) \, \big(E^{-1} (X - V_0)\big)_d,
$$

Here $E$ has columns $V_{d+1}-V_0$. The implementation uses a derived index:

```python
def affine_map(F, points):
    l = Index(("library.affine_map", F, points, 0))
    k, V = F.dimension, F.vertices
    if k == 0:
        return component_tensor(points[0, l], (l,))
    E = tensor(tuple(tensor(tuple(V[d + 1, r] - V[0, r] for d in range(k))) for r in range(k)))
    Y = inverse(E) @ (coordinate(F) - V[0])
    total = points[0, l]
    for d in range(k):
        total = total + (points[d + 1, l] - points[0, l]) * Y[d]
    return component_tensor(total, (l,))
```

Those $k+1$ vertices must be affinely independent; later rows of `points`
are ignored. A sub-entity with vertex numbers $(w_0,w_1,\dots)$ in `K` uses
row `K.vertices[w_v]` for point $v$. For example:

| Parent | Sub-entity vertices | Reference cell | Inclusion |
| ------ | ------------------- | -------------- | --------- |
| Unit triangle | `(1, 2)` | Unit interval | $t\mapsto(1-t,t)$ |
| Unit interval | `(1,)` | Reference point | The endpoint $1$ |

For group $i$ with reference cell `F_i`, let `e` be its cell index and
`positions[e, v]` the position of vertex `v` in the side's reference cell
`K_j`. With derived indices `v` and `k`, its inclusion is:

```python
points = component_tensor(positions[e, v, k], (v, k))
inclusion = domain_map(F_i, affine_map(F_i, points), K_j)
```

For a group of facets, the reference normal is the cofactor vector of the
Jacobian of the inclusion, turned away from the centroid of $K_j$:

```python
J_A = jacobian(inclusion)  # value shape (k, k - 1)
C_r = tensor(tuple(
    tensor((Integer(1 if q == r else 0), J_A[q, 0], ..., J_A[q, k - 2]))
    for q in range(k)
))
m0 = tensor((determinant(C_0), ..., determinant(C_(k - 1))))
centroid = (V[0] + V[1] + ... + V[m_j - 1]) / Integer(m_j)
away = m0 @ (points[0] - centroid)
reference_normal = away / absolute(away) * m0
```

Here `V = K_j.vertices`, $m_j$ is its vertex count, and the sums nest to the
left. All sides preserve the vertex order of cell `e`, so their inclusions
parameterise it in the same way.

### The unit simplex as a reference cell

`ReferenceSimplex(n)` is the cell type `Simplex(n)` as a reference cell, the
unit simplex in $R^n$. Vertex 0 is the origin, and vertex $v \ge 1$ is the unit
vector $e_{v-1}$. `ReferenceSimplex(0)` is the reference point.

```python
class ReferenceSimplex(Simplex, AbstractReferenceCell):
    identifier = "library.ReferenceSimplex"

    @property
    def vertices(self) -> AbstractExpression:
        n = self.dimension
        if n == 0:
            return Zero((1, 0), (), ())
        return tensor(
            tuple(
                tensor(tuple(Integer(1 if k == v - 1 else 0) for k in range(n)))
                for v in range(n + 1)
            )
        )
```

`vertices` is a property of the library, not of `AbstractReferenceCell`. The
side tables of a `StaticTopology` read it (see [Side tables of a static
topology](#side-tables-of-a-static-topology)).

`ReferenceSimplex(2)` and `Simplex(2)` are different cell types, so a topology
of `Simplex(2)` cells has no reference cells. `ReferenceSimplex(2)` and
`UnitSimplex(2)` are the same set but different coordinate domains, so an
expression on one is not an expression on the other. Likewise,
`ReferenceSimplex(0)` and `RD(0)` are different coordinate domains.

### Implementation of the FEniCS layer

The [Poisson example](#a-classic-fenics-layer-for-the-poisson-problem) uses
these wrappers. `Lagrange(topology, degree, shape)`, `Constants(shape)` and
`Fields(domain, shape)` stand for library spaces: finite element functions,
domainless constants and fields without an element, respectively. The plain
Python wrapper classes are outside the core.

```python
# uflx/legacy.py
from uflx import (
    AbstractExpression,
    Base,
    Function,
    conjugate,
    curl as core_curl,
    divergence,
    form,
    function,
    gradient,
    index,
    index_sum,
    integral,
    normal,
    pull_back,
    push_forward,
    side_domain,
    side_map,
)


class FunctionSpace:
    # A space on the topology of the mesh, with the mesh that its functions are pushed to.
    def __init__(self, mesh, family, degree, shape=()):
        self.mesh = mesh
        self.space = Lagrange(mesh.domain.topology, degree, shape)


def VectorFunctionSpace(mesh, family, degree):
    return FunctionSpace(mesh, family, degree, (mesh.domain.ambient_domain.dimension,))


def TestFunction(V):
    return push_forward(function(V.space, ("argument", 0)), V.mesh.domain)


def TrialFunction(V):
    return push_forward(function(V.space, ("argument", 1)), V.mesh.domain)


def Coefficient(V, name):
    return push_forward(function(V.space, name), V.mesh.domain)


def Constant(name, shape=()):
    return function(Constants(shape), name)


def FacetNormal(mesh):
    # A placeholder on the mesh, which restrict replaces by the normal of the facets.
    n = mesh.domain.ambient_domain.dimension
    return function(Fields(mesh.domain, (n,)), "facet normal")


def is_facet_normal(obj):
    return isinstance(obj, Function) and obj.label == "facet normal"


def grad(f):
    return gradient(f, f.support)


def div(w):
    return divergence(w, w.support)


def curl(w):
    return core_curl(w, w.support)


def inner(a, b):
    # For scalars and vectors, conjugating the second operand as legacy UFL does.
    if a.value_shape == ():
        return a * conjugate(b)
    i = index("inner")
    return index_sum(a[i] * conjugate(b[i]), (i,))


def contains(obj, predicate):
    if predicate(obj):
        return True
    if isinstance(obj, tuple):
        return any(contains(item, predicate) for item in obj)
    return isinstance(obj, Base) and any(
        contains(item, predicate) for item in obj.construction_arguments
    )


def arguments(expression):
    # The legacy extract_arguments: the functions labelled ("argument", n), by n,
    # which the walk finds inside the push-forwards.
    found = {}

    def visit(obj):
        if isinstance(obj, Function) and isinstance(obj.label, tuple):
            if obj.label[0] == "argument":
                found[obj.label[1]] = obj
        if isinstance(obj, tuple):
            for item in obj:
                visit(item)
        elif isinstance(obj, Base):
            for item in obj.construction_arguments:
                visit(item)

    visit(expression)
    return tuple(found[n] for n in sorted(found))


def restrict(obj, E):
    # The legacy apply_restrictions, at construction: the largest subtrees on the
    # mesh are pulled back to E, and the facet normal of the mesh becomes the
    # normal of E.
    if is_facet_normal(obj):
        return normal(side_map(obj.support, E, "inside"))
    if isinstance(obj, tuple):
        return tuple(restrict(item, E) for item in obj)
    if not isinstance(obj, AbstractExpression) or obj.support is None:
        return obj
    if not contains(obj, is_facet_normal):
        return pull_back(obj, side_map(obj.support, E, "inside"))
    return obj.__class__(*restrict(obj.construction_arguments, E))


class Mesh:
    def __init__(self, domain, boundary):
        self.domain, self.boundary = domain, boundary


class Measure:
    def __init__(self, integral_type, domain=None, subdomain_id="everywhere", subdomain_data=None):
        self.integral_type = integral_type
        self.domain = domain
        self.subdomain_id = subdomain_id
        self.subdomain_data = subdomain_data

    def __call__(self, subdomain_id="everywhere", domain=None, subdomain_data=None):
        # dx(1), or dx(domain=mesh, subdomain_data=markers), fills the slots.
        return Measure(
            self.integral_type,
            self.domain if domain is None else domain,
            subdomain_id,
            self.subdomain_data if subdomain_data is None else subdomain_data,
        )

    def __rmul__(self, integrand):
        # expr * dx is a form in legacy UFL, so the measure binds the arguments.
        mesh = self.domain
        if self.subdomain_id != "everywhere":
            # The marker value names a topology with a side "inside" into the mesh.
            E = side_domain(mesh.domain, self.subdomain_data[self.subdomain_id])
        elif self.integral_type == "ds":
            E = side_domain(mesh.domain, mesh.boundary)
        else:
            E = integrand.support if mesh is None else mesh.domain
        if E != integrand.support:
            integrand = restrict(integrand, E)
        expression = integral(integrand, E)
        args = arguments(expression)
        if not args:
            return expression
        # Legacy forms are antilinear in the test function.
        return form(expression, (conjugate(args[0]),) + args[1:])


dx = Measure("dx")
ds = Measure("ds")
```

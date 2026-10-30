# Facets and coupled meshes

These examples use the core names from [language.md](../language.md) without
imports. Shared helper definitions are in the
[reference library](library.md).

## An interior facet integral

The interior-facet terms of the symmetric interior penalty method combine
traces from two cells on their shared facet. Here the cells are two triangles
forming a unit square, and the facet is their diagonal. Gradients are taken
on the cells before restriction to the facet.

### Mesh and fields

Use the library classes [ReferenceSimplex](library.md#the-unit-simplex-as-a-reference-cell)
and [StaticTopology](library.md#a-static-topology):

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
[kernel inputs](integration.md#integrating-over-a-cut-cell)).

```python
x = function(X_space, "x")
D = image(x, RD(2))
E = side_domain(D, diagonal)
e = cell_index(diagonal, 0)

u, v = function(V, "u"), function(V, "v")
u_D, v_D = push_forward(u, D), push_forward(v, D)
grad_u, grad_v = gradient(u_D, D), gradient(v_D, D)
```

### Traces and the facet form

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

### From the facet integral to a kernel

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
[affine_map](library.md#side-tables-of-a-static-topology). For example, the `"+"` trace
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

### Variations

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

## A fluid and a solid

Couple Stokes flow to linear elasticity through a shared interface. The
fluid velocity and solid displacement live on separate meshes; a multiplier
enforces equality of their velocities over a time step.

### Domains

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

### Fields and bulk terms

The library supplies these real Lagrange spaces and the
[runtime-array inputs](integration.md#integrating-over-a-cut-cell):

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

### Interface coupling

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

### Variations

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

## A 1D network in a 3D mesh

Couple pressure in a vessel network to pressure in the surrounding tissue.
Each field lives on its own topology; an exchange term compares them on the
network.

### A conforming network

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

### Averaging around the network

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

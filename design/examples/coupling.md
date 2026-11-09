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
square = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)), ())
diagonal = StaticTopology((I,), ((1, 2),), ((("+", square, 0), ("-", square, 1)),))
```

The library supplies `X_space`, a vector Lagrange space of degree 1 and value
shape `(2,)`, and `V`, a real scalar discontinuous Lagrange space, both on `square`.
The solver gives the coordinate function `x` the vertex values
`((0, 0), (1, 0), (0, 1), (1, 1))`, in vertex-number order.
`runtime_array` supplies an array of solver data with no domain (see
[kernel inputs](integration.md#integrating-over-a-cut-cell)).

```python
x = function(X_space, "x")
D = mesh(x, RD(2))
E = side_domain(D, diagonal, "+")

u, v = function(V, "u"), function(V, "v")
u_D, v_D = push_forward(u, D, transformation_law(V)), push_forward(v, D, transformation_law(V))
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
`alpha` is the supplied scalar penalty for this local facet, including the
usual cell-size factor $1/h$.

```python
i = index("i")
n = normal(plus)
half = Integer(1) / Integer(2)
jump_u, jump_v = u_p - u_m, v_p - v_m
avg_dn_u = half * index_sum((du_p[i] + du_m[i]) * n[i], (i,))
avg_dn_v = half * index_sum((dv_p[i] + dv_m[i]) * n[i], (i,))
alpha = runtime_array("alpha", ())

f = -avg_dn_u * jump_v - jump_u * avg_dn_v + alpha * jump_u * jump_v
a_facet = form(integral(f, E), (v, u))
```

`a_facet` binds the functions `v` and `u` on `square`. Its integrand has
domain `E` and no free indices. The package supplies the two cells' local
inputs and the facet penalty for each kernel call. Swapping the sides negates the normal and both jumps,
leaving every term unchanged.

For this square, the reference-facet inclusions are
$A_+(t) = (1-t, t)$ and $A_-(t) = (0, t)$. Both select the physical point
$(1-t, t)$. The cell Jacobians are $J_+ = I_2$ and
$J_- = \begin{pmatrix}0 & -1 \\ 1 & 1\end{pmatrix}$, so the normal is
$(1,1)/\sqrt{2}$ and the facet density is $\sqrt{2}$.

### Evaluating a trace

Use the topology `square`, diagonal and reference coordinate `X = coordinate(K)`
above. For a scalar field `u` on `square`, set
`plus_T = side_map(square, diagonal, "+")` and
`minus_T = side_map(square, diagonal, "-")`. The traces are
`pull_back(u, plus_T)` and `pull_back(u, minus_T)`.

From `"+"`, the parent is slot 0 and the inclusion is
$A_+(t) = (1 - t, t)$, through the local vertices 1 and 2 of `K`. From
`"-"`, the parent is slot 1 and the inclusion is
$A_-(t) = (0, t)$, through the local vertices 0 and 2. The coordinates of
`E` are $(1 - t, t)$.

Suppose the local fields on the two slots have affine reference values:

```python
u0 = a0[0] + a0[1] * X[0] + a0[2] * X[1]
u1 = a1[0] + a1[1] * X[0] + a1[2] * X[1]
```

Here `a0` and `a1` are separate local coefficient arrays and
`X = coordinate(K)`. Composing with the two inclusions gives:

$$
\begin{aligned}
u_+(t) &= a_{0,0} + a_{0,1}(1 - t) + a_{0,2}t, \\
u_-(t) &= a_{1,0} + a_{1,2}t.
\end{aligned}
$$

The two are mathematically equal for all $t$ if `u` is continuous at vertices 1 and 2.

Listing the diagonal as `(2, 1)` changes the inclusion from `"+"` to
$A_+(t) = (t, 1 - t)$. Both restrictions then run the other way along the
diagonal, leaving the jump at each point unchanged.

### From the facet integral to a kernel

`expand_integrals` expresses this single-facet contribution on the reference
interval, without a cell-index sum:

```python
facet_cell = cell_map(E, 0)
form(integral(pull_back(f, facet_cell) * density(facet_cell), I), (v, u))
```

The side declarations identify the local parents. Side data supply the maps
and normals. For example, the `"+"` trace lowers through:

```python
C_p = cell_map(square, 0)
A_p, m_p = diagonal.side_data(0, "+")
pull_back(pull_back(u, C_p), A_p)
```

The coordinate trace is `pull_back(pull_back(x, C_p), A_p)`; its Jacobian
gives the facet density. At a reference-facet point $t$, let $J_+$ be the
parent Jacobian evaluated at $A_+(t)$ and $\hat m_+$ the reference normal.

| Expression on `E` | Value at the facet point |
| ----------------- | ------------------------ |
| `u_p` | $\hat u_+(A_+(t))$ |
| `du_p` | $J_+^{-T}\widehat\nabla\hat u_+(A_+(t))$ |
| `n` | $J_+^{-T}\hat m_+ / \lvert J_+^{-T}\hat m_+\rvert$ |
| `alpha` | Supplied scalar local penalty |

The minus trace uses slot 1 and its own inclusion. The resulting form is
`form(integral(g, I), (v, u))`. The package selects a concrete interface,
packs both parents' local geometry and coefficient data, and calls this
kernel. Neighbour lookup and mesh iteration are outside the expression.

### Variations

- A computed cell size can use `B_p = cell_map(D, 0)` and
  `integral(density(B_p), K)`. Alternatively supply the parent's size as
  a scalar local input.
- Restricting values alone needs no physical mesh:
  `pull_back(u, side_map(square, diagonal, "+"))` is a trace on the topology.
- A triangle–quadrilateral interface uses parent types `(Kt, Kq)` and an
  interval integration cell. The same facet formula applies with the
  library's corresponding local bases and reference inclusions.

```{dropdown} Legacy UFL
This is the `dS` form
`-avg(dot(grad(u), n)) * jump(v) * dS - jump(u) * avg(dot(grad(v), n)) * dS + alpha / h * jump(u) * jump(v) * dS`.
Here the side maps make restriction explicit; jumps and averages are ordinary
arithmetic on the traces.
```

## A triangle–quadrilateral interface

Let `K` be the unit reference triangle, `I` the unit interval, and `Q` a
library reference quadrilateral with vertices `(0,0)`, `(1,0)`, `(1,1)`,
`(0,1)` in that order. The package supplies the local data for one parent
of each type. Their local configuration and common integration facet are:

```python
mixed = StaticTopology((K, Q), ((0, 1, 2), (1, 3, 4, 2)), ())
edge = StaticTopology((I,), ((1, 2),), ((("+", mixed, 0), ("-", mixed, 1)),))
A_plus, m_plus = edge.side_data(0, "+")
A_minus, m_minus = edge.side_data(0, "-")
```

The inclusions are $A_+(t)=(1-t,t)$ into `K` and $A_-(t)=(0,t)$ into `Q`.
For a library field `u` on `mixed`, its two traces are:

```python
u_plus = pull_back(pull_back(u, cell_map(mixed, 0)), A_plus)
u_minus = pull_back(pull_back(u, cell_map(mixed, 1)), A_minus)
jump = u_plus - u_minus
```

Both have domain `I`, although their parent bases and local coefficient
counts can differ. For instance, parent coordinate maps $x_+(r,s)=(r,s)$
and $x_-(r,s)=(1+r-s,r+s)$ pair both sides at $(1-t,t)$.
The SIPG construction above uses the corresponding physical local mesh and
side domain unchanged. A library constructs the required configurations from
one form recipe; the package dispatches concrete interfaces to their kernels.

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
square = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)), ())
D = mesh(x, RD(2))
```

Each local material cell has a side into its parent slot. Its side domain
uses the parent coordinate map:

```python
fluid_cells = StaticTopology((K,), ((0, 1, 2),), ((("inside", square, 0),),))
solid_cells = StaticTopology((K,), ((1, 3, 2),), ((("inside", square, 1),),))

F = side_domain(D, fluid_cells, "inside")
S = side_domain(D, solid_cells, "inside")

interface = StaticTopology(
    (I,), ((1, 2),), ((("fluid", fluid_cells, 0), ("solid", solid_cells, 0)),)
)
G = side_domain(F, interface, "fluid")

inlet = StaticTopology((I,), ((2, 0),), ((("fluid", fluid_cells, 0),),))
G_in = side_domain(F, inlet, "fluid")
```

`interface` names the fluid and solid topologies. Restricting a field from
`D` would require a third side `("parent", square, 0)`. The package uses markers and entity maps to select the concrete cells,
then supplies their data to this local configuration.

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

u_F, v_F, p_F, q_F = (push_forward(f, F, transformation_law(f.space)) for f in (u, v, p, q))
d_S, w_S, d_old_S = (push_forward(f, S, transformation_law(f.space)) for f in (d, w, d_old))
lam_G, mu_G = (
    push_forward(lam, G, transformation_law(L)),
    push_forward(mu, G, transformation_law(L)),
)

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
fluid = integral(nu * index_sum(grad_u[i, j] * grad_v[i, j], (i, j)) - p_F * div_v - q_F * div_u, F)

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

The fields on `F` and `S` cannot combine directly: their domains differ.
Their traces share domain `G`. The integrals have no domain and can be
added; tuple arguments bind functions on all three topologies. No field is
extended by zero onto the parent mesh.

The interface is the image of the coordinate trace from `F`.
`normal(from_fluid)` points out of the fluid, and `normal(from_solid)` is
its negative. Lowering follows the [facet example](#from-the-facet-integral-to-a-kernel):
the coupling kernel reads the fluid and solid reference bases on their
respective cells and the multiplier basis on the facet. The package chooses the cells to pack; the side declarations distinguish
their local roles.

### Variations

- Two marker values with the same cell type and sides use separate topologies;
  legacy `ds((3, 4))` is a sum of two integrals. Their kernels agree after
  renaming the topology placeholders. Different cell types use separate
  local configurations.
- A material coefficient varying by cell is supplied as local field data;
  a single-cell kernel can receive a scalar `nu` directly.
- For a moving fluid mesh, use `mesh(F.coordinate_function + m, RD(2))`,
  where the mesh displacement `m` on `fluid_cells` extends the solid motion.
  Functions on that topology push forward to the moved mesh unchanged.

```{dropdown} Legacy UFL
`dx(1)` and `ds(3)` become integrals over side domains of marked topologies.
The package uses entity maps to supply local data; side inclusions express
the reference-coordinate relation.
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
The library supplies `runtime_array` for solver data with no domain and
`runtime_extent` for sizes known at run time.

```python
D3 = mesh(x3, RD(3))
D1 = side_domain(D3, T1, "tissue")
identity = transformation_law((), 0)

p3, q3 = function(V3, "p3"), function(V3, "q3")
p1, q1 = function(V1, "p1"), function(V1, "q1")
beta = runtime_array("beta", ())
m = index("m")

grad_p3 = gradient(push_forward(p3, D3, transformation_law(V3)), D3)
grad_q3 = gradient(push_forward(q3, D3, transformation_law(V3)), D3)
grad_p1 = gradient(push_forward(p1, D1, transformation_law(V1)), D1)
grad_q1 = gradient(push_forward(q1, D1, transformation_law(V1)), D1)
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
combine on their different domains. Restricting `p3` to `T1` makes the
exchange possible; pushing it forward then gives a field on `D1`.
The three integrals have no domain and add into one form.

For one segment kernel, the package supplies a containing tetrahedron and
its local data. The side declaration names that parent slot, and its side
data give the reference-edge inclusion. Continuity of `p3` makes the trace
independent of which containing tetrahedron the package chooses.

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
average tissue pressure over a ring around each network point. The package
constructs the ring quadrature, locates its points in concrete tissue cells,
and supplies their reference positions. A local sampling kernel evaluates:

```python
R = domain_map(RD(0), positions[l], K3)
p_samples = pull_back(pull_back(p3, cell_map(T3, 0)), R)
weighted_samples = weights[l] * p_samples
```

Here `T3` has one tetrahedron slot with reference cell `K3`; `positions` and
`weights` are local arrays for points assigned to that cell. The package
combines the sampled contributions into the network average. If an entire
fixed neighbourhood is represented in one local configuration, its separate
parent slots can be read explicitly instead. A physical point alone does
not select a cell, and core lowering performs no point-location search.

## Side topology examples

Use reference cells `K = ReferenceSimplex(2)`, `I = ReferenceSimplex(1)`
and `P = ReferenceSimplex(0)`. Each example describes a local configuration.

- A boundary facet has types `(I,)` and sides `((("inside", T, 0),),)`
  into a one-triangle configuration `T`.
- An interior facet has sides `((("+", T, 0), ("-", T, 1)),)` into
  a two-triangle configuration `T` with types `(K, K)`.
- A corner shared by two boundary edges has types `(P,)` and sides
  `((("+", B, 0), ("-", B, 1)),)`, where `B` contains two interval slots.
- A local boundary submesh `M` can be paired with a parent through
  `((("mesh", T, 0), ("submesh", M, 0)),)`.
- A periodic interface names two parent slots and supplies the two reference
  inclusions. Their physical images may be separated by a period.
- A triangle–quadrilateral pair `M` has types `(K, Q)`. Its interface has
  types `(I,)` and sides `((("+", M, 0), ("-", M, 1)),)`. A library
  builds the other required pair configurations from the same form recipe.
- A prism configuration `U3` with types `(W,)` can have a local boundary
  containing one triangular and one quadrilateral facet, with types `(K, Q)`
  and sides `((("inside", U3, 0),), (("inside", U3, 0),))`. Each slot has
  its own inclusion; integrating this configuration adds both contributions.

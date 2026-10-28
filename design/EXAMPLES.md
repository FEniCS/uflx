# Worked examples for the UFLx core

These examples accompany [LANGUAGE.md](LANGUAGE.md) and use its names. The
first part defines the library classes that the examples use; none of them
is part of the core. A reader of LANGUAGE.md
needs none of this to follow the rules.

## Reference library

### The unit simplex

Requires: [`AbstractCoordinateDomain`](LANGUAGE.md#abstractcoordinatedomain).

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

> LLM generated content.

Requires: [`AbstractCellType`](LANGUAGE.md#abstractcelltype).

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

> LLM generated content.

Requires: [`AbstractTopology`](LANGUAGE.md#abstracttopology),
[The simplex as a cell type](#the-simplex-as-a-cell-type),
[Sides of a topology](LANGUAGE.md#sides-of-a-topology).

`StaticTopology(cell_types, cells, sides)` is a topology whose cells are known
during construction. Entry $j$ of `cells` is a `tuple` with the vertices of
each cell of type $j$, as vertex numbers of the topology. Local vertex $i$ of
cell $c$ of type $j$ is the vertex `cells[j][c][i]` of the topology. A
`StaticTopology` with sides has one topology on every side, uses the vertex
numbers of that topology, pairs each vertex with the same vertex, and computes
the side tables from the vertex numbers: the cell on the side, and the
positions of the vertices of each cell as the coordinates `vertices` of the
paired local vertices of the reference cell on the side. So a `StaticTopology`
covers sub-entities only, not pieces of a cell. For a cell seen from two cells of the
same group, the cell with the lower number is on the first side, and otherwise
the group of each side decides. For a sub-entity of lower dimension than a
facet, the cell with the lowest number that contains the sub-entity is chosen.

> `StaticTopology` is mainly for testing UFLx ideas, but also covers Alfeld
> splits.

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
# A square divided into two triangles, with the numbers of cells (2,).
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())

# The Alfeld split of a triangle at vertex 3, with the numbers of cells (3,).
alfeld = StaticTopology((K,), (((0, 1, 3), (1, 2, 3), (2, 0, 3)),), ())

# The facets of the boundary of the square, each seen from one triangle.
boundary = StaticTopology(
    (I,), (((0, 1), (1, 3), (3, 2), (2, 0)),), ((("inside", square, 0),),)
)

# The facet of the diagonal of the square, seen from two triangles.
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
```

A facet of a cell is a set of vertex numbers of the topology, found through the
facets of `K`. In `square`, the facet $\{1, 2\}$ is a facet of two cells and
the other four facets are facets of one cell. In `alfeld`, the three facets
that contain vertex 3 are facets of two cells, and the three other facets are
facets of one cell. `boundary` and `diagonal` have the numbers of cells `(4,)`
and `(1,)`.

`StaticTopology` does not check the sides. A `StaticTopology` has no
coordinates. The positions of the vertices come from the
[coordinates of a decomposed domain](LANGUAGE.md#decomposeddomain).

### The unit simplex as a reference cell

> LLM generated content.

Requires: [The unit simplex](#the-unit-simplex),
[The simplex as a cell type](#the-simplex-as-a-cell-type),
[Literals](LANGUAGE.md#literals), [`Tensor`](LANGUAGE.md#tensor), [Reference cell](LANGUAGE.md#reference-cell).

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

`ReferenceSimplex(2)` and `Simplex(2)` are different cell types, so a topology
of `Simplex(2)` cells has no reference cells. `ReferenceSimplex(2)` and
`UnitSimplex(2)` are the same set but different coordinate domains, so an
expression on one is not an expression on the other. Likewise,
`ReferenceSimplex(0)` and `RD(0)` are different coordinate domains.

## Examples

### Transpose

Requires: [`Index`](LANGUAGE.md#index), [`Indexed`](LANGUAGE.md#indexed), [`ComponentTensor` and
`IndexSum`](LANGUAGE.md#componenttensor-and-indexsum).

See also: [Tensor arithmetic](LANGUAGE.md#tensor-arithmetic).

`transpose` is a free function of the core (see [Tensor
arithmetic](LANGUAGE.md#tensor-arithmetic)), built with derived indices (see
[Derived indices](LANGUAGE.md#derived-indices)). `transpose(A)` is
`component_tensor(A[k, l], (l, k))`, with
$k$ = `Index(("uflx.transpose", A, 0))` and
$l$ = `Index(("uflx.transpose", A, 1))`.

An operand that does not have rank 2, including a scalar operand, is rejected
with `IndexError` by the constructor of `Indexed`. The matrix product is the
free function `matrix_product` of the core, with the overload `@` (see [Tensor
arithmetic](LANGUAGE.md#tensor-arithmetic)).

Examples, with a matrix `A` of value shape `(2, 2)`:

- `transpose(A)` is not equal to `component_tensor(A[i, j], (j, i))` for user
  indices `i` and `j`, because the index labels differ.

### Tensor algebra with `EllipsisIndex`

Requires: [`EllipsisIndex`](LANGUAGE.md#ellipsisindex), [`Indexed`](LANGUAGE.md#indexed), [`ComponentTensor`
and `IndexSum`](LANGUAGE.md#componenttensor-and-indexsum), [`Sum`](LANGUAGE.md#sum), [`Product`](LANGUAGE.md#product),
[`Conjugate`](LANGUAGE.md#conjugate).

See also: [Derived indices](LANGUAGE.md#derived-indices), [New
operations](LANGUAGE.md#new-operations), [Transpose](#transpose), [Tensor
arithmetic](LANGUAGE.md#tensor-arithmetic).

An `EllipsisIndex` writes an operation on tensors of every rank as one fixed
composition. The scaling `a * T` and the sum `A + B` of tensors are
overloads of the core, built with one index for each axis (see [Tensor
arithmetic](LANGUAGE.md#tensor-arithmetic)). An operation of the table does
not reduce when the `EllipsisIndex` meets a binder with one index for each
axis, such as `inner(transpose(A), B)` (see [Derived
indices](LANGUAGE.md#derived-indices)). The operations in the table below
are not part of the core. A library may provide the operations as free
functions (see [New operations](LANGUAGE.md#new-operations)), as for
`transpose` above. With derived ellipsis indices `r` and `s` and an ordinary
index `k`:

| Operation          | Fixed composition                                              |
| ------------------ | -------------------------------------------------------------- |
| `outer(u, v)`      | `component_tensor(u[r] * v[s], (r, s))`                        |
| `dot(u, v)`        | `component_tensor(index_sum(u[r, k] * v[k, s], (k,)), (r, s))` |
| `inner(A, B)`      | `index_sum(A[r] * conjugate(B[r]), (r,))`                      |

For operands of rank 1 and 2, `dot(u, v)` contracts the same axes as `u @ v`
with a different composition, so `dot(u, v)` and `u @ v` are not equal.

### An affinely mapped triangle

Requires: [The unit simplex](#the-unit-simplex),
[`Coordinate`](LANGUAGE.md#coordinate), [`Tensor`](LANGUAGE.md#tensor),
[`DomainMap`](LANGUAGE.md#domainmap), [Jacobian and
density](LANGUAGE.md#jacobian-and-density).

The physical triangle with vertices $v_0 = (1, 0)$, $v_1 = (3, 1)$ and
$v_2 = (2, 3)$ is the image of the unit triangle under an affine map.

The affine map from `UnitSimplex(2)` onto the physical triangle is

$$\psi(X) = v_0 + X_0 (v_1 - v_0) + X_1 (v_2 - v_0)
= (1 + 2 X_0 + X_1, X_0 + 3 X_1).$$

```python
# AbstractCoordinateDomain: UnitSimplex(2), no edges in the graph.
T = UnitSimplex(2)

# AbstractExpression: Coordinate(UnitSimplex(2)), value shape (2,).
X = coordinate(T)

# AbstractExpression: Tensor(..., ...), value shape (2,).
one, two, three = Integer(1), Integer(2), Integer(3)
psi = tensor((one + two * X[0] + X[1], X[0] + three * X[1]))

# DomainMap(UnitSimplex(2), psi, RD(2)), k = n = 2.
K = domain_map(T, psi, RD(2))
```

The Jacobian matrix of $\psi$ is the constant matrix

$$J = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix},$$

so $\det J = 5$, and the area of the physical triangle is $5$ times the area
$1/2$ of the unit triangle. `jacobian(K)` evaluates to $J$,
`determinant(jacobian(K))` to $5$, and `density(K)` to $|\det J| = 5$.
Swapping $v_1$ and $v_2$ gives $\det J = -5$ and the same density.

### Circulation along an edge

Requires: [Concept: integral over a domain](LANGUAGE.md#concept-integral-over-a-domain),
[Jacobian and density](LANGUAGE.md#jacobian-and-density), [An affinely
mapped triangle](#an-affinely-mapped-triangle).

The examples use `T`, `X`, `one`, `S`, `t` and `e` of the
[`PullBack`](LANGUAGE.md#pullback) example.

The edge `e`, with $\psi(t) = (1 - t, t)$, runs from vertex 1 to vertex 2 of
the unit triangle. The edge with $\psi(t) = (t, 1 - t)$ is the same set, run
the other way. The circulation of a vector field along an edge depends on the
direction, and the integral of a function does not:

```python
x, y = X[0], X[1]
l = index("l")

# The vector field F = (0, x) on T.
F = tensor((Integer(0), x))

# The edge e run the other way.
e_reversed = domain_map(S, tensor((t, one - t)), T)

def circulation(edge):
    return integral(index_sum(pull_back(F, edge)[l] * jacobian(edge)[l, 0], (l,)), S)

# A later pass evaluates the two circulations to 1/2 and -1/2.
circulation(e)
circulation(e_reversed)

# A later pass evaluates both integrals of x to sqrt(2) / 2.
integral(pull_back(x, e) * density(e), S)
integral(pull_back(x, e_reversed) * density(e_reversed), S)
```

On `e`, $F \circ \psi = (0, 1 - t)$ and $\psi' = (-1, 1)$, so the
circulation is $\int_0^1 (1 - t) \, \mathrm{d}t = 1/2$. On `e_reversed`,
$F \circ \psi = (0, t)$ and $\psi' = (1, -1)$, so the circulation is
$\int_0^1 -t \, \mathrm{d}t = -1/2$. The density is $|\psi'| = \sqrt 2$ on
both.

### Integrating over an edge

Requires: [`PullBack`](LANGUAGE.md#pullback), [`Integral`](LANGUAGE.md#integral),
[Jacobian and density](LANGUAGE.md#jacobian-and-density), [Circulation along
an edge](#circulation-along-an-edge).

The integral over the edge `e` of the [`PullBack`](LANGUAGE.md#pullback)
example is the integral over `S` of the pull-back by `e`, times the density.
The examples use `T`, `X`, `one`, `S`, `t` and `e` of that example, and `x`
and `y` of the [circulation example](#circulation-along-an-edge).

The flux of a vector field $u$ through the edge is the integral of $u$ against
the normal $(\psi'_1, -\psi'_0)$, whose length is the density. For $u = X$:

```python
J = jacobian(e)

# A later pass evaluates the flux to 1.
integral(pull_back(x, e) * J[1, 0] - pull_back(y, e) * J[0, 0], S)
```

The flux is 1. The divergence of $X$ is 2 and the area of `T` is $1/2$, so
the flux of $X$ out of `T` is 1. $X$ is tangent to the other two edges of
`T`, so the whole flux passes through `e`. The equality of the two fluxes is
the divergence theorem.

The integral of a scalar function with respect to arc length uses the
density $|\psi'|$:

```python
# A later pass evaluates the integral to the length sqrt(2) of e.
integral(density(e), S)

# A later pass evaluates the integral to sqrt(2) / 6.
integral(pull_back(x * y, e) * density(e), S)
```

`density(e)` is `square_root(determinant(G))`, with `G` the $1 \times 1$
matrix $J^T J$ (see [Jacobian and density](LANGUAGE.md#jacobian-and-density)).
For a decomposed domain, `expand_integrals` builds the same density (see
[`Integral` over a decomposed
domain](LANGUAGE.md#integral-over-a-decomposed-domain)).

### Integrating over the boundary of a triangle

Requires: [`DomainMap`](LANGUAGE.md#domainmap),
[`ComponentTensor` and `IndexSum`](LANGUAGE.md#componenttensor-and-indexsum),
[Circulation along an edge](#circulation-along-an-edge).

The three edges of `T` are one family of domains, with a free index. The
example uses the names of the [circulation
example](#circulation-along-an-edge).

```python
zero = Integer(0)
i, k = index("i"), index("k")

# The start and the end vertex of each edge, counterclockwise.
starts = tensor((tensor((zero, zero)), tensor((one, zero)), tensor((zero, one))))
ends = tensor((tensor((one, zero)), tensor((zero, one)), tensor((zero, zero))))

# The family of the three edges of T, with the free index i.
psi_edges = component_tensor(starts[i, k] + t * (ends[i, k] - starts[i, k]), (k,))
edges = domain_map(S, psi_edges, T)

# A later pass evaluates the integral to 1/2.
index_sum(circulation(edges), (i,))
```

The result is the circulation of $F = (0, x)$ around the boundary of `T`,
which is $\int_T (\partial_x F_1 - \partial_y F_0) = 1/2$ by Green's
theorem. Without the `IndexSum`, the integral has the free index `i`, and the
value for each `i` is the circulation along one edge.

The direction of each edge comes from the order of `starts` and `ends`. Every
domain of a family has the same parameter domain, here `S`.

### Quadrature

Requires: [Points and point evaluation](LANGUAGE.md#points-and-point-evaluation),
[`DomainMap`](LANGUAGE.md#domainmap), [`PullBack`](LANGUAGE.md#pullback),
[`Integral`](LANGUAGE.md#integral), [`Tensor`](LANGUAGE.md#tensor),
[`ComponentTensor` and `IndexSum`](LANGUAGE.md#componenttensor-and-indexsum),
[Extents](LANGUAGE.md#extents).

A quadrature rule replaces the integral by the inner product of the weights
with the values at the points. With `T`, `X` and `f = X[0] * X[1]` of
[Points and point evaluation](LANGUAGE.md#points-and-point-evaluation), the edge midpoint
rule is:

```python
zero, half = Integer(0), Integer(1) / Integer(2)
sixth = Integer(1) / Integer(6)
i = index("i")

points = tensor((tensor((half, zero)), tensor((half, half)), tensor((zero, half))))
weights = tensor((sixth, sixth, sixth))

# A family of three points in T, with the free index i.
xs = domain_map(RD(0), points[i], T)

q = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

The points are one domain with the free index `i`, so `f` occurs once in `q`,
and a form compiler can evaluate `q` as a loop over the points.

A later pass evaluates `q` to $1/24$, which is `integral(f, T)`, because the
rule is exact for $f$. The two expressions are not equal, so quadrature is a
pass, and the rule is an argument of the pass.

After the passes, `q` is

```python
index_sum(weights[i] * (points[i, 0] * points[i, 1]), (i,))
```

`expand_pull_backs` replaces `X` by the parameterisation of `xs`. A pass
replaces the integral over `RD(0)` of an expression by the expression, and
`reduce_component_tensors` replaces `points[i][0]`, the entry 0 of a
`ComponentTensor` (see [Partial indexing](LANGUAGE.md#partial-indexing)), by
`points[i, 0]`. A form compiler could generate:

```python
weights = (1 / 6, 1 / 6, 1 / 6)
points = ((1 / 2, 0.0), (1 / 2, 1 / 2), (0.0, 1 / 2))

def kernel() -> float:
    q = 0.0
    for i in range(3):
        q += weights[i] * (points[i][0] * points[i][1])
    return q
```

The `IndexSum` over `i` is the loop over the points, and the literals are
tables.

With `points` and `weights` given by leaves of a
[language extension](LANGUAGE.md#language-extensions), of value shapes `(q, 2)` and
`(q,)` for an `AbstractExtent` `q`, the same expression is runtime quadrature
(see [Integrating over a cut
cell](#integrating-over-a-cut-cell)).
The number of points is known only at assembly, and the loop bound `q` is an
input of the kernel.

### Integrating over a cut cell

Requires: [Quadrature](#quadrature), [Extents](LANGUAGE.md#extents).

In CutFEM, a level set $\phi$ cuts a cell, and an integral covers only the
part of the cell where $\phi < 0$. The solver computes a quadrature rule for
that part of each cut cell, so the number of points differs from cell to cell.

> The points and the weights are leaves of a small
> [language extension](LANGUAGE.md#language-extensions) for kernel inputs:
> `RuntimeExtent(label)`, an `AbstractExtent`,
> `RuntimeArray(label, value_shape)`, an array such as quadrature data, the
> coefficients of a function on a cell, or a constant, and
> `RuntimeIndex(label, extent)`, an `AbstractIndexValue` of the shape `()`,
> such as the local number of a facet, read as
> `gathered_index(RuntimeIndex(label, extent), ())`.

With `T`, `f` and `i` of the quadrature example:

```python
# These are language extensions, not the core.
q = runtime_extent("q")

points = runtime_array("points", (q, 2))
weights = runtime_array("weights", (q,))

# The family of q points in the part of T where phi < 0. The level set phi
# is not written here, because the core has no conditional.
xs = domain_map(RD(0), points[i], T)

a = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

For all the cells of a topology, see the open item on ragged indexing in
[Extents](LANGUAGE.md#extents).

### The centroid of a triangle

> LLM came up with this example.

Requires: [An affinely mapped
triangle](#an-affinely-mapped-triangle), [`PullBack`](LANGUAGE.md#pullback),
[`Integral`](LANGUAGE.md#integral), [`Division`](LANGUAGE.md#division),
[`ComponentTensor` and `IndexSum`](LANGUAGE.md#componenttensor-and-indexsum),
[Points and point evaluation](LANGUAGE.md#points-and-point-evaluation).

An integral is an expression, so the centroid of `K` is a ratio of integrals.
The example uses `T`, `psi` and `K` of the affinely mapped triangle.

```python
Y = coordinate(RD(2))
k = index("k")

# A later pass evaluates the area of K to 5/2.
area = integral(density(K), T)

# A later pass evaluates the centroid of K to (2, 4/3).
centroid = component_tensor(integral(pull_back(Y[k], K) * density(K), T) / area, (k,))
```

`density(K)` is $|\det J| = 5$ (see [Jacobian and
density](LANGUAGE.md#jacobian-and-density)). The centroid is also `pull_back(psi, c)`, the image of the
centroid `c` of `T` of [Points and point
evaluation](LANGUAGE.md#points-and-point-evaluation). The two expressions are not
equal, and a later pass evaluates both to $(2, 4/3)$.

Over a decomposed domain `D`, the centroid is one integral of a vector:

```python
centroid_D = integral(coordinate(D), D) / integral(Integer(1), D)
```

[`expand_integrals`](LANGUAGE.md#expanding-integrals) splits the integrand
into its components, and each component is a scalar integral over the cells.

### An interior facet integral

Requires: [Sides of a topology](LANGUAGE.md#sides-of-a-topology),
[Restriction](LANGUAGE.md#restriction), [`Integral` over a decomposed
domain](LANGUAGE.md#integral-over-a-decomposed-domain), [`Normal`](LANGUAGE.md#normal),
[`Form`](LANGUAGE.md#form), [Evaluation](LANGUAGE.md#evaluation).

The interior facet terms of the symmetric interior penalty method for the
Poisson problem, the legacy UFL `dS` integral, on the square divided into two
triangles. The example uses `D`, `square`, `K`, `X`, `y` and `cell_vertices`
of the [`DecomposedDomain`](LANGUAGE.md#decomposeddomain) example, `I = ReferenceSimplex(1)`
and the `diagonal` topology of the [Restriction](LANGUAGE.md#restriction) example, a
space `V` on `D` of the value shape `()`, and functions `u` and `v` in
`V`. `runtime_array` is the language extension of the [cut cell example](#integrating-over-a-cut-cell).

```python
# The interior facets, seen from the sides "+" and "-", with the free index e.
E = side_domain(D, diagonal)
e = cell_index(diagonal, 0)

# The gradients on D, in the frame of RD(2).
i = index("i")
grad_u = gradient(u, D)
grad_v = gradient(v, D)

# The restrictions from each side, expressions on E.
u_p, u_m = restriction(u, E, "+"), restriction(u, E, "-")
v_p, v_m = restriction(v, E, "+"), restriction(v, E, "-")
du_p, du_m = restriction(grad_u, E, "+"), restriction(grad_u, E, "-")
dv_p, dv_m = restriction(grad_v, E, "+"), restriction(grad_v, E, "-")

# The normal out of the cell on the side "+", and the jumps and averages.
n = normal(D, E, "+")
half = Integer(1) / Integer(2)
jump_u, jump_v = u_p - u_m, v_p - v_m
avg_dn_u = half * index_sum((du_p[i] + du_m[i]) * n[i], (i,))
avg_dn_v = half * index_sum((dv_p[i] + dv_m[i]) * n[i], (i,))

# The penalty parameter of each facet, an array of the solver.
alpha = runtime_array("alpha", (diagonal.num_cells[0],))

# The facet terms of the symmetric interior penalty method.
f = -avg_dn_u * jump_v - jump_u * avg_dn_v + alpha[e] * jump_u * jump_v
a_facet = form(integral(f, E), (v, u))
```

- `restriction(u, E, "+")` is `pull_back(u, side_map(D, E, "+"))`, an
  expression with the support `E` and no free index. `jump_u` is the
  difference of two expressions on `E`, and `f` has the support `E` and the
  free index `e`, from `alpha[e]`.
- `n` points out of the cell on the side `"+"`, the triangle `(0, 1, 2)`, so
  `jump_u * n` is the jump of `u` with a fixed sign convention. Swapping the
  two sides in `diagonal` negates `n` and every jump, and leaves `f`
  unchanged, because every term of `f` has two jumps or a jump and a normal.
- `alpha[e]` stands for the legacy $\alpha / h$. The array is a choice, not a
  limitation. A size built from the two cells is a composition with the free
  index `e`: the measure of the cell on the side `"+"` is
  `integral(density(B_p), K)`, with `B_p` the map $B_+$ of
  [Restriction](LANGUAGE.md#restriction) (see [Jacobian and
  density](LANGUAGE.md#jacobian-and-density)), and a per-cell array `h` of the
  solver is read as `h[gathered_index(cells_p, (e,))]`.

The draft passes rewrite `a_facet` as follows. With $B_\pm$ =
`gathered_cell_domain(D, 0, gathered_index(cells_pm, (e,)))` the cell on each
side and $A_\pm$ = `domain_map(I, affine_map(I, points_pm), K)` the inclusion
of the facet into that cell, both of [Restriction](LANGUAGE.md#restriction), with the side
tables `(cells_pm, positions_pm)` of each side:

1. `expand_integrals` rewrites `integral(f, E)` as
   `index_sum(integral(pull_back(f, cell_domain(E, 0)) * density(cell_domain(E, 0)), I), (e,))`,
   with the density the length of the facet, from the Jacobian of the
   parameterisation `pull_back(pull_back(coordinate(D), B_+), A_+)` of the
   piece.
2. `expand_pull_backs` moves the pull-back by `cell_domain(E, 0)` into `f`.
   On `u_p` the rule for `SideMap` gives `pull_back(pull_back(u, B_+), A_+)`.
   The inner pull-back remains, for the library of `V`, as `u` on the cell
   `cells_p[e]` (see [`AbstractFunctionSpace`](LANGUAGE.md#abstractfunctionspace)), and
   the outer pull-back remains outside the inner one, by the rule that
   pull-backs nest.
   On `du_p` the same rule gives the pull-back of `gradient(u, D)` by
   $B_+$, which the rule for `Gradient` rewrites as
   $K^T \hat g$, with $K$ = `inverse(J_+)` the inverse Jacobian of the cell and
   $\hat g$ = `gradient(pull_back(u, B_+), K)` the
   reference gradient, and then pulled back by $A_+$, which remains on
   $\hat g$ and on $J_+$, because `I` has lower dimension than `K`: the
   reference gradient and the cell Jacobian are read at the points of the
   facet. On `n` the rule for `Normal` gives $K^T \hat m / \lvert K^T \hat m \rvert$,
   with $\hat m$ the normal of the edge `(1, 2)` of `K` pointing away from the
   centroid of `K`. On `alpha[e]`, an expression with no support whose free
   index `e` is free in `cell_domain(E, 0)`, the pull-back is `alpha[e]`.
3. The gradients pass leaves `gradient(pull_back(u, B_+), K)` on
   the function, and `reduce_component_tensors` reduces the gathers of the
   vertex positions inside $J_+$ to `Indexed`s with gathered indices, nested
   through the cell on the side, and `points_pm` to `Indexed`s of
   `positions_pm`.

After the passes, `a_facet` is `form(index_sum(integral(g, I), (e,)), (v, u))`
with `g` an expression on `I` with the free index `e`. The `IndexSum` over `e`
is the loop of the assembler over the interior facets. The gathered indices in
`g` read, for the facet `e`, the two cells `cells_p[e]` and `cells_m[e]`, the
vertex positions of each through `cell_vertices`, the positions
`positions_p[e, v, k]` and `positions_m[e, v, k]`, and `alpha[e]`, which the assembler
packs for a kernel. The kernel is `g` without `e`: an expression on the
reference facet `I` that reads the reference basis of `u` and `v` on the two
cells, through the library of `V`, at the points of the facet in each cell's
reference coordinates $A_\pm(t)$, and the two cell Jacobians at those points.
This is the interface of a legacy interior facet kernel, with the local facet
numbers replaced by the positions of the facet vertices, which also fix the
direction of the facet from each side.

- In the example, $A_+(t) = (1 - t, t)$ and $A_-(t) = (0, t)$ (see
  [Restriction](LANGUAGE.md#restriction)), $J_+$ is the identity and $J_-$ is
  $\begin{pmatrix} 0 & -1 \\ 1 & 1 \end{pmatrix}$, and $\hat m = (1, 1) / \sqrt 2$,
  so `n` evaluates to $(1, 1) / \sqrt 2$ on the diagonal, pointing from
  triangle 0 into triangle 1, and the density to $\sqrt 2$.
- `restriction(u, E, "-")` on the same `t` reads `u` on triangle 1 at
  $(0, t)$, the same physical point $(1 - t, t)$ as from `"+"`, so `jump_u`
  is a difference at one point.

```{dropdown} Legacy UFL
The legacy UFL form is `-avg(dot(grad(u), n)) * jump(v) * dS - jump(u) * avg(dot(grad(v), n)) * dS + alpha / h * jump(u) * jump(v) * dS`, with `n` the normal from `'+'` and `h` a cell size. In UFLx the side names, the normal and the facet domain are explicit, `avg` and `jump` are ordinary arithmetic on restrictions, and the size is an array of the solver or a composition of integrals over the cells.
```

### A Bank-Weiser error estimator

Requires: [An interior facet
integral](#an-interior-facet-integral), [Restriction](LANGUAGE.md#restriction),
[`Normal`](LANGUAGE.md#normal), [`Form`](LANGUAGE.md#form).

The local problems of the Bank-Weiser estimator for the Poisson problem, as
in Bulle et al. On each cell, with `u_h` the Galerkin solution, the local
problem seeks `e` in a space of quadratic functions with zero vertex values
such that

$$
\begin{aligned}
\int_K \nabla e \cdot \nabla v ={} & \int_K (f + \Delta u_h) \, v \\
& - \frac{1}{2} \int_{\partial K \cap \Omega} n_K \cdot (\nabla u_h|_K - \nabla u_h|_{K'}) \, v \\
& + \int_{\partial K \cap \Gamma_N} (g - \partial_n u_h) \, v,
\end{aligned}
$$

and the estimator on the cell is the energy norm of `e`. The facet terms need
the solution on every cell attached to `K`. The example writes the local
problems as one global form on a discontinuous space, so the local matrix and
vector of each cell are the blocks of the global matrix and vector, and an
assembler may form them cell by cell without a global array.

The example uses `D`, `square`, `K`, `I` and the interior facet domain `E`
with the sides `"+"` and `"-"` of the [interior facet
example](#an-interior-facet-integral), the boundary facet domain `B`
with the side `"inside"`, a Lagrange space `V` of polynomial degree 1 on `D`,
a space `F` on `D` for the data, and the Bank-Weiser space `W` on `D`, a
discontinuous space whose pull-back by each cell is the space of quadratic
functions on `K` with zero vertex values, a representation of the library:

```python
u_h, f, g = function(V, "u_h"), function(F, "f"), function(F, "g")
e, v = function(W, "e"), function(W, "v")
i = index("i")
half = Integer(1) / Integer(2)

# The local bilinear form, one block for each cell because W is discontinuous.
grad_e = gradient(e, D)
grad_v = gradient(v, D)
a = form(integral(index_sum(grad_e[i] * grad_v[i], (i,)), D), (v, e))

# The interior residual. The Laplacian of u_h is zero for degree 1, and is kept
# for higher degrees.
grad_u = gradient(u_h, D)
lap_u = index_sum(gradient(grad_u, D)[i, i], (i,))
r_cell = integral((f + lap_u) * v, D)

# The jump of the normal derivative with the normal of "+". From "-" the normal
# and the jump both change sign, so the same scalar serves both cells.
n = normal(D, E, "+")
jump_n = index_sum(
    (restriction(grad_u, E, "+")[i] - restriction(grad_u, E, "-")[i]) * n[i], (i,)
)
r_facets = integral(
    -half * jump_n * (restriction(v, E, "+") + restriction(v, E, "-")), E
)

# The Neumann residual on the boundary.
n_B = normal(D, B, "inside")
dn_u = index_sum(restriction(grad_u, B, "inside")[i] * n_B[i], (i,))
r_boundary = integral((g - dn_u) * restriction(v, B, "inside"), B)

r = form(r_cell + r_facets + r_boundary, (v,))
```

- `W` is discontinuous, so the matrix of `a` is block diagonal with one block
  for each cell, the local matrix of the cell. The vector of `r` holds, in the
  block of each cell, the interior residual, one half-jump term for each
  interior facet of the cell and the Neumann residual for each boundary facet
  of the cell. These are the local matrix and the local vector of the
  Bank-Weiser problem of the cell.
- `jump_n` is the same expression from both sides. For the cell on `"-"`, the
  outward normal is `-n` and the jump of the gradient changes sign, so the
  product does not. One `Integral` over `E` therefore writes the facet term of
  both cells, and `v` is restricted from both sides so that each cell receives
  its half.
- A Dirichlet facet contributes nothing, so `B` here is the topology of the
  Neumann facets only. The library decides which facets are in `B`.
- The core does not represent the basis of `W`. The library may represent a
  function in `W` on `K` through a quadratic Lagrange basis and a change of
  basis to the functions with zero vertex values, as Bulle et al. do. That
  choice stays behind the pull-back and is invisible to the form.

The draft passes rewrite `a` and `r` as in the [interior facet
example](#an-interior-facet-integral). After the passes `a` is one
`IndexSum` over the cells, `r_cell` and `r_boundary` are an `IndexSum` over
the cells and over the boundary facets, and `r_facets` is an `IndexSum` over
the interior facets of an expression on `I` that reads `u_h` on the two cells
through the side tables. The kernels are:

- the cell kernel of `a`, which returns the local matrix;
- the cell kernel of `r_cell`, which returns the interior residual of the cell;
- the facet kernel of `r_facets`, which reads the two cells, evaluates
  `jump_n` once and returns the contribution to the rows of both cells,
  because `v` appears from both sides;
- the facet kernel of `r_boundary`, which returns the rows of the one cell.

Where the results go is outside the language. A global assembler scatters
them into a global matrix and vector, and the blocks are the local systems. A
local assembler instead loops over the cells and, for each cell, calls the
cell kernels, calls the facet kernel for each facet of the cell and keeps the
rows of the cell, solves the local problem and computes the estimator, with
no global array. That is the loop of the custom assembler of Bulle et al.,
with the two kernels of the global form in place of hand-written code.

- With a topology of the pairs (cell, local facet) in place of `E`, with the
  sides `"self"` and `"other"` and `v` restricted from `"self"` only, the
  facet kernel returns the rows of the one cell and the pairs come packed by
  cell. The form is the same with `restriction(v, E, "self")` in place of the
  sum of the two restrictions.
- The per-cell estimator $\eta_K^2 = \int_K \lvert \nabla e \rvert^2$ is
  `integral(pull_back(index_sum(grad_e[i] * grad_e[i], (i,)), cell_domain(D, 0)) * density(cell_domain(D, 0)), K)`,
  with the density of the cell (see [Jacobian and
  density](LANGUAGE.md#jacobian-and-density)), an expression with no support
  and the cell index free: one value for each cell.

### A fluid and a solid

Requires: [An interior facet
integral](#an-interior-facet-integral), [Sides of a
topology](LANGUAGE.md#sides-of-a-topology), [Restriction](LANGUAGE.md#restriction),
[`Form`](LANGUAGE.md#form).

A fluid and a solid that share an interface, with the fields of each part in
spaces on that part only: the legacy UFL `dx(1)`, `dx(2)` and `ds(3)` with
markers, and the submeshes and entity maps of a solver such as DOLFINx. The
square of the [`DecomposedDomain`](LANGUAGE.md#decomposeddomain) example is the
whole, triangle `(0, 1, 2)` is the fluid, triangle `(1, 3, 2)` is the solid,
and the diagonal is the interface, the smallest mesh with every part of the
problem. The example uses `D`, `square`, `K` and
`I` of the [interior facet example](#an-interior-facet-integral), and
`runtime_array` of the [cut cell
example](#integrating-over-a-cut-cell).

A marked set of entities is a topology with sides, and a marker value is the
identity of the topology. A submesh is the decomposed domain of such a
topology, and the entity map of the submesh is the side table into the
parent:

```python
# dx(1) and dx(2): the fluid cells and the solid cells, each seen as itself in square.
fluid_cells = StaticTopology((K,), (((0, 1, 2),),), ((("inside", square, 0),),))
solid_cells = StaticTopology((K,), (((1, 3, 2),),), ((("inside", square, 0),),))

# The fluid and the solid, each a decomposed domain with its own topology.
F = side_domain(D, fluid_cells)
S = side_domain(D, solid_cells)

# The interface, seen from a fluid cell and from a solid cell.
interface = StaticTopology(
    (I,), (((1, 2),),), ((("fluid", fluid_cells, 0), ("solid", solid_cells, 0)),)
)
G = side_domain(F, interface)

# ds(3) on the fluid: the left edge, seen from a fluid cell.
inlet = StaticTopology((I,), (((2, 0),),), ((("fluid", fluid_cells, 0),),))
G_in = side_domain(F, inlet)
```

The sides of `interface` name the two submesh topologies, not `square`. A
side is a name and two tables, so a facet may have more sides than
neighbours: a field on `D` restricted to the interface needs a third side
`("parent", square, 0)`, because a restriction needs a side with the
topology of the support of the field.

On the mesh of a solver, the four `StaticTopology`s are the only lines that
change. The library provides `fluid_cells`, `solid_cells`, `interface` and
`inlet` as topologies with `AbstractExtent`s as the numbers of cells and the
side tables as arrays of the solver, from the cell markers, the facet markers
and the entity maps of the submeshes. The domains, the spaces and the forms
below are written in the same way for two triangles and for a million cells,
and a solver mesh with several cell types adds a group to each topology, as
in the [interior facet example](#an-interior-facet-integral).

The problem is Stokes flow in the fluid and linear elasticity in the solid,
coupled by a multiplier on the interface that equates the fluid velocity to
the velocity of the solid over a time step. With spaces of a library, `V` and
`Q` on `F`, `W` on `S` and `L` on `G`, with `V`, `W` and `L` of the value
shape `(2,)` and `Q` of the value shape `()`:

```python
u, v = function(V, "u"), function(V, "v")
p, q = function(Q, "p"), function(Q, "q")
d, w, d_old = function(W, "d"), function(W, "w"), function(W, "d_old")
lam, mu = function(L, "lam"), function(L, "mu")

# Material constants, the time step and the inlet traction, arrays of the solver.
nu, mu_s, lam_s = runtime_array("nu", ()), runtime_array("mu_s", ()), runtime_array("lam_s", ())
dt, g = runtime_array("dt", ()), runtime_array("g", (2,))
i, j = index("i"), index("j")
half, two = Integer(1) / Integer(2), Integer(2)

# Stokes flow on F.
grad_u = gradient(u, F)
grad_v = gradient(v, F)
div_u, div_v = index_sum(grad_u[i, i], (i,)), index_sum(grad_v[i, i], (i,))
fluid = integral(
    nu * index_sum(grad_u[i, j] * grad_v[i, j], (i, j)) - p * div_v - q * div_u, F
)

# Linear elasticity on S.
grad_d = gradient(d, S)
grad_w = gradient(w, S)
eps_d, eps_w = half * (grad_d[i, j] + grad_d[j, i]), half * (grad_w[i, j] + grad_w[j, i])
div_d, div_w = index_sum(grad_d[i, i], (i,)), index_sum(grad_w[i, i], (i,))
solid = integral(two * mu_s * index_sum(eps_d * eps_w, (i, j)) + lam_s * div_d * div_w, S)

# The restrictions to the interface, expressions on G.
u_G, v_G = restriction(u, G, "fluid"), restriction(v, G, "fluid")
d_G, w_G = restriction(d, G, "solid"), restriction(w, G, "solid")
d_old_G = restriction(d_old, G, "solid")

# The coupling: lam is the traction on both parts, and mu enforces u = (d - d_old) / dt.
coupling = integral(
    index_sum(lam[i] * (w_G[i] / dt - v_G[i]) + mu[i] * (d_G[i] / dt - u_G[i]), (i,)), G
)
a = form(fluid + solid + coupling, ((v, q, w, mu), (u, p, d, lam)))

# The right-hand side: the inlet traction and the old displacement.
traction = integral(index_sum(g[i] * restriction(v, G_in, "fluid")[i], (i,)), G_in)
l = form(traction + integral(index_sum(mu[i] * d_old_G[i], (i,)) / dt, G), ((v, q, w, mu),))
```

- `u` has the support `F` and `d` the support `S`, two decomposed domains with
  different topologies, so `u` and `d` do not combine. They meet only on `G`,
  through the restrictions, as two expressions on `G` at the same points of
  the reference facet `I`. Every integral has no support, so `fluid`, `solid`
  and `coupling` add into one expression, and the arguments of `a` are tuples
  of functions on three domains (see [`Form`](LANGUAGE.md#form)).
- `side_domain(F, interface)` takes the geometry of the first side, so the
  first argument is `F`, whose topology is `fluid_cells`. `normal(F, G, "fluid")`
  is the normal out of the fluid, and `normal(S, G, "solid")` is the negative
  of it.
- An integral is over the smallest domain, and the fields of the larger domains
  are restricted into it. A field on a submesh is never extended by zero onto
  the parent. This is the rule of a mixed-domain assembler, which fixes the
  mesh of the integration domain and reads the other fields through entity
  maps.
- Two groups of one topology differ in cell type or in sides, so a second
  marked set of inlet facets with the same side is a second topology, and the
  legacy `ds((3, 4))` is the sum of two integrals. The two kernels are one
  kernel after renaming, because the two topologies become the same
  placeholder.
- A material coefficient that differs between the fluid and the solid needs no
  submeshes when the fields live on `D`. With `c = cell_index(square, 0)`, an
  array `nu[c]` under `integral(..., D)` takes the value of the cell (see
  [`Integral` over a decomposed domain](LANGUAGE.md#integral-over-a-decomposed-domain)).
- A fluid domain that moves with the solid is not covered. The moved domain is
  a decomposed domain with the topology `fluid_cells` and new coordinates, so a
  function in a space on `F` has a different support from an expression on the
  moved domain, and the two do not combine. A map between two decomposed
  domains with the same topology, the identity in the reference coordinates of
  each cell, is a fixed composition that a later version may add.

The draft passes rewrite `a` as in the [interior facet
example](#an-interior-facet-integral). `expand_integrals` gives one
`IndexSum` for each of `F`, `S` and `G`, over the cell index of its topology.
On `G`, with `e = cell_index(interface, 0)`, the restriction `u_G` lowers to
the pull-back of `u` by the fluid cell `cells_fluid[e]` and then by the
inclusion of the facet into that cell, and `d_G` to the pull-back of `d` by
the solid cell `cells_solid[e]`. The kernel of `coupling` is an expression on
`I` that reads the reference basis of `V` on one fluid cell, of `W` on one
solid cell and of `L` on the facet, and the assembler packs the two cells
through the two side tables, the entity maps of the two submeshes.

```{dropdown} Legacy UFL
The legacy UFL `dx(1)`, `ds(3)` and `ds((3, 4))` are integrals over the decomposed domains of marked topologies, one topology for each marker value. A DOLFINx submesh is the decomposed domain of a marked topology, and its entity map is the side table into the parent.
```

### Constants on a domain

Requires: [A fluid and a solid](#a-fluid-and-a-solid), [Spaces of
constants](LANGUAGE.md#spaces-of-constants).

With the names of the [fluid and solid example](#a-fluid-and-a-solid), `R` a
space of constants of value shape `()`, and `R_S` a space of a library on `S`
whose functions are constant on `S`:

```python
# Multipliers for the mean pressure of the fluid and the mean volume change of the solid.
r, s = function(R, "r"), function(R, "s")
k, z = function(R_S, "k"), function(R_S, "z")
means = integral(p * s + q * r, F) + integral(k * div_w + z * div_d, S)
```

- `r` and `s` have no support, so the integrals decide where `r` and `s` act.
  `k` and `z` have the support `S`, so `k * q` raises `ValueError`.
- The multipliers join the `tuple` arguments, as in
  `((v, q, w, s, z), (u, p, d, r, k))`.
- A value for each cell marker is a function in a space of a library that is
  constant on each cell.

Known data is a constant that no form binds. With `nu` and `mu_s` as functions
in place of the arrays of the solver:

```python
nu, h = function(R, "nu"), function(R, "h")
mu_s = function(R_S, "mu_s")

# The sensitivity of the solid form to mu_s.
e = integral(two * mu_s * index_sum(eps_d * eps_w, (i, j)), S)
de_dmu = form(directional_derivative(e, mu_s, Integer(1)), (w, d))

# The gradient of the dissipation with respect to nu.
dissipation = integral(nu * index_sum(grad_u[i, j] * grad_u[i, j], (i, j)), F)
g_nu = form(directional_derivative(dissipation, nu, h), (h,))
```

`g_nu` is a linear form on `R`, a function in `dual_space(R)`. The variable
of a derivative is a function, so a constant whose sensitivity is wanted is a
function in `R`, not a `runtime_array`.

### A 1D network in a 3D mesh

Requires: [`DecomposedDomain`](LANGUAGE.md#decomposeddomain),
[Sides of a topology](LANGUAGE.md#sides-of-a-topology),
[Restriction](LANGUAGE.md#restriction), [A fluid and a
solid](#a-fluid-and-a-solid), [Points and point
evaluation](LANGUAGE.md#points-and-point-evaluation), [`Form`](LANGUAGE.md#form).

A network of vessels, a mesh of segments in $R^3$, coupled to the tissue
around it, a mesh of tetrahedra, with the pressure of each in a space on its
own mesh. The network is conforming: each segment is an edge of a
tetrahedron. The topology of the segments then has a side into the tissue,
as the interface of the [fluid and solid example](#a-fluid-and-a-solid)
has a side into each part, and the restriction of the tissue pressure to the
network is a pull-back by that side.

With `K3 = ReferenceSimplex(3)` and `I = ReferenceSimplex(1)`, the
coordinates `x3` of the tetrahedra with the free index `cell_index(T3, 0)`,
built as in the [`DecomposedDomain`](LANGUAGE.md#decomposeddomain) example, the
topology `T1` of the segments with the cell type `I` and the side
`(("tissue", T3, 0),)`, spaces of a library `V3` on `D3` and `V1` on `D1` of
the value shape `()`, and `runtime_extent` and
`runtime_array` of the [cut cell
example](#integrating-over-a-cut-cell):

```python
# The tissue, of dimension 3, and the network, of dimension 1, both in RD(3).
D3 = decomposed_domain((domain_map(K3, x3, RD(3)),))
D1 = side_domain(D3, T1)

p3, q3 = function(V3, "p3"), function(V3, "q3")
p1, q1 = function(V1, "p1"), function(V1, "q1")
beta = runtime_array("beta", ())
m = index("m")

# The Laplacian of each part on its own mesh. On D1 the gradient is the
# derivative along the segment, in the frame of RD(3).
grad_p3 = gradient(p3, D3)
grad_q3 = gradient(q3, D3)
grad_p1 = gradient(p1, D1)
grad_q1 = gradient(q1, D1)
tissue = integral(index_sum(grad_p3[m] * grad_q3[m], (m,)), D3)
vessel = integral(index_sum(grad_p1[m] * grad_q1[m], (m,)), D1)

# The tissue pressure on the network, and the exchange between the two.
p3_1, q3_1 = restriction(p3, D1, "tissue"), restriction(q3, D1, "tissue")
coupling = integral(beta * (p1 - p3_1) * (q1 - q3_1), D1)
a = form(tissue + vessel + coupling, ((q3, q1), (p3, p1)))
```

- The side tables of `T1` hold, for each segment `e`, one tetrahedron
  `cells_t[e]` that has the segment as an edge, chosen by the library among
  those that share the edge, and the positions `positions_t[e, v, k]` of the
  two ends in `K3`, the coordinates of two vertices of `K3`. The tables come
  from the connectivity of the mesh, and a segment and its tetrahedron share
  vertices, so the partition of `T3` gives the partition of `T1`.
- `restriction(p3, D1, "tissue")` is `pull_back(pull_back(p3, B_t), A_t)` with
  $B_t$ the tetrahedron of each segment and $A_t$ =
  `domain_map(I, affine_map(I, points_t), K3)` the edge, by
  [Restriction](LANGUAGE.md#restriction), an expression on `D1` that combines
  with `p1`. `p3` is continuous, so the restriction is the same from every
  tetrahedron that shares the edge, and the choice of the library does not
  change the form.
- `side_domain(D3, T1)` takes its geometry from the tissue side, the edge of
  the tetrahedron, so the integrals over `D1` use the length element of the
  edge.
- `p3` and `p1` have the supports `D3` and `D1` and never combine directly.
  Every integral has no support, so the three terms add, and the arguments of
  `a` are tuples of functions on the two meshes.
- The positions take one of twelve values, an edge of `K3` in one of two
  directions. A library may write the table as a gather from the twelve by
  a `RuntimeIndex` of the [cut cell
  example](#integrating-over-a-cut-cell), as `StaticTopology` computes
  the positions from vertex numbers, and a form compiler may then tabulate
  the basis of `V3` on each edge of `K3` once and select by the index, as a
  legacy facet kernel selects by the local facet number.
- For a network that is not conforming, the library provides a topology of
  the parts of the segments that lie in one tetrahedron each, with a side into
  the network and a side into the tissue, and the positions from the tissue
  two points of `K3` (see [Sides of a
  topology](LANGUAGE.md#sides-of-a-topology)). The same lines apply with
  `side_domain(D1, parts)` in place of `D1` and `p1` restricted from the
  network side, and the kernel evaluates the basis of `V3` at run time, as for
  a cut cell.

The draft passes rewrite `a` as in the [interior facet
example](#an-interior-facet-integral). `expand_integrals` gives one
`IndexSum` for `tissue`, over `cell_index(T3, 0)`, and one for `vessel` and
one for `coupling`, over `cell_index(T1, 0)`. On `coupling`,
`expand_pull_backs` lowers the two restrictions to pull-backs by one
tetrahedron, selected through the side table, and then by the affine map
through the positions. The kernel of `coupling` for a segment reads the
reference basis of `V1` on the segment and of `V3` on one tetrahedron, at the
points of the edge in the reference coordinates of each: the interface of a
facet kernel, with the cells of two meshes and of two dimensions. The local
matrix couples the degrees of freedom of `V1` on the segment with those of
`V3` on the tetrahedron, and the assembler scatters the matrix through the
maps of degrees of freedom of both.

The trace of a function in $H^1(\Omega)$ on a curve is not defined, so a
model may replace `p3` on the network by its average over a ring of radius
$R$ around each point of the network, in the plane of two unit normals `n1`
and `n2` of the segment, functions in a space `N` on `D1` of the value shape
`(3,)` that the solver provides. The ring points are not on the edges, so
the average is a family of maps from `D1` into `D3` and point evaluation,
with the cosines and sines of the angles and the weights as tables, because
the core has no `cos`:

```python
L = runtime_extent("L")
l = index("l")
R = runtime_array("R", ())
n1, n2 = function(N, "n1"), function(N, "n2")
cs, w = runtime_array("cs", (L, 2)), runtime_array("w", (L,))
X1 = coordinate(D1)

# The ring of L points around each point of the network, a family of maps into D3.
ring = domain_map(
    D1, component_tensor(X1[m] + R * (cs[l, 0] * n1[m] + cs[l, 1] * n2[m]), (m,)), D3
)
p3_avg = index_sum(w[l] * pull_back(p3, ring), (l,))
q3_avg = index_sum(w[l] * pull_back(q3, ring), (l,))

coupling_avg = integral(beta * (p1 - p3_avg) * (q1 - q3_avg), D1)
```

- `ring` is a `DomainMap` with the parameter domain `D1`, a family with the
  free index `l`, and the ambient domain `D3`, a decomposed domain.
  `pull_back(p3, ring)` has the support `D1` and the free index `l`, so
  `p3_avg` is an expression on `D1`. With `R = 0`, `pull_back(p3, ring)` is
  the trace of `p3` on the network, `pull_back(p3, domain_map(D1, X1, D3))`,
  equal to `p3_1` for a conforming network.
- The pull-back of `p3` by `ring` remains after `expand_pull_backs`, and the
  pass of [Points and point evaluation](LANGUAGE.md#points-and-point-evaluation)
  lowers it after a quadrature pass on `I` as in the [quadrature example](#quadrature). The maps below `p3` are then `ring`,
  `cell_domain(D1, 0)` and the quadrature points, a chain with the free
  indices of the segment, the quadrature point and the ring point, which the
  pass reads as one family of points. The side tables of that family hold
  the tetrahedron that contains each point and its position, and the kernel
  for a segment reads `V3` on one tetrahedron for each point, a list of
  cells rather than one cell. A kernel allows that, because the interface of
  a kernel follows from the expression.

### The Laplace-Beltrami operator

Requires: [`DecomposedDomain`](LANGUAGE.md#decomposeddomain),
[`Gradient`](LANGUAGE.md#gradient), [`Integral` over a decomposed
domain](LANGUAGE.md#integral-over-a-decomposed-domain), [`Form`](LANGUAGE.md#form).

See also: [`Inverse`](LANGUAGE.md#inverse).

The Poisson forms on a surface in $R^3$ are the forms of the Laplace–Beltrami
operator. The mesh is built as the mesh of the
[`DecomposedDomain`](LANGUAGE.md#decomposeddomain) example, with the vertex
positions in $R^3$. A mesh library gives a topology `T` of $N$ triangles with
the cell type `K = ReferenceSimplex(2)`, the vertex positions `y` of value
shape `(m, 3)` and the table `cell_vertices` of the shape `(N, 3)` and the
extent $m$, the vertices of each triangle. The ambient coordinate domain is `RD(3)`, and `V`
is a space of a library on `M`:

```python
# The affine map x of each triangle, gathered from the vertex positions y.
X = coordinate(K)
c = cell_index(T, 0)
k, w = index("k"), index("w")
p = component_tensor(y[gathered_index(cell_vertices, (c, w)), k], (c, w, k))
x = component_tensor(
    p[c, 0, k] + (p[c, 1, k] - p[c, 0, k]) * X[0] + (p[c, 2, k] - p[c, 0, k]) * X[1],
    (k,),
)
M = decomposed_domain((domain_map(K, x, RD(3)),))

# Functions in a space of a library on M.
u, v, f = function(V, "u"), function(V, "v"), function(V, "f")

# The surface gradients, of value shape (3,) in the frame of RD(3).
i = index("i")
grad_u = gradient(u, M)
grad_v = gradient(v, M)

a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), M), (v, u))
L = form(integral(f * v, M), (v,))
```

- `M` has the dimension 2 and the ambient coordinate domain `RD(3)`, so
  [`Gradient`](LANGUAGE.md#gradient) on `M` is the tangential gradient, and
  `grad_u` is the surface gradient of `u`.
- `integral(..., M)` multiplies the integrand on `K` by the density
  $\rho = \sqrt{\det(J^T J)}$, with $J$ the Jacobian matrix of `x`, of value
  shape `(3, 2)`. The forms need no orientation of `M`.
- `expand_pull_backs` rewrites the pull-back of `grad_u` by the cells as
  $K^T \hat g$, with $K$ = `inverse_jacobian` of the cells, here `inverse(J)`, the pseudo-inverse, of value shape
  `(2, 3)`, and $\hat g$ the gradient of the pull-back of `u` on `K`. The integrand of `a` on `K` is then
  $\hat g_u^T K K^T \hat g_v \, \rho$, and $K K^T = (J^T J)^{-1}$ is the
  inverse of the metric that `x` induces.
- With `P = gradient(coordinate(M), M)`, the projection onto the tangent
  plane, `gradient(grad_u, M)` is the tangential
  derivative $H$ of the surface gradient, which is not symmetric in general.
  The covariant Hessian of `u` is $P H P$.

```{dropdown} Legacy UFL
Legacy UFL writes `inner(grad(u), grad(v)) * dx` on a mesh of topological dimension 2 in $R^3$, and `grad` is the surface gradient because of the dimensions of the mesh. In UFLx, the same rule of `Gradient` gives the gradient in every dimension.
```

### A deformed domain

Requires: [`DecomposedDomain`](LANGUAGE.md#decomposeddomain),
[`DomainMap`](LANGUAGE.md#domainmap), [Jacobian and
density](LANGUAGE.md#jacobian-and-density), [`Gradient`](LANGUAGE.md#gradient),
[`Form`](LANGUAGE.md#form), [The Laplace-Beltrami
operator](#the-laplace-beltrami-operator).

Shape optimisation and updated Lagrangian mechanics write a form on the image
of the mesh under a deformation, and assemble on the mesh. The deformation is
a domain map with the mesh as parameter domain. With `D` the square mesh of
the [`DecomposedDomain`](LANGUAGE.md#decomposeddomain) examples, `W` a vector
space of a library on `D`, and `V_s` a space of a library on the deformed
mesh:

```python
k, i = index("k"), index("i")
s = function(W, "s")
phi = component_tensor(coordinate(D)[k] + s[k], (k,))

# The deformed square, the image of the mesh under phi.
D_s = decomposed_domain((domain_map(D, phi, RD(2)),))

# The volume of the deformed square, a functional of s.
volume = integral(Integer(1), D_s)

# The Poisson stiffness form on the deformed square.
u, v = function(V_s, "u"), function(V_s, "v")
a = form(integral(index_sum(gradient(u, D_s)[i] * gradient(v, D_s)[i], (i,)), D_s), (v, u))
```

- `expand_integrals` rewrites `volume` as the integral over `D` of
  `density(T)`, with `T = cell_domain(D_s, 0)`, and then, because `D` is
  itself decomposed, over the cells of `D`. `jacobian(T)` is
  `gradient(phi, D)`, which the gradients pass rewrites
  as $I + \nabla s$, so `density(T)` is $|\det(I + \nabla s)|$, the Jacobian
  determinant of shape calculus.
- `expand_pull_backs` rewrites the pull-back of `gradient(u, D_s)` by `T` as
  $K^T \nabla(u \circ \phi)$ with $K$ = `inverse_jacobian(T)`, here `inverse(jacobian(T))`, the
  $(I + \nabla s)^{-T}$ of shape calculus. The library represents `V_s` on
  `D`, so `u` has the degrees of freedom of a space on `D` and the geometry of
  `D_s`.
- The shape derivative of `volume` is the Gateaux derivative with respect to
  `s` (see [Gateaux derivative](LANGUAGE.md#gateaux-derivative)).

A surface deforms in the same way. With `M` the surface mesh of the
[Laplace-Beltrami example](#the-laplace-beltrami-operator) and `w` a
displacement in a vector space on `M`:

```python
w = function(W_M, "w")
phi_M = component_tensor(coordinate(M)[k] + w[k], (k,))
M_w = decomposed_domain((domain_map(M, phi_M, RD(3)),))

# The area of the deformed surface.
area = integral(Integer(1), M_w)
```

`M` has dimension 2 in `RD(3)`, so `jacobian(cell_domain(M_w, 0))` is
$\Pi + \nabla_M w$, with $\Pi$ the projection onto the tangent plane and
$\nabla_M w$ the surface gradient, a `(3, 3)` matrix of rank 2, and the density
is $\sqrt{\det(J^T J + I - \Pi)}$ (see [Jacobian and
density](LANGUAGE.md#jacobian-and-density)). The cells of `M` then contribute
their own density $\sqrt{\det(J_c^T J_c)}$ as in the Laplace-Beltrami example.

```{dropdown} Isogeometric analysis
Isogeometric analysis by extraction, as in tIGAr, builds a mesh of the
parametric knot spans, stores the NURBS map $F$ as a function in the spline
space on that mesh, and inserts $\det DF$ and $DF^{-T}$ into every form by
hand. Here the physical patch is `domain_map(parametric_mesh, F, RD(3))`, a
multi-patch geometry is a decomposed domain with one piece per patch, and the
passes insert the factors.
```

```{dropdown} Shells
A shell code such as FEniCS-Shells defines the midsurface as a map $\phi_0$
from a flat parametric mesh into $R^3$ and forms the fundamental forms from
$\nabla \phi_0$ by hand. Here the midsurface is `domain_map(flat_mesh, phi_0, RD(3))`
with $k = m = 2 < n = 3$, `jacobian(T)` is $\nabla \phi_0$, its Gram matrix
$J^T J$ is the first fundamental form, and `density(T)` is
$\sqrt{\det(J^T J)}$, the area element of the midsurface.
```

### The flux through a surface

Requires: [The Laplace-Beltrami
operator](#the-laplace-beltrami-operator), [Concept: integral over a
domain](LANGUAGE.md#concept-integral-over-a-domain), [Jacobian and
density](LANGUAGE.md#jacobian-and-density).

The flux of a vector field through a surface depends on which side of the
surface is out, so the flux needs a consistent orientation of the cells. With
`K`, `T`, `c` and `x` of the [Laplace-Beltrami
example](#the-laplace-beltrami-operator), where `M` is a closed surface such
as a triangulated sphere:

```python
# The triangles of M, a family with the free index c.
cells = domain_map(K, x, RD(3))
J = jacobian(cells)
r, l = index("r"), index("l")

def cross(a, b):
    return tensor((a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]))

# The normal of each triangle, with the length of the density.
nu = cross(component_tensor(J[r, 0], (r,)), component_tensor(J[r, 1], (r,)))

# The sign of each triangle from the mesh library, +1 where the order of the
# vertices of the triangle agrees with the outward orientation of M.
N = T.num_cells[0]
o = runtime_array("o", (N,))

# The outward flux of the position field, 3 times the volume inside M.
Y = coordinate(RD(3))
flux = index_sum(
    integral(o[c] * index_sum(pull_back(Y, cells)[l] * nu[l], (l,)), K), (c,)
)

# The area of M, with no orientation.
area = integral(Integer(1), M)
```

- The flux is an integral over the reference triangle `K` of each cell, not
  over `M`, so no density enters. `nu` is the normal of the triangle with the
  length of the density, and the direction of `nu` follows the order of the
  vertices of the triangle in `T`. If the mesh library numbers some triangles
  clockwise and others anticlockwise, seen from outside, the contributions of
  the two groups have opposite signs, and the sum is not the flux.
- `o[c]` has the free index `c` of the triangles, so each triangle takes the
  sign of that triangle. With the signs, every contribution uses the outward
  normal. A mesh library may instead number the vertices of every triangle
  consistently, and then the flux needs no `o`.
- `area` integrates against the density $\sqrt{\det(J^T J)}$, which does not
  depend on the order of the vertices, so `area` needs no signs.
- A surface with no consistent orientation, such as a Möbius strip, has no
  signs `o`. The area of the surface is defined, and the flux is not.

### Helmholtz with a perfectly matched layer

Requires: [Literals](LANGUAGE.md#literals), [`Conjugate`](LANGUAGE.md#conjugate),
[`ConjugateSpace`](LANGUAGE.md#conjugatespace), [`Form`](LANGUAGE.md#form),
[A classic FEniCS layer for the Poisson
problem](#a-classic-fenics-layer-for-the-poisson-problem).

The Helmholtz equation $-\Delta u - k_0^2 u = f$ on the square mesh `D` of
the [`DecomposedDomain`](LANGUAGE.md#decomposeddomain) example, with a
perfectly matched layer that absorbs outgoing waves. The layer stretches each
coordinate into the complex plane by a factor $s_x = 1 + i \sigma_x / k_0$,
and the stretched equation is
$-\nabla \cdot (A \nabla u) - k_0^2 s_x s_y u = f$ with
$A = \mathrm{diag}(s_y / s_x, s_x / s_y)$. `V` is a complex Lagrange space on
`D`, `F` a real space on `D` for the data, and `C` a real space of
constants. The solver sets `sigma_x` and `sigma_y` to zero outside the layer:

```python
k0 = function(C, "k0")
sigma_x, sigma_y = function(F, "sigma_x"), function(F, "sigma_y")
one, i_unit = Integer(1), ImaginaryUnit()
s_x = one + i_unit * sigma_x / k0
s_y = one + i_unit * sigma_y / k0

u, v, f = function(V, "u"), function(V, "v"), function(F, "f")
grad_u, grad_v = gradient(u, D), gradient(v, D)
j = index("j")

# The stretched Laplacian and the stretched mass, sesquilinear in u and v.
A = tensor((s_y / s_x, s_x / s_y))
stiffness = index_sum(A[j] * grad_u[j] * conjugate(grad_v[j]), (j,))
mass = k0 * k0 * s_x * s_y * u * conjugate(v)
a = form(integral(stiffness - mass, D), (conjugate(v), u))
L = form(integral(f * conjugate(v), D), (conjugate(v),))
```

- `a` is antilinear in `v` and linear in `u`, which the entries state, so
  the argument spaces are `(conjugate_space(V), V)`, and `L` has
  `(conjugate_space(V),)`. `form(integral(stiffness - mass, D), (v, u))`
  states a form that is linear in `v`, which construction accepts and a
  linearity check rejects, because the expression conjugates `v`.
- Every expression here is complex, because `ImaginaryUnit()` occurs in
  `s_x`. With `s_x` and `s_y` replaced by `one`, the same lines are the
  real Helmholtz forms when `V` is real, and the conjugates change nothing.
- The layer is often derived as the pull-back of the Helmholtz equation by the
  complex map $\tilde x = x + (i / k_0) \int \sigma_x$. A `DomainMap` is real,
  so the pull-back is done once by hand and enters as the tensor `A` and the
  factor $s_x s_y$, which is the form every code assembles.
- The adjoint of `a` is `form(conjugate(stiffness - mass), (conjugate(u), v))`
  (see [`Action`](LANGUAGE.md#action)). The first-order absorbing condition
  $\int_{\partial D} i k_0 u \bar v$ is one more term, written with
  restrictions to the boundary as in the [fluid and solid
  example](#a-fluid-and-a-solid).

### Mixed Poisson with Raviart-Thomas elements

Requires: [Representation on the reference
cell](LANGUAGE.md#representation-on-the-reference-cell),
[`Form`](LANGUAGE.md#form).

See also: [A Piola extension](LANGUAGE.md#a-piola-extension).

The mixed form of the Poisson problem seeks a flux $\sigma$ and a potential
$u$ with $\sigma = \nabla u$ and $\mathrm{div}\, \sigma = -f$. `RT` is a
Raviart-Thomas space on the square mesh `D` of the
[`DecomposedDomain`](LANGUAGE.md#decomposeddomain) example, of value shape
`(2,)`, and `DG` a discontinuous space of degree 0:

```python
sigma, tau = function(RT, "sigma"), function(RT, "tau")
u, v, f = function(DG, "u"), function(DG, "v"), function(DG, "f")
i = index("i")

def divergence(w):
    return index_sum(gradient(w, D)[i, i], (i,))

e = index_sum(sigma[i] * tau[i], (i,)) + u * divergence(tau) + divergence(sigma) * v
a = form(integral(e, D), ((tau, v), (sigma, u)))
L = form(integral(-f * v, D), ((tau, v),))
```

- There is no mixed space. Each argument is a `tuple` of a function in `RT`
  and a function in `DG`, so the argument spaces are `((RT, DG), (RT, DG))`
  and `((RT, DG),)`. The block of `a` for `tau` and `u` is `a` with `v` and
  `sigma` replaced by `Zero`, the term `u * divergence(tau)`. `L` depends on
  `v` only, which a `tuple` argument allows.
- The library of `RT` represents a function on the reference cell with the
  contravariant Piola map, so `pull_back(sigma, cell_domain(D, 0))` is
  $J \hat\sigma / \det J$ with $\hat\sigma$ the function on `K`, and the
  library of `DG` uses the identity. The core sees only the compositions.
- `expand_pull_backs` rewrites the pull-back of `gradient(sigma, D)` as
  $K^T \nabla(J \hat\sigma / \det J)$, and `divergence` is its trace. On the
  affine cells of `D`, the gradients pass finds the derivatives of $J$ and
  $\det J$ to be zero, and the trace reduces to
  $\widehat{\mathrm{div}}\, \hat\sigma / \det J$. On a curved cell the extra
  terms cancel only in value, which [A Piola
  extension](LANGUAGE.md#a-piola-extension) makes structural.
- A Nédélec space differs only in the map, $K^T \hat u$, with the curl in
  place of the divergence as a free function.

### Dual-weighted residual error estimation

Requires: [`Form`](LANGUAGE.md#form), [`Action`](LANGUAGE.md#action),
[`DualSpace`](LANGUAGE.md#dualspace), [A Bank-Weiser error
estimator](#a-bank-weiser-error-estimator).

The dual-weighted residual method estimates the error in a goal functional
$J(u)$ by the residual of the Galerkin solution weighted by the solution $z$
of the adjoint problem, $J(u) - J(u_h) = R(z) = R(z - z_h)$ for any $z_h$ in
the trial space. The example takes an advection-diffusion problem on the
square mesh `D`, with `V` a Lagrange space of degree 1, `W` a Lagrange space
of degree 2 in which the adjoint solution is computed, `F` a space for the
data and `B` a vector space for the advection field:

```python
u, v, u_h, z_h = function(V, "u"), function(V, "v"), function(V, "u_h"), function(V, "z_h")
z = function(W, "z")
f, psi, b = function(F, "f"), function(F, "psi"), function(B, "b")
i = index("i")

def a_expr(w, y):
    grad_w, grad_y = gradient(w, D), gradient(y, D)
    return index_sum(grad_w[i] * grad_y[i], (i,)) + index_sum(b[i] * grad_w[i], (i,)) * y

# The primal problem, the goal, and the adjoint problem with the arguments swapped.
a = form(integral(a_expr(u, v), D), (v, u))
L = form(integral(f * v, D), (v,))
goal = form(integral(psi * v, D), (v,))
a_adjoint = form(integral(a_expr(u, v), D), (u, v))

# The residual at u_h, a linear form, and the error in the goal.
R = L - action(a, 1, u_h)
eta = action(R, 0, z - z_h)

# The same with the assembled load, a function in the dual space.
l = function(dual_space(V), "l")
R_h = l - action(a, 1, u_h)
```

- `a_adjoint` is `a` with the entries swapped, so the function at position 0
  is `u` and the adjoint problem seeks `z` at position 1 with `goal` as the
  load. For a complex `V` the adjoint also conjugates the expression (see
  [`Action`](LANGUAGE.md#action)).
- `R` has the argument spaces `(V,)`, because `action(a, 1, u_h)` removes the
  entry at position 1 and the difference with `L` is a form. `eta` is a
  functional, with the argument spaces `()` and no support, one number after
  assembly.
- The argument of an action needs only the static type of a function in `V`,
  so `z` from the richer space `W` is accepted, and so is `z - z_h`. Galerkin
  orthogonality gives `action(R, 0, z_h) == 0` in value, which is why the
  interpolant is subtracted.
- `R_h` is the residual as the solver holds it, with `l` the assembled load,
  a function in `dual_space(V)`. `action(R_h, 0, w)` pairs it with any `w` of
  the static type of `V`, and the library evaluates `action(l, 0, w)` as the
  product of the load vector with the coefficients of `w`.
- The cell-wise estimator integrates `R` by parts into the interior residual
  and the facet jumps of the [Bank-Weiser example](#a-bank-weiser-error-estimator),
  with `z - z_h` in place of `v`, so that each cell receives its own weighted
  residual.

### Hyperelasticity and the stress

Requires: [`DirectionalDerivative`](LANGUAGE.md#directionalderivative),
[Gateaux derivative](LANGUAGE.md#gateaux-derivative), [Tensor
arithmetic](LANGUAGE.md#tensor-arithmetic), [Transpose](#transpose).

See also: [Symmetric tensors](#symmetric-tensors).

A Saint Venant–Kirchhoff material on the square mesh `D`, with a real vector
space `V` and a real space `M` for the Lamé parameters, both of a library on
`D`, and a real space `R` of constants of value shape `(2, 2)`. The stress is
$S = 2 \operatorname{sym}(\partial \psi / \partial C)$:

```python
one, zero, two = Integer(1), Integer(0), Integer(2)
I = tensor((tensor((one, zero)), tensor((zero, one))))
i, j, k, l, a = index("i"), index("j"), index("k"), index("l"), index("a")

u, v, w = function(V, "u"), function(V, "v"), function(V, "w")
mu, lmbda = function(M, "mu"), function(M, "lambda")

def strain_energy(C):
    E = (C - I) / two
    return lmbda / two * index_sum(E[a, a], (a,)) ** two + mu * index_sum(E[i, j] * E[i, j], (i, j))

F = I + gradient(u, D)
C = transpose(F) @ F
psi = strain_energy(C)

# The derivative with respect to a symbol at C, in the symmetric unit directions.
C_ = function(R, "C")
N = component_tensor((I[i, k] * I[j, l] + I[i, l] * I[j, k]) / two, (i, j))
S = two * component_tensor(directional_derivative(strain_energy(C_), C_, N, C), (k, l))
sigma = F @ S @ transpose(F) / determinant(F)

Pi = integral(psi, D)
residual = gateaux_derivative(Pi, u, v)
tangent = gateaux_derivative(residual, u, w)
```

- `S` and the Cauchy stress `sigma` are expressions on `D`, which a solver
  interpolates or evaluates for output. `apply_derivatives` evaluates `S` to
  an expression equal to $\lambda \operatorname{tr}(E) I + 2 \mu E$.
- `V` is real, so `residual` is the form of `directional_derivative(Pi, u, v)` that binds `v`, and `tangent` binds `w` after `v`.

### Symmetric tensors

Requires: [`GatheredIndex`](LANGUAGE.md#gatheredindex), [Hyperelasticity and
the stress](#hyperelasticity-and-the-stress).

The core has no symmetry of a tensor. A symmetric array of data is read from
its stored values through a table, and a derivative with respect to a
symmetric tensor states the symmetry with the direction.

A stiffness tensor $\mathbb{C}_{ijkl}$ with the minor and the major
symmetries has 6 values in 2D. With the names of the [hyperelasticity
example](#hyperelasticity-and-the-stress), a real space `K` of constants of
value shape `(6,)`, and a table `packed` of the shape `(2, 2, 2, 2)` and the
extent 6 with the position of each entry in `c`:

```python
c = function(K, "c")
stiffness = component_tensor(c[gathered_index(packed, (i, j, k, l))], (i, j, k, l))

eps_u = (gradient(u, D) + transpose(gradient(u, D))) / two
eps_v = (gradient(v, D) + transpose(gradient(v, D))) / two
elastic = form(integral(index_sum(stiffness[i, j, k, l] * eps_u[k, l] * eps_v[i, j], (i, j, k, l)), D), (v, u))
```

`stiffness` has 16 entries, and the kernel receives the 6 values of `c`.

```{dropdown} Directions and Voigt notation
For `psi_2 = C_[0, 1] ** two`, the direction
`component_tensor(I[i, k] * I[j, l], (i, j))` gives a derivative with the
entry $2 C_{01}$ at `(0, 1)` and 0 at `(1, 0)`, and the symmetric direction
`N` gives $C_{01}$ at both. The two agree in every symmetric direction, and
only the second is symmetric.

Voigt and Mandel notation are representations of the full arrays. A Voigt
strain or a Mandel vector scales the entries off the diagonal, by 2 or by
$\sqrt{2}$, so a library reads the vector with a table and an array of the
factors. The derivative with respect to the stored values of a Voigt vector
differs from the derivative with respect to the full tensor by a factor 2 off
the diagonal.
```

### A classic FEniCS layer for the Poisson problem

Requires: [The base object `Base`](LANGUAGE.md#the-base-object-base),
[`Function`](LANGUAGE.md#function), [`Form`](LANGUAGE.md#form),
[Restriction](LANGUAGE.md#restriction), [`Normal`](LANGUAGE.md#normal),
[`Integral` over a decomposed
domain](LANGUAGE.md#integral-over-a-decomposed-domain).

The Poisson problem, written as classic FEniCS writes it:
`inner(grad(u), grad(v)) * dx`. The layer is a module `uflx.legacy`, which
does not exist, of a page of free functions and two plain Python classes over
the core. Nothing in the module subclasses a core class, so the forms that the
module produces are ordinary core graphs, the same graphs as the core free
functions build by hand. The example shows what the core leaves to a layer: the argument
numbers of test and trial functions, measures that find or hold their domain,
restrict the integrand and bind the arguments, a facet normal that belongs to
the mesh, and `grad`, `div` and `curl` that read the domain from the function.
`Lagrange`, `Constants` and `Fields` stand for spaces of a library: Lagrange
elements of a value shape, constants with no domain, and fields of a value
shape with no element.

```python
# uflx/legacy.py
from uflx import (
    AbstractExpression,
    Base,
    Function,
    conjugate,
    form,
    function,
    gradient,
    index,
    index_sum,
    integral,
    normal,
    restriction,
    side_domain,
    tensor,
)


def TestFunction(V):
    return function(V, ("argument", 0))


def TrialFunction(V):
    return function(V, ("argument", 1))


def Coefficient(V, name):
    return function(V, name)


def Constant(name, shape=()):
    return function(Constants(shape), name)


def FunctionSpace(mesh, family, degree, shape=()):
    return Lagrange(mesh.domain, degree, shape)


def VectorFunctionSpace(mesh, family, degree):
    return FunctionSpace(mesh, family, degree, (mesh.domain.ambient_domain.dimension,))


def FacetNormal(mesh):
    # A placeholder on the mesh, which restrict replaces by the normal of the facets.
    n = mesh.domain.ambient_domain.dimension
    return function(Fields(mesh.domain, (n,)), "facet normal")


def is_facet_normal(obj):
    return isinstance(obj, Function) and obj.label == "facet normal"


def grad(f):
    return gradient(f, f.support)


def div(w):
    i = index("div")
    return index_sum(grad(w)[i, i], (i,))


def curl(w):
    # grad(w)[i, j] is the derivative of w[i] along x[j], as in legacy UFL.
    G = grad(w)
    if w.value_shape == (2,):
        return G[1, 0] - G[0, 1]
    return tensor((G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]))


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
    # The legacy extract_arguments: the functions labelled ("argument", n), by n.
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
        return normal(obj.support, E, "inside")
    if isinstance(obj, tuple):
        return tuple(restrict(item, E) for item in obj)
    if not isinstance(obj, AbstractExpression) or obj.support is None:
        return obj
    if not contains(obj, is_facet_normal):
        return restriction(obj, E, "inside")
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

With the module, a Poisson problem with a conductivity, a source on a marked
subdomain, a second source from a flux field and a Neumann condition from the
same field on a marked part of the boundary reads as a classic FEniCS demo.
`D` is the mesh of the [`DecomposedDomain`](LANGUAGE.md#decomposeddomain)
example and `boundary` the topology of its boundary facets of the [static
topology example](#a-static-topology). `lower` and `left` are marker
data that a mesh library supplies: the topology of the lower triangle and the
topology of the left edge, each with a side `"inside"` into the mesh, built as
`fluid_cells` and `inlet` of the [fluid and solid
example](#a-fluid-and-a-solid) are. A user script of the module sees
none of the topologies, as a legacy script sees only a `MeshFunction`.

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

- For a real `V`, `a` is `kappa` times the Poisson bilinear form
  `form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), D), (v, u))`, with
  `grad_u = gradient(u, D)`, `grad_v = gradient(v, D)` and `i = index("i")`,
  up to the labels of `u`, `v` and `i` and the conjugates, which change nothing in
  real arithmetic: `conjugate_space(V)` is `V`, so the argument spaces are
  `(V, V)`. `L` is a sum of three forms with the argument spaces `(V,)`.
  `f * dx` with no argument stays an `Integral`, a functional.
- For a complex `V`, the same lines give the sesquilinear forms of a
  complex-mode legacy demo. `inner` conjugates the second operand, as the
  legacy `inner` does, and the measure binds `conjugate(v)` at position 0, so
  `a` has the argument spaces `(conjugate_space(V), V)` and is antilinear in
  `v`, and `L` has `(conjugate_space(V),)`. The convention is the module's,
  not the core's: the core records the linearity that the form states (see
  [`Form`](LANGUAGE.md#form)). `f * v * dx` in complex mode is linear in `v`
  but bound as antilinear, which a pass that checks the linearity rejects, as
  legacy UFL rejects an unconjugated test function in complex mode.
- The argument numbers of legacy UFL live in the labels of the module. The
  core has none: `u` and `v` are functions, and the measure states the
  arguments. `arguments` is the legacy `extract_arguments`, a walk of the
  construction arguments, which rule 5 of [`Base`](LANGUAGE.md#the-base-object-base)
  makes possible for every expression. A module with `MixedElement` and
  `split` would collect the members of a mixed space into a `tuple` argument.
- `Measure` has the four slots of the legacy `Measure`: the integral type,
  the domain, the subdomain id and the subdomain data, and calling a measure
  fills slots, so `dx(1)` and `Measure("dx", domain=mesh, subdomain_data=...)`
  read as in a DOLFINx script. The module-level `dx` has no domain and finds
  the domain in the support of the integrand, as the legacy `dx` finds the
  mesh in the integrand; `f * dx` with a constant `f` raises in `integral`,
  the legacy "missing an integration domain". `ds` needs the mesh, because
  the boundary facets are mesh data. The legacy subdomain data is opaque to
  UFL and read by the assembler; here it is a `dict` from marker value to
  topology, and `dx(1)` and `ds(1)` are the measures over `lower` and
  `left`.
- A measure other than `dx` restricts the integrand to its domain from the
  side `"inside"`, which legacy UFL does without saying so for exterior facets
  and subdomains. `restrict` is the legacy `apply_restrictions`, run when the
  form is built: the largest subtrees on the mesh become one restriction each,
  so `grad(u)` becomes `restriction(grad(u), E, "inside")` and not the
  derivative of a restriction, and the rebuilding of the other nodes uses
  rule 4 of [`Base`](LANGUAGE.md#the-base-object-base). `dx` with the mesh
  as domain restricts nothing, since the domain of the measure is the support
  of the integrand. A `dS` would be a
  measure over the interior facets with a side for each restriction, the
  legacy `u('+')` would be `restriction(u, E, "+")`, and `avg` and `jump`
  would be arithmetic, as in the [interior facet
  example](#an-interior-facet-integral).
- The legacy `FacetNormal` belongs to the mesh and is resolved per integral.
  `n` is a placeholder function of the module in a space of fields on the
  mesh, known by its label, so that `inner(q, n)` is an expression on the
  mesh like any other, and `restrict` replaces `n` by `normal(D, E, "inside")`
  on the domain `E` of the measure. A subtree that contains `n` is rebuilt
  rather than restricted, because the normal has no meaning on the mesh.
- `grad` reads the domain from the support of the function, as the legacy
  `grad` reads the mesh from the coefficient, and `grad(w)[i, j]` is the
  derivative of `w[i]` along `x[j]`, as in legacy UFL, so `div` is the trace
  and `curl` reads off the antisymmetric part. On the surface mesh of
  the [Laplace-Beltrami example](#the-laplace-beltrami-operator)
  the same `grad` is the surface gradient. The legacy `grad` of a `Constant`
  is the zero of the shape of the mesh. `kappa` has no support, so `grad(kappa)`
  raises in the module; a module that wants the legacy behaviour gives
  `Constant` the mesh and returns a `Zero`.
- A form compiler sees only core graphs, so the module costs nothing at
  compile time, and the passes apply unchanged.

Not in the module: `derivative`, which is the Gateaux derivative, not yet
covered, and `u('+')`, which calls a core `Function`, so a module writes a
free function for it. Boundary conditions, assembly and solving are outside
legacy UFL too.

```{dropdown} Legacy UFL
The legacy spelling `inner(grad(u), grad(v)) * dx` is a module `uflx.legacy` of free functions over the core, whose forms are ordinary core graphs.
```


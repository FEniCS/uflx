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
the side tables from the vertex numbers: the cell on the side, the inclusion of
each cell as the affine map through the coordinates `vertices` of the paired
local vertices of the reference cell on the side, and the reference normal of
a facet (see [Side tables of a static topology](#side-tables-of-a-static-topology)).
So a `StaticTopology` covers sub-entities only, not pieces of a cell. For a cell seen from two cells of the
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
coordinates. The positions of the vertices come from the coordinate function
of a [mesh](LANGUAGE.md#meshes).

#### Side tables of a static topology

The [side tables](LANGUAGE.md#side-tables) of a `StaticTopology` are built from
the vertices of the reference cells, which the core does not read. The
library function `affine_map(F, points)` returns the affine map on a reference
cell $F$ of dimension $k$ with the vertices $V_v$ that sends vertex $v$ of $F$
to row $p_v$ of `points` for $v \le k$:

$$
X \mapsto p_0 + \sum_{d < k} (p_{d+1} - p_0) \, \big(E^{-1} (X - V_0)\big)_d,
$$

with $E$ the matrix with the columns $V_{d+1} - V_0$. With the derived index
`l = Index(("library.affine_map", F, points, 0))`, the map is

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

Only the first $k + 1$ rows of `points` enter the map, and the first $k + 1$
vertices of $F$ must not lie in a hyperplane of $F$. A sub-entity of $K$ with
the vertex numbers $(w_0, w_1, \dots)$ and the reference cell $F$ is
`domain_map(F, affine_map(F, points), K)`, with row $v$ of `points` the
coordinates `K.vertices[w_v]`. For the unit triangle, the sub-entity with the
vertex numbers $(1, 2)$ and the unit interval as $F$ is the hypotenuse
$t \mapsto (1 - t, t)$, and for the unit interval, the sub-entity with the
vertex number $(1)$ and the reference point as $F$ is the end point.

For the group $i$ and a side with the reference cell $K_j$, with the cell
index `e` of the group, the positions `positions[e, v]` of the vertices of
each cell $e$ in $K_j$, and derived indices `v` and `k` of the library:

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

with `V = K_j.vertices` and $m_j$ vertices. The sums are nested to the left.
Every side sends the vertices of $F_i$ to the positions of the same vertices
of the cell $e$, in the order of the cell, so the inclusions of all sides
parameterise the cell $e$ in the same way.

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

`vertices` is a property of the library, not of `AbstractReferenceCell`. The
side tables of a `StaticTopology` read it (see [Side tables of a static
topology](#side-tables-of-a-static-topology)).

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
with the side `"inside"`, a Lagrange space `V` of polynomial degree 1 on
`square`, a space `F` on `square` for the data, and the Bank-Weiser space `W`
on `square`, a discontinuous space whose pull-back by each cell is the space
of quadratic functions on `K` with zero vertex values, a representation of the
library:

```python
u_h, f, g = function(V, "u_h"), function(F, "f"), function(F, "g")
e, v = function(W, "e"), function(W, "v")
f_D, g_D, v_D = push_forward(f, D), push_forward(g, D), push_forward(v, D)
i = index("i")
half = Integer(1) / Integer(2)

# The local bilinear form, one block for each cell because W is discontinuous.
grad_e = gradient(push_forward(e, D), D)
grad_v = gradient(v_D, D)
a = form(integral(index_sum(grad_e[i] * grad_v[i], (i,)), D), (v, e))

# The interior residual. The Laplacian of u_h is zero for degree 1, and is kept
# for higher degrees.
grad_u = gradient(push_forward(u_h, D), D)
lap_u = divergence(grad_u, D)
r_cell = integral((f_D + lap_u) * v_D, D)

# The jump of the normal derivative with the normal of "+". From "-" the normal
# and the jump both change sign, so the same scalar serves both cells.
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

The passes rewrite `a` and `r` as in the [interior facet
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
  cell. The form is the same with `pull_back(v_D, side_map(D, E, "self"))` in place of the
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
identity of the topology. A submesh is the side domain of such a topology, and
the entity map of the submesh is the side table into the parent:

```python
# dx(1) and dx(2): the fluid cells and the solid cells, each seen as itself in square.
fluid_cells = StaticTopology((K,), (((0, 1, 2),),), ((("inside", square, 0),),))
solid_cells = StaticTopology((K,), (((1, 3, 2),),), ((("inside", square, 0),),))

# The fluid and the solid, each a mesh with its own topology.
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
the velocity of the solid over a time step. With Lagrange spaces of a library,
`V` and `Q` on `fluid_cells`, `W` on `solid_cells` and `L` on `interface`,
with `V`, `W` and `L` of the value shape `(2,)` and `Q` of the value shape
`()`:

```python
u, v = function(V, "u"), function(V, "v")
p, q = function(Q, "p"), function(Q, "q")
d, w, d_old = function(W, "d"), function(W, "w"), function(W, "d_old")
lam, mu = function(L, "lam"), function(L, "mu")

# The fields on the meshes.
u_F, v_F, p_F, q_F = (push_forward(f, F) for f in (u, v, p, q))
d_S, w_S, d_old_S = (push_forward(f, S) for f in (d, w, d_old))
lam_G, mu_G = push_forward(lam, G), push_forward(mu, G)

# Material constants, the time step and the inlet traction, arrays of the solver.
nu, mu_s, lam_s = runtime_array("nu", ()), runtime_array("mu_s", ()), runtime_array("lam_s", ())
dt, g = runtime_array("dt", ()), runtime_array("g", (2,))
i, j = index("i"), index("j")
half, two = Integer(1) / Integer(2), Integer(2)

# Stokes flow on F.
grad_u = gradient(u_F, F)
grad_v = gradient(v_F, F)
div_u, div_v = divergence(u_F, F), divergence(v_F, F)
fluid = integral(
    nu * index_sum(grad_u[i, j] * grad_v[i, j], (i, j)) - p_F * div_v - q_F * div_u, F
)

# Linear elasticity on S.
grad_d = gradient(d_S, S)
grad_w = gradient(w_S, S)
eps_d, eps_w = half * (grad_d[i, j] + grad_d[j, i]), half * (grad_w[i, j] + grad_w[j, i])
div_d, div_w = divergence(d_S, S), divergence(w_S, S)
solid = integral(two * mu_s * index_sum(eps_d * eps_w, (i, j)) + lam_s * div_d * div_w, S)

# The restrictions to the interface, pull-backs by the side maps, expressions on G.
from_fluid, from_solid = side_map(F, G, "fluid"), side_map(S, G, "solid")
u_G, v_G = pull_back(u_F, from_fluid), pull_back(v_F, from_fluid)
d_G, w_G = pull_back(d_S, from_solid), pull_back(w_S, from_solid)
d_old_G = pull_back(d_old_S, from_solid)

# The coupling: lam is the traction on both parts, and mu enforces u = (d - d_old) / dt.
coupling = integral(
    index_sum(lam_G[i] * (w_G[i] / dt - v_G[i]) + mu_G[i] * (d_G[i] / dt - u_G[i]), (i,)), G
)
a = form(fluid + solid + coupling, ((v, q, w, mu), (u, p, d, lam)))

# The right-hand side: the inlet traction and the old displacement.
traction = integral(index_sum(g[i] * pull_back(v_F, side_map(F, G_in, "fluid"))[i], (i,)), G_in)
l = form(traction + integral(index_sum(mu_G[i] * d_old_G[i], (i,)) / dt, G), ((v, q, w, mu),))
```

- `u` has the support `fluid_cells` and `d` the support `solid_cells`, so `u`
  and `d` do not combine, and neither do `u_F` and `d_S`. They meet only on
  `G`, through the restrictions, as two expressions on `G` at the same points
  of the reference facet `I`. Every integral has no support, so `fluid`,
  `solid` and `coupling` add into one expression, and the arguments of `a` are
  tuples of functions on three topologies (see [`Form`](LANGUAGE.md#form)).
- `side_domain(F, interface)` is the image of the trace of the coordinate
  function of `F` on `interface`. `normal(from_fluid)`
  is the normal out of the fluid, and `normal(from_solid)` is the negative
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
  submeshes when the fields live on `square`. With `c = cell_index(square, 0)`, an
  array `nu[c]` under `integral(..., D)` takes the value of the cell (see
  [`Integral` over a decomposed domain](LANGUAGE.md#integral-over-a-decomposed-domain)).
- A fluid domain that moves with the solid is the ALE mesh
  `image(x_F + m, RD(2))`, with `x_F = F.coordinate_function` and a mesh
  displacement `m` on `fluid_cells` that extends `d` across the interface. The
  functions on `fluid_cells` push forward to the moved mesh unchanged.

The passes rewrite `a` as in the [interior facet
example](#an-interior-facet-integral). `expand_integrals` gives one
`IndexSum` for each of `F`, `S` and `G`, over the cell index of its topology.
On `G`, with `e = cell_index(interface, 0)`, the restriction `u_G` lowers to
the pull-back of `u` by the cell map of the fluid cell `cells_fluid[e]` and
then by the inclusion of the facet into that cell, and `d_G` to the pull-back
of `d` by the cell map of the solid cell `cells_solid[e]`. The kernel of `coupling` is an expression on
`I` that reads the reference basis of `V` on one fluid cell, of `W` on one
solid cell and of `L` on the facet, and the assembler packs the two cells
through the two side tables, the entity maps of the two submeshes.

```{dropdown} Legacy UFL
The legacy UFL `dx(1)`, `ds(3)` and `ds((3, 4))` are integrals over the side domains of marked topologies, one topology for each marker value. A DOLFINx submesh is the side domain of a marked topology, and its entity map is the side table into the parent.
```

### Constants on a domain

Requires: [A fluid and a solid](#a-fluid-and-a-solid), [Spaces of
constants](LANGUAGE.md#spaces-of-constants).

With the names of the [fluid and solid example](#a-fluid-and-a-solid), `R` a
space of constants of value shape `()`, and `R_S` a space of a library on
`solid_cells` whose functions are constant on `solid_cells`:

```python
# Multipliers for the mean pressure of the fluid and the mean volume change of the solid.
r, s = function(R, "r"), function(R, "s")
k, z = function(R_S, "k"), function(R_S, "z")
k_S, z_S = push_forward(k, S), push_forward(z, S)
means = integral(p_F * s + q_F * r, F) + integral(k_S * div_w + z_S * div_d, S)
```

- `r` and `s` have no support, so the integrals decide where `r` and `s` act.
  `k_S` and `z_S` have the support `S`, so `k_S * q_F` raises `ValueError`.
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
e = integral(two * push_forward(mu_s, S) * index_sum(eps_d * eps_w, (i, j)), S)
de_dmu = form(directional_derivative(e, mu_s, mu_s, Integer(1)), (w, d))

# The gradient of the dissipation with respect to nu.
dissipation = integral(nu * index_sum(grad_u[i, j] * grad_u[i, j], (i, j)), F)
g_nu = form(directional_derivative(dissipation, nu, nu, h), (h,))
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

With `K3 = ReferenceSimplex(3)` and `I = ReferenceSimplex(1)`, a topology
`T3` of tetrahedra and its coordinate function `x3`, the topology `T1` of the
segments with the cell type `I` and the side `(("tissue", T3, 0),)`, Lagrange
spaces of a library `V3` on `T3` and `V1` on `T1` of the value shape `()`,
and `runtime_extent` and `runtime_array` of the [cut cell
example](#integrating-over-a-cut-cell):

```python
# The tissue, of dimension 3, and the network, of dimension 1, both in RD(3).
D3 = image(x3, RD(3))
D1 = side_domain(D3, T1)
identity = transformation_law((), 0)

p3, q3 = function(V3, "p3"), function(V3, "q3")
p1, q1 = function(V1, "p1"), function(V1, "q1")
beta = runtime_array("beta", ())
m = index("m")

# The Laplacian of each part on its own mesh. On D1 the gradient is the
# derivative along the segment, in the frame of RD(3).
grad_p3 = gradient(push_forward(p3, D3), D3)
grad_q3 = gradient(push_forward(q3, D3), D3)
grad_p1 = gradient(push_forward(p1, D1), D1)
grad_q1 = gradient(push_forward(q1, D1), D1)
tissue = integral(index_sum(grad_p3[m] * grad_q3[m], (m,)), D3)
vessel = integral(index_sum(grad_p1[m] * grad_q1[m], (m,)), D1)

# The tissue pressure on the network, a restriction on the topologies, and the
# exchange between the two, pushed forward to D1.
from_tissue = side_map(T3, T1, "tissue")
p3_1, q3_1 = pull_back(p3, from_tissue), pull_back(q3, from_tissue)
exchange = push_forward(beta * (p1 - p3_1) * (q1 - q3_1), D1, identity)
coupling = integral(exchange, D1)
a = form(tissue + vessel + coupling, ((q3, q1), (p3, p1)))
```

- The side tables of `T1` hold, for each segment `e`, one tetrahedron
  `cells_t[e]` that has the segment as an edge, chosen by the library among
  those that share the edge, and the inclusion of the segment into it, the
  affine map through the positions `positions_t[e, v, k]` of the two ends in
  `K3`, the coordinates of two vertices of `K3`. The tables come
  from the connectivity of the mesh, and a segment and its tetrahedron share
  vertices, so the partition of `T3` gives the partition of `T1`.
- `p3_1` is `pull_back(pull_back(p3, B_t), A_t)` on each segment, with
  $B_t$ the cell map of the tetrahedron of each segment and $A_t$ the
  inclusion of the side tables, the edge, by
  [Restriction](LANGUAGE.md#restriction), an expression on `T1` that combines
  with `p1`. The exchange is an expression on `T1`, which pushes forward to
  `D1` with the identity law. `p3` is continuous, so the restriction is the same from every
  tetrahedron that shares the edge, and the choice of the library does not
  change the form.
- `side_domain(D3, T1)` is the image of the trace of `x3` on `T1`, the edges
  of the tetrahedra, so the integrals over `D1` use the length element of the
  edge.
- `p3` and `p1` have the supports `T3` and `T1` and never combine directly.
  Every integral has no support, so the three terms add, and the arguments of
  `a` are tuples of functions on the two topologies.
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

The passes rewrite `a` as in the [interior facet
example](#an-interior-facet-integral). `expand_integrals` gives one
`IndexSum` for `tissue`, over `cell_index(T3, 0)`, and one for `vessel` and
one for `coupling`, over `cell_index(T1, 0)`. On `coupling`,
`expand_pull_backs` lowers the two restrictions to pull-backs by the cell map
of one tetrahedron, selected through the side table, and then by the affine
map through the positions. The kernel of `coupling` for a segment reads the
reference basis of `V1` on the segment and of `V3` on one tetrahedron, at the
points of the edge in the reference coordinates of each: the interface of a
facet kernel, with the cells of two meshes and of two dimensions. The local
matrix couples the degrees of freedom of `V1` on the segment with those of
`V3` on the tetrahedron, and the assembler scatters the matrix through the
maps of degrees of freedom of both.

The trace of a function in $H^1(\Omega)$ on a curve is not defined, so a
model may replace `p3` on the network by its average over a ring of radius
$R$ around each point of the network, in the plane of two unit normals `n1`
and `n2` of the segment, functions in a space `N` on `T1` of the value shape
`(3,)` that the solver provides. The ring points are not on the edges, so
the average is a family of maps from `T1` into `D3` and point evaluation,
with the cosines and sines of the angles and the weights as tables, because
the core has no `cos`:

```python
L = runtime_extent("L")
l = index("l")
R = runtime_array("R", ())
n1, n2 = function(N, "n1"), function(N, "n2")
cs, w = runtime_array("cs", (L, 2)), runtime_array("w", (L,))
x1 = D1.coordinate_function

# The ring of L points around each point of the network, a family of maps into D3.
ring = domain_map(
    T1, component_tensor(x1[m] + R * (cs[l, 0] * n1[m] + cs[l, 1] * n2[m]), (m,)), D3
)
p3_avg = index_sum(w[l] * pull_back(push_forward(p3, D3), ring), (l,))
q3_avg = index_sum(w[l] * pull_back(push_forward(q3, D3), ring), (l,))

exchange_avg = push_forward(beta * (p1 - p3_avg) * (q1 - q3_avg), D1, identity)
coupling_avg = integral(exchange_avg, D1)
```

- `ring` is a `DomainMap` with the parameter domain `T1`, a family with the
  free index `l`, and the ambient domain `D3`, a mesh.
  `pull_back(push_forward(p3, D3), ring)` has the support `T1` and the free
  index `l`, so `p3_avg` is an expression on `T1`. With `R = 0`, the
  pull-back is the trace of `p3` on the network, equal to `p3_1` for a
  conforming network.
- The pull-back by `ring` remains after the lowering stage, and the pass of
  [Points and point evaluation](LANGUAGE.md#points-and-point-evaluation)
  lowers it after a quadrature pass on `I` as in the [quadrature
  example](#quadrature). The maps below the push-forward are then `ring`,
  `cell_map(T1, 0)` and the quadrature points, a chain with the free
  indices of the segment, the quadrature point and the ring point, which the
  pass reads as one family of points. The side tables of that family hold
  the tetrahedron that contains each point and its position, and the kernel
  for a segment reads `V3` on one tetrahedron for each point, a list of
  cells rather than one cell. A kernel allows that, because the interface of
  a kernel follows from the expression.

### The Laplace-Beltrami operator

Requires: [`DecomposedDomain`](LANGUAGE.md#decomposeddomain),
[`Gradient`](LANGUAGE.md#gradient), [`PushForward`](LANGUAGE.md#pushforward),
[`Integral` over a decomposed domain](LANGUAGE.md#integral-over-a-decomposed-domain),
[`Form`](LANGUAGE.md#form).

See also: [`Inverse`](LANGUAGE.md#inverse).

The Poisson forms on a surface in $R^3$ are the forms of the Laplace–Beltrami
operator. A mesh library gives a topology `T` of $N$ triangles with the cell
type `K = ReferenceSimplex(2)`, and the coordinate function `x` in a vector
Lagrange space of degree 1 on `T` of the value shape `(3,)`. `V` is a Lagrange
space of a library on `T`:

```python
M = image(x, RD(3))

# Functions on T, and their fields on M.
u, v, f = function(V, "u"), function(V, "v"), function(V, "f")
u_M, v_M, f_M = push_forward(u, M), push_forward(v, M), push_forward(f, M)

# The surface gradients, of value shape (3,) in the frame of RD(3).
i = index("i")
grad_u = gradient(u_M, M)
grad_v = gradient(v_M, M)

a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), M), (v, u))
L = form(integral(f_M * v_M, M), (v,))
```

- `M` has the dimension 2 and the ambient coordinate domain `RD(3)`, so
  [`Gradient`](LANGUAGE.md#gradient) on `M` is the tangential gradient, and
  `grad_u` is the surface gradient of `u_M`.
- `integral(..., M)` multiplies the integrand on `K` by the density
  $\rho = \sqrt{\det(J^T J)}$, with $J$ the reference gradient of `x` on each
  cell, of value shape `(3, 2)`. The forms need no orientation of `M`.
- `expand_pull_backs` rewrites the pull-back of `grad_u` by the cells as
  $K^T \hat g$, with $K$ = `inverse_jacobian` of the cells, the pseudo-inverse
  `inverse(J)` of value shape `(2, 3)`, and $\hat g$ the gradient of the
  pull-back of `u` on `K`. The integrand of `a` on `K` is then
  $\hat g_u^T K K^T \hat g_v \, \rho$, and $K K^T = (J^T J)^{-1}$ is the
  inverse of the metric that `x` induces.
- With `P = gradient(coordinate(M), M)`, the projection onto the tangent
  plane, `gradient(grad_u, M)` is the tangential derivative $H$ of the surface
  gradient, which is not symmetric in general. The covariant Hessian of `u_M`
  is $P H P$.

```{dropdown} Legacy UFL
Legacy UFL writes `inner(grad(u), grad(v)) * dx` on a mesh of topological dimension 2 in $R^3$, and `grad` is the surface gradient because of the dimensions of the mesh. In UFLx, the same rule of `Gradient` gives the gradient in every dimension.
```

### A deformed domain

Requires: [`DecomposedDomain`](LANGUAGE.md#decomposeddomain),
[`PushForward`](LANGUAGE.md#pushforward), [Jacobian and
density](LANGUAGE.md#jacobian-and-density), [`Gradient`](LANGUAGE.md#gradient),
[`Form`](LANGUAGE.md#form), [Shape derivatives](LANGUAGE.md#shape-derivatives),
[The Laplace-Beltrami operator](#the-laplace-beltrami-operator).

Shape optimisation and updated Lagrangian mechanics write a form on the image
of the mesh under a deformation, and assemble on the mesh. The deformed mesh
is the image of moved coordinates. With `D`, `x` and `square` of the
[`DecomposedDomain`](LANGUAGE.md#decomposeddomain) example, a displacement `s`
in a vector space `W` of a library on `square`, and a Lagrange space `V` on
`square`:

```python
i = index("i")
s = function(W, "s")

# The deformed square, the image of the topology under the moved coordinates.
D_s = image(x + s, RD(2))

# The volume of the deformed square, a functional of s.
volume = integral(Integer(1), D_s)

# The Poisson stiffness form on the deformed square.
u, v = function(V, "u"), function(V, "v")
grad_u, grad_v = gradient(push_forward(u, D_s), D_s), gradient(push_forward(v, D_s), D_s)
a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), D_s), (v, u))
```

- `expand_integrals` rewrites `volume` as the sum over the cells of
  `integral(density(C), K)`, with `C = cell_domain(D_s, 0)`. `jacobian(C)` is
  the reference gradient of `pull_back(x + s, cell_map(square, 0))`, $J + \hat\nabla s$,
  so `density(C)` is $|\det(J + \hat\nabla s)|$, the Jacobian determinant of
  shape calculus with respect to the reference cell.
- `expand_pull_backs` rewrites the pull-back of `grad_u` by `C` as
  $K^T \hat\nabla \hat u$ with $K$ = `inverse(jacobian(C))`, the
  $(J + \hat\nabla s)^{-T}$ of shape calculus. `u` is a function on `square`,
  so `u` has the degrees of freedom of the topology and moves with the mesh.
- The derivatives of the core hold the domains fixed, so the shape derivative
  of `volume` in a direction `h` on `square` is `ShapeDerivative(volume, s, h)`
  of the [extension](LANGUAGE.md#shape-derivatives) (see [A shape
  derivative](#a-shape-derivative)). After a lowering stage,
  `s` occurs explicitly in the Jacobians and the densities, and the node
  lowers to `DirectionalDerivative(volume, s, s, h)`. With an argument as `h`,
  the Gateaux derivative gives the shape gradient.

A surface deforms in the same way. With `M`, `x` and `T` of the
[Laplace-Beltrami example](#the-laplace-beltrami-operator) and `w` a
displacement in a vector space on `T` of the value shape `(3,)`:

```python
w = function(W_T, "w")
M_w = image(x + w, RD(3))

# The area of the deformed surface.
area = integral(Integer(1), M_w)
```

`jacobian(cell_domain(M_w, 0))` is the reference gradient of `x + w`, of value
shape `(3, 2)`, and the density is $\sqrt{\det(J^T J)}$ (see [Jacobian and
density](LANGUAGE.md#jacobian-and-density)).

```{dropdown} Isogeometric analysis
Isogeometric analysis by extraction, as in tIGAr, builds a mesh of the
parametric knot spans, stores the NURBS map $F$ as a function in the spline
space on that mesh, and inserts $\det DF$ and $DF^{-T}$ into every form by
hand. Here the physical patch is `image(F, RD(3))`, with $F$ a function on
the topology of the knot spans, a multi-patch geometry is the image of the
topology of the knot spans of all the patches, and the passes insert the
factors.
```

```{dropdown} Shells
A shell code such as FEniCS-Shells defines the midsurface as a map $\phi_0$
from a flat parametric mesh into $R^3$ and forms the fundamental forms from
$\nabla \phi_0$ by hand. Here the midsurface is `image(phi_0, RD(3))`, with
$\phi_0$ a function on the topology of the flat mesh, $k = 2 < n = 3$,
`jacobian(C)` is the reference gradient of $\phi_0$, its Gram matrix $J^T J$
is the first fundamental form in the reference coordinates, and `density(C)`
is $\sqrt{\det(J^T J)}$, the area element of the midsurface.
```

### A shape derivative

Requires: [Shape derivatives](LANGUAGE.md#shape-derivatives),
[`DirectionalDerivative`](LANGUAGE.md#directionalderivative), [A deformed
domain](#a-deformed-domain).

With `square` of the [static topology example](#a-static-topology), a vector
Lagrange space `X` of degree 1 and a Lagrange space `V`, both on `square`:

```python
x, h, v = function(X, "x"), function(X, "h"), function(X, "v")
D = image(x, RD(2))
grad_u = gradient(push_forward(function(V, "u"), D), D)
i = index("i")
energy = integral(index_sum(grad_u[i] * grad_u[i], (i,)), D)

d_energy = shape_derivative(energy, x, h)
g_energy = form(shape_derivative(energy, x, v), (v,))
```

`g_energy` is the shape gradient. The pass of the extension rewrites
`ShapeDerivative(g, x, h)` as `directional_derivative(g, x, x, h)`, after a
lowering stage and before `apply_derivatives` and a second lowering stage.
Written before the lowering stage, `directional_derivative(energy, x, x, h)`
is `Zero`, because `x` occurs only inside `D` and the substitution does not
enter domains. After the lowering stage, `x` occurs in the Jacobians and the
densities.

The whole extension is short. In a prototype it is about 60 lines, almost
all of them the node class, and the pass is the one rule above, because the
directional derivative of the core does the work.

### The flux through a surface

Requires: [The Laplace-Beltrami
operator](#the-laplace-beltrami-operator), [Concept: integral over a
domain](LANGUAGE.md#concept-integral-over-a-domain), [Jacobian and
density](LANGUAGE.md#jacobian-and-density),
[`AbstractTopology`](LANGUAGE.md#abstracttopology).

The flux of a vector field through a surface depends on which side of the
surface is out, so the flux needs a consistent orientation of the cells. With
`K`, `T`, `x` and `M` of the [Laplace-Beltrami
example](#the-laplace-beltrami-operator), where `M` is a closed surface such
as a triangulated sphere:

```python
# The triangles of M, a family with the free index cell_index(T, 0).
cells = cell_domain(M, 0)
J = jacobian(cells)
r, l = index("r"), index("l")

def cross(a, b):
    return tensor((a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]))

# The normal of each triangle, with the length of the density.
nu = cross(component_tensor(J[r, 0], (r,)), component_tensor(J[r, 1], (r,)))

# The sign of each triangle, +1 where the order of the vertices of the
# triangle agrees with the outward orientation of M.
o = pull_back(T.orientation, cell_map(T, 0))

# The outward flux of the position field, 3 times the volume inside M.
flux = index_sum(
    integral(o * index_sum(pull_back(coordinate(M), cells)[l] * nu[l], (l,)), K),
    (cell_index(T, 0),),
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
- `o` is the orientation of `T` on each triangle (see
  [`AbstractTopology`](LANGUAGE.md#abstracttopology)). With the signs, every
  contribution uses the outward normal. A mesh library may instead number the
  vertices of every triangle consistently, and then the flux needs no `o`.
- `area` integrates against the density $\sqrt{\det(J^T J)}$, which does not
  depend on the order of the vertices, so `area` needs no signs.
- A surface with no consistent orientation, such as a Möbius strip, has the
  orientation `None`. The area of the surface is defined, and the flux is not.

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
`square`, `F` a real Lagrange space on `square` for the data, and `C` a real
space of constants. The solver sets `sigma_x` and `sigma_y` to zero outside
the layer:

```python
k0 = function(C, "k0")
sigma_x, sigma_y = push_forward(function(F, "sigma_x"), D), push_forward(function(F, "sigma_y"), D)
one, i_unit = Integer(1), ImaginaryUnit()
s_x = one + i_unit * sigma_x / k0
s_y = one + i_unit * sigma_y / k0

u, v, f = function(V, "u"), function(V, "v"), function(F, "f")
u_D, v_D, f_D = push_forward(u, D), push_forward(v, D), push_forward(f, D)
grad_u, grad_v = gradient(u_D, D), gradient(v_D, D)
j = index("j")

# The stretched Laplacian and the stretched mass, sesquilinear in u and v.
A = tensor((s_y / s_x, s_x / s_y))
stiffness = index_sum(A[j] * grad_u[j] * conjugate(grad_v[j]), (j,))
mass = k0 * k0 * s_x * s_y * u_D * conjugate(v_D)
a = form(integral(stiffness - mass, D), (conjugate(v), u))
L = form(integral(f_D * conjugate(v_D), D), (conjugate(v),))
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

Requires: [`TransformationLaw`](LANGUAGE.md#transformationlaw),
[`PushForward`](LANGUAGE.md#pushforward), [`Divergence` and
`Curl`](LANGUAGE.md#divergence-and-curl), [`Form`](LANGUAGE.md#form).

The mixed form of the Poisson problem seeks a flux $\sigma$ and a potential
$u$ with $\sigma = \nabla u$ and $\mathrm{div}\, \sigma = -f$. `RT` is a
Raviart-Thomas space on the topology `square` of the
[`DecomposedDomain`](LANGUAGE.md#decomposeddomain) example, of value shape
`(2,)`, and `DG` a discontinuous Lagrange space of degree 0 on `square`:

```python
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

- There is no mixed space. Each argument is a `tuple` of a function in `RT`
  and a function in `DG`, so the argument spaces are `((RT, DG), (RT, DG))`
  and `((RT, DG),)`. The block of `a` for `tau` and `u` is `a` with `v` and
  `sigma` replaced by `Zero`, the term `u_D * divergence(tau_D, D)`. `L`
  depends on `v` only, which a `tuple` argument allows.
- The law of `RT` is `transformation_law(("contravariant",), 1)`, so
  `pull_back(sigma_D, cell_domain(D, 0))` is $J \hat\sigma / \det J$ with
  $\hat\sigma$ = `pull_back(sigma, cell_map(square, 0))`. The law of `DG` is
  `transformation_law((), 0)`.
- `expand_pull_backs` applies the identity of Piola to
  `divergence(sigma_D, D)`, which gives
  $\widehat{\mathrm{div}}\, \hat\sigma / \det J$ on every cell, affine or
  curved, with no derivative of $J$.
- A Nédélec space has the law `transformation_law(("covariant",), 0)`, and
  `curl` of its push-forward is the reference curl carried by the
  contravariant law of weight 1.

### Dual-weighted residual error estimation

Requires: [`Form`](LANGUAGE.md#form), [`Action`](LANGUAGE.md#action),
[`DualSpace`](LANGUAGE.md#dualspace), [A Bank-Weiser error
estimator](#a-bank-weiser-error-estimator).

The dual-weighted residual method estimates the error in a goal functional
$J(u)$ by the residual of the Galerkin solution weighted by the solution $z$
of the adjoint problem, $J(u) - J(u_h) = R(z) = R(z - z_h)$ for any $z_h$ in
the trial space. The example takes an advection-diffusion problem on the
square mesh `D`, with `V` a Lagrange space of degree 1 on `square`, `W` a
Lagrange space of degree 2 on `square` in which the adjoint solution is
computed, `F` a space on `square` for the data and `B` a vector space on
`square` for the advection field:

```python
u, v, u_h, z_h = function(V, "u"), function(V, "v"), function(V, "u_h"), function(V, "z_h")
z = function(W, "z")
f, psi = push_forward(function(F, "f"), D), push_forward(function(F, "psi"), D)
b = push_forward(function(B, "b"), D)
i = index("i")

def a_expr(w, y):
    w_D, y_D = push_forward(w, D), push_forward(y, D)
    grad_w, grad_y = gradient(w_D, D), gradient(y_D, D)
    return index_sum(grad_w[i] * grad_y[i], (i,)) + index_sum(b[i] * grad_w[i], (i,)) * y_D

# The primal problem, the goal, and the adjoint problem with the arguments swapped.
v_D = push_forward(v, D)
a = form(integral(a_expr(u, v), D), (v, u))
L = form(integral(f * v_D, D), (v,))
goal = form(integral(psi * v_D, D), (v,))
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
  so `z` from the richer space `W` is accepted, and so is `z - z_h`. The
  action substitutes `z - z_h` for `v` in `push_forward(v, D)`, which gives
  the push-forward of an expression by the law of `V`. Galerkin
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
Lagrange space `V` and a real space `M` for the Lamé parameters, both of a
library on `square`, and a real space `R` of constants of value shape
`(2, 2)`. The stress is
$S = 2 \operatorname{sym}(\partial \psi / \partial C)$:

```python
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

# The derivative with respect to a symbol at C, in the symmetric unit directions.
C_ = function(R, "C")
N = component_tensor((I[i, k] * I[j, l] + I[i, l] * I[j, k]) / two, (i, j))
S = two * component_tensor(directional_derivative(strain_energy(C_), C_, C, N), (k, l))
sigma = F @ S @ transpose(F) / determinant(F)

Pi = integral(psi, D)
residual = gateaux_derivative(Pi, u, u, v)
tangent = gateaux_derivative(residual, u, u, w)
```

- `S` and the Cauchy stress `sigma` are expressions on `D`, which a solver
  interpolates or evaluates for output. `apply_derivatives` evaluates `S` to
  an expression equal to $\lambda \operatorname{tr}(E) I + 2 \mu E$.
- `V` is real, so `residual` is the form of `directional_derivative(Pi, u, u, v)` that binds `v`, and `tangent` binds `w` after `v`.

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

grad_u, grad_v = gradient(push_forward(u, D), D), gradient(push_forward(v, D), D)
eps_u = (grad_u + transpose(grad_u)) / two
eps_v = (grad_v + transpose(grad_v)) / two
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
[`Function`](LANGUAGE.md#function), [`PushForward`](LANGUAGE.md#pushforward),
[`Form`](LANGUAGE.md#form),
[Restriction](LANGUAGE.md#restriction), [`Normal`](LANGUAGE.md#normal),
[`Integral` over a decomposed
domain](LANGUAGE.md#integral-over-a-decomposed-domain).

The Poisson problem, written as classic FEniCS writes it:
`inner(grad(u), grad(v)) * dx`. The layer is a module `uflx.legacy`, which
does not exist, of a page of free functions and three plain Python classes
over the core. Nothing in the module subclasses a core class, so the forms
that the module produces are ordinary core graphs, the same graphs as the core
free functions build by hand. The example shows what the core leaves to a
layer: the argument numbers of test and trial functions, one name for a
function and its push-forward to the mesh, measures that find or hold their
domain, restrict the integrand and bind the arguments, a facet normal that
belongs to the mesh, and `grad`, `div` and `curl` that read the domain from
the function. `Lagrange`, `Constants` and `Fields` stand for spaces of a
library: Lagrange elements of a value shape on a topology, constants with no
domain, and fields of a value shape on a mesh with no element.

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

- `u` and `v` of the script are push-forwards to `D` of the functions `u_T`
  and `v_T` on the topology of `D`, so a user script has one name for each.
  For a real `V`, `a` is `kappa` times the Poisson bilinear form
  `form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), D), (v_T, u_T))`,
  with `grad_u = gradient(push_forward(u_T, D), D)`, `grad_v` likewise and
  `i = index("i")`, up to the labels of `u_T`, `v_T` and `i` and the
  conjugates, which change nothing in real arithmetic: with `S = V.space`,
  `conjugate_space(S)` is `S`, so the argument spaces are `(S, S)`. `L` is a
  sum of three forms with the argument spaces `(S,)`.
  `f * dx` with no argument stays an `Integral`, a functional.
- For a complex `V`, the same lines give the sesquilinear forms of a
  complex-mode legacy demo. `inner` conjugates the second operand, as the
  legacy `inner` does, and the measure binds `conjugate(v_T)` at position 0,
  so `a` has the argument spaces `(conjugate_space(S), S)` and is antilinear
  in `v_T`, and `L` has `(conjugate_space(S),)`. The convention is the module's,
  not the core's: the core records the linearity that the form states (see
  [`Form`](LANGUAGE.md#form)). `f * v * dx` in complex mode is linear in `v`
  but bound as antilinear, which a pass that checks the linearity rejects, as
  legacy UFL rejects an unconjugated test function in complex mode.
- The argument numbers of legacy UFL live in the labels of the module. The
  core has none: `u_T` and `v_T` are functions, and the measure states the
  arguments, the functions on the topology. `arguments` is the legacy `extract_arguments`, a walk of the
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
  so `grad(u)` becomes the pull-back of `grad(u)` by the side map and not the
  derivative of a restriction, and the rebuilding of the other nodes uses
  rule 4 of [`Base`](LANGUAGE.md#the-base-object-base). `dx` with the mesh
  as domain restricts nothing, since the domain of the measure is the support
  of the integrand. A `dS` would be a
  measure over the interior facets with a side for each restriction, the
  legacy `u('+')` would be `pull_back(u, side_map(D, E, "+"))`, and `avg` and `jump`
  would be arithmetic, as in the [interior facet
  example](#an-interior-facet-integral).
- The legacy `FacetNormal` belongs to the mesh and is resolved per integral.
  `n` is a placeholder function of the module in a space of fields on the
  mesh, known by its label, so that `inner(q, n)` is an expression on the
  mesh like any other, and `restrict` replaces `n` by `normal(side_map(D, E, "inside"))`
  on the domain `E` of the measure. A subtree that contains `n` is rebuilt
  rather than restricted, because the normal has no meaning on the mesh.
- `grad` reads the domain from the support of the function, as the legacy
  `grad` reads the mesh from the coefficient, and `grad(w)[i, j]` is the
  derivative of `w[i]` along `x[j]`, as in legacy UFL. `div` and `curl` are
  `Divergence` and `Curl` of the core, so a Raviart–Thomas or a Nédélec
  coefficient gets the identity of Piola. On the surface mesh of
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


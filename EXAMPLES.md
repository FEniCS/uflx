# Worked examples for the UFLx core

These examples accompany [LANGUAGE.md](LANGUAGE.md) and use its names. The
first part defines the library classes that the examples and the tour of
LANGUAGE.md use; none of them is part of the core. The second part exercises
the core on small problems with known answers, and then on five larger
problems: an interior facet integral, a fluid and a solid, a 1D network in a
3D mesh, the Laplace-Beltrami operator and the flux through a surface. A
reader of LANGUAGE.md needs none of this to follow the rules.

## Reference library

### Example: the unit simplex

Requires: [AbstractCoordinateDomain](LANGUAGE.md#abstractcoordinatedomain).

`UnitSimplex` is not part of the core. The examples use `UnitSimplex` as a
concrete `AbstractCoordinateDomain`.

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

### Example: the simplex as a cell type

> LLM generated content.

Requires: [AbstractCellType](LANGUAGE.md#abstractcelltype).

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

### Example: a static topology

> LLM generated content.

Requires: [AbstractTopology](LANGUAGE.md#abstracttopology),
[Example: the simplex as a cell type](#example-the-simplex-as-a-cell-type),
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

### Example: the unit simplex as a reference cell

> LLM generated content.

Requires: [Example: the unit simplex](#example-the-unit-simplex),
[Example: the simplex as a cell type](#example-the-simplex-as-a-cell-type),
[Literals](LANGUAGE.md#literals), [Tensor](LANGUAGE.md#tensor), [Reference cell](LANGUAGE.md#reference-cell).

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
            return Zero((1, 0), (), (), ())
        return tensor(
            tuple(
                tensor(tuple(Integer(1 if k == v - 1 else 0) for k in range(n)))
                for v in range(n + 1)
            )
        )
```

`ReferenceSimplex(2)` and `ReferenceSimplex(1)` have the vertices of the unit
triangle and the unit interval in the [Reference cell](LANGUAGE.md#reference-cell)
examples.

`ReferenceSimplex(2)` and `Simplex(2)` are different cell types, so a topology
of `Simplex(2)` cells has no reference cells. `ReferenceSimplex(2)` and
`UnitSimplex(2)` are the same set but different coordinate domains, so an
expression on one is not an expression on the other. Likewise,
`ReferenceSimplex(0)` and `RD(0)` are different coordinate domains.

## Examples

### Example: transpose and matrix product

Requires: [Index](LANGUAGE.md#index), [Indexed](LANGUAGE.md#indexed), [ComponentTensor and
IndexSum](LANGUAGE.md#componenttensor-and-indexsum), [Product](LANGUAGE.md#product).

`transpose` and `matrix_product` are not part of the core. A library provides
the two as free functions built with derived indices (see [New
operations](LANGUAGE.md#new-operations)), and the examples of this note use the two:

- `transpose(A)` is `component_tensor(A[k, l], (l, k))`, with
  $k$ = `Index(("library.transpose", A, 0))` and
  $l$ = `Index(("library.transpose", A, 1))`.
- `matrix_product(A, B)` is
  `component_tensor(index_sum(A[i, k] * B[k, j], (k,)), (i, j))`, with $i$, $k$
  and $j$ the derived indices `Index(("library.matrix_product", A, B, m))` for
  $m = 0, 1, 2$.

An operand that does not have rank 2, including a scalar operand, is rejected
with `IndexError` by the constructor of `Indexed`. Unequal inner extents in
`matrix_product` are rejected with `ValueError` by the constructor of `Product`,
because the index $k$ has two extents. So `matrix_product` accepts only rank 2
operands. Contractions of other ranks, such as a matrix-vector product, are
written in index notation. The core defines no overload for the matrix
product, because the core has no free function for the matrix product (see
[Operator overloads](LANGUAGE.md#operator-overloads)).

Examples, with a matrix `A` of value shape `(2, 2)`:

- `transpose(A)` is not equal to `component_tensor(A[i, j], (j, i))` for user
  indices `i` and `j`, because the index labels differ.

### Example: a paraboloid

Requires: [RD](LANGUAGE.md#rd), [Coordinate](LANGUAGE.md#coordinate), [Tensor](LANGUAGE.md#tensor),
[Indexed](LANGUAGE.md#indexed), [DomainMap](LANGUAGE.md#domainmap).

The paraboloid over `RD(2)` collects the examples of the sections above.

```python
# AbstractCoordinateDomain: RD(2), no edges in the graph.
P = RD(2)

# AbstractExpression: Coordinate(RD(2)), value shape (2,).
X = coordinate(P)

# AbstractExpression: Tensor((Indexed(X, (0,)), Indexed(X, (1,)), Sum(...))),
# value shape (3,).
psi = tensor((X[0], X[1], X[0]**Integer(2) + X[1]**Integer(2)))

# DomainMap(RD(2), psi, RD(3)), k = 2, n = 3.
D = domain_map(P, psi, RD(3))
```

### Example: an affinely mapped triangle

Requires: [Example: the unit simplex](#example-the-unit-simplex),
[Coordinate](LANGUAGE.md#coordinate), [Tensor](LANGUAGE.md#tensor), [Components](LANGUAGE.md#components),
[DomainMap](LANGUAGE.md#domainmap).

The physical triangle with vertices $v_0 = (1, 0)$, $v_1 = (3, 1)$ and
$v_2 = (2, 3)$ is the image of the unit triangle under an affine map.

There is no topology here, only geometry.

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
$1/2$ of the unit triangle.

The physical triangle is the image of a domain map with $k = n$, not a
coordinate domain, although the physical triangle has non-empty interior in
$R^2$. The class records how the domain is built.

### Example: orientation of an edge and of a triangle

Requires: [Concept: orientation of a domain](LANGUAGE.md#concept-orientation-of-a-domain),
[Concept: integral over a domain](LANGUAGE.md#concept-integral-over-a-domain),
[Components](LANGUAGE.md#components), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

The examples use `T`, `X`, `one`, `S`, `t` and `e` of the
[PullBack](LANGUAGE.md#pullback) example, and `K` of the
[affinely mapped triangle example](#example-an-affinely-mapped-triangle).

The edge `e`, with $\psi(t) = (1 - t, t)$, runs from vertex 1 to vertex 2 of
the unit triangle. The edge with $\psi(t) = (t, 1 - t)$ is the same set with
the opposite orientation, so the integral of $x \, \mathrm{d}y$ changes sign:

```python
x, y = X[0], X[1]
dy = exterior_derivative(y, T)

# The edge e reversed, the same set with the opposite orientation.
e_reversed = domain_map(S, tensor((t, one - t)), T)

# A later pass evaluates the two integrals to 1/2 and -1/2.
integral(pull_back(x * dy, e), S)
integral(pull_back(x * dy, e_reversed), S)
```

On `e_reversed`, $x = t$ and $\mathrm{d}y$ pulls back to $-\mathrm{d}t$, so
the integral is $\int_0^1 -t \, \mathrm{d}t = -1/2$.

The triangle `K` has $\det J = 5$, so the orientation of `K` agrees with the
orientation of `RD(2)`:

```python
# The pull-back of the volume form of RD(2), det J dX_0 ^ dX_1, a 2-form on T.
w = pull_back(volume_form(RD(2)), K)

# The 2-form as the antisymmetric matrix [[0, det J], [-det J, 0]] of value shape (2, 2).
W = components(components(w, T), T)

# Entry [0, 1] follows the order of dX_0 ^ dX_1, so the entry is det J and not -det J.
# A later pass evaluates det J to 5.
det_J = W[0, 1]

# The sign of det J, 1 here.
s = det_J / absolute(det_J)

# The area of K. A later pass evaluates the integral to 5/2, whatever the sign of det J.
integral(s * w, T)
```

Swapping $v_1$ and $v_2$ gives $\det J = -5$ and $s = -1$. The integral of `w`
alone is then $-5/2$, and the integral of `s * w` is still $5/2$.

### Example: integrating over an edge

Requires: [PullBack](LANGUAGE.md#pullback), [Integral](LANGUAGE.md#integral),
[Components](LANGUAGE.md#components), [Elementary functions (very
incomplete)](LANGUAGE.md#elementary-functions-very-incomplete).

The integral over the edge `e` of the [PullBack](LANGUAGE.md#pullback) example is the
integral over `S` of the pull-back by `e`. The examples use `T`, `X`, `one`,
`S`, `t`, `e` and `dx` of that example.

```python
x, y = X[0], X[1]
dy = exterior_derivative(y, T)
```

A 1-form integrates over the edge without a metric. On the edge,
$x = 1 - t$ and $\mathrm{d}y$ pulls back to $\mathrm{d}t$, so

```python
# A later pass evaluates the integral to 1/2.
integral(pull_back(x * dy, e), S)
```

is $\int_0^1 (1 - t) \, \mathrm{d}t = 1/2$.

The flux of a vector field $u$ through the edge is the integral of the 1-form
$u_0 \, \mathrm{d}y - u_1 \, \mathrm{d}x$, which also needs no metric. For
$u = X$, the 1-form pulls back to $(1 - t) \, \mathrm{d}t + t \, \mathrm{d}t$:

```python
# A later pass evaluates the integral to 1.
integral(pull_back(x * dy - y * dx, e), S)
```

The flux is 1. The divergence of $X$ is 2 and the area of `T` is $1/2$, so
the flux of $X$ out of `T` is 1. The 1-form $x \, \mathrm{d}y - y \,
\mathrm{d}x$ is zero on the other two edges of `T`, so the whole flux passes
through `e`. This is the divergence theorem, without a normal vector and
without a metric.

The integral of a scalar function with respect to arc length needs the length
element $|\psi'| \, \mathrm{d}t$, which comes from the Euclidean metric of
`T`. The metric enters as an expression that the user writes:

> This example writes the length element by hand, for a curve in the
> Euclidean plane, because `Integral` over a coordinate domain takes the
> density from the user. The general concept is the volume form of the metric
> that $\psi$ induces, $\sqrt{\det(J^T J)} \, \mathrm{d}t$. For an `Integral`
> over a decomposed domain, `expand_integrals` builds that density for each
> piece (see [Integral over a decomposed
> domain](LANGUAGE.md#integral-over-a-decomposed-domain)).

```python
# The Jacobian matrix of the edge, with value shape (2, 1).
J = components(exterior_derivative(e.parameterisation, S), S)

# The length element of the edge.
ds = square_root(J[0, 0] * J[0, 0] + J[1, 0] * J[1, 0]) * volume_form(S)

# A later pass evaluates the integral to the length sqrt(2) of e.
integral(ds, S)

# A later pass evaluates the integral to sqrt(2) / 6.
integral(pull_back(x * y, e) * ds, S)
```

### Example: integrating over the boundary of a triangle

Requires: [DomainMap](LANGUAGE.md#domainmap),
[ComponentTensor and IndexSum](LANGUAGE.md#componenttensor-and-indexsum),
[Example: integrating over an edge](#example-integrating-over-an-edge).

The three edges of `T` are one family of domains, with a free index. The
example uses the names of
[Example: integrating over an edge](#example-integrating-over-an-edge).

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
index_sum(integral(pull_back(x * dy, edges), S), (i,))
```

The result is the integral of $x \, \mathrm{d}y$ over the boundary of `T`,
which is $\int_T \mathrm{d}x \wedge \mathrm{d}y = 1/2$ by Stokes' theorem.
Without the `IndexSum`, the integral has the free index `i`, and the value
for each `i` is the integral over one edge.

The orientation of each edge comes from the order of `starts` and `ends`. Every
domain of a family has the same parameter domain, here `S`.

### Example: quadrature

> This example is inspired by GEM and TSFC, but is written in the primitives
> of the core. Combining these ideas gives UFLx additional expressivity.
> Quadrature in UFLx may make TSFC simpler, if TSFC eventually targets UFLx.

Requires: [Points and point evaluation](LANGUAGE.md#points-and-point-evaluation),
[DomainMap](LANGUAGE.md#domainmap), [PullBack](LANGUAGE.md#pullback),
[Integral](LANGUAGE.md#integral), [Components](LANGUAGE.md#components), [Tensor](LANGUAGE.md#tensor),
[ComponentTensor and IndexSum](LANGUAGE.md#componenttensor-and-indexsum),
[Extents](LANGUAGE.md#extents).

This example looks ahead to quadrature in the core language.

A quadrature rule replaces the integral by the inner product of the weights
with the values at the points. With `T`, `X` and `f = X[0] * X[1]` of
[Points and point evaluation](LANGUAGE.md#points-and-point-evaluation), the edge midpoint
rule is:

```python
zero, half = Integer(0), Integer(1) / Integer(2)
sixth = Integer(1) / Integer(6)
i, k = index("i"), index("k")

points = tensor((tensor((half, zero)), tensor((half, half)), tensor((zero, half))))
weights = tensor((sixth, sixth, sixth))

# A family of three points in T, with the free index i.
xs = domain_map(RD(0), component_tensor(points[i, k], (k,)), T)

q = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

The points are one domain with the free index `i`, so `f` occurs once in `q`,
and a form compiler can evaluate `q` as a loop over the points.

A later pass evaluates `q` to $1/24$, which is
`integral(f * volume_form(T), T)`, because the rule is exact for $f$. The two
expressions are not equal, so quadrature is a pass, and the rule is an
argument of the pass. The pass finds the coefficient $c$ of a top form
$c \, \mathrm{d}X_0 \wedge \mathrm{d}X_1$ on `T` as
`components(components(c * volume_form(T), T), T)[0, 1]`. The result is an
expression of the core.

After the passes, `q` is

```python
index_sum(weights[i] * (points[i, 0] * points[i, 1]), (i,))
```

`expand_pull_backs` replaces `X` by the parameterisation of `xs`. A pass
replaces the integral over `RD(0)` of a 0-form by the 0-form, and
`reduce_component_tensors` replaces `component_tensor(points[i, k], (k,))[0]`
by `points[i, 0]` (see [Reducing component
tensors](LANGUAGE.md#reducing-component-tensors)). A form compiler could generate:

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
(see [Example: integrating over a cut
cell](#example-integrating-over-a-cut-cell)).
The number of points is known only at assembly, and the loop bound `q` is an
input of the kernel.

### Example: integrating over a cut cell

Requires: [Example: quadrature](#example-quadrature), [Extents](LANGUAGE.md#extents).

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

With `T`, `f`, `i` and `k` of the quadrature example:

```python
# These are language extensions, not the core.
q = runtime_extent("q")

points = runtime_array("points", (q, 2))
weights = runtime_array("weights", (q,))

# The family of q points in the part of T where phi < 0. The level set phi
# is not written here, because the core has no conditional.
xs = domain_map(RD(0), component_tensor(points[i, k], (k,)), T)

a = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

For all the cells of a topology, see the open item on ragged indexing in
[Extents](LANGUAGE.md#extents).

### Example: the centroid of a triangle

> LLM came up with this example.

Requires: [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle), [PullBack](LANGUAGE.md#pullback),
[Integral](LANGUAGE.md#integral), [Division](LANGUAGE.md#division),
[ComponentTensor and IndexSum](LANGUAGE.md#componenttensor-and-indexsum),
[Points and point evaluation](LANGUAGE.md#points-and-point-evaluation).

An integral is an expression, so the centroid of `K` is a ratio of integrals.
The example uses `T` and `K` of the affinely mapped triangle.

```python
Y = coordinate(RD(2))
vol = volume_form(RD(2))
k = index("k")

# A later pass evaluates the area of K to 5/2.
area = integral(pull_back(vol, K), T)

# A later pass evaluates the centroid of K to (2, 4/3).
centroid = component_tensor(integral(pull_back(Y[k] * vol, K), T) / area, (k,))
```

The pull-back of `vol` by `K` gives the factor $\det J = 5$, without a
determinant. The centroid is the parameterisation of the point `C` of
[Points and point evaluation](LANGUAGE.md#points-and-point-evaluation), which is the
image of the centroid of `T`. The two expressions are not equal, and a later
pass evaluates both to $(2, 4/3)$.

> Geometry stops here until topology (cells, facets, families, orientations)
> and functions are defined. Later sections return to geometry.

### Example: an interior facet integral

Requires: [Sides of a topology](LANGUAGE.md#sides-of-a-topology),
[Restriction](LANGUAGE.md#restriction), [Integral over a decomposed
domain](LANGUAGE.md#integral-over-a-decomposed-domain), [Normal](LANGUAGE.md#normal),
[Form](LANGUAGE.md#form), [Core transformations](LANGUAGE.md#core-transformations).

The interior facet terms of the symmetric interior penalty method for the
Poisson problem, the legacy UFL `dS` integral, on the square divided into two
triangles. The example uses `D`, `square`, `K`, `X`, `y` and `cell_vertices`
of the [DecomposedDomain](LANGUAGE.md#decomposeddomain) example, `I = ReferenceSimplex(1)`
and the `diagonal` topology of the [Restriction](LANGUAGE.md#restriction) example, a
space `V` on `D` of the degree `()` and the value shape `()`, and functions
`u` and `v` in `V`. `runtime_array` is the language extension of [Example:
integrating over a cut cell](#example-integrating-over-a-cut-cell).

```python
# The interior facets, seen from the sides "+" and "-", with the free index e.
E = side_domain(D, diagonal)
e = cell_index(diagonal, 0)

# The gradients on D, in the frame of RD(2).
i = index("i")
grad_u = components(exterior_derivative(u, D), D)
grad_v = components(exterior_derivative(v, D), D)

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
  `integral(rho_p * volume_form(K), K)`, with `rho_p` the density of $B_+$ of
  [Restriction](LANGUAGE.md#restriction) built as in [Integral over a decomposed
  domain](LANGUAGE.md#integral-over-a-decomposed-domain), and a per-cell array `h` of the
  solver is read as `h[gathered_index(cells_p, (e,))]`.
- `integral(f, E)` is one `Integral` over the decomposed domain `E`. On a mesh
  of triangles and quadrilaterals, the facet topology has three groups (see
  [Sides of a topology](LANGUAGE.md#sides-of-a-topology)), `E` has three pieces, and the
  same lines give the integral over all three.

The core transformations rewrite `a_facet` as follows. With $B_\pm$ =
`gathered_cell_domain(D, 0, gathered_index(cells_pm, (e,)))` the cell on each
side and $A_\pm$ = `domain_map(I, affine_map(I, points_pm), K)` the inclusion
of the facet into that cell, both of [Restriction](LANGUAGE.md#restriction), with the side
tables `(cells_pm, positions_pm)` of each side:

1. `expand_integrals` rewrites `integral(f, E)` as
   `index_sum(integral(pull_back(f, cell_domain(E, 0)) * rho * volume_form(I), I), (e,))`,
   with `rho` the length of the facet, from the Jacobian of the parameterisation
   `pull_back(pull_back(coordinate(D), B_+), A_+)` of the piece.
2. `expand_pull_backs` moves the pull-back by `cell_domain(E, 0)` into `f`.
   On `u_p` the rule for `SideMap` gives `pull_back(pull_back(u, B_+), A_+)`.
   The inner pull-back remains, for the library of `V`, as `u` on the cell
   `cells_p[e]` (see [AbstractFunctionSpace](LANGUAGE.md#abstractfunctionspace)), and
   the outer pull-back remains outside the inner one, by the rule that
   pull-backs nest (see [Expanding pull-backs](LANGUAGE.md#expanding-pull-backs)).
   On `du_p` the same rule gives the pull-back of `components(du, D)` by
   $B_+$, which the rule for `Components` rewrites as
   $K^T \hat g$, with $K$ = `inverse(J_+)` the inverse Jacobian of the cell and
   $\hat g$ = `components(exterior_derivative(pull_back(u, B_+), K), K)` the
   reference gradient, and then pulled back by $A_+$, which remains on
   $\hat g$ and on $J_+$, because `I` has lower dimension than `K`: the
   reference gradient and the cell Jacobian are read at the points of the
   facet. On `n` the rule for `Normal` gives $K^T \hat m / \lvert K^T \hat m \rvert$,
   with $\hat m$ the normal of the edge `(1, 2)` of `K` pointing away from the
   centroid of `K`. On `alpha[e]`, an expression with no support whose free
   index `e` is free in `cell_domain(E, 0)`, the pull-back is `alpha[e]`.
3. `apply_derivatives` leaves `exterior_derivative(pull_back(u, B_+), K)` on
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
numbers replaced by the positions of the facet vertices, which also fix the orientation of the
facet from each side.

- In the example, $A_+(t) = (1 - t, t)$ and $A_-(t) = (0, t)$ (see
  [Restriction](LANGUAGE.md#restriction)), $J_+$ is the identity and $J_-$ is
  $\begin{pmatrix} 0 & -1 \\ 1 & 1 \end{pmatrix}$, and $\hat m = (1, 1) / \sqrt 2$,
  so `n` evaluates to $(1, 1) / \sqrt 2$ on the diagonal, pointing from
  triangle 0 into triangle 1, and `rho` to $\sqrt 2$.
- `restriction(u, E, "-")` on the same `t` reads `u` on triangle 1 at
  $(0, t)$, the same physical point $(1 - t, t)$ as from `"+"`, so `jump_u`
  is a difference at one point.

### Example: a fluid and a solid

Requires: [Example: an interior facet
integral](#example-an-interior-facet-integral), [Sides of a
topology](LANGUAGE.md#sides-of-a-topology), [Restriction](LANGUAGE.md#restriction),
[Form](LANGUAGE.md#form), [Concept: kernels](LANGUAGE.md#concept-kernels).

A fluid and a solid that share an interface, with the fields of each part in
spaces on that part only: the legacy UFL `dx(1)`, `dx(2)` and `ds(3)` with
markers, and the submeshes and entity maps of a solver such as DOLFINx. The
square of the [DecomposedDomain](LANGUAGE.md#decomposeddomain) example is the
whole, triangle `(0, 1, 2)` is the fluid, triangle `(1, 3, 2)` is the solid,
and the diagonal is the interface, the smallest mesh with every part of the
problem. The example uses `D`, `square`, `K` and
`I` of the [interior facet example](#example-an-interior-facet-integral), and
`runtime_array` of [Example: integrating over a cut
cell](#example-integrating-over-a-cut-cell).

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
in the [interior facet example](#example-an-interior-facet-integral).

The problem is Stokes flow in the fluid and linear elasticity in the solid,
coupled by a multiplier on the interface that equates the fluid velocity to
the velocity of the solid over a time step. With spaces of a library, `V` and
`Q` on `F`, `W` on `S` and `L` on `G`, all of the degree `()`, with `V`, `W`
and `L` of the value shape `(2,)` and `Q` of the value shape `()`:

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
grad_u = components(exterior_derivative(u, F), F)
grad_v = components(exterior_derivative(v, F), F)
div_u, div_v = index_sum(grad_u[i, i], (i,)), index_sum(grad_v[i, i], (i,))
fluid = integral(
    nu * index_sum(grad_u[i, j] * grad_v[i, j], (i, j)) - p * div_v - q * div_u, F
)

# Linear elasticity on S.
grad_d = components(exterior_derivative(d, S), S)
grad_w = components(exterior_derivative(w, S), S)
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
  of functions on three domains (see [Form](LANGUAGE.md#form)).
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
  kernel after the renaming of [Concept: kernels](LANGUAGE.md#concept-kernels),
  because the two topologies become the same placeholder.
- A material coefficient that differs between the fluid and the solid needs no
  submeshes when the fields live on `D`. With `c = cell_index(square, 0)`, an
  array `nu[c]` under `integral(..., D)` takes the value of the cell (see
  [Integral over a decomposed domain](LANGUAGE.md#integral-over-a-decomposed-domain)).
- A solver that builds the fluid submesh from the coordinates of the submesh
  gives a decomposed domain with the points of `F` but not equal to `F`, and a
  space on one does not combine with a space on the other. The restriction
  rules hold for either, because the inclusion of a cell into itself is the
  identity and only the side tables enter. A library picks one construction.
- Dirichlet conditions, such as a clamped solid and a no-slip wall, belong to
  the solver.
- A fluid domain that moves with the solid is not covered. The moved domain is
  a decomposed domain with the topology `fluid_cells` and new coordinates, so a
  function in a space on `F` has a different support from an expression on the
  moved domain, and the two do not combine. A map between two decomposed
  domains with the same topology, the identity in the reference coordinates of
  each cell, is a fixed composition that a later version may add.

The core transformations rewrite `a` as in the [interior facet
example](#example-an-interior-facet-integral). `expand_integrals` gives one
`IndexSum` for each of `F`, `S` and `G`, over the cell index of its topology.
On `G`, with `e = cell_index(interface, 0)`, the restriction `u_G` lowers to
the pull-back of `u` by the fluid cell `cells_fluid[e]` and then by the
inclusion of the facet into that cell, and `d_G` to the pull-back of `d` by
the solid cell `cells_solid[e]`. The kernel of `coupling` is an expression on
`I` that reads the reference basis of `V` on one fluid cell, of `W` on one
solid cell and of `L` on the facet, and the assembler packs the two cells
through the two side tables, the entity maps of the two submeshes.

### Example: a 1D network in a 3D mesh

Requires: [DecomposedDomain](LANGUAGE.md#decomposeddomain),
[Sides of a topology](LANGUAGE.md#sides-of-a-topology),
[Restriction](LANGUAGE.md#restriction), [Example: a fluid and a
solid](#example-a-fluid-and-a-solid), [Points and point
evaluation](LANGUAGE.md#points-and-point-evaluation), [Form](LANGUAGE.md#form).

A network of vessels, a mesh of segments in $R^3$, coupled to the tissue
around it, a mesh of tetrahedra, with the pressure of each in a space on its
own mesh. The segments pass through the tetrahedra and are not sub-entities
of them. The mesh library cuts each segment where it crosses a face, so that
each piece lies in one tetrahedron, and the pieces are a topology with two
sides, into the network and into the tissue, as the interface of [Example: a
fluid and a solid](#example-a-fluid-and-a-solid) has a side into each part.
The positions in the side tables are not vertices: from the network, a piece
is a sub-interval of the reference interval, and from the tissue, a segment
between two points of the reference tetrahedron.

With `K3 = ReferenceSimplex(3)` and `I = ReferenceSimplex(1)`, the
coordinates `x3` of the tetrahedra with the free index `cell_index(T3, 0)`
and `x1` of the segments with the free index `cell_index(T1, 0)`, both built
as in the [DecomposedDomain](LANGUAGE.md#decomposeddomain) example, a
topology `pieces` of the library, spaces of a library `V3` on `D3` and `V1`
on `D1` of the degree `()` and the value shape `()`, and `runtime_extent` and
`runtime_array` of [Example: integrating over a cut
cell](#example-integrating-over-a-cut-cell):

```python
# The tissue, of dimension 3, and the network, of dimension 1, both in RD(3).
D3 = decomposed_domain((domain_map(K3, x3, RD(3)),))
D1 = decomposed_domain((domain_map(I, x1, RD(3)),))

# The pieces of the segments, each inside one tetrahedron, with the cell type I
# and the sides (("vessel", T1, 0), ("tissue", T3, 0)), and their domain.
P = side_domain(D1, pieces)

p3, q3 = function(V3, "p3"), function(V3, "q3")
p1, q1 = function(V1, "p1"), function(V1, "q1")
beta = runtime_array("beta", ())
m = index("m")

# The Laplacian of each part on its own mesh. On D1 the gradient is the
# derivative along the segment, in the frame of RD(3).
grad_p3 = components(exterior_derivative(p3, D3), D3)
grad_q3 = components(exterior_derivative(q3, D3), D3)
grad_p1 = components(exterior_derivative(p1, D1), D1)
grad_q1 = components(exterior_derivative(q1, D1), D1)
tissue = integral(index_sum(grad_p3[m] * grad_q3[m], (m,)), D3)
vessel = integral(index_sum(grad_p1[m] * grad_q1[m], (m,)), D1)

# The two pressures on the pieces, each from its own mesh.
p1_P, q1_P = restriction(p1, P, "vessel"), restriction(q1, P, "vessel")
p3_P, q3_P = restriction(p3, P, "tissue"), restriction(q3, P, "tissue")

# The exchange between the vessel and the tissue around it.
coupling = integral(beta * (p1_P - p3_P) * (q1_P - q3_P), P)
a = form(tissue + vessel + coupling, ((q3, q1), (p3, p1)))
```

- `pieces` has the cell type `I` and the sides
  `(("vessel", T1, 0), ("tissue", T3, 0))`. For a piece `e`, the side table
  from `"vessel"` holds the segment `cells_v[e]` and the positions
  `positions_v[e, v, 0]`, the two parameters in `I` of the ends of the piece
  along the segment, and the side table from `"tissue"` holds the tetrahedron
  `cells_t[e]` and the positions `positions_t[e, v, k]`, the two ends in
  `K3`. The library computes the pieces and both tables, outside the
  language, once for the pair of meshes.
- `restriction(p3, P, "tissue")` is `pull_back(pull_back(p3, B_t), A_t)` with
  $B_t$ the tetrahedron of each piece and $A_t$ =
  `domain_map(I, affine_map(I, points_t), K3)` the segment between the two
  ends, by [Restriction](LANGUAGE.md#restriction). `restriction(p1, P, "vessel")`
  is the same with the sub-interval of `I`. Both are expressions on `P` at
  the same points of each piece, so `p1_P - p3_P` is a difference at one
  point.
- `side_domain(D1, pieces)` takes its geometry from the first side, the
  segment, so the first argument is `D1`, whose topology is `T1`, and the
  integral over `P` uses the length element along the segment. The pieces
  cover the segments, so `integral(f, P)` and `integral(f, D1)` agree for an
  `f` on `D1`, and only the coupling, which needs `p3`, is written over `P`.
- `p3` and `p1` have the supports `D3` and `D1` and never combine directly.
  Every integral has no support, so the three terms add, and the arguments of
  `a` are tuples of functions on the two meshes.
- In a curved tetrahedron a piece is a curve in `K3`, and the straight segment
  between its ends is an approximation, as a quadrature rule is.

The core transformations rewrite `a` as in the [interior facet
example](#example-an-interior-facet-integral). `expand_integrals` gives one
`IndexSum` for `tissue`, over `cell_index(T3, 0)`, one for `vessel`, over
`cell_index(T1, 0)`, and one for `coupling`, over `cell_index(pieces, 0)`.
On `coupling`, `expand_pull_backs` lowers the four restrictions to pull-backs
by one segment and one tetrahedron, selected through the two side tables, and
then by the affine maps through the positions. The kernel of `coupling` for a
piece reads the reference basis of `V1` on one segment and of `V3` on one
tetrahedron, at the points of the piece in the reference coordinates of each:
the interface of a two-sided facet kernel, with the cells of two meshes and
of two dimensions. The local matrix couples the degrees of freedom of `V1` on
the segment with those of `V3` on the tetrahedron, and the assembler scatters
the matrix through the maps of degrees of freedom of both.

For a network whose segments are edges of the tetrahedra, the pieces are the
segments themselves, so the library gives `T1` the tissue side directly,
`(("tissue", T3, 0),)`, with the tetrahedron of each segment chosen by the
library among those that share the edge, and the positions the coordinates of
two vertices of `K3`. `D1` is then `side_domain(D3, T1)`, `p1` needs no
restriction, and the coupling is an integral over `D1`:

```python
D1 = side_domain(D3, T1)
p3_1, q3_1 = restriction(p3, D1, "tissue"), restriction(q3, D1, "tissue")
coupling = integral(beta * (p1 - p3_1) * (q1 - q3_1), D1)
```

- The side tables come from the connectivity of the mesh, with no cutting,
  and a segment and its tetrahedron share vertices, so the partition of `T3`
  gives the partition of `T1`.
- `p3` is continuous, so the restriction is the same from every tetrahedron
  that shares the edge, and the choice of the library does not change the
  form.
- The positions take one of twelve values, an edge of `K3` in one of two
  orientations. A library may write the table as a gather from the twelve by
  a `RuntimeIndex` of [Example: integrating over a cut
  cell](#example-integrating-over-a-cut-cell), as `StaticTopology` computes
  the positions from vertex numbers, and a form compiler may then tabulate
  the basis of `V3` on each edge of `K3` once and select by the index, as a
  legacy facet kernel selects by the local facet number. For the pieces of a
  non-conforming network the positions are arbitrary points of `K3`, and the
  kernel evaluates the basis of `V3` at run time, as for a cut cell.

The trace of a function in $H^1(\Omega)$ on a curve is not defined, so a
model may replace `p3` on the network by its average over a ring of radius
$R$ around each point of the network, in the plane of two unit normals `n1`
and `n2` of the segment, functions in a space `N` on `D1` of the value shape
`(3,)` that the solver provides. The ring points are not pieces of the
network, so the average is a family of maps from `D1` into `D3` and point
evaluation, with the cosines and sines of the angles and the weights as
tables, because the core has no `cos`:

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
  `p3_avg` is an expression on `D1`, and `coupling_avg` needs no pieces. With
  `R = 0`, `pull_back(p3, ring)` is the trace of `p3` on the network,
  `pull_back(p3, domain_map(D1, X1, D3))`.
- The pull-back of `p3` by `ring` remains after `expand_pull_backs`, and the
  pass of [Points and point evaluation](LANGUAGE.md#points-and-point-evaluation)
  lowers it after a quadrature pass on `I` as in [Example:
  quadrature](#example-quadrature). The maps below `p3` are then `ring`,
  `cell_domain(D1, 0)` and the quadrature points, a chain with the free
  indices of the segment, the quadrature point and the ring point, which the
  pass reads as one family of points. The side tables of that family hold
  the tetrahedron that contains each point and its position, and the kernel
  for a segment reads `V3` on one tetrahedron for each point, a list of
  cells rather than one cell, which [Concept:
  kernels](LANGUAGE.md#concept-kernels) allows because the interface of a
  kernel follows from the expression.
- A ring point outside the tissue, near its boundary, has no containing cell,
  and what the table holds for the point is the choice of the solver.

### Example: the Laplace-Beltrami operator

Requires: [A tour: the Poisson problem](LANGUAGE.md#a-tour-the-poisson-problem),
[Components](LANGUAGE.md#components), [Integral over a decomposed
domain](LANGUAGE.md#integral-over-a-decomposed-domain), [Form](LANGUAGE.md#form).

See also: [Expanding pull-backs](LANGUAGE.md#expanding-pull-backs), [Inverse](LANGUAGE.md#inverse).

The Poisson forms of the [tour](LANGUAGE.md#a-tour-the-poisson-problem) on a surface in
$R^3$ are the forms of the Laplace–Beltrami operator. Only the geometry
changes. The vertex positions `y` have value shape `(m, 3)`, and the ambient
coordinate domain is `RD(3)`:

```python
# A triangulated surface in RD(3), with K, T, c and p of the tour.
x = component_tensor(
    p[c, 0, k] + (p[c, 1, k] - p[c, 0, k]) * X[0] + (p[c, 2, k] - p[c, 0, k]) * X[1],
    (k,),
)
M = decomposed_domain((domain_map(K, x, RD(3)),))

# Functions in a space of a library on M.
u, v, f = function(V, "u"), function(V, "v"), function(V, "f")

# The surface gradients, of value shape (3,) in the frame of RD(3).
i = index("i")
grad_u = components(exterior_derivative(u, M), M)
grad_v = components(exterior_derivative(v, M), M)

a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), M), (v, u))
L = form(integral(f * v, M), (v,))
```

- `M` has the dimension 2 and the ambient coordinate domain `RD(3)`, so
  [`Components`](LANGUAGE.md#components) projects the frame of `RD(3)` onto the tangent
  plane of each cell, and `grad_u` is the surface gradient of `u`.
- `integral(..., M)` multiplies the integrand on `K` by the density
  $\rho = \sqrt{\det(J^T J)}$, with $J$ the Jacobian matrix of `x`, of value
  shape `(3, 2)`. The integrands are 0-forms, so the forms need no orientation
  of `M`.
- `expand_pull_backs` rewrites the pull-back of `grad_u` by the cells as
  $K^T \hat g$, with $K$ = `inverse(J)` the pseudo-inverse, of value shape
  `(2, 3)`, and $\hat g$ the components of the exterior derivative of the
  pull-back of `u` on `K`. The integrand of `a` on `K` is then
  $\hat g_u^T K K^T \hat g_v \, \rho$, and $K K^T = (J^T J)^{-1}$ is the
  inverse of the metric that `x` induces.
- With a function `w` in a space on `M` of the value shape `(3,)`,
  `index_sum(components(exterior_derivative(w, M), M)[i, i], (i,))` is the
  surface divergence of `w`.
- With `P = components(exterior_derivative(coordinate(M), M), M)`, the
  projection onto the tangent plane,
  `components(exterior_derivative(grad_u, M), M)` is the tangential
  derivative $H$ of the surface gradient, which is not symmetric in general.
  The covariant Hessian of `u` is $P H P$.

### Example: the flux through a surface

Requires: [Example: the Laplace-Beltrami
operator](#example-the-laplace-beltrami-operator), [Concept: orientation of a
domain](LANGUAGE.md#concept-orientation-of-a-domain), [Integral over a decomposed
domain](LANGUAGE.md#integral-over-a-decomposed-domain).

The flux of a vector field through a surface is the integral of a 2-form, so
the flux needs an orientation of the surface. With `M`, `T` and `c` of the
[Laplace-Beltrami example](#example-the-laplace-beltrami-operator), where `M` is
a closed surface such as a triangulated sphere:

```python
# The 2-form of the flux of the position field through a surface in RD(3).
Y = coordinate(RD(3))
dY = exterior_derivative(Y, RD(3))
omega = Y[0] * dY[1] * dY[2] + Y[1] * dY[2] * dY[0] + Y[2] * dY[0] * dY[1]

# The sign of each triangle from the mesh library, +1 where the order of the
# vertices of the triangle agrees with the outward orientation of M.
N = T.num_cells[0]
o = runtime_array("o", (N,))

# The outward flux, 3 times the volume inside M.
flux = integral(o[c] * omega, M)

# The area of M, with no orientation.
area = integral(Integer(1), M)
```

- `omega` has the support `RD(3)`, so `integral(omega, M)` pulls `omega` back
  by the triangles of `M` (see [Integral over a decomposed
  domain](LANGUAGE.md#integral-over-a-decomposed-domain)). Each triangle contributes the
  flux with the orientation of the parameterisation of the triangle, which the
  order of the vertices of the triangle in `T` fixes. If the mesh library
  numbers some triangles clockwise and others anticlockwise, seen from
  outside, the contributions of the two groups have opposite signs, and the
  sum is not the flux.
- `o[c]` has the free index `c` of the triangles, so each triangle takes the
  sign of that triangle. With the signs, every contribution uses the outward
  orientation. A mesh library may instead number the vertices of every
  triangle consistently, and then `integral(omega, M)` is the flux without
  `o`.
- `area` integrates a 0-form against the density $\sqrt{\det(J^T J)}$, which
  does not depend on the order of the vertices, so `area` needs no signs.
- A surface with no consistent orientation, such as a Möbius strip, has no
  signs `o`. The area of the surface is defined, and the flux is not.

`runtime_extent` and `runtime_array` are not in the core (see [Example:
integrating over a cut cell](#example-integrating-over-a-cut-cell)). The same
signs enter a term that turns a 2-form on `M` into a 0-form, such as the
scalar curl of a 1-form on `M`.

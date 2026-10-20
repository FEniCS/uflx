# Worked examples for the UFLx core

These examples accompany [LANGUAGE.md](LANGUAGE.md) and use its names. The
first part defines the library classes that the examples and the tour of
LANGUAGE.md use; none of them is part of the core. The second part exercises
the core on small problems with known answers, and then on three larger
problems: an interior facet integral, the Laplace-Beltrami operator and the
flux through a surface. A reader of LANGUAGE.md needs none of this to follow
the rules.

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
the side tables from the vertex numbers. For a cell seen from two cells of the
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
tables `(cells_pm, vertices_pm)` of each side:

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
   vertex positions inside $J_+$ and `points_pm` to `Indexed`s with gathered
   indices, nested through the cell on the side.

After the passes, `a_facet` is `form(index_sum(integral(g, I), (e,)), (v, u))`
with `g` an expression on `I` with the free index `e`. The `IndexSum` over `e`
is the loop of the assembler over the interior facets. The gathered indices in
`g` read, for the facet `e`, the two cells `cells_p[e]` and `cells_m[e]`, the
vertex positions of each through `cell_vertices`, the paired vertices
`vertices_p[e, v]` and `vertices_m[e, v]`, and `alpha[e]`, which the assembler
packs for a kernel. The kernel is `g` without `e`: an expression on the
reference facet `I` that reads the reference basis of `u` and `v` on the two
cells, through the library of `V`, at the points of the facet in each cell's
reference coordinates $A_\pm(t)$, and the two cell Jacobians at those points.
This is the interface of a legacy interior facet kernel, with the local facet
numbers replaced by the paired vertices, which also fix the orientation of the
facet from each side.

- In the example, $A_+(t) = (1 - t, t)$ and $A_-(t) = (0, t)$ (see
  [Restriction](LANGUAGE.md#restriction)), $J_+$ is the identity and $J_-$ is
  $\begin{pmatrix} 0 & -1 \\ 1 & 1 \end{pmatrix}$, and $\hat m = (1, 1) / \sqrt 2$,
  so `n` evaluates to $(1, 1) / \sqrt 2$ on the diagonal, pointing from
  triangle 0 into triangle 1, and `rho` to $\sqrt 2$.
- `restriction(u, E, "-")` on the same `t` reads `u` on triangle 1 at
  $(0, t)$, the same physical point $(1 - t, t)$ as from `"+"`, so `jump_u`
  is a difference at one point.

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

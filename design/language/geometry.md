# Geometry

## Maps between domains

### `DomainMap`

`DomainMap` is a map $\psi$ from a parameter domain $P$ into an ambient coordinate domain
$Q$, with geometric image $\psi(P)$. The sections below write $D$ for the map
or its parameterised image, retaining multiplicity. A `DomainMap` is not an
`AbstractDomain` and is not the domain of an expression.

| `DomainMap(P, psi, Q)` |                                                                                                                                                                                                             |
| -------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                          | `Base` |
| Free function                    | `domain_map(P, psi, Q)`                                                                                                                                                                           |
| Takes                            | `P: AbstractCoordinateDomain \| AbstractTopology`, `psi: AbstractExpression`, `Q: AbstractCoordinateDomain` |
| Properties                       | `parameter_domain` $P$, an [`AbstractCoordinateDomain`](rules.md#abstractcoordinatedomain) or an [`AbstractTopology`](#abstracttopology) of dimension $k$, `parameterisation` $\psi$, and `ambient_domain` $Q$, an `AbstractCoordinateDomain`                                                   |
| Requires                         | $\psi$ with the domain $P$ or no domain, $\psi$ real with value shape `(n,)` for the dimension $n$ of the ambient coordinate domain of $Q$ and the argument spaces `()`, and $k \le n$ |
| Derived                          | `free_indices`, the free indices of $\psi$ |
| Raises                           | `TypeError` for $P$ or $Q$, and `ValueError` if another requirement fails                                                                                                                                   |

The coordinates of $Q$ locate the image of
$D$. The 'topological dimension' of $D$ is $k$,
`D.parameter_domain.dimension`, and the 'geometric dimension' is $n$, the
dimension of the ambient coordinate domain of $Q$.

Maps need not be injective.

$P$ is a coordinate domain, such as a reference cell, or a topology. The
target $Q$ is a coordinate domain, never a topology or a mesh. Coordinates
alone do not select a cell or a reference preimage in a mesh; evaluation of
mesh fields uses explicit cell or side maps (see [A family of points in a
mesh](#a-family-of-points-in-a-mesh)).

A $\psi$ with free indices makes $D$ a family of domains, with one domain for
each value of the free indices. For example, with an array `x` of value shape
`(m, n)` and an index `i`, `domain_map(RD(0), x[i], RD(n))` is a family of $m$
points in `RD(n)`, retaining distinct members even at identical coordinates.

$Q$ is usually `RD(n)`. $Q$ may be smaller, such as the unit triangle for an
edge of the unit triangle.

Examples:

- The hypotenuse of the unit triangle is the image of a domain map with $k = 1$
  and $n = 2$. $P$ is the unit interval in $R^1$, and
  $\psi(X) = (1 - X_0, X_0)$. $Q$ is `RD(2)` for the hypotenuse in the plane,
  and `UnitSimplex(2)` for the hypotenuse in the unit triangle.

A curved example is the paraboloid, with $k = 2$ and $n = 3$:

```python
P = RD(2)
X = coordinate(P)
h = X[0] ** Integer(2) + X[1] ** Integer(2)
psi = tensor((X[0], X[1], h))
paraboloid = domain_map(P, psi, RD(3))
```

### `PullBack`

`PullBack` is the pull-back $\psi^* f$ of an expression $f$ on $Q$ by a map $D$
with parameter domain $P$ and ambient domain $Q$. A `DomainMap` additionally
specifies its coordinate parameterisation $\psi$.

| `PullBack(f, D)` |                                                                                                            |
| ---------------- | ---------------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `pull_back(f, D)`                                                                                          |
| Takes            | `f: AbstractExpression`, `D: DomainMap \| SideMap \| CellMap` |
| Properties       | `expression` $f$, and `map` $D$, a `DomainMap`, a [`SideMap`](#sidemap) or a [`CellMap`](#cellmap) |
| Requires         | $f$ with the domain $Q$ or no domain, and an equal extent for an `Index` free in both $f$ and $D$          |
| Static type      | The value shape of $f$, and the free indices of $f$ followed by those of $D$, with each `Index` taken once |
| Raises           | `TypeError` for $D$, and `ValueError` if another requirement fails                                         |

The second operand is a map that identifies its source $P$ and target $Q$.

For a `DomainMap`, the pull-back is the composition $f \circ \psi$, the
expression on $P$ with the value $f(\psi(p))$ at each point $p$ of $P$.
For a `CellMap` or `SideMap`, the map selects cells and reference positions;
its physical coordinates alone do not define the evaluation. The axes of the
value shape are not transformed. An axis of an expression on $Q$ is an array in the
coordinate frame of $Q$, so the pull-back keeps the frame of $Q$.

[`PushForward`](#pushforward) gives the corresponding forward value
transformation from $T$ to the mesh $D$. With the identity law it reverses
the pull-back cellwise; other laws also transform the values. It does not
denote a dual or adjoint operator.

Pull-backs of terminals remain for their supplying libraries to evaluate
(see [Functions in a space](forms.md#functions-in-a-space)).

For coordinate maps $D$ and $E$, with $E$ mapping into the parameter domain
of $D$, their composition is the domain map with the parameterisation
`pull_back(D.parameterisation, E)` on the parameter domain of $E$, in the
ambient domain of $D$, such as `E` in the example below. Composition through
a cell map instead keeps nested pull-backs, preserving the cell selection.

This example maps an edge of the unit triangle `T` into a triangle `K` in
`RD(2)`:

```python
T = UnitSimplex(2)
X = coordinate(T)
one, two = Integer(1), Integer(2)
psi = tensor((one + two * X[0] + X[1], one + two * X[1]))

# The unit interval, the parameter domain of the edge, with coordinate t.
S = UnitSimplex(1)
t = coordinate(S)[0]

# The edge of T from (1, 0) to (0, 1), as a domain in T.
e = domain_map(S, tensor((one - t, t)), T)

# PullBack(f, e), with the domain UnitSimplex(1) and value shape ().
f = X[0] * X[1]
fe = pull_back(f, e)

# The edge of K from (3, 1) to (2, 3), as a domain in RD(2).
E = domain_map(S, pull_back(psi, e), RD(2))
```

- A later pass evaluates `fe` to $(1 - t) t$, and the parameterisation of `E`
  to $(3 - t, 1 + 2 t)$.

Changing the ambient domain makes this pull-back invalid: `f` has domain
`T`, not `RD(2)`.

```python
edge_in_plane = domain_map(S, tensor((one - t, t)), RD(2))
pull_back(f, edge_in_plane)  # Raises ValueError.
```

### Jacobian and density

The free functions `jacobian(D)`, `inverse_jacobian(D)` and `density(D)` are
the Jacobian matrix, the left inverse of the Jacobian matrix and the density
of a `DomainMap` or a `CellMap` of a mesh, with coordinate parameter domain
$P$ of dimension $k$ and parameterisation $\psi$. $n$ is the physical
coordinate dimension of the target. A `CellMap` of a topology has no
parameterisation and is rejected with `TypeError`.

Each helper checks these input requirements before building its composition.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `jacobian(D)` | | `D: DomainMap \| CellMap` | `gradient(D.parameterisation, D.parameter_domain)` |
| `inverse_jacobian(D)` | | `D: DomainMap \| CellMap` | `inverse(jacobian(D))` |
| `density(D)` | | `D: DomainMap \| CellMap` | The density $\rho$ below |

The Jacobian matrix $J$ has value shape `(n, k)`, and $K$ =
`inverse_jacobian(D)`, of value shape `(k, n)`, is the inverse of $J$ for
$k = n$ and the pseudo-inverse $(J^T J)^{-1} J^T$ for $k < n$. The density
$\rho$ is the factor by which $\psi$ scales the length, the area or the volume
of a small region of $P$:

| Case | `density(D)` |
| ---- | ------------ |
| $k = 0$ | `Integer(1)` |
| $0 < k = n$ | `absolute(determinant(J))` |
| $0 < k < n$ | `square_root(determinant(G))`, with $G = J^T J$ |

The row $k = 0$ gives the value of the other rows, and keeps the Jacobian, of
value shape `(n, 0)`, out of an integral over points.

$J^T J$ is `component_tensor(index_sum(J[l, a] * J[l, b], (l,)), (a, b))`, with
derived indices.

An inverse Jacobian requires $\operatorname{rank} J = k$ wherever evaluated.
The density alone does not require full rank and may vanish.

The core, passes and libraries use `jacobian(D)`, `determinant(jacobian(D))`
and `inverse_jacobian(D)` for $J$, $\det J$ and $K$. These expressions let a
pass identify and evaluate each geometric quantity once per $D$.

### Concept: integral over a domain

Integration of an expression $f$ on $Q$ over $D$ = `DomainMap(P, psi, Q)`
is defined by the parameterisation:

$$
\int_D f = \int_P (f \circ \psi) \, \rho \, \mathrm{d}X,
$$

which is `integral(pull_back(f, D) * density(D), P)` for a scalar $f$. The
integral does not depend on the direction in which $\psi$ traverses $D$.

Coverage is counted with multiplicity. The result equals the ordinary integral
over the geometric image when the parameterisation is injective almost everywhere.

The core has no orientation of a domain. An integral that depends on a
direction, such as a circulation or a flux, is written with the Jacobian
matrix, which carries the direction of $\psi$. For example, the circulation
of a vector field $F$ along an edge $e$ = `DomainMap(S, psi, Q)` with $k = 1$
is

$$
\int_e F \cdot \mathrm{d}r = \int_S (F \circ \psi) \cdot \psi' \, \mathrm{d}t,
$$

which is `integral(index_sum(pull_back(F, e)[i] * jacobian(e)[i, 0], (i,)), S)`,
and the circulation changes sign when $\psi$ runs the other way. Where the
cells of a mesh must be oriented consistently, such as for the
flux through a surface in $R^3$, a separately supplied scalar field gives
a sign $+1$ or $-1$ for each cell (see
[The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface)).

### Mapped integration families

An indexed `DomainMap` describes a family of maps and retains its indices.
For a scalar field $f$ on $Q$ and `M = domain_map(P, psi[i, k], Q)`, with a
coordinate parameter domain $P$:

```python
per_member = integral(pull_back(f, M) * density(M), P)
per_configuration = index_sum(per_member, (i,))
```

The first expression retains `i` and `k`. The second sums only over `i` and
retains `k`. For several unrelated maps, add their reference-domain integrals
explicitly. No field domain is created by this construction. Integration counts
each map member separately, including coincident members.

For example, three segment maps of lengths 2, 3 and 5 give an indexed integral
of `Integer(1)` with values 2, 3 and 5. An explicit sum over their member index
gives 10. A topology-backed mesh of those segments also has total length 10,
but additionally identifies a domain on which fields and spaces can be defined.

### Points and point evaluation

A point is the image of a domain map. A point in a domain $Q$ is
`domain_map(RD(0), x, Q)`, with $x$ a constant expression of value shape
`(n,)` for the dimension $n$ of the ambient coordinate domain of $Q$. The core
needs no node class and no free function for a point.

#### Point evaluation

The value of an expression `f` at a point or indexed family of points `p` is:

```python
integral(pull_back(f, p), RD(0))
```

`p` has parameter domain `RD(0)`, and `f` is defined at its points. The
result has no domain, preserves the value shape, scalar type and argument
spaces of `f`, and carries the free indices of `f` followed by those of `p`,
each taken once. Point indices are retained, not summed.

Evaluation preserves symbolic dependence on the expression and the points.
Derivatives and value transformations are applied to the expression before
sampling; sampling `gradient(f, Q)` gives the gradient on `Q` at the points.
Basis tabulation and point observations use the same evaluation. The sampled
values alone define no field between the points.

For points in a mesh, see [A family of points in a
mesh](#a-family-of-points-in-a-mesh).

#### Quadrature

A quadrature rule on a coordinate domain is a weighted sum of point
evaluations. For a family `p` with index `i` and weights `w`:

```python
index_sum(w[i] * integral(pull_back(f, p), RD(0)), (i,))
```

`integral(f, Q)` denotes the integral; this sum approximates it and is
mathematically equal when the rule is exact for `f`. A quadrature pass
applies a supplied rule. On a mapped domain, the integrand includes the density
from [integral expansion](#integral-over-a-mesh).

#### Examples

The examples use the unit triangle `T = UnitSimplex(2)`, `X = coordinate(T)`,
and `f = X[0] * X[1]`:

```python
third = Integer(1) / Integer(3)

# DomainMap: the centroid of T, as a point in T.
c = domain_map(RD(0), tensor((third, third)), T)

# The value at the centroid: 1/9, with no domain and value shape ().
fc = integral(pull_back(f, c), RD(0))
```

Worked examples for this chapter: [an affinely mapped triangle](../examples/integration.md#an-affinely-mapped-triangle), [circulation along an edge](../examples/integration.md#circulation-along-an-edge), [integrating over an edge](../examples/integration.md#integrating-over-an-edge), [integrating over the boundary of a triangle](../examples/integration.md#integrating-over-the-boundary-of-a-triangle), [quadrature](../examples/integration.md#quadrature), [integrating over a cut cell](../examples/integration.md#integrating-over-a-cut-cell), [the centroid of a triangle](../examples/integration.md#the-centroid-of-a-triangle).

## Topology

### `AbstractCellType`

`AbstractCellType` describes a cell type with numbered vertices.

| `AbstractCellType` |                                                                         |
| ------------------ | ----------------------------------------------------------------------- |
| Extends            | `Base` |
| Provides           | `dimension`, an `int` $n \ge 0$, and `num_vertices`, an `int` $m \ge 1$ |
| Defined by         | A subclass outside the core, such as `Simplex`                          |

The core does not number the facets or the other sub-entities of a cell type.
A sub-entity of a cell is named by the vertex numbers of the cell that the
sub-entity contains, such as the vertex numbers of a facet, and the mesh or a
library such as Basix decides which sets of vertices are sub-entities (see
[Sides of a topology](#sides-of-a-topology)).

A topology of polygons, as for the virtual element method, has a cell type for
each number of vertices. A number of vertices known only at run time is not yet
covered.

### `AbstractTopology`

`AbstractTopology` describes a set of cells of given cell types, without
coordinates.

| `AbstractTopology` |                                                                                             |
| ------------------ | ------------------------------------------------------------------------------------------- |
| Extends            | `Base` |
| Provides           | `cell_types`, `num_cells`, `sides`, and `side_tables(i, side)` for a topology with sides |
| Derived            | `dimension`, the dimension of the cell types                                                |
| Defined by         | A subclass outside the core, such as the topology of the mesh of a solver                   |

- `cell_types` is a nonempty `tuple` of `AbstractCellType`s $C_0, C_1, \dots$ of the
  same dimension, one for each group of cells. A cell of group $j$ has type
  $C_j$ and its vertices are numbered as in $C_j$.
  The cell types are distinct, except that two groups with different sides may
  have the same cell type (see [Sides of a topology](#sides-of-a-topology)).
- `num_cells` is a `tuple` of the same length as `cell_types`, with the number
  of cells in each group: a nonnegative `int`, such as for a macro element,
  or an `AbstractExtent`, such as for the mesh of a solver.
- `sides` is `()`, or the sides of [Sides of a topology](#sides-of-a-topology).

Zero-cell groups are valid, including when every group is empty. An empty
topology retains its cell types, which determine its dimension. An empty
boundary selection likewise retains its facet types and side declarations;
its side tables have zero entries along the cell axis. The supplying library
must satisfy these contracts, including for counts supplied at evaluation.

A cell-group number `j` for a topology `T` is an `int` satisfying
`0 <= j < len(T.cell_types)`. This applies to cell maps, side triples and
`side_tables`; negative numbers do not select groups from the end.

A topology whose cell types are reference cells may be the domain of an
expression and the domain of a space, as the cell complex of its cells (see
[`CellMap`](#cellmap)). A topology is not an `AbstractDomain`, so
`Coordinate`, `Gradient` and `Integral` reject a topology.

Topology identity is structural: a concrete `AbstractTopology` stores its
identifying arguments, such as a label.

The free function `cell_index(T, j)` returns the derived index
`Index(("uflx.cell_index", T, j, 0))` for the cells of group $j$, numbered from
0. The index differs between groups and between topologies, even if the numbers of cells are equal.

`cell_index` checks that $T$ is a topology and $j$ is an `int` with
`0 <= j < len(T.cell_types)`.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_index(T, j)` | | `T: AbstractTopology`, `j: int` | `Index` |

An expression with the free index `cell_index(T, j)` has one value per cell of
group $j$.

Examples:

- The Alfeld split of a triangle, a macro element, has the numbers of cells
  `(3,)`.
- A mesh of a solver with triangles and quadrilaterals has two cell types, and
  two `AbstractExtent`s as the numbers of cells.

```{dropdown} Legacy UFL
In legacy UFL, a `Mesh` holds both the topology and the coordinate element, and the measure `dx` hides the sum over the cells. In UFLx, the topology, the parameterisations of the pieces of a mesh and the `IndexSum` over `cell_index(T, j)` are separate.
```

### Sides of a topology

A topology $S$ whose cells are sub-entities of cells of other topologies has
sides. A cell of $S$ is seen from one cell of another topology on each side,
and [Restriction](#restriction) picks the value of an expression on a cell of
$S$ by the side name.

#### `sides`

`S.sides` describes the sides of $S$:

- `S.sides` has one entry for each cell group of $S$.
- Entry $i$ is a `tuple` of one or more sides with distinct names.
- Each side is a triple `(name, T, j)` of a side name, a `str`, a topology
  $T$ and a cell-group number $j$ of $T$.
- Entries for groups with the same cell type must differ. A group of cells of $S$ is
  therefore fixed by the cell type and the sides together.
- A topology whose cells are not seen from other cells has the sides `()`.

Notation, for the group $i$ of $S$ with a side `(name, T, j)`:

- $K_j$ is the reference cell of group $j$ of $T$, and $k_j$ is the
  dimension of $K_j$.
- $m$ is the number of vertices of the cell type of group $i$ of $S$.

#### The image of a cell on a side

On each side, a cell of group $i$ of $S$ is the image of the cell type $F_i$
of group $i$ under the inclusion of the side tables, a map into $K_j$ (see
[Restriction](#restriction)). The image is one of the following:

- A sub-entity of the cell on the side. The sub-entity is the cell itself, with
  the vertices in the order of the cell, a facet, or a sub-entity of lower
  dimension, such as an edge of a tetrahedron or a vertex.
- A piece of the cell on the side that is not a sub-entity, such as half of a
  facet at a hanging node or a segment of a line through the cell.

The concrete mesh gives the inclusion of each cell of $S$ into the cell on
each side. A mesh library may build the inclusion as the affine map through
the positions of the vertices of the cell in $K_j$ (see [A static
topology](../examples/library.md#a-static-topology)):

- For a matching facet, the position is that of the same vertex.
- For a periodic boundary, the position is that of the vertex one period away.
- For a piece, the position is that of a point of the cell.

#### Side tables

A concrete topology with sides also provides `side_tables(i, side)`, the
triple `(cells, inclusion, reference_normal)` for the side named `side` of
group $i$, with the side `(side, T, j)`, the cell type $F_i$ of group $i$ and
`e = cell_index(S, i)`:

- `cells` is an [`AbstractIndexTable`](expressions.md#abstractindextable) of the shape
  `(S.num_cells[i],)` and the extent `T.num_cells[j]`. Entry $e$ selects a cell of
  group $j$ of $T$ on the side of the cell $e$ of $S$.
- `inclusion` is a [`DomainMap`](#domainmap) from $F_i$ into $K_j$ with the
  free index `e`. Member $e$ is the inclusion of the cell $e$ into the cell on
  the side.
- `reference_normal` is an expression with the domain $F_i$ or no domain,
  the value shape `(k_j,)` and the free index `e`, for a group of cells of
  dimension $k_j - 1$, and `None` for any other group. A normal that varies
  along the cell, such as on a curved piece, has the domain $F_i$. Member $e$ is a normal in $K_j$ of
  the image of the cell $e$. The normal points out of $K_j$ for an image on the
  boundary of $K_j$. For an image inside $K_j$, the supplying library chooses
  its orientation.

Side inclusions pair points across sides. Physical coincidence is a
library-supplied assumption for ordinary interfaces, not a core requirement;
periodic and mortar couplings may pair distinct physical points.
[Restriction](#restriction) and [`Normal`](#normal) read the tables.

#### What the concrete mesh decides

The concrete mesh decides the following:

- Which sub-entities or pieces are the cells of $S$. For example, the cells of
  $S$ may be the facets of a marked boundary or the cells of a marked
  subdomain.
- From which cells on each side each cell of $S$ is seen. For example, a facet
  of two cells may be seen from one cell only.
- For a sub-entity of lower dimension than a facet, such as a vertex of a
  triangle mesh, which one of the several cells that contain the sub-entity is
  on a side. A function that is continuous across the cells has the same
  restriction from each of the cells. The restriction of any other function
  depends on the choice.
- The names of the sides, and which cell is on which side. For example, the
  facets of an interface between fluid and solid triangles of the one cell
  type 0 of `T` have the sides `(("fluid", T, 0), ("solid", T, 0))`. The cell
  on the side `"fluid"` is the fluid triangle.

The sides may have different topologies, such as a mesh and a submesh, or two
meshes that meet. The topology of a side may itself have sides, so the cells
of $S$ may be facets of facets.

#### Cells that are pieces

A cell of $S$ that is a piece of the cell on a side covers the cases where one
entity meets several cells:

- At a hanging node, the facet of the small cell is half of a facet of the
  large cell. The large cell sees the half facet through a vertex and an edge
  midpoint of $K_j$.
- At a mortar interface between two meshes, the cells of $S$ are the pieces of
  the intersections of the facets. Each piece is seen from one cell of each
  mesh.
- For a line through the cells of a 3D mesh, the cells of $S$ are the segments
  of the line inside one cell each. Each segment is seen from the cell of the
  line and from the cell of the mesh (see [A 1D network in a 3D
  mesh](../examples/coupling.md#a-1d-network-in-a-3d-mesh)).

The concrete mesh computes which parts of an entity meet which cells, outside
the language, as the concrete mesh locates the cell that contains a point (see
[A family of points in a mesh](#a-family-of-points-in-a-mesh)). The inclusion of
a piece of a curved cell may be an approximation, as a quadrature rule is.

#### Examples

The examples use a triangle topology `T` with cell types `(K,)`,
`I = Simplex(1)` and `P = Simplex(0)`. For restriction, use reference cells
instead, such as `ReferenceSimplex(1)` for `I`.

- A boundary topology `B` has cell types `(I,)` and one side,
  `((("inside", T, 0),),)`. Its facet index is `cell_index(B, 0)`.
- Interior facets have two sides: `((("+", T, 0), ("-", T, 0)),)`.
- Corners of the boundary have cell types `(P,)` and sides
  `((("+", B, 0), ("-", B, 0)),)`: a side topology can itself have sides.
- For a boundary submesh `M` with cell types `(I,)`, the sides
  `((("mesh", T, 0), ("submesh", M, 0)),)` refer to different topologies.
- A periodic left boundary has sides `((("left", T, 0), ("right", T, 0)),)`.
  From `"right"`, each vertex is located one period to the right.
- For a topology `M` of triangles and quadrilaterals, interior facets have
  cell types `(I, I, I)` and sides
  `((("+", M, 0), ("-", M, 0)), (("+", M, 0), ("-", M, 1)), (("+", M, 1), ("-", M, 1)))`.
  The same facet type occurs in three groups, distinguished by the adjacent
  cell types. An integral has one term per group.
- For a prism topology `U3` with cell types `(W,)`, the boundary has two
  different facet types `(K, Q)`, triangular and quadrilateral, with sides
  `((("inside", U3, 0),), (("inside", U3, 0),))`.

```{dropdown} Legacy UFL
The legacy UFL measures `ds` and `dS` integrate over a topology of facets with one and two sides, and legacy UFL names the two sides `'+'` and `'-'`. The legacy UFL measure `dx(1)` integrates over a topology of the cells of a subdomain, and the legacy UFL measure `dP` integrates over a topology of vertices.
```

Worked examples for this chapter: [the simplex as a cell type](../examples/library.md#the-simplex-as-a-cell-type), [a static topology](../examples/library.md#a-static-topology).

## Reference cells

### Reference cell

An abstract subclass of `AbstractCoordinateDomain` and of `AbstractCellType`
for a reference cell `AbstractReferenceCell`, a cell type that is also a
coordinate domain. The core reads no vertices of a reference cell. A library
that builds inclusions, such as from vertex positions, keeps the vertices of
its reference cells (see [The unit simplex as a reference
cell](../examples/library.md#the-unit-simplex-as-a-reference-cell)).

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `AbstractCoordinateDomain`, `AbstractCellType` |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](../examples/library.md#the-unit-simplex-as-a-reference-cell)                                                   |

### `CellMap`

`CellMap(X, j, e)` selects cells of group $j$ of a topology or mesh $X$.
Let $T$ be $X$ when $X$ is a topology, and `X.topology` when $X$ is a mesh.
The selected cells retain their identities even when their physical
coordinates coincide.

| `CellMap(X, j, e)` | |
| ------------------ | - |
| Extends | `Base` |
| Takes | `X: AbstractTopology \| Mesh`, `j: int`, `e: int \| Index \| GatheredIndex` |
| Properties | `ambient_domain` $X$, `cell_group` $j$, and `entry` $e$ |
| Requires | $j$ a valid cell-group number of $T$, `T.cell_types[j]` a reference cell $K_j$; an integer $e$ satisfies `0 <= e < T.num_cells[j]`, and a `GatheredIndex` has extent `T.num_cells[j]` |
| Derived | `topology` $T$, `parameter_domain` $K_j$, `free_indices` and, for a mesh, `parameterisation`, as below |
| Raises | `TypeError` for an invalid argument type, and `ValueError` if another requirement fails |

For an ordinary `Index` $e$, `free_indices` contains $e$ with extent
`T.num_cells[j]`; the index itself stores no extent. For a `GatheredIndex`,
it is the free indices of $e$. For an integer entry it is `()`.

A cell map of a topology has no parameterisation. For a mesh $D$ with
coordinate expression $x$, `CellMap(D, j, e).parameterisation` is
`pull_back(x, cell_map(T, j, e))`. This expression is derived, not a
construction argument. Substitution reconstructs the map from $D$, $j$ and
$e$, preserving the selected cell independently of its coordinates.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_map(X, j, e)` | | `X: AbstractTopology \| Mesh`, `j: int`, `e: int \| Index \| GatheredIndex` | `CellMap(X, j, e)` |
| `cell_map(X, j)` | | `X: AbstractTopology \| Mesh`, `j: int` | `CellMap(X, j, cell_index(T, j))`, with $T$ the topology of $X$ |

`pull_back(f, cell_map(X, j, e))` is the value of an expression $f$ on $X$
on the reference cell $K_j$ of each selected cell. For a function on a
topology, the library of its space gives the value (see [Functions in a
space](forms.md#functions-in-a-space)). For a mesh field, lowering uses the selected
cell and its parameterisation (see [Expanding pull-backs](evaluation.md#expanding-pull-backs)).

Worked examples for this chapter: [the unit simplex as a reference cell](../examples/library.md#the-unit-simplex-as-a-reference-cell).

## Meshes

### `Mesh`

`Mesh(x, Q)` is one physical field domain, with coordinate expression $x$ on
an `AbstractTopology` $T$ and ambient domain $Q$. Cells keep their identities
when their physical coordinates coincide. The mesh stores no connectivity or
coordinate values; those remain with the topology and coordinate-function library.

| `Mesh(x, Q)` | |
| ------------ | - |
| Extends | `AbstractDomain` |
| Free function | `image(x, Q)` |
| Takes | `x: AbstractExpression`, `Q: AbstractCoordinateDomain` |
| Properties | `coordinate_function` $x$, and `ambient_domain` $Q$ |
| Requires | $x$ with domain $T$, a topology whose cell types are reference cells; $x$ real, with value shape `(n,)` for the ambient coordinate dimension of $Q$, free indices `()` and argument spaces `()`; `T.dimension <= n` |
| Derived | `topology` $T$ = `x.domain`, `dimension` = `T.dimension`, and `pieces`, below |
| Raises | `TypeError` for an invalid coordinate domain or ambient-domain type, and `ValueError` if another requirement fails |

`image(x, Q)` returns `Mesh(x, Q)`. The topology is determined by $x$ and is
not a separate construction argument. Mesh identity depends on $x$ and $Q$;
replacing the coordinates reconstructs the mesh from the new coordinates.

A mesh describes one geometry, not an indexed family of geometries. Coordinates
may depend on symbolic coefficients, deformation fields or scalar parameters,
but have no free indices. Construction and reconstruction enforce this rule.
A solver may evaluate a form for different coordinate values or manage several
meshes outside the language. Shape differentiation does not require an indexed
mesh (see [Shape derivatives](extensions.md#shape-derivatives)).

Deformed surfaces, isogeometric patches and analytic maps are meshes of
coordinate functions. Continuity across cells is a property of the space of
$x$. Independent patches or points that support fields use a topology with
independent cells; no sides are required. An indexed family used only for
integration or sampling needs a `DomainMap`, not a mesh (see
[Mapped integration families](#mapped-integration-families)).

#### Cell maps and pieces

For group $j$ of $T$, let $K_j$ be its reference cell. `pieces` is the derived
`tuple` of maps
`domain_map(K_j, pull_back(x, cell_map(T, j)), Q)`.
These maps are views of the mesh, not construction arguments. They carry the
cell index `cell_index(T, j)` and do not introduce extra geometry parameters.

The ambient domain $Q$ is a coordinate domain. A mesh embedded in another
mesh has its own physical coordinates; cell or side maps specify its
relation to the other mesh's cells.

#### Mesh cell maps

`cell_map(D, j, e)` returns `CellMap(D, j, e)`, whose target is the mesh $D$.
Its parameterisation is `pull_back(x, cell_map(T, j, e))`. Piece $j$ instead
has target $Q$ and records only physical coordinates. These maps have the
same parameterisation but different meanings: only the cell map selects a
cell of $D$ for evaluating a mesh field.

For a field $f$ on $D$, `pull_back(f, cell_map(D, j))` is its local
expression on $K_j$, indexed by the cells of group $j$. Pulling back
`coordinate(D)` gives `pull_back(x, cell_map(T, j))`.

#### Examples

Use the two-triangle topology of the [static topology
example](../examples/library.md#a-static-topology), on reference cells:

```python
K = ReferenceSimplex(2)
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
```

Let `x` be a vector Lagrange function of degree 1 on `square`, of value shape
`(2,)`, with vertex values $(0, 0)$, $(1, 0)$, $(0, 1)$ and $(1, 1)$.
Its image is the square:

```python
D = image(x, RD(2))
```

- `D.topology` is `square`. `cell_map(D, 0)` is the family of the two
  triangles, with the free index `cell_index(square, 0)`.
- `pull_back(coordinate(D), cell_map(D, 0))` is
  `pull_back(x, cell_map(square, 0))`, which the library of the space of `x`
  evaluates on triangle 1 to $(1 - X_1, X_0 + X_1)$.

For a topology `M` of triangles and quadrilaterals, with the cell types
`(Kt, Kq)`, `image(x, RD(2))` has one piece for each cell group. A function `u`
on `M` is one expression, pulled back through each cell group:

```python
# On Kt with the free index cell_index(M, 0), and on Kq with cell_index(M, 1).
u_tri = pull_back(u, cell_map(M, 0))
u_quad = pull_back(u, cell_map(M, 1))
```

```{dropdown} Legacy UFL
In legacy UFL, a `Mesh` is `image(x, Q)` for the coordinate element of `x`, and `SpatialCoordinate(mesh)` is `coordinate(D)`.
```

### `TransformationLaw`

`TransformationLaw` specifies how function values on a reference cell transform
to a cell of a mesh.

| `TransformationLaw(variances, weight)` |                                                                       |
| -------------------------------------- | --------------------------------------------------------------------- |
| Extends                                | `Base` |
| Free function                          | `transformation_law(variances, weight)`                               |
| Takes                                  | `variances: tuple[str, ...]`, `weight: int` |
| Properties                             | `variances`, a `tuple` of `"covariant"` and `"contravariant"`, and `weight` $w$ |
| Raises                                 | `TypeError` for an argument of another type, and `ValueError` for another string |

There is one variance for each of the last axes of the value shape, and the
leading axes are unchanged. A library may type the variances as a `Literal`.

For a cell $\psi$ of a mesh, with $J$, $K$ and $\rho$ of [Jacobian and
density](#jacobian-and-density), $k$ the dimension of the topology and $n$ the
dimension of the ambient coordinate domain, the law $m$ maps a value $\hat f$
on the reference cell to $m_\psi(\hat f)$. $m_\psi$ multiplies each covariant
axis by $K^T$ and each contravariant axis by $J$, and divides by
$(\det J)^w$ for $k = n$ and by $\rho^w$ for $k < n$. Each transformed axis
has the extent $k$ in $\hat f$ and $n$ in $m_\psi(\hat f)$.

The pull-back of a gradient by a cell and the normal of a side use the
covariant law of weight 0 (see [Expanding pull-backs](evaluation.md#expanding-pull-backs)).
The density of an integral is not an instance, because a measure uses
$\lvert \det J \rvert$ and the law uses $\det J$.

The maps of Basix are instances. The core names no instances:

| Space                             | Law                                       | Legacy UFL map         |
| --------------------------------- | ----------------------------------------- | ---------------------- |
| Lagrange, vector Lagrange         | `((), 0)`                                 | identity               |
| Nédélec, $H(\mathrm{curl})$       | `(("covariant",), 0)`                     | covariant Piola        |
| Raviart–Thomas, $H(\mathrm{div})$ | `(("contravariant",), 1)`                 | contravariant Piola    |
| Discontinuous densities           | `((), 1)`                                 | $L^2$ Piola            |
| Regge                             | `(("covariant", "covariant"), 0)`         | double covariant       |
| Hellan–Herrmann–Johnson           | `(("contravariant", "contravariant"), 2)` | double contravariant   |

The law has no orientation. A space of odd weight on a topology of lower
dimension than the ambient coordinate domain gives its own signs for each pair
of cell and facet, as the space gives the transformations of its degrees of
freedom.

### `PushForward`

`PushForward` gives the field on a mesh $D$ of an expression $f$ on its topology
$T$.

| `PushForward(f, D, m)` |                                                                     |
| ---------------------- | ------------------------------------------------------------------- |
| Extends                | `AbstractExpression` |
| Takes                  | `f: AbstractExpression`, `D: Mesh`, `m: TransformationLaw` |
| Properties             | `expression` $f$, `mesh` $D$ and `law` $m$                        |
| Requires               | $D$ a mesh with the topology $T$, $f$ with the domain $T$ or no domain, and the last axes of $f$, one for each variance of $m$, of the extent $k$ |
| Static type            | The value shape of $f$ with the extent $n$ for each transformed axis, the domain $D$ |
| Raises                 | `TypeError` for $D$ or $m$, and `ValueError` if another requirement fails |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `push_forward(f, D, m)` | | `f: AbstractExpression`, `D: Mesh`, `m: TransformationLaw` | `PushForward(f, D, m)` |

The pull-back of `PushForward(f, D, m)` by `cell_map(D, j, e)` is
$m_\psi(\hat f)$, with $\psi$ the cell and $\hat f$ =
`pull_back(f, cell_map(T, j, e))`. The push-forward of the coordinate function
$x$ of $D$ by `transformation_law((), 0)` is `coordinate(D)`.

The transformation law is always explicit, for example
`push_forward(u, D, transformation_law(V))` for `u` in `V`, where the
one-argument `transformation_law` is supplied by the library of `V`.
The arguments of a form, the variable and the direction of a derivative,
and the keys of the data are the
functions on $T$. A library may wrap `push_forward`, as the
[classic FEniCS layer](../examples/forms.md#a-classic-fenics-layer-for-the-poisson-problem) does.

The identity of Piola holds on each cell of $D$, with no derivative of $J$:

- The divergence of a push-forward by `(("contravariant",), 1)` is the
  divergence on $K_j$ of $\hat f$, carried by `((), 1)`.
- For $k = n$, the curl of a push-forward by `(("covariant",), 0)` is the curl
  on $K_j$ of $\hat f$, carried by `(("contravariant",), 1)` for $n = 3$ and
  by `((), 1)` for $n = 2$.

[`expand_pull_backs`](evaluation.md#expanding-pull-backs) applies the identity, and every
other `Divergence` and `Curl` lowers to its composition of `Gradient`.

### Restriction

The restriction of an expression $f$ on a topology $T$ or a mesh $D$ to the
cells of a topology $S$ with sides, from the side $s$, is the pull-back of $f$
along the inclusion of each cell of $S$ into the cell on the side $s$. The
side tables of $S$ give the cell on the side, through a
[`GatheredIndex`](expressions.md#gatheredindex) entry, and the inclusion of each cell of $S$
into that cell. The restriction is `pull_back(f, side_map(T, S, s))` for $f$
on $T$, and `pull_back(f, side_map(D, E, s))` for $f$ on $D$ and a side
domain $E$ (see [`SideMap`](#sidemap)).

Notation, for the group $i$ of $S$ with the side `(s, T, j)`:

- `(cells, inclusion, reference_normal)` are the side tables of the side $s$
  of group $i$.
- $F_i$ is the cell type of group $i$, and $K_j$ the reference cell of
  group $j$ of $T$.
- `e = cell_index(S, i)` and `g = gathered_index(cells, (e,))`.

#### The inclusion

Two families of maps with the free index `e` make the inclusion of the cells
of group $i$:

- $B_i$, the cell on the side, is `cell_map(T, j, g)` on $T$, and
  `cell_map(D, j, g)` on a mesh $D$ with the topology $T$.
- $A_i$, the inclusion into the cell on the side, is `inclusion` of the side
  tables, a domain map from $F_i$ into $K_j$. $F_i$ must be a reference cell.

The restriction of $f$ to group $i$ from the side $s$ is
`pull_back(pull_back(f, B_i), A_i)`, an expression on $F_i$ with the free index
`e`.

For an explicit cell selection `q` (an `int`, `Index` or `GatheredIndex`),
use the same side data with their cell index replaced:

```python
B_selected = replace(B_i, {e: q})
A_selected = replace(A_i, {e: q})
```

The selected restriction is `pull_back(pull_back(f, B_selected), A_selected)`.
Substitute only into the side data, not into `f`; its free indices remain
independent. The default selection `q = e` leaves the maps unchanged.

Restrictions are compared at paired points; their difference is an ordinary
interface jump when those points physically coincide (see [Side tables](#side-tables)).
A cell with more than two sides has no jump; a form uses each side's
restriction separately.

#### `SideMap`

`SideMap` is the lazy map of the cells of $S$ into $T$, or of the cells of a
side domain $E$ into a mesh $D$, through the side $s$.

| `SideMap(X, P, s)` |                                                                                                         |
| ------------------ | ------------------------------------------------------------------------------------------------------- |
| Extends            | `Base` |
| Takes              | `X, P: AbstractTopology \| Mesh`, `s: str` |
| Properties         | `ambient_domain` $X$, `parameter_domain` $P$, and `side` $s$ |
| Derived            | `free_indices`, always `()` |
| Requires           | Conditions 1 to 3 below                                                                                 |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_map(X, P, s)` | | `X, P: AbstractTopology \| Mesh`, `s: str` | `SideMap(X, P, s)` |

1. $X$ is a topology $T$ and $P$ a topology $S$, or $X$ is a mesh with the
   topology $T$ and $P$ a side domain with the topology $S$.
2. Every cell type of $S$ is an `AbstractReferenceCell`.
3. Every group of $S$ has a side named $s$ with the topology $T$.

A side map of a mesh takes a side domain, so that a coupling of two meshes
shares one side domain. `pull_back(f, side_map(X, P, s))` is an expression
with the domain $P$, the restriction of $f$ to all the groups at once. The
pull-back of that expression by `cell_map(P, i, q)` uses
`pull_back(pull_back(f, B_selected), A_selected)`, with the selected maps
defined above. This covers both topology and mesh targets.

#### Side domains

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_domain(D, S, s)` | | `D: Mesh`, `S: AbstractTopology`, `s: str` | `image(pull_back(x, side_map(T, S, s)), Q)` |

For a mesh $D$ = `image(x, Q)` with the topology $T$, the side domain of $S$
is the image of the trace of $x$ from side $s$, so the side domain moves with
$x$. The result has topology $S$ and geometry supplied by that side.

The chosen side supplies the integration measure. Any measure conversion
belongs in the integrand.

#### Traces and normals

The restriction of an expression is the value of the expression on the cell
of $S$, from the side:

- The restriction of `gradient(f, D)` is the whole gradient of `f` on the cell
  on the side, as for a normal derivative. The last axis is in the frame of
  the ambient coordinate domain from every side.
- `gradient(pull_back(f, side_map(D, E, s)), E)` is the tangential gradient of
  the restriction along the cells of $E$.

The unit side conormal is [`normal(side_map(D, E, s))`](#normal), outward on
boundary pieces.

#### Example

The example uses `D` and `square` of the
[`Mesh`](#mesh) example, `I = ReferenceSimplex(1)`, a
function `u` on `square`, and the diagonal of `square`:

```python
diagonal = StaticTopology((I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),))
plus = side_map(square, diagonal, "+")
minus = side_map(square, diagonal, "-")
jump = pull_back(u, plus) - pull_back(u, minus)
E = side_domain(D, diagonal, "+")
```

From `"+"`, the side tables are the cells `(0,)` and the inclusion
$A_+(t) = (1 - t, t)$, through the local vertices 1 and 2 of `K`. From
`"-"`, the side tables are the cells `(1,)` and the inclusion
$A_-(t) = (0, t)$, through the local vertices 0 and 2. The coordinates of
`E` are $(1 - t, t)$.

Suppose the reference-cell value `pull_back(u, cell_map(square, 0))` is:

```python
a[c, 0] + a[c, 1] * X[0] + a[c, 2] * X[1]
```

Here `c = cell_index(square, 0)` and `X = coordinate(K)`. Substituting the cell
number and inclusion from each side gives:

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

```{dropdown} Legacy UFL
The legacy UFL restriction `f('+')` is `pull_back(f, side_map(T, S, "+"))` on the topology $S$ of the interior facets, whose sides are named `'+'` and `'-'`.
```

### A family of points in a mesh

The following constructions use cells and reference positions supplied by
the solver (see [Functions in a space](forms.md#functions-in-a-space)).

For a mesh $D$ with topology $T$, suppose a family indexed by $i$ has selected
cells in group $j$, with reference cell $K_j$. The solver supplies a table
`cells` with extent `T.num_cells[j]` and an expression `positions[i]` for the
reference positions. Evaluation is the nested pull-back:

```python
g = gathered_index(cells, (i,))
R = domain_map(RD(0), positions[i], K_j)
values = pull_back(pull_back(f, cell_map(D, j, g)), R)
```

Families spanning several cell groups are split by group. The inner
pull-back selects the cell; the outer pull-back evaluates within it. The
cell selection is retained when maps are substituted or expressions lowered.

#### Adjacent and coincident cells

Let group 0 contain two interval cells with reference cell $K = [0, 1]$.
On a mesh $D$ with cells $[0, 1]$ and $[1, 2]$, let $f$ equal 2 on cell 0
and 7 on cell 1. These snippets describe symbolic values, not numerical
execution by the core:

```python
f0 = pull_back(f, cell_map(D, 0, 0))
f1 = pull_back(f, cell_map(D, 0, 1))
midpoint = domain_map(RD(0), tensor((Float(0.5),)), K)
left_endpoint = domain_map(RD(0), tensor((Integer(0),)), K)
right_endpoint = domain_map(RD(0), tensor((Integer(1),)), K)

pull_back(f0, midpoint)  # Physical coordinate 0.5; value 2.
pull_back(f0, right_endpoint)  # Physical coordinate 1; value 2.
pull_back(f1, left_endpoint)  # Physical coordinate 1; value 7.
```

The last two expressions are the two traces at the interface. Physical
coordinate 1 alone does not choose between them.

For a mesh `D_overlap` whose two cells both occupy $[0, 1]$, let $g$ again
have values 2 and 7 on the respective cells:

```python
g0 = pull_back(g, cell_map(D_overlap, 0, 0))
g1 = pull_back(g, cell_map(D_overlap, 0, 1))
pull_back(g0, midpoint)  # Physical coordinate 0.5; value 2.
pull_back(g1, midpoint)  # Physical coordinate 0.5; value 7.
```

The physical and reference coordinates coincide, but the selected cells differ.

#### Points forming a domain

A reusable domain of points may instead have a topology $S$ with one
zero-dimensional cell per point and a side `("inside", T, j)`. Its side
tables store the selected cells and reference positions supplied by the
solver. Set `P = side_domain(D, S, "inside")`:

```python
inside = side_map(D, P, "inside")
values = pull_back(pull_back(f, inside), cell_map(P, 0))
total = integral(pull_back(f, inside), P)
```

This representation supports fields, restrictions and integration on the
points. A single evaluation or indexed family needs only the nested
pull-backs above, without constructing a point topology.

### `Integral` over a mesh

`Integral(f, D)` sums over all cells of mesh $D$, counting coincident cells
separately. The integrand has domain $D$, domain `D.ambient_domain`, or no
domain. Every free index of $f$ remains free in the result. The cell indices
used by integration are internal to the expansion; a matching name in $f$
does not make that index a cell selector.

If all cell counts are zero, the mesh is empty and the integral has value
zero with its declared static type. Each zero-cell group contributes an empty
index sum, without evaluating the integrand or geometry for any cell. This
also applies when symbolic cell counts are supplied as zero.

#### Expansion over cells

For group $j$ of `D.topology`, let $K_j$ be its reference cell and $e_j$ an
internal cell index fresh for $f$. The deterministic choice of $e_j$ is stated
in [Expanding integrals](evaluation.md#expanding-integrals). Set
`D_j = cell_map(D, j, e_j)`. For a scalar field on $D$, the contribution is:

```python
index_sum(integral(pull_back(f, D_j) * density(D_j), K_j), (e_j,))
```

For a field on `D.ambient_domain`, use the map with the same parameterisation
and that ambient domain in the pull-back. A domain-free integrand is also
pulled back, so the internal cell index is present even for zero-dimensional
cells, whose density is `Integer(1)`. Array-valued integrands are integrated
component by component.

The density does not depend on orientation. For ambient fields, the result
agrees with integration over the geometric union when each cell map is injective
almost everywhere and distinct images overlap only on sets of zero measure.

#### Cell data and indexed loads

Cell-dependent data are fields on the topology or mesh. A field equal to 2
on one unit cell and 7 on another has integral 9. In contrast, a domain-free
indexed load with values `[2, 7]` is integrated over the whole two-cell mesh
for each index value, giving `[4, 14]`. Its index remains free, even when it
has the name `cell_index(D.topology, j)`.

To pair a raw indexed datum with explicit cell maps, pull back to reference
cells and use `IndexSum` explicitly. Mesh integration does not infer that
pairing from the names. Substitution into a mesh integral therefore cannot
capture an integrand index; the expansion enforces the same rule.

#### Examples

The examples use the following objects:

- `D` and `square` of the [`Mesh`](#mesh) example.
- A field `f` on `D`, such as `push_forward(u, D, m)` for a function `u` on
  `square` (see [`PushForward`](#pushforward)).
- The mesh `B` of the boundary facets of `D`, with the side
  `"inside"`.

- `integral(Integer(1), D)` is the area 1 of the square, the sum of two
  integrals over `K`.
- `integral(pull_back(f, side_map(D, B, "inside")), B)` restricts `f` to the
  boundary before integrating it.

```{dropdown} Legacy UFL
`f * dx`, `f * ds` and `f * dS` are `integral(f, D)` over the mesh of the cells, of the boundary facets and of the interior facets, with the fields restricted to the facets in `f`. The three are one integral over three domains.
```

### `Normal`

`Normal` is the unit side conormal: tangent to the parent cell and orthogonal
to the side cell, with the orientation supplied by `reference_normal`.
It points outward on boundary pieces.

| `Normal(sm)` |                                                                                                       |
| ------------ | ----------------------------------------------------------------------------------------------------- |
| Extends      | `AbstractExpression` |
| Free function | `normal(sm)`                                                                                         |
| Takes        | `sm: SideMap` |
| Properties   | `side_map`, a [`SideMap`](#sidemap) `SideMap(D, E, s)`                                               |
| Requires     | $D$ a mesh, and $E$ of the dimension $k - 1$ for the dimension $k$ of $D$                             |
| Static type  | Value shape `(n,)` for the dimension $n$ of the ambient coordinate domain, the domain $E$, real      |
| Raises       | `TypeError` for an `sm` that is not a `SideMap`, and `ValueError` if the requirement fails           |

The normal of `sm` belongs to the same side as the restriction
`pull_back(f, sm)`, so a form writes both from one side map.

Notation, for `sm` = `SideMap(D, E, s)` and the piece $i$ of $E$:

- $B_i$ and $A_i$ are the maps of [Restriction](#restriction) for the side $s$.
- $K_j$ is the reference cell of the cell on the side.
- $k$ is the dimension of $D$, and $n$ is the dimension of the ambient
  coordinate domain.

#### The pull-back of `Normal`

The pull-back of `Normal(SideMap(D, E, s))` by `cell_map(E, i, q)` is
$\nu = w / \lvert w \rvert$, with $w = m_{B_{\mathrm{selected}}}(\hat m)$,
using the selected maps of [Restriction](#restriction), where:

- $m$ is `transformation_law(("covariant",), 0)`, applied with the Jacobian
  of `B_selected` on the side at the points of `A_selected` (see
  [`TransformationLaw`](#transformationlaw)).
- $\hat m$ is `replace(reference_normal, {cell_index(E.topology, i): q})`,
  using the normal from the side tables of side $s$ of group $i$ (see
  [Side tables](#side-tables)).

The composition is `w / square_root(w @ w)`. It requires full column rank
of the cell Jacobian and $w \ne 0$ wherever evaluated.

The core needs no numbering of the facets and no vertices of $K_j$, because
the library of the topology gives $\hat m$.

#### The direction of the normal

- For a cell of $E$ on the boundary of $K_j$, such as a facet or a piece of a
  facet, $\hat m$ points out of $K_j$ (see [Side tables](#side-tables)).
- For a cell of $E$ inside $K_j$, such as a cut through the cell, $\hat m$
  has the orientation chosen by the supplying library.
- For $k = n$, outward normals are opposite when the two regions adjoin on
  opposite sides of the interface. Coincident geometry or side names alone
  do not imply this relation.
- For $k < n$, the same conclusion also requires matching tangent spaces.
  On a surface, a jump therefore uses each side's conormal, as in
  `u_p * n_p + u_m * n_m`.

`Normal` is the normal of a cell on a side, not the unit normal field of a
domain of dimension $n - 1$, such as of a surface in $R^3$. The unit normal
field uses a separately supplied orientation field (see
[The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface)).

#### Examples

The examples use the following objects:

- `B` of the [`Integral` over a mesh](#integral-over-a-mesh) examples.
- The mesh `E` of the interior facets of `D`, with the sides
  `"+"` and `"-"`.
- `grad_f = gradient(f, D)` and `i = index("i")`.

```python
inside, plus, minus = side_map(D, B, "inside"), side_map(D, E, "+"), side_map(D, E, "-")
n_B = normal(inside)

# The flux of grad f out of the square, the legacy UFL dot(grad(f), n) * ds.
flux = integral(index_sum(pull_back(grad_f, inside)[i] * n_B[i], (i,)), B)

# The jump of the normal derivative across the diagonal, from the side "+".
n_E = normal(plus)
jump = pull_back(grad_f, plus)[i] - pull_back(grad_f, minus)[i]
penalty = integral(index_sum(jump * n_E[i], (i,)) ** Integer(2), E)
```

Worked examples for this chapter: [an interior facet integral](../examples/coupling.md#an-interior-facet-integral), [a Bank-Weiser error estimator](../examples/estimation.md#a-bank-weiser-error-estimator), [a fluid and a solid](../examples/coupling.md#a-fluid-and-a-solid), [a 1D network in a 3D mesh](../examples/coupling.md#a-1d-network-in-a-3d-mesh), [a deformed domain](../examples/surfaces.md#a-deformed-domain).

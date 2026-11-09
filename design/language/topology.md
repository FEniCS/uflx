# Topological meshes

UFLx makes the topological mesh an explicit object on which fields and function
spaces with geometry-independent reference-cell representations can be defined
before physical coordinates are supplied. A coordinate expression on that
topology then defines a physical mesh, represented by
[`Mesh`](geometry.md#mesh) in the next chapter.

```{dropdown} Legacy UFL and Firedrake
Legacy UFL places functions on a mesh domain, with the mapping supplied by
the element. UFLx places suitable functions on a topology and represents them
on a physical mesh through an explicit [`PushForward`](geometry.md#pushforward).
Firedrake already makes this separation at the solver level.
```

```{dropdown} Spaces defined directly on a mesh
Spaces such as Hermite, Argyris and Morley have geometry-dependent
reference-cell representations and are defined directly on a physical mesh.
The element library constructs the local basis from the cell's coordinate map
and evaluates cell pull-backs and reference derivatives (see
[Functions in a space](forms.md#functions-in-a-space)).
```

A topology identifies cells and their side relations. It supports fields and
function spaces when its cell types are reference cells.

Reference-cell inclusions and reference normals describe relations between
cells. These use reference coordinates, independently of physical geometry.
Distinct cells retain their identities even when a later geometric realisation
places them at the same physical coordinates.

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

`AbstractTopology` represents a topological mesh: a set of cells of given
cell types, without physical coordinates.

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

## Domains

Reference cells have coordinates of their own, independently of any physical
mesh. The following domain classes provide the foundation for these reference
coordinates and for the physical domains of the next chapter.

### `AbstractDomain`

`AbstractDomain` represents a domain, the set on which an expression is defined.
The core provides coordinate domains and `Mesh`. Extensions may subclass
`AbstractDomain`, supplying their own geometric semantics and rewrite rules.

| `AbstractDomain` |                                                                             |
| ---------------- | --------------------------------------------------------------------------- |
| Extends          | `Base` |
| Provides         | `dimension`, an `int` $k \ge 0$, and `ambient_domain`, an `AbstractDomain` |
| Defined by       | `AbstractCoordinateDomain`, `Mesh`, and extension subclasses |

The ambient domain of $D$ supplies the coordinates naming its points. The
first coordinate domain reached by following `ambient_domain` repeatedly is
the *ambient coordinate domain*, of dimension $n \ge k$. A coordinate domain
is its own ambient domain.

### `AbstractCoordinateDomain`

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$ whose
points are named by the coordinates of the points themselves. The core cannot
build the set from smaller objects. A library defines a new coordinate domain,
such as a reference cell, by extending this abstract class. On a coordinate
domain the coordinates are canonical: the coordinate function is the identity
(see [`Coordinate`](geometry.md#coordinate)), and the integral is the Lebesgue integral in
the coordinates (see [`Integral`](geometry.md#integral)). On a core `Mesh`,
coordinates and integration are defined through cell parameterisations.

| `AbstractCoordinateDomain` |                                              |
| -------------------------- | -------------------------------------------- |
| Extends                    | `AbstractDomain` |
| Provides                   | `dimension`, an `int` $n \ge 0$              |
| Derived                    | `ambient_domain`, the coordinate domain itself |
| Defined by                 | A subclass, such as `RD` or a reference cell |

Examples:

- The unit simplex in $R^n$ is a coordinate domain. A codimension-1 face is
  not: it has empty interior in $R^n$ and is instead the image of a domain map
  from the unit simplex in $R^{n-1}$.
- A mesh of triangles is *not* a coordinate domain, because the core builds the
  mesh from a topology and a coordinate expression; its cells are parameterised
  by reference triangles.

### `RD`

`RD(n)` is the whole of $R^n$.

| `RD(n)`    |                                                                 |
| ---------- | --------------------------------------------------------------- |
| Extends    | `AbstractCoordinateDomain` |
| Properties | `dimension`, an `int` $n \ge 0$                                 |
| Raises     | `ValueError` if $n < 0$ |

Worked examples for this chapter: [the unit simplex](../examples/library.md#the-unit-simplex).

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
| Requires | $j$ a valid cell-group number of $T$, `T.cell_types[j]` a reference cell $K_j$; an integer $e$ satisfies the bounds below, and a `GatheredIndex` has extent `T.num_cells[j]` |
| Derived | `topology` $T$, `parameter_domain` $K_j$, `free_indices` and, for a mesh, `parameterisation`, as below |
| Raises | `ValueError` if another requirement fails |

For a symbolic cell count, construction cannot check the upper bound. An
integer entry $e$ must be nonnegative. If `T.num_cells[j]` is an `int`,
construction also checks `e < T.num_cells[j]`; otherwise the caller must ensure
that $e$ is smaller than the supplied cell count.

For an ordinary `Index` $e$, `free_indices` contains $e$ with extent
`T.num_cells[j]`; the index itself stores no extent. For a `GatheredIndex`,
it is the free indices of $e$. For an integer entry it is `()`.

A cell map of a topology has no parameterisation. The same API selects cells
of a physical mesh; its derived parameterisation and reconstruction are
defined under [Mesh cell maps](geometry.md#mesh-cell-maps).

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_map(X, j, e)` | | `X: AbstractTopology \| Mesh`, `j: int`, `e: int \| Index \| GatheredIndex` | `CellMap(X, j, e)` |
| `cell_map(X, j)` | | `X: AbstractTopology \| Mesh`, `j: int` | `CellMap(X, j, cell_index(T, j))`, with $T$ the topology of $X$ |

[`pull_back(f, cell_map(X, j, e))`](geometry.md#pullback) is the value of an expression $f$ on $X$
on the reference cell $K_j$ of each selected cell. For a function on a
topology, the library of its space gives the value (see [Functions in a
space](forms.md#functions-in-a-space)). For a mesh field, lowering uses the selected
cell and its parameterisation (see [Expanding pull-backs](evaluation.md#expanding-pull-backs)).

Worked examples for this chapter: [the unit simplex as a reference cell](../examples/library.md#the-unit-simplex-as-a-reference-cell).

## Sides of a topology

A topology $S$ whose cells are sub-entities of cells of other topologies has
sides. A cell of $S$ is seen from one cell of another topology on each side,
and [Restriction](#restriction) picks the value of an expression on a cell of
$S$ by the side name.

### `sides`

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

### The image of a cell on a side

Each side inclusion maps the cell type $F_i$ into the parent reference cell
$K_j$. Its image may be the whole cell, a sub-entity, or a piece such as half
a facet or a cut through the cell. A whole-cell inclusion preserves vertex
order. The mesh library supplies the map; an affine inclusion can be built
from the vertex positions in $K_j$ (see [A static topology](../examples/library.md#a-static-topology)).
For periodic sides, those positions refer to the corresponding periodic image.

### Side tables

A concrete topology with sides also provides `side_tables(i, side)`, the
triple `(cells, inclusion, reference_normal)` for the side named `side` of
group $i$, with the side `(side, T, j)`, the cell type $F_i$ of group $i$ and
`e = cell_index(S, i)`:

- `cells` is an [`AbstractIndexTable`](expressions.md#abstractindextable) of the shape
  `(S.num_cells[i],)` and the extent `T.num_cells[j]`. Entry $e$ selects a cell of
  group $j$ of $T$ on the side of the cell $e$ of $S$.
- `inclusion` is a [`DomainMap`](geometry.md#domainmap) from $F_i$ into $K_j$ with the
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
[Restriction](#restriction) and [`Normal`](geometry.md#normal) read the tables.

### What the concrete mesh decides

The mesh library chooses the cells or pieces of $S$, the side names, and
one parent cell on each side. For a lower-dimensional entity shared by several
cells, this choice affects the restriction of a discontinuous field.
Sides may refer to different topologies, and a side topology may itself have
sides, allowing facets of facets.

### Cells that are pieces

Pieces cover hanging facets, mortar intersections and segments of a network
inside mesh cells. The mesh library computes these intersections and their
inclusions; a curved inclusion may be approximate. See the
[network example](../examples/coupling.md#a-1d-network-in-a-3d-mesh).

### Examples

For a triangle topology `T` and facet cell type `I`, a boundary topology has
cell types `(I,)` and sides `((("inside", T, 0),),)`. Interior facets instead
have sides `((("+", T, 0), ("-", T, 0)),)`.

See [Side topology examples](../examples/coupling.md#side-topology-examples)
for corners, submeshes, periodic boundaries and mixed cell types.

```{dropdown} Legacy UFL
The legacy UFL measures `ds` and `dS` integrate over a topology of facets with one and two sides, and legacy UFL names the two sides `'+'` and `'-'`. The legacy UFL measure `dx(1)` integrates over a topology of the cells of a subdomain, and the legacy UFL measure `dP` integrates over a topology of vertices.
```

Worked examples for this chapter: [the simplex as a cell type](../examples/library.md#the-simplex-as-a-cell-type), [a static topology](../examples/library.md#a-static-topology).

## Restriction

The restriction of an expression $f$ on a topology $T$ or a mesh $D$ to the
cells of a topology $S$ with sides, from the side $s$, is the pull-back of $f$
along the inclusion of each cell of $S$ into the cell on the side $s$. The
side tables of $S$ give the cell on the side, through a
[`GatheredIndex`](expressions.md#gatheredindex) entry, and the inclusion of each cell of $S$
into that cell. The restriction is `pull_back(f, side_map(T, S, s))` for $f$
on $T$, and `pull_back(f, side_map(D, E, s))` for $f$ on $D$ and a side
domain $E$ (see [`SideMap`](#sidemap) and [Side domains](geometry.md#side-domains)).

Notation, for the group $i$ of $S$ with the side `(s, T, j)`:

- `(cells, inclusion, reference_normal)` are the side tables of the side $s$
  of group $i$.
- $F_i$ is the cell type of group $i$, and $K_j$ the reference cell of
  group $j$ of $T$.
- `e = cell_index(S, i)` and `g = gathered_index(cells, (e,))`.

### The inclusion

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

### `SideMap`

`SideMap` is the lazy map of the cells of $S$ into $T$, or of the cells of a
side domain $E$ into a mesh $D$, through the side $s$.

| `SideMap(X, P, s)` |                                                                                                         |
| ------------------ | ------------------------------------------------------------------------------------------------------- |
| Extends            | `Base` |
| Takes              | `X, P: AbstractTopology \| Mesh`, `s: str` |
| Properties         | `ambient_domain` $X$, `parameter_domain` $P$, and `side` $s$ |
| Derived            | `free_indices`, always `()` |
| Requires           | The domain conditions below |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_map(X, P, s)` | | `X, P: AbstractTopology \| Mesh`, `s: str` | `SideMap(X, P, s)` |

| Domain condition | Requirement |
| ---------------- | ----------- |
| $X$, $P$ | Topologies $T$, $S$, or a mesh on $T$ and a side domain on $S$ |
| Cell types of $S$ | All are `AbstractReferenceCell` |
| Groups of $S$ | Each has a side named $s$ with topology $T$ |

A side map of a mesh takes a side domain, so that a coupling of two meshes
shares one side domain. `pull_back(f, side_map(X, P, s))` is an expression
with the domain $P$, the restriction of $f$ to all the groups at once. The
pull-back of that expression by `cell_map(P, i, q)` uses
`pull_back(pull_back(f, B_selected), A_selected)`, with the selected maps
defined above. This covers both topology and mesh targets.

## A field on a topological mesh

Use the library's [static topology](../examples/library.md#a-static-topology)
and reference triangle:

```python
K = ReferenceSimplex(2)
T = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
```

Let `u` be a scalar function in a space on `T`. Its library supplies a
reference-cell expression for each cell (see [Functions in a
space](forms.md#functions-in-a-space)):

```python
u0 = pull_back(u, cell_map(T, 0, 0))
u1 = pull_back(u, cell_map(T, 0, 1))
```

Both expressions have domain `K`, but select different cells of `T`.
No physical vertex positions are needed. Reference derivatives such as
`gradient(u0, K)` are defined on `K`; `gradient(u, T)` is not defined.
The next chapter supplies physical coordinates with `D = mesh(x, Q)`.


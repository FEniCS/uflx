# Topological meshes

UFLx describes local kernels as mathematics. A topology describes the local
cell configuration used by a kernel: its cell types and named side relations.
The finite element package selects concrete mesh cells and supplies their
local coefficient and coordinate data. Mesh-wide cell counts, global cell
numbers, neighbour lookup, iteration over a mesh and global assembly are
outside the core. Runtime-sized local arrays and index maps can still describe
a complete patch calculation.

A local configuration may contain one cell, the two parents of an interface,
or the cells of a fixed macro element. Distinct local cells retain their
identities even when they have the same type or physical coordinates. A mixed
mesh uses configurations for its cell and interface types; a library can
specialise one mathematical recipe for these configurations.

Fields with geometry-independent reference-cell representations can be defined
on a topology before coordinates are supplied. A coordinate expression defines
a physical local mesh, represented by [`Mesh`](geometry.md#mesh).

```{dropdown} Spaces defined directly on a mesh
Spaces such as Hermite, Argyris and Morley have geometry-dependent
reference-cell representations and are defined directly on a physical local
mesh. The element library constructs the local basis from the coordinate map
and evaluates cell pull-backs and reference derivatives (see
[Functions in a space](forms.md#functions-in-a-space)).
```

## Topology

### `AbstractCellType`

`AbstractCellType` describes a cell type with numbered vertices.

| `AbstractCellType` | |
| ------------------ | - |
| Extends | `Base` |
| Provides | `dimension`, the nonnegative topological dimension as an `int`, and `num_vertices`, the positive number of vertices as an `int` |
| Defined by | A subclass outside the core, such as `Simplex` |

The core does not number facets or other sub-entities. A library decides which
sets of local vertices describe sub-entities and supplies their reference
maps. A vertex set identifies an entity; an ordered correspondence determines
its parameterisation. Pieces such as cuts need not be sub-entities.

A polygon's vertex count belongs to its cell type. A vertex count known only
at evaluation is not yet covered; a library specialises kernels by cell type.

### `AbstractTopology`

`AbstractTopology` describes a finite local cell configuration without physical
coordinates. Each position in `cell_types` is a distinct local cell, called a
*cell slot*. Slots identify the mathematical inputs of a kernel, not concrete
cells in the package's mesh.

| `AbstractTopology` | |
| ------------------ | - |
| Extends | `Base` |
| Provides | `cell_types`, `sides`, and `side_data(cell_slot, side)` for declared sides |
| Derived | `dimension`, the dimension of the cell types |
| Defined by | A subclass outside the core |

- `cell_types` is a nonempty `tuple` of `AbstractCellType`s of the same
  dimension. Repeated types are allowed: `(K, K)` describes two distinct
  local cells of type `K`.
- `sides` is `()`, or one tuple of side declarations per slot, as defined
  under [Sides of a topology](#sides-of-a-topology).
- For a topology `topology`, a `cell_slot` is an `int` satisfying
  `0 <= cell_slot < len(topology.cell_types)`.
  Negative numbers do not select slots from the end.

For example, `(K,)` describes a single-cell kernel, `(K, K)` the parents of
an interior-facet kernel, and `(Kt, Kq)` a triangle–quadrilateral pair.
An Alfeld macro element has three triangle slots `(K, K, K)`. A list of
alternative kernel types is managed by the library; it is not a configuration
whose cells should all be integrated together.

The number of slots is fixed when the expression is constructed. An empty
selection in the package causes no kernel calls; it does not require an empty
symbolic topology. A runtime-sized patch can instead be expressed through
packed local arrays, index tables and sums, as in the
[patch problem](../examples/patches.md#a-mean-constrained-dg-patch-correction).
The package supplies the packing and local-to-patch degree-of-freedom maps;
no runtime topology is required.

A topology whose cell types are reference cells may be the domain of an
expression or space. Its points are cell-labelled reference points; continuity
between cells is a property of the space. A topology is not an `AbstractDomain`,
so `Coordinate`, `Gradient` and `Integral` reject it.

Identity is structural. A concrete topology stores its identifying arguments,
such as a label and the local configuration. The supplying library associates
local coefficient and coordinate inputs with each slot. The same configuration
can be evaluated repeatedly with data from different concrete cells.

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
| Provides         | `dimension`, the nonnegative topological dimension as an `int`, and `ambient_domain`, the `AbstractDomain` supplying coordinates for its points |
| Defined by       | `AbstractCoordinateDomain`, `Mesh`, and extension subclasses |

For a domain `domain`, `domain.ambient_domain` supplies the coordinates
naming its points. The first coordinate domain reached by following
`ambient_domain` repeatedly is the *ambient coordinate domain*. Its dimension
is at least `domain.dimension`. A coordinate domain
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
| Provides                   | `dimension`, the nonnegative number of coordinate axes as an `int`              |
| Derived                    | `ambient_domain`, the coordinate domain itself |
| Defined by                 | A subclass, such as `EuclideanSpace` or a reference cell |

Examples:

- The unit simplex in $R^n$ is a coordinate domain. A codimension-1 face is
  not: it has empty interior in $R^n$ and is instead the image of a domain map
  from the unit simplex in $R^{n-1}$.
- A mesh of triangles is *not* a coordinate domain, because the core builds the
  mesh from a topology and a coordinate expression; its cells are parameterised
  by reference triangles.

### `EuclideanSpace`

`EuclideanSpace(dimension)` is the whole coordinate space with `dimension` axes.
Writing $n$ for `dimension`, this is $\mathbb{R}^n$.

| `EuclideanSpace(dimension)`    |                                                                 |
| ---------- | --------------------------------------------------------------- |
| Extends    | `AbstractCoordinateDomain` |
| Properties | `dimension`, the nonnegative number of coordinate axes as an `int`                                 |
| Raises     | `ValueError` if `dimension` is negative |

Worked examples for this chapter: [the unit simplex](../examples/library.md#the-unit-simplex).

## Reference cells

### `AbstractReferenceCell`

`AbstractReferenceCell` is a cell type that is also a coordinate domain.
It extends `AbstractCoordinateDomain` and `AbstractCellType`. The core reads no vertices of a reference cell. A library
that builds inclusions, such as from vertex positions, keeps the vertices of
its reference cells (see [The unit simplex as a reference
cell](../examples/library.md#the-unit-simplex-as-a-reference-cell)).

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `AbstractCoordinateDomain`, `AbstractCellType` |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](../examples/library.md#the-unit-simplex-as-a-reference-cell)                                                   |

### `CellMap`

`CellMap(parent, cell_slot)` selects one local cell of `parent`, which is a
topology or a physical local mesh. The map goes from that cell's reference
cell into the labelled local cell of `parent`. Pulling a field back through
it gives the field's reference-cell representation on the selected slot.
The slot is part of the kernel's local configuration, not a runtime cell
number in the package's mesh.

| `CellMap(parent, cell_slot)` | |
| --------------------------- | - |
| Extends | `Base` |
| Takes | `parent: AbstractTopology \| Mesh`, `cell_slot: int` |
| Properties | `parent`, the topology or mesh carrying the field, and `cell_slot`, the local cell to read |
| Requires | `cell_slot` is a valid slot of `topology`, and `topology.cell_types[cell_slot]` is a reference cell |
| Derived | `topology`, `reference_cell`, `free_indices`, and, for a mesh, `parameterisation`, as described below |
| Raises | `ValueError` if a requirement fails |

| Derived property | Meaning |
| ---------------- | ------- |
| `topology` | `parent` itself when it is a topology; otherwise `parent.topology`. |
| `reference_cell` | `topology.cell_types[cell_slot]`, the domain of the field after pull-back. |
| `free_indices` | Always `()`; selecting a cell introduces no array index. |
| `parameterisation` | For a mesh, its coordinate expression pulled back to the selected reference cell. |

A cell map of a topology has no coordinate parameterisation. For a physical
mesh its parameterisation and reconstruction are defined under
[Mesh cell maps](geometry.md#mesh-cell-maps).

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_map(parent, cell_slot)` | | `parent: AbstractTopology \| Mesh`, `cell_slot: int` | `CellMap(parent, cell_slot)` |

```python
local_cell = cell_map(parent=volume_mesh, cell_slot=0)
local_field = pull_back(field, local_cell)
# local_field has domain local_cell.reference_cell.
```

The library of the field's space supplies its local representation and data
(see [Functions in a space](forms.md#functions-in-a-space)). Ordinary indices
still describe basis functions, components, quadrature points and other local
arrays.

## Sides of a topology

A topology `S` can describe integration cells or pieces viewed from named
parent cells. Side declarations identify parent slots; local side data supply
the mathematical inclusions and normals. The package supplies the concrete
neighbouring cells' data before kernel evaluation.

### `sides`

`S.sides` is `()` when there are no sides. Otherwise it has one tuple per slot
of `S`. A tuple may be empty; its entries have distinct names. Each entry is
`(side, parent_topology, parent_cell_slot)`:

| Entry | Meaning |
| ----- | ------- |
| `side` | A `str` naming this side, such as `"+"` or `"inside"`. |
| `parent_topology` | The topology containing the parent cell from which a field is read. |
| `parent_cell_slot` | A valid integer cell slot in `parent_topology`. |

Both the side cell and the selected parent cell must be reference
cells. Side order has no meaning; names select sides. Parent dependencies
obey the acyclic construction rules of `Base`.

For a parent configuration `T` with cell types `(K, K)` and an interval `F`,
a single interior facet has cell types `(F,)` and sides
`((("+", T, 0), ("-", T, 1)),)`. The different parent slots preserve the two
local fields. A boundary facet needs only `((("inside", T, 0),),)`.
Parents may have different cell types or belong to different topologies.

### The image of a cell on a side

For slot `i` of `S` with reference cell $F_i$ and side `(s, T, j)`, an
inclusion $A_i:F_i\to K_j$ describes the points at which the parent field
is evaluated. It can describe a whole cell, a facet, part of a facet or an
interior cut. Vertex permutations are represented by the map, including for
whole-cell inclusions.

All sides of one integration cell use the same parameter point in $F_i$.
For an ordinary interface the parent coordinate maps composed with their
inclusions agree. Physical coincidence is a library-supplied assumption, not
a core check; periodic and mortar couplings may pair distinct physical points.

### Side data

`topology.side_data(cell_slot, side)` supplies the geometry of a declared
side. `cell_slot` selects a cell of `topology`; `side` is the name of one of
that cell's side declarations. The method returns
`(inclusion, reference_normal)`.

The returned objects have the following roles:

- `inclusion` is a [`DomainMap`](geometry.md#domainmap) from the selected
  reference cell into the parent reference cell, with no free indices. Its
  real parameterisation has no argument spaces.
- When the selected cell has dimension one less than the parent cell,
  `reference_normal` is a real vector expression on the selected reference
  cell or with no domain. It has one component per parent reference
  coordinate, no free indices and no argument spaces. It is nonzero and
  orthogonal to the image tangents wherever evaluated. Its length need not
  be one. It points out of the parent reference cell on a boundary piece;
  the library chooses its orientation on an interior cut.
- In other dimensions, `reference_normal` is `None`.

The inclusion maps into the parent reference cell. Its differential has full
column rank wherever geometry formulas require it; integration uses the given
parameterisation and counts multiplicity. The core assumes these mathematical
conditions rather than validating concrete mesh geometry.

The data may depend on supplied local coefficients. Expressing an inclusion
in reference coordinates does not imply independence from physical geometry:
a cut or mortar intersection can move. Libraries expose this dependence as
symbolic operands or explicit reconstruction and derivative rules (see
[Geometry dependence](geometry.md#geometry-dependence)). The package computes
intersections and supplies data; UFLx describes the resulting mathematics.

### Mixed cells and facet types

A cell-type pair alone does not determine an interface. The configuration also
specifies its reference integration cell and one inclusion per side. For a
triangle and quadrilateral, take $F=[0,1]$ and, for example,
$A_+(t)=(1-t,t)$ and $A_-(t)=(0,t)$. In three dimensions $F$ can be a reference
triangle or quadrilateral. A library specialises the same form for the needed
cell and interface configurations and supplies orientation-dependent local
maps; the core does not enumerate combinations or dispatch kernels.

### Cells that are pieces

Hanging facets, mortar intersections and embedded network segments use the
same side data. The package supplies each piece and its inclusions. A fixed
local collection can be described by several slots; a runtime number of
pieces can instead be processed by repeated kernel calls or local quadrature
arrays. See the [network example](../examples/coupling.md#a-1d-network-in-a-3d-mesh).

## Restriction

Restriction is evaluation through a named side's inclusion. For a field `f`
on `T`, use `pull_back(f, side_map(T, S, s))`. For a field on a physical local
mesh `D`, use `pull_back(f, side_map(D, E, s))`, where `E` is a mesh on `S`.
Its geometry determines the integration measure, while the chosen side
specifies where the field is read.

### The inclusion

For slot `i` with side `(s, T, j)`, let `(A_i, reference_normal)` be
`S.side_data(i, s)`. Let `B_i = cell_map(T, j)` for a topology field, or
`B_i = cell_map(D, j)` for a field on a mesh with topology `T`.
The local restriction is

```python
pull_back(pull_back(f, B_i), A_i)
```

It has domain $F_i$ and the free indices of `f`. There is no neighbour-table
lookup, cell-index substitution or implicit summation. Restrictions can be
combined at the paired points. For more than two sides the core defines no
canonical jump; the form combines the named restrictions explicitly.

### `SideMap`

`SideMap(parent, restriction_domain, side)` describes how to read a field on
`parent` from a named side of `restriction_domain`. Each cell of
`restriction_domain` uses its side declaration to select a parent cell and
an inclusion into that cell.

The map goes from `restriction_domain` into `parent`. Pulling a field back
through it produces a field on `restriction_domain`. The restriction remains
lazy until a cell of that domain is read.

| `SideMap(parent, restriction_domain, side)` | |
| ----------------------------------------- | - |
| Extends | `Base` |
| Takes | `parent, restriction_domain: AbstractTopology \| Mesh`, `side: str` |
| Properties | `parent`, `restriction_domain`, and `side`, as described below |
| Derived | `free_indices`, always `()` |
| Requires | The domain conditions below |
| Raises | `ValueError` if a domain condition fails |

| Property | Meaning |
| -------- | ------- |
| `parent` | The topology or mesh carrying the field to restrict. |
| `restriction_domain` | The topology or mesh on which the restricted field is defined. |
| `side` | The name selecting a parent cell and inclusion for each cell of `restriction_domain`. |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_map(parent, restriction_domain, side)` | | `parent, restriction_domain: AbstractTopology \| Mesh`, `side: str` | `SideMap(parent, restriction_domain, side)` |

For the conditions below, the parent topology is `parent` itself when it is
a topology, and `parent.topology` otherwise. The restriction topology is
obtained from `restriction_domain` in the same way.

| Domain condition | Requirement |
| ---------------- | ----------- |
| Kind of inputs | `parent` and `restriction_domain` are both topologies or both meshes. |
| Side declarations | Every slot of the restriction topology has a side named `side`, whose declared parent topology is the parent topology. |
| Reference cells | Every slot of the restriction topology and every selected parent slot is a reference cell. |

The geometry of `restriction_domain` need not equal the selected parent's
trace. The helper [side_domain](geometry.md#side-domains) constructs a mesh
from a chosen side. Coupled fields share the restriction domain even when
their parent geometries differ.

```python
interface_side = side_map(
    parent=volume_mesh,
    restriction_domain=interface_mesh,
    side="+",
)
restricted_field = pull_back(field, interface_side)
# restricted_field has domain interface_side.restriction_domain.
```

Pulling `restricted_field` back by `cell_map(interface_mesh, i)` gives
`pull_back(pull_back(field, B_i), A_i)`, with the selected parent cell map
`B_i` and inclusion `A_i` defined under [The inclusion](#the-inclusion).

## A field on a topological mesh

Use the library's [static topology](../examples/library.md#a-static-topology)
to describe two distinct triangle slots:

```python
K = ReferenceSimplex(2)
T = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)), ())
u0 = pull_back(u, cell_map(T, 0))
u1 = pull_back(u, cell_map(T, 1))
```

Here `u` is a scalar function in a library space on `T`. Both expressions have
domain `K`, but use different local field data. The vertex labels describe
this fixed configuration, not global mesh connectivity. The package can reuse
the kernel for many pairs of triangles.

Reference derivatives such as `gradient(u0, K)` are defined on `K`;
`gradient(u, T)` is not defined. The next chapter supplies physical local
coordinates with `D = mesh(x, Q)`.

# Topological meshes

A topology describes the local cell configuration of a kernel: one cell, the
parents of an interface, or the cells of a macro element. Side maps relate
these configurations. The library selects concrete mesh cells, supplies their
coefficient and coordinate data, and assembles the results.

Reference cells have **reference geometry**, including their own coordinates.
A topology labels cell slots; when their types are reference cells, each slot
carries that reference geometry. A coordinate field maps the reference cells
into physical space, supplying **physical geometry** through
[`Mesh`](geometry.md#mesh). “Local” describes the scope of the cell configuration;
both its reference geometry and its physical geometry can be local.

Fields whose reference representations are independent of physical geometry
can be defined on a topology. Spaces whose bases depend on physical geometry
are defined directly on a mesh (see [Functions in a space](forms.md#functions-in-a-space)).

## Topology

### `AbstractCellType`

`AbstractCellType` describes a cell type with numbered vertices.

| `AbstractCellType` | |
| ------------------ | - |
| Extends | `Base` |
| Provides | `dimension`, the nonnegative topological dimension as an `int`, and `num_vertices`, the positive number of vertices as an `int` |
| Defined by | A subclass outside the core, such as `Simplex` |

The library defines facets and other sub-entities and supplies their reference
maps. General inclusions, such as cuts, need not follow sub-entities.

### `AbstractTopology`

`AbstractTopology` describes a finite local cell configuration without specifying
its physical geometry. Each position in `cell_types` is a distinct local cell,
called a *cell slot*. Slots identify the mathematical inputs of a kernel.
For reference cells, the topology is a labelled disjoint union: each slot has
reference geometry and retains its identity even when reference cells or
physical coordinates coincide.

| `AbstractTopology` | |
| ------------------ | - |
| Extends | `Base` |
| Provides | `cell_types` |
| Derived | `dimension`, the dimension of the cell types |
| Defined by | A subclass outside the core |

- `cell_types` is a nonempty `tuple` of `AbstractCellType`s of the same
  dimension. Repeated types are allowed: `(K, K)` describes two distinct
  local cells of type `K`.
- For a topology `topology`, a `cell_slot` is an `int` satisfying
  `0 <= cell_slot < len(topology.cell_types)`.
  Negative numbers do not select slots from the end.

For example, `(K,)` describes a single-cell kernel, `(K, K)` the parents of
an interior-facet kernel, and `(Kt, Kq)` a triangle–quadrilateral pair.
The number of slots is fixed when the expression is constructed.

A topology of reference cells may be the domain of an expression or space.
Its points are cell-labelled reference points. The library supplies field
continuity through its [spaces and bindings](forms.md#abstractfunctionspace).
A topology is not an `AbstractDomain`: coordinates, gradients and integrals
are defined on its reference cells or on a physical mesh.

Identity is structural, determined by the concrete topology's construction
arguments, such as a label and its cell configuration. A one-slot topology is
distinct from its reference cell. The library may evaluate the same
configuration repeatedly with data from different concrete cells.

## Domains

Reference coordinates describe points within a reference cell. A physical
coordinate field maps those points into physical space. The following classes
describe the domains used for both kinds of geometry.

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

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$,
defined by a library. Its [`Coordinate`](geometry.md#coordinate) is the identity,
and its [`Integral`](geometry.md#integral) is the Lebesgue integral in those
coordinates. On a `Mesh`, coordinates and integration are instead defined
through cell parameterisations.

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
Libraries may supply vertex coordinates for constructing inclusions; these
are not required by the core.

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `AbstractCoordinateDomain`, `AbstractCellType` |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](../examples/library.md#the-unit-simplex-as-a-reference-cell)                                                   |

### `CellMap`

`CellMap(parent, cell_slot)` selects one local cell of `parent`, which is a
topology or a mesh carrying physical geometry. The map goes from that cell's reference
cell into the labelled local cell of `parent`. Pulling a field back through
it gives the field's reference-cell representation on the selected slot.

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

For a physical mesh, see [Mesh cell maps](geometry.md#mesh-cell-maps) for the
parameterisation and reconstruction rules.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_map(parent, cell_slot)` | | `parent: AbstractTopology \| Mesh`, `cell_slot: int` | `CellMap(parent, cell_slot)` |

The examples below build one configuration step by step. They use the library's
[`ReferenceSimplex`](../examples/library.md#the-unit-simplex-as-a-reference-cell),
[`StaticTopology`](../examples/library.md#a-static-topology), and scalar
`Lagrange(topology, degree, shape)` space. All maps and pull-backs use the core
language directly.

First declare two triangle slots and a scalar field on them:

```{figure} ../figures/reference-slots.*
:label: reference-cell-slots
:alt: Two copies of the unit reference triangle. Slot 0 has vertex labels 0, 1, 2 and its diagonal edge highlighted in blue. Slot 1 has labels 1, 3, 2 and its left edge highlighted in amber. Both selected edges run from label 1 to label 2.
:width: 700px

These are reference triangles, not placements of cells in physical space.
Each slot has its own copy. Bold numbers are library vertex labels. The
highlighted edges show the library-supplied inclusions $A_0$ and $A_1$;
the arrows follow the source ordering from label 1 to label 2.
```

```python
K = ReferenceSimplex(2)                 # Unit reference triangle.
T = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)))
V = Lagrange(T, 1, ())                   # Scalar, degree-one local representation.
u = function(V, "u")                    # Field on the two-slot topology T.

c0 = cell_map(T, 0)                     # K -> slot 0 of T.
c1 = cell_map(T, 1)                     # K -> slot 1 of T.
u0 = pull_back(u, c0)                   # Slot 0's field, expressed on K.
u1 = pull_back(u, c1)                   # Slot 1's field, expressed on K.
```

Both `u0` and `u1` have domain `K`; their cell maps select different local
field data. The vertex labels describe the correspondence along edge `(1, 2)`
in the two reference cells. The library supplies the reference inclusions for
this correspondence; physical geometry has not yet been specified.

## Restriction

A side map is a relation from a source configuration to a parent configuration.
Each source slot selects a parent slot and a reference inclusion.

### `SideMap`

| `SideMap(parent, source, entries)` | |
| -------------------------------- | - |
| Extends | `Base` |
| Takes | `parent, source: AbstractTopology \| Mesh`, `entries: tuple[tuple[int, DomainMap], ...]` |
| Properties | `parent`, `source`, `entries` |
| Derived | `free_indices`, always `()` |
| Requires | One entry per source slot; valid parent slots; inclusions without free indices and with endpoints as below |
| Raises | `ValueError` if an entry or domain condition is invalid |

For source reference cell $F_i$, entry `(j, A)` requires
`A.parameter_domain == F_i` and `A.ambient_domain == K_j`, the selected parent
reference cell. Both must be reference cells. Parent slots may repeat.
Each endpoint may independently be a topology or a mesh; for a mesh, its
topology determines its reference cells.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_map(parent, source, *entries)` | | `parent, source: AbstractTopology \| Mesh`, `entries: tuple[int, DomainMap]` | `SideMap(parent, source, entries)` |
| `compose(outer, inner)` | | `outer, inner: SideMap` | The composed `SideMap`, requiring `inner.parent == outer.source` |

Continue with one interval slot for that edge, ordered from vertex label 1
to vertex label 2. Let `A0` and `A1` be library-supplied `DomainMap`s from
its reference interval `I` into `K`, selecting the indicated edges. Their
[construction](../examples/library.md#worked-edge-inclusions) is separate from
the core side-map declarations:

```{mermaid}
flowchart LR
    accTitle: Two inclusions of the source interval
    accDescr: A0 maps the source interval to the diagonal of the reference triangle for slot 0. A1 maps it to the left edge of the reference triangle for slot 1.
    I["Source interval I<br/>Vertex order 1 → 2"]
    K0["K for slot 0<br/>Diagonal edge"]
    K1["K for slot 1<br/>Left edge"]
    I -->|A0| K0
    I -->|A1| K1
    classDef source fill:#f1f5f9,stroke:#64748b,color:#1e293b
    classDef first fill:#eff6ff,stroke:#2563a6,color:#1e293b
    classDef second fill:#fffbeb,stroke:#b45309,color:#1e293b
    class I source
    class K0 first
    class K1 second
    linkStyle 0 stroke:#2563a6,stroke-width:2px
    linkStyle 1 stroke:#b45309,stroke-width:2px
```

```python
I = ReferenceSimplex(1)                 # Unit reference interval.
S = StaticTopology((I,), ((1, 2),))      # One edge slot.

# A0, A1: library-supplied reference inclusions I -> K.

plus = side_map(T, S, (0, A0))           # Read slot 0 through A0.
minus = side_map(T, S, (1, A1))          # Read slot 1 through A1.
u_plus = pull_back(u, plus)             # Both traces have domain S.
u_minus = pull_back(u, minus)
jump = u_plus - u_minus                 # A field on the edge topology S.
```

There is one entry in each side map because `S` has one slot. In `(1, A1)`,
the integer selects parent slot 1, while `A1` specifies the reference points
at which its field is read. Both inclusions respect the source vertex order,
so the two traces read corresponding points along the edge.

`plus` and `minus` are names for these maps. Neither chooses a normal;
orientation is supplied separately to [Normal](geometry.md#normal).
The same construction accepts meshes carrying physical geometry as endpoints.
A source mesh then has its own physical geometry and measure, independently
of either parent trace. For the two traces to describe a matching physical
interface, the parent coordinate fields must satisfy
$x_0(A_0(t)) = x_1(A_1(t))$. The reference inclusions alone do not impose this.

### The inclusion

To express the edge trace on the reference interval, select the sole
cell of `S`. Reading through the side map is equivalent to selecting the
parent cell first, then applying its inclusion:

```python
edge = cell_map(S, 0)                   # I -> the edge slot of S.
trace0 = pull_back(u_plus, edge)        # Domain I.

# The same field, with the parent read written out:
trace0_explicit = pull_back(u0, A0)     # Slot 0's trace on I.
trace1_explicit = pull_back(u1, A1)     # Slot 1's trace on I.
jump_on_I = pull_back(jump, edge)       # Their difference, on I.
```

In general, for `sm.entries[i] = (j, A)`, restriction satisfies

```python
pull_back(pull_back(f, sm), cell_map(sm.source, i))
# becomes
pull_back(pull_back(f, cell_map(sm.parent, j)), A)
```

The result has domain `A.parameter_domain` and the free indices of `f`.

### Composition

If `inner.entries[i] = (j, A)` and `outer.entries[j] = (k, B)`, the composite
entry is `(k, C)`, where

```python
C = domain_map(A.parameter_domain, pull_back(B.parameterisation, A), B.ambient_domain)
```

The composite is a `SideMap` from `inner.source` to `outer.parent`, requiring
`inner.parent == outer.source`. This covers permutations, repeated selections
and nested inclusions.

For a concrete composition, reverse the edge ordering from the preceding
example. Let `B` be the library-supplied reference map that reverses `I`:

```{mermaid}
flowchart LR
    accTitle: Composition with a reversed edge
    accDescr: Reverse maps R to S, then plus maps S to triangle slot 0 of T. The direct map reversed_plus is their composite.
    R["R · edge slot<br/>Vertex order 2 → 1"]
    S["S · edge slot<br/>Vertex order 1 → 2"]
    T["T · slot 0"]
    R -->|reverse| S
    S -->|plus| T
    R -->|reversed_plus| T
    classDef edge fill:#f1f5f9,stroke:#64748b,color:#1e293b
    classDef parent fill:#eff6ff,stroke:#2563a6,color:#1e293b
    class R,S edge
    class T parent
    linkStyle 2 stroke:#2563a6,stroke-width:2px
```

```python
R = StaticTopology((I,), ((2, 1),))      # Edge ordered from label 2 to label 1.
# B: library-supplied reversal I -> I.
reverse = side_map(S, R, (0, B))        # R -> S, reversing the edge ordering.
reversed_plus = compose(plus, reverse)  # R -> T, reading parent slot 0.

reversed_trace = pull_back(u, reversed_plus)
# Equivalent to:
reversed_trace_nested = pull_back(pull_back(u, plus), reverse)
```

Both expressions have domain `R`. Composition selects the same parent slot
as `plus`, with its inclusion composed with `B`.

## A field on a topological mesh

The domains of the resulting expressions are:

| Expression | Domain | Meaning |
| ---------- | ------ | ------- |
| `u` | `T` | Field on two labelled triangle slots |
| `u0`, `u1` | `K` | One selected triangle's reference representation |
| `u_plus`, `u_minus`, `jump` | `S` | Traces and their difference on one labelled edge slot |
| `trace0`, `jump_on_I` | `I` | Edge fields expressed in the interval coordinate |

Reference derivatives such as `gradient(u0, K)` are defined on `K`.
Selecting a cell is what exposes those coordinates; `gradient(u, T)` is
not defined. A coordinate field `x` on `T` supplies physical geometry through
`D = mesh(x, Q)`, where `Q` is the physical coordinate space. On slot 0, for
example, the reference edge point `A0(t)` has physical position `x0(A0(t))`,
where `x0 = pull_back(x, c0)`. See [Geometry](geometry.md) for this construction.

# Core UFLx language definition

This note describes the smallest core language behind UFLx. The note leaves
out implementation details, compiler algorithms and broad redesign ideas.
The argument behind the design, with a reading guide, is in DESIGN.md. The
library classes of the examples and further worked examples are in
EXAMPLES.md. Implementation principles and the enforcement of the rules of
the language are in DESIGN_CHOICES.md. The note assumes knowledge of
differential geometry and functional analysis, and explains some advanced
concepts of differential geometry in more detail.

This document is written for users and developers who know domain-specific
languages for variational forms, such as UFL and ngsolve.

## Conventions

The note uses "must" for a rule that an implementation or a library has to
satisfy, "should" for a rule that an implementation or a library satisfies
unless the note states a reason not to, and "may" for a choice that the note
leaves open. A statement in the present tense, such as "the degree is a
`tuple`", is a rule of the language. "Not yet covered" marks a decided gap, a
part of the language that the authors intend to specify later. "Open" marks a
question that the authors have not decided. A "Legacy UFL" subsection relates
a section to legacy UFL and is not a rule.

## Scope

### What UFLx is and isn't

In one sentence:

> UFLx is a typed symbolic language for writing variational forms.

The UFLx core language, henceforth the core, should not use finite element or
numerical analysis terminology. Differential geometry and functional analysis
should be all that the core requires. This draft requires calculus on
manifolds in coordinate-free form. The coordinate-free form gives a compact,
dimension-free statement of many concepts already in UFL.

Where finite element or numerical analysis concepts are needed, a layer above
the core or an external library provides them. Basix is the in-house example.
Neither Basix nor DOLFINx, nor their assumptions, should be visible in the
core.

UFLx is *not* a symbolic language for reasoning about problems in functional
analysis or differential geometry. The core should therefore contain only the
concepts and relations needed to write variational forms as code for later
algorithmic transformation.

A UFLx program is ordinary Python code that constructs symbolic objects and
connects them through a graph. The result is usually a variational form, such
as a bilinear form, a linear form or a functional.

UFLx does not handle concrete data or numerical algorithms, such as meshes,
boundary markers, values of coefficients or assembly procedures. Those belong
to the surrounding problem-solving environment, such as DOLFINx, Firedrake,
FIAT and Basix. UFLx specifies a set of leaf nodes, and the external library
provides the concrete data behind each leaf node.

### Basic principles

- Each core object definition should follow the mathematical definition of
  the object closely.
- The core must not add definitions, properties or terminology outside
  differential geometry and functional analysis.
- The core values strictness and expressivity over brevity. Higher-level APIs
  can add brevity later.
- The mathematical fields state variational problems uniformly across
  geometric dimensions. Special cases for particular dimensions do not exist
  in the mathematics, so they should not exist in the core.
- The core does not contain algorithms for the convenience of a form
  compiler. A form compiler states the forms that the form compiler accepts.
  A pass outside the core must lower a form to those forms and must produce
  only node classes of the core and of the
  [language extensions](#language-extensions) that the form compiler accepts.
  For example, a basic form compiler might only accept integrals over a
  reference cell, such as a unit simplex, as in each cell term of a legacy UFL
  `dx` integral.

#### Legacy UFL

TSFC lowers a UFL form to GEM, which repeats much of the algebra of UFL, partly because UFL simplifies during construction. The UFLx core does not simplify during construction, so a form compiler may be able to use a lowered UFLx form directly, adding only leaves for tables and kernel inputs. This has not been confirmed with the TSFC or FFCx developers.

## Rules of the language

### The base object `Base`

The base object in UFLx must be extended by all objects that are included in
the graph. The graph is a directed acyclic graph (DAG). An edge of the graph
means "is defined from". The note writes "edge" for an edge of the graph, and
"edge of a triangle" or "edge of a cell" for the geometric edge.

Everything needed to define the variational problem must be recoverable from
the graph. No information is held outside of the graph.

| `Base`     |                                                                    |
| ---------- | ------------------------------------------------------------------ |
| Properties | `construction_arguments`, the construction arguments, as in rule 4 |

Objects must satisfy the following rules:

1. A construction argument must be one of a small curated set of types: a base
   object, or exactly an `int`, a finite `float`, a `str`, or a `tuple` whose
   items are of these types. A subclass of these Python types is not accepted,
   so a `bool` is not an `int` and a named tuple is not a `tuple`. Construction
   arguments are therefore deeply immutable, have a hash that is deterministic
   across processes on the same machine, and are compared by type and value.
   `1`, `1.0` and `True` are distinct, `0.0` and `-0.0` are distinct, and every
   value is equal to itself. Construction raises `TypeError` for an argument of
   any other type, and `ValueError` for a `float` that is not finite.
2. Construction must depend only on the construction arguments, not on global
   state.
3. The state of an object is immutable after construction. Every result that
   the object produces after construction depends only on the class and the
   construction arguments.
4. An object must be reconstructible from its class and construction arguments,
   and the reconstructed object must be equal to the original. The base object
   records the construction arguments, and the property
   `construction_arguments` returns the construction arguments, in order, as a
   `tuple`. The constructor of a subclass receives the construction arguments,
   validates the construction arguments, and may store data derived from the
   construction arguments.
5. The edges of an object must be exactly the base objects among its
   construction arguments, including the base objects nested in `tuple`s, in
   order.
6. Any order that the core or a pass relies on must be deterministic across
   processes. A pass, and any other walk of the graph such as hashing,
   serialisation or printing, follows the edges of each object in order (rule
   5). Where unrelated objects must be sorted, the sort uses a total order
   defined by the class and the construction arguments alone.
7. Equality (`==`) must be structural and return a `bool`. Whether `is` is
   structural is implementation dependent. An implementation in which `is` is
   structural must intern objects, so that equal objects are the same object.
   Equality does not decide mathematical equivalence, so $a + b$ and $b + a$
   are not equal. Mathematical equivalencies are applied in later graph
   algorithms.
8. All validation of an object happens in the constructor of the class of the
   object. This rule is a hard rule. A free function never validates and never
   converts arguments. A free function only reads derived properties of
   arguments and constructs objects. An invalid argument is rejected by the
   constructor of the first object that the argument makes invalid. Unless the
   section of a class states otherwise, the constructor raises `TypeError` for
   an argument of a type or a class that the constructor does not accept, and
   `ValueError` for an argument of an accepted type or class with a value that
   the constructor does not accept.
9. A constructor has no default values, so every construction argument is
   written at construction. A free function may supply a default value.
10. A concrete class states an `identifier`, a `str` with the public name of
    the class qualified by the package name, such as `"uflx.Sum"` for the core
    or `"library.UnitSimplex"` for a library. The hash of rule 1, the total
    order of rule 6 and serialisation use the identifier and the construction
    arguments, and never the location of the class in the code. Construction
    raises `TypeError` for a class with no identifier. The package name keeps
    the identifiers of different libraries apart, and a library keeps the
    identifiers of the classes of the library distinct.

Hence two objects are equal if and only if the two objects have the same class
and equal construction arguments.

Because objects are immutable, an object need only be validated once, at
construction. Every object that exists is therefore valid.

### Leaves and subclassing

Requires: [The base object `Base`](#the-base-object-base).

See also: [`RD`](#rd), [`DomainMap`](#domainmap), [Language
extensions](#language-extensions).

Only leaves are subclassed. A leaf is an object that the core cannot build
from other objects, such as the coordinate domain `RD(n)` or a function whose
values an external library supplies. A terminal is a leaf that is an
expression, such as `Coordinate`, `Function` and `Normal` of the core, or a
leaf expression of a downstream library. An object built from other objects is a
single class of the core or a free function that returns objects of the core,
never a subclass. For example, every domain is a `DomainMap`, and an
affinely mapped triangle is a free function that returns a
`DomainMap`. A language extension is the one exception (see
[Language extensions](#language-extensions)).

A concrete class of the core is final. A library subclasses only an abstract
class of the core, such as `AbstractCoordinateDomain`, so a pass that
dispatches on a concrete class of the core never meets a subclass of the
class.

A parent class never constructs a subclass. A property of a parent class may
return an object whose class is a subclass, because the property only
returns an object that already exists. A subclass may narrow the return
type of a property.

### Properties and free functions

Requires: [The base object `Base`](#the-base-object-base).

See also: [`AbstractExpression`](#abstractexpression),
[`DomainMap`](#domainmap), [`Coordinate`](#coordinate),
[Transpose and matrix product](EXAMPLES.md#transpose-and-matrix-product).

A property reads and a free function builds. Data that an object stores, or
derives from the class and the construction arguments of the object, is a
property of the object, such as `value_shape` of an expression or
`parameter_domain` of a domain. A property never constructs an object. Any
operation that constructs objects is a pure free function, such as
`coordinate(P)` or `volume_form(P)`. A class method, where
needed, immediately calls a free function.

A constructor takes the construction arguments by position only, in order, and
raises `TypeError` for an argument by name. So an object has one `tuple` of
construction arguments. A free function may take arguments by name, such as
`domain_map(P=T, psi=psi, Q=RD(2))`.

A stored property returns a construction argument. A derived property returns
data derived from the class and the construction arguments. A concrete class
has the derived properties of the abstract classes of the concrete class.

### No simplification during construction

Requires: [The base object `Base`](#the-base-object-base).

Simplification and rewriting should be avoided entirely while the user
constructs the graph, so the graph keeps the intent of the user. Only later
passes, chosen by the user, should rewrite. For example, in `(a + b) + c` in
floating point arithmetic, where `a`, `b` and `c` are terminals supplied
later, the user may need `(a + b)` evaluated before `c` is added.

### Passes

Requires: [The base object `Base`](#the-base-object-base).

See also: [`Index`](#index), [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum).

A pass is a pure free function that walks a graph after construction of the
graph. The result of a pass depends only on the arguments of the pass. A pass
returns a value, such as a unit, or a new graph. A pass builds a new graph by
constructing an object for each node that the pass replaces, and by
reconstructing every object that depends on a replaced node, up to the root
(see [The base object](#the-base-object-base)). An object depends on a
replaced node through any edge, so a pass that replaces an object also
reconstructs an [`Index`](#index) whose label contains the object. Hence a
replacement commutes with the free functions: replacing `A` by `B` in
`transpose(A)` gives `transpose(B)`, with the derived indices of `B`. A pass
may cache results for each object, because objects are immutable.

Substitution avoids capture. A `ComponentTensor` or an `IndexSum` binds
indices (see [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum)), and a [`Form`](#form) binds
functions. If a pass substitutes an expression
into the expression of a binding node, and a bound index of the binding node is
free in the substituted expression, the pass renames the bound index. The new
name is the derived index `Index(("uflx.rename", node, m))`, or the derived
`EllipsisIndex` with the same label, or the function
`Function(V, ("uflx.rename", node, m))` in the space $V$ of a bound function,
with `node` the binding node before the substitution and `m` the position of
the bound index or function. For example, replacing
`c` by `b[i]` in `index_sum(a[i] * c, (i,))` gives
`index_sum(a[j] * b[i], (j,))`, with `j` the renamed index, and the result has
the free index `i`. Without the renaming, the result would be
`index_sum(a[i] * b[i], (i,))`, with no free index.

A pass does one thing, such as differentiating or lowering. Passes compose,
and the user of the core chooses which passes run and in what order.

### Kinds of code that produce objects

Requires: [The base object `Base`](#the-base-object-base), [Properties and free
functions](#properties-and-free-functions), [Passes](#passes).

See also: [Transpose and matrix product](EXAMPLES.md#transpose-and-matrix-product),
[`Indexed`](#indexed), [New questions](#new-questions).

Three kinds of code produce objects:

| Kind                                               | When the code runs              | Inspects arguments to choose what to build                                | Validates                     |
| -------------------------------------------------- | ------------------------------- | ------------------------------------------------------------------------- | ----------------------------- |
| Constructor                                        | Whenever an object is made      | No, a constructor inspects arguments only to validate                     | Yes, and only constructors do |
| Construction free function, such as `volume_form`  | While the user builds the graph | Only the static type, to build one fixed composition for each static type | No                            |
| Pass, such as differentiation, lowering or `units` | After construction              | Yes, and a pass may rewrite the graph                                     | No                            |

A construction free function fails when a constructor rejects the
arguments, for example `restriction(f, E, s)` raises `ValueError` from the
constructor of `SideMap` when no group of the topology of `E` has a side named
`s`. A pass fails with `TypeError` for a node class that
the pass has no rule for, or with an error of the pass, such as inconsistent
units.

## Domains

The function spaces and operators of UFLx are defined on a domain and need the
topological structure of the domain, which is the dimension, boundary,
orientation and decomposition into cells of the domain. The core needs no
metric beyond the Euclidean metric of the ambient coordinate domain, which
gives the components of a differential form on a domain of lower dimension
(see [`Components`](#components)) and the unit normal (see [`Normal`](#normal)).

### `AbstractDomain`

Requires: [The base object `Base`](#the-base-object-base).

See also: [Support](#support), [`Coordinate`](#coordinate),
[`AbstractCoordinateDomain`](#abstractcoordinatedomain),
[`DecomposedDomain`](#decomposeddomain).

An abstract base class for domains `AbstractDomain`. A domain is the set on
which an expression is defined. A domain is a coordinate domain, a leaf that a
library supplies (see [`AbstractCoordinateDomain`](#abstractcoordinatedomain)),
or a decomposed domain, the union of the images of families of domain maps
that the core builds (see [`DecomposedDomain`](#decomposeddomain)). Every node class that
takes a domain accepts either kind, and a core transformation gives the meaning
of the node on a decomposed domain through the pieces of the decomposed domain.

| `AbstractDomain` |                                                                             |
| ---------------- | --------------------------------------------------------------------------- |
| Extends          | `Base` |
| Provides         | `dimension`, an `int` $k \ge 0$, and `ambient_domain`, an `AbstractDomain` |
| Defined by       | `AbstractCoordinateDomain` and `DecomposedDomain`                           |

The ambient domain of a domain $D$ is the domain whose coordinates name the
points of $D$. The ambient domain of a coordinate domain is the coordinate
domain itself. The ambient domains of $D$ are `D.ambient_domain`, the ambient
domain of `D.ambient_domain`, and so on. The *ambient coordinate domain* of
$D$ is the first coordinate domain among the ambient domains of $D$, with the
dimension $n \ge k$.

### `AbstractCoordinateDomain`

Requires: [`AbstractDomain`](#abstractdomain).

See also: [Support](#support), [`Coordinate`](#coordinate),
[`DomainMap`](#domainmap).

An abstract subclass of `AbstractDomain` for a coordinate domain, a leaf
domain.

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$ whose
points are named by the coordinates of the points themselves. The core cannot
build the set from smaller objects, so a new coordinate domain, such as a
reference cell, is a subclass of `AbstractCoordinateDomain`. On a coordinate
domain the coordinates are canonical: the coordinate function is the identity
(see [`Coordinate`](#coordinate)), the volume form is the coordinate volume form
(see [Volume form](#volume-form)), and the integral is the Lebesgue integral in
the coordinates (see [`Integral`](#integral)). Every other domain gets these by
pull-back.

| `AbstractCoordinateDomain` |                                              |
| -------------------------- | -------------------------------------------- |
| Extends                    | `AbstractDomain` |
| Provides                   | `dimension`, an `int` $n \ge 0$              |
| Derived                    | `ambient_domain`, the coordinate domain itself |
| Defined by                 | A subclass, such as `RD` or a reference cell |

Examples:

- A face of codimension 1 of the unit simplex in $R^n$ is *not* a coordinate
  domain, because the face has empty interior in $R^n$. The face is the image
  of the unit simplex in $R^{n-1}$ under a map, so the face is the image of a
  domain map (defined below).
- A mesh of triangles is *not* a coordinate domain, because the core builds the
  mesh from the reference triangle and the maps of the cells. The mesh is a
  decomposed domain.

### `RD`

Requires: [`AbstractCoordinateDomain`](#abstractcoordinatedomain).

A concrete subclass of `AbstractCoordinateDomain` for the whole of $R^n$.

| `RD(n)`    |                                                                 |
| ---------- | --------------------------------------------------------------- |
| Extends    | `AbstractCoordinateDomain` |
| Properties | `dimension`, an `int` $n \ge 0$                                 |
| Raises     | `TypeError` if $n$ is not an `int`, and `ValueError` if $n < 0$ |

Worked examples for this chapter in EXAMPLES.md: [the unit simplex](EXAMPLES.md#the-unit-simplex).

## Expressions and static type

### `AbstractExpression`

Requires: [The base object `Base`](#the-base-object-base).

| `AbstractExpression` |                                                                                      |
| -------------------- | ------------------------------------------------------------------------------------ |
| Extends              | `Base`                                                                               |
| Provides             | `value_shape`, `free_indices`, `scalar_type`, `support`, `degree`, `argument_spaces` |

Every `AbstractExpression` derives the value shape, the free indices, the
scalar type, the support, the degree and the [argument
spaces](#argument-spaces) of the expression. The value shape is
a `tuple` of [extents](#extents), one for each axis. An expression with value
shape `()` is a scalar.

### `Index`

Requires: [The base object `Base`](#the-base-object-base), [Properties and free
functions](#properties-and-free-functions).

See also: [`Indexed`](#indexed), [Transpose and matrix
product](EXAMPLES.md#transpose-and-matrix-product).

> LLM generated content.

A concrete subclass of the base object for a free index. An `Index` is not an
expression. Index notation writes tensor algebra as fixed compositions, so
that operators such as `transpose` and `matrix_product` are free functions
rather than node classes.

| `Index(label)` |                                                                                                 |
| -------------- | ----------------------------------------------------------------------------------------------- |
| Extends        | `Base` |
| Free function  | `index(label)`                                                                                  |
| Takes          | `label: str \| tuple` |
| Properties     | `label`, a `str` for a user index, or a `tuple` that contains graph objects for a derived index |
| Raises         | `TypeError` if the label is neither a `str` nor a `tuple`                                       |

An `Index` has no extent of its own, and takes the extent of the axis that the
`Index` indexes.

The user names user indices. Two free occurrences of one `Index` in one scope
are the same index. An index that the user means to be new needs a new name.
The scope of an index is given in [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum).

A free function that needs fresh indices derives the indices from the
operands. The label of a derived index is the identifier of the free function,
followed by the arguments of the free function and a number. The identifier is
the public name of the free function, qualified by the package name, such as
`"uflx.cell_index"` for the core or `"library.transpose"` for a library. For
example, `transpose(A)` of the [transpose
example](EXAMPLES.md#transpose-and-matrix-product) uses
`Index(("library.transpose", A, 0))` and `Index(("library.transpose", A, 1))`.
An index whose label contains `A` has an
edge to `A`. Every node that uses an `Index` has an edge to the `Index`. So an index
whose label contains `A` cannot occur inside `A`, because `A` would then
contain `A`, and the graph is acyclic. Hence
a derived index is fresh with respect to the operands, and the free function
does not inspect the operands to find a fresh index. A global counter is not
allowed, because a counter makes construction depend on global state.

A pass that needs fresh indices derives the indices in the same way, with the
identifier of the pass, the arguments of the rule and a number. For example,
`expand_components` writes the new axis of `Components(g, P)` with the index
`Index(("uflx.expand_components", g, P, 0))`. The renamed index of
[Passes](#passes) follows the same rule.

The edges of an index to the objects in the label of the index identify the
index and are not operands. A pass that derives a property from the value of
an expression, such as the free indices or the support, reads the edges of
each node that are expressions, and an `Index` is not an expression, so such a
pass never meets the label of an index. A pass that replaces an object
reconstructs every index whose label contains the object (see
[Passes](#passes)), so the label of a derived index always names the current
operands, and the derived index stays fresh with respect to the current
operands by the argument above.

### `EllipsisIndex`

Requires: [`Index`](#index).

See also: [`Indexed`](#indexed), [Tensor algebra with
`EllipsisIndex`](#tensor-algebra-with-ellipsisindex).

> LLM generated content.

A concrete subclass of the base object for an index that covers any number of
axes. An `EllipsisIndex` is not an expression.

| `EllipsisIndex(label)` |                                                           |
| ---------------------- | --------------------------------------------------------- |
| Extends                | `Base` |
| Free function          | `ellipsis_index(label)`                                   |
| Takes                  | `label: str \| tuple` |
| Properties             | `label`, with the rules of the label of an `Index`        |
| Raises                 | `TypeError` if the label is neither a `str` nor a `tuple` |

An `EllipsisIndex` covers every axis of an `Indexed` that the other entries do not
cover, so an `EllipsisIndex` may cover no axis. The extent of an `EllipsisIndex` is the
`tuple` of the extents of the covered axes. In every other rule of index
notation, an `EllipsisIndex` behaves as an `Index`.

### Extents

Requires: [The base object `Base`](#the-base-object-base).

See also: [Free indices](#free-indices), [`Indexed`](#indexed),
[Quadrature](EXAMPLES.md#quadrature), [Integrating over a cut
cell](EXAMPLES.md#integrating-over-a-cut-cell).

An extent is the length of an axis of a value shape, or the extent of a free
index. An extent is an `int` $n \ge 0$ or an `AbstractExtent`.

An abstract subclass of the base object `AbstractExtent` for an extent whose
value is not known during construction. An `AbstractExtent` is not an
expression. A concrete `AbstractExtent` is a subclass outside the core, such as
`RuntimeExtent` in the
[cut cell example](EXAMPLES.md#integrating-over-a-cut-cell).

| `AbstractExtent` |                                      |
| ---------------- | ------------------------------------ |
| Extends          | `Base`                               |
| Defined by       | A subclass outside the core          |

An `AbstractExtent` has no arithmetic and no ordering, so the rules of the
core, which compare extents by equality, accept an `AbstractExtent` unchanged.
A concrete `AbstractExtent` stores the arguments that identify the extent, such
as a label.

Example:

- With an `AbstractExtent` `q`, an expression `points` of value shape `(q, 2)`
  and an expression `weights` of value shape `(q,)` share the extent `q`, so
  `weights[i] * points[i, 0]` is accepted. With a second `AbstractExtent` `r`,
  an expression `weights` of value shape `(r,)` raises `ValueError` in
  `weights[i] * points[i, 0]`, because the index `i` has two extents.

> Open: ragged indexing. An axis whose length depends on another index, such
> as the quadrature points of each cut cell, cannot be written. A dependent
> extent, an `AbstractExtent` whose label contains the other index, would cover
> the case. An index value read from a table is a
> [`GatheredIndex`](#gatheredindex).

### `AbstractIndexValue`

Requires: [The base object `Base`](#the-base-object-base), [Extents](#extents).

See also: [`Indexed`](#indexed), [Sides of a topology](#sides-of-a-topology).

An abstract subclass of the base object `AbstractIndexValue` for an index whose
value is fixed but not known during construction, such as the local number of
the facet of a cell in a facet integral, which the solver gives for each facet.
An `AbstractIndexValue` is not an expression and not a free index.

| `AbstractIndexValue` |                                                        |
| -------------------- | ------------------------------------------------------ |
| Extends              | `Base` |
| Provides             | `extent`, an extent, and `shape`, a `tuple` of extents |
| Defined by           | A subclass outside the core                            |

An `AbstractIndexValue` is a table of values, each an `int` $v$ with
$0 \le v <$ `extent`, such as the cell on one side of each facet, which the
mesh gives. An `AbstractIndexValue` with the shape `()` is one value. A table
is read only through a [`GatheredIndex`](#gatheredindex), so one value is read
as `gathered_index(v, ())`.

Example:

- With an `AbstractIndexValue` `l` of the shape `()` and the extent 3, and the
  table `starts` of value shape `(3, 2)` of the
  [boundary example](EXAMPLES.md#integrating-over-the-boundary-of-a-triangle),
  `starts[gathered_index(l, ()), k]` is the start of edge `l`, with the free
  index `k`. With an `AbstractIndexValue` with the extent 4, the same
  expression raises `IndexError`.

### `GatheredIndex`

Requires: [`AbstractIndexValue`](#abstractindexvalue), [`Index`](#index),
[Extents](#extents).

See also: [`Indexed`](#indexed), [Restriction](#restriction).

A concrete subclass of the base object for an index value read from a table.
A `GatheredIndex` is not an expression and not a free index.

| `GatheredIndex(table, indices)` |                                                                                                                                    |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| Extends                         | `Base` |
| Free function                   | `gathered_index(table, indices)`                                                                                                   |
| Takes                           | `table: AbstractIndexValue`, `indices: tuple` |
| Properties                      | `table`, an `AbstractIndexValue`, and `indices`, an `int`, an `Index` or a `GatheredIndex` for each axis of the shape of the table |
| Derived                         | `extent` and `free_indices`, below |
| Raises                          | The errors below                                                                                                                   |

The value of a `GatheredIndex` is the entry of the table at the indices.

`extent` and `free_indices`:

- The extent is the extent of the table.
- The free indices are the free indices of the `Index` entries, each with the
  extent of the axis, and the free indices of the `GatheredIndex` entries, in
  order. Each `Index` is taken once.

Raises:

- `TypeError` for an entry that is neither exactly `int`, `Index` nor
  `GatheredIndex`.
- `IndexError` if the number of entries differs from the number of axes.
- `IndexError` for an `int` entry out of range.
- `IndexError` for a `GatheredIndex` entry whose extent differs from the
  extent of the axis.

#### Where a `GatheredIndex` appears

- A `GatheredIndex` is an entry of [`Indexed`](#indexed), and covers one axis
  whose extent equals the extent of the `GatheredIndex`. The free indices of
  the `GatheredIndex` are free in the `Indexed`.
- A `GatheredIndex` is also an entry of a `GatheredIndex`, with the same rule,
  so gathers nest. For example, the vertices of the cell on the side of each
  facet are read through the cell on the side of the facet.

#### Examples

The examples use the following objects:

- A table `cells` of the shape `(4,)`, the extent 2 and the values
  `(0, 1, 1, 0)`, such as the cell of each of the four boundary edges of the
  square divided into two triangles.
- An array `a` of value shape `(2,)`, with one value for each triangle.

```python
e = index("e")
g = gathered_index(cells, (e,))

# The value of a on the triangle of each edge, with the free index e of extent 4.
a_e = a[g]

# The array (a[0], a[1], a[1], a[0]) of value shape (4,).
b = component_tensor(a_e, (e,))

# The sum over the edges, 2 * (a[0] + a[1]) after evaluation.
total = index_sum(a_e, (e,))
```

- Let `cell_vertices` be a table of the vertices of each triangle, of the shape
  `(2, 3)` and the extent 4. `gathered_index(cell_vertices, (g, w))` is the
  vertex `w` of the triangle of edge `e`, with the free indices `e` and `w`.
- For an array `c` of value shape `(3,)`, `c[g]` raises `IndexError`, because
  the extent 2 of `g` differs from the extent 3 of the axis.

### Free indices

Requires: [`AbstractExpression`](#abstractexpression), [`Index`](#index),
[`EllipsisIndex`](#ellipsisindex), [Extents](#extents).

See also: [`Sum`](#sum), [`Tensor`](#tensor), [`Indexed`](#indexed), [`ComponentTensor`
and `IndexSum`](#componenttensor-and-indexsum).

Every `AbstractExpression` derives the free indices, a `tuple` of pairs of an
`Index` and an extent, or of an `EllipsisIndex` and a `tuple` of extents. The
order is the order of first occurrence, never the order of a `set`.

Unless a node class states otherwise, the free indices of a node are the free
indices of the edges of the node that are expressions, in edge order, with each
`Index` taken once at the first occurrence. An `Index` that occurs in two edges
must have the same extent in both edges, or construction raises `ValueError`.

Summation is never implicit. `a[i] * b[i]` has the free index `i` and denotes
the component-wise product. The sum over `i` is `index_sum(a[i] * b[i], (i,))`.
The operands of a node may have different free indices, so `a[i] + b[j]` has
the free indices `i` and `j` and denotes the outer sum.

### Scalar type

Requires: [`AbstractExpression`](#abstractexpression).

See also: [Literals](#literals), [`Absolute`](#absolute).

The scalar type is real or complex. Unless the [static type
table](#appendix-static-type-table) states otherwise, an expression is complex
if and only if at least one edge is a complex expression.

The scalar type is the arithmetic in which the expression is evaluated, not a
statement about the value. A real `SquareRoot` or a real `Power` whose value is
not a real number, such as `square_root(Integer(-1))`, is undefined.

### Support

Requires: [`AbstractExpression`](#abstractexpression),
[`AbstractDomain`](#abstractdomain).

The support is an [`AbstractDomain`](#abstractdomain). The support is absent
for an expression with no support, such as a literal. Unless the [static type
table](#appendix-static-type-table) states otherwise, the support
of a node is the support of the edges of the node that are expressions with a support.
All such edges must have the same support, or construction raises `ValueError`. A
node with no such edge has no support.

### Degree

Requires: [`AbstractExpression`](#abstractexpression).

See also: [Appendix: Static type table](#appendix-static-type-table), [New
subclasses](#new-subclasses), [`ExteriorDerivative`](#exteriorderivative).

The degree is a `tuple` with the form degree of each slot of the expression,
and every entry is at least 1. An expression with the degree `()` has no
slots and is a 0-form, such as a scalar or an array. The axes of the value
shape are not slots. The degree has no default. Every concrete class states
the degree of the class, and a leaf class of a downstream library states the
degree as part of the contract of the class (see
[New subclasses](#new-subclasses)).

The degree of every node class of the core is in the
[static type table](#appendix-static-type-table). If the
degree of an edge breaks the rule of the node, construction raises
`ValueError`.

Examples of the degree, on $R^3$ with the coordinate function $X$, the
coordinates $x, y, z$ and a scalar function $f$:

| Object                                                                     | Value shape | Degree   |
| -------------------------------------------------------------------------- | ----------- | -------- |
| The scalar $f$, or the literal `Integer(1)`                                | `()`        | `()`     |
| The coordinate function `X`                                                | `(3,)`      | `()`     |
| The 1-form $\mathrm{d}f$                                                   | `()`        | `(1,)`   |
| The 2-form $\mathrm{d}x \wedge \mathrm{d}y$                                | `()`        | `(2,)`   |
| The volume form $\mathrm{d}x \wedge \mathrm{d}y \wedge \mathrm{d}z$        | `()`        | `(3,)`   |
| The vector-valued 1-form $\mathrm{d}X$, with one 1-form for each component | `(3,)`      | `(1,)`   |
| The double form $\mathrm{d}x \otimes \mathrm{d}y$, a covariant 2-tensor    | `()`        | `(1, 1)` |

### Expressions as differential forms

Requires: [Support](#support), [Degree](#degree).

See also: [`Zero`](#zero).

An expression denotes an array-valued differential form, a multi-form, on the
support of the expression, with the degree of the expression. A 0-form is an
array-valued function. An expression with no support and the degree `()` denotes a
single array. An expression with argument spaces other than `()` denotes a
multilinear map from the argument spaces to such differential forms (see
[Variational forms](#variational-forms)). A leaf with a non-empty degree has a
support, except `Zero`. So an
expression with no support and a non-empty degree is built from `Zero`s, such as
`Zero((), (), (1,), ()) + Zero((), (), (1,), ())`, and denotes the zero form of that
degree on every support.

### Node class or free function

Requires: [Properties and free functions](#properties-and-free-functions),
[`AbstractExpression`](#abstractexpression), [Free indices](#free-indices),
[Scalar type](#scalar-type), [Support](#support), [Degree](#degree),
[`AbstractDomain`](#abstractdomain).

See also: [Transpose and matrix product](EXAMPLES.md#transpose-and-matrix-product), [Volume
form](#volume-form), [`ExteriorDerivative`](#exteriorderivative),
[`Division`](#division), [Elementary functions (very
incomplete)](#elementary-functions-very-incomplete),
[Tensor algebra with
`EllipsisIndex`](#tensor-algebra-with-ellipsisindex),
[Language extensions](#language-extensions).

A concept must be a free function, not a node class, if the free function can
build the concept as one fixed composition of existing node classes for each
static type of the arguments. The static type of an expression is the value
shape, the free indices, the scalar type, the support, the degree and the
argument spaces. The static type of a domain is the stored properties of the
domain, such as the dimension of a coordinate domain or the pieces of a
decomposed domain. A free function may read the static type to choose the
composition, and never reads the structure of an argument, such as the class,
the edges or the values. Otherwise the concept gets its own node class. A
concept that is evaluated by a numerical algorithm that the user supplies also
gets its own node class, even if a fixed composition exists, because the
composition would fix the algorithm. A concept that a core transformation must
meet as the user wrote the concept also gets its own node class, even if a
fixed composition exists, and the core transformation states the composition.
`Integral` over a decomposed domain is such a concept (see
[`Integral`](#integral)).

For example, the transpose is always `component_tensor(A[k, l], (l, k))`, and
`volume_form(P)` is the product of the $n$ coordinate 1-forms of a coordinate
domain of dimension $n$, so both are free functions.

A free function is part of the core only if a definition or a rule of the
core uses the free function, such as `volume_form` in [`Integral` over a
decomposed domain](#integral-over-a-decomposed-domain) and `cell_domain` in
[Restriction](#restriction). Every other free function belongs to a library,
such as `transpose` (see [Transpose and matrix
product](EXAMPLES.md#transpose-and-matrix-product)) and the operations of
[Tensor algebra with `EllipsisIndex`](#tensor-algebra-with-ellipsisindex). The
operator overloads of the core are those of the free functions of the core,
so `A + B` on tensors is not an overload (see
[Operator overloads](#operator-overloads)).

A free function that builds a fixed composition is a definition, not a
simplification. The free function does not read the structure of the
operands, so the free function never removes structure that the user built.

A node class is part of the core only if a definition or a rule of the core
needs the node class, or if the node class is an operation of the scalar
field, such as conjugation. A rule for the node class itself does not count.
Otherwise the node class belongs to a
[language extension](#language-extensions). For example, the core writes a
rational number as the `Division` of two `Integer`s, so `Division` is part of
the core.

The node classes of the core, and why each is a node class and not a
composition:

| Node class                                       | Why a node class                                                                                              |
| ------------------------------------------------ | ------------------------------------------------------------------------------------------------------------- |
| `Integer`, `Float`, `ImaginaryUnit`              | A number written down is a leaf                                                                               |
| `Zero`                                           | No composition gives a zero with free indices or a non-empty degree without edges to other expressions        |
| `Sum`, `Product`, `Power`, `Absolute`, `Conjugate` | Operations of the scalar field, which nothing else composes                                                 |
| `Division`                                       | `multiply(f, power(g, Integer(-1)))` rounds twice, and the quotient the user wrote rounds once                |
| `SquareRoot`                                     | IEEE 754 rounds the square root correctly and the power not, so `power(f, Integer(1) / Integer(2))` differs   |
| `TensorProduct`                                  | Nothing else adds a slot                                                                                      |
| `Inverse`                                        | No one composition covers every $k \le n$, and a pass chooses the algorithm, as for `Division`                |
| `Tensor`, `Indexed`, `ComponentTensor`, `IndexSum` | The primitives of index notation                                                                            |
| `Coordinate`, `Function`, `Normal`               | Terminals: the coordinate function of a domain, a member of a space, the outward normal of a cell on a side   |
| `ExteriorDerivative`                             | Nothing composes to differentiate                                                                             |
| `Components`                                     | Nothing composes to turn a slot into an axis                                                                  |
| `Integral`                                       | Nothing composes to integrate over a coordinate domain, and over a decomposed domain a pass must meet the integral as written |
| `PullBack`                                       | Nothing composes to move an expression between domains                                                        |
| `Form`, `Action`                                 | Binding arguments and giving an expression to an argument are not compositions                                |

## Algebra

### Literals

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices), [Scalar type](#scalar-type), [Support](#support),
[Degree](#degree).

See also: [`Division`](#division).

A literal is a number written down in an expression, such as
`Integer(2) * u`. Every literal has value shape
`()`, no free indices, no support and the degree `()`.

The literals are concrete subclasses of `AbstractExpression`:

| `Integer(v)` |                       |
| ------------ | --------------------- |
| Extends      | `AbstractExpression`  |
| Properties   | `value` $v$, an `int` |

| `Float(v)` |                               |
| ---------- | ----------------------------- |
| Extends    | `AbstractExpression`          |
| Properties | `value` $v$, a finite `float` |

| `ImaginaryUnit()` |                      |
| ----------------- | -------------------- |
| Extends           | `AbstractExpression` |

Each literal constructor accepts only exact types, as in rule 1 of [the base
object](#the-base-object-base). So `Integer(True)` raises `TypeError`, and
`Float(-0.0) != Float(0.0)`. `ImaginaryUnit()` denotes $i$.

A rational number is the `Division` of two `Integer`s, such as
`Integer(1) / Integer(3)`. In floating-point arithmetic, the division of two
integers rounds once, the same as the conversion of a rational number, so the
core has no rational literal. A complex number is a composition such as
`a + ImaginaryUnit() * b`, so the core has no complex literal.

### `Zero`

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices), [Support](#support), [Degree](#degree).

See also: [`Indexed`](#indexed), [Applying derivatives](#applying-derivatives).

A concrete subclass of `AbstractExpression` for the zero of any static type.

| `Zero(value_shape, free_indices, degree, argument_spaces)` |                                                                                                                                                                                    |
| ---------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                                                    | `AbstractExpression` |
| Properties                                                 | `value_shape`, `free_indices`, `degree` and `argument_spaces`, as in [Extents](#extents), [Free indices](#free-indices), [Degree](#degree) and [Argument spaces](#argument-spaces) |
| Static type                                                | The stored properties, with no support                                                                                                                                             |
| Raises                                                     | `TypeError` for an argument of another form, and `ValueError` for a negative extent, an index that occurs twice, or a slot degree less than 1                                      |

`apply_derivatives` produces `Zero` (see [Applying
derivatives](#applying-derivatives)). For example,
`Zero((), (), (1,), ())` is the zero 1-form, `Zero((3,), (), (), ())` is the
zero array of value shape `(3,)`, and `Zero((), (), (), (V,))` is the zero
linear form on a space `V`.

### Arithmetic

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices).

See also: [`ComponentTensor` and `IndexSum`](#componenttensor-and-indexsum), [`Tensor`
algebra with `EllipsisIndex`](#tensor-algebra-with-ellipsisindex).

`Sum`, `Product`, `Division` and `Power` are concrete subclasses of
`AbstractExpression` that each store two expressions. Arithmetic is scalar.
Both operands must have value shape `()`, or construction raises
`ValueError`, and the result has value shape `()`. Free indices pass through
arithmetic (see [Free indices](#free-indices)), so tensor arithmetic is
written in index notation. For example, a scalar `a` times a tensor `T` of any
rank is `component_tensor(a * T[r], (r,))`, with `r` an
[`EllipsisIndex`](#ellipsisindex).

| Class      | Free function | Overload |
| ---------- | ------------- | -------- |
| `Sum`      | `add`         | `+`      |
| `Product`  | `multiply`    | `*`      |
| `Division` | `divide`      | `/`      |
| `Power`    | `power`       | `**`     |

### `Sum`

Requires: [Arithmetic](#arithmetic), [Degree](#degree).

See also: [`ExteriorDerivative`](#exteriorderivative).

A concrete subclass of `AbstractExpression` for the sum of two differential
forms of the same degree.

| `Sum(first, second)` |                                                                 |
| -------------------- | --------------------------------------------------------------- |
| Extends              | `AbstractExpression` |
| Free function        | `add(first, second)`, with the overload `+`                     |
| Takes                | `first, second: AbstractExpression` |
| Properties           | `first` and `second` |
| Requires             | The operands of [Arithmetic](#arithmetic), with the same degree |
| Raises               | `ValueError` if a requirement fails                             |

For example, $\mathrm{d}x + \mathrm{d}y$ has the degree `(1,)`, and
$f + \mathrm{d}x$ raises `ValueError`.

### `Product`

Requires: [Arithmetic](#arithmetic), [Free indices](#free-indices),
[Support](#support), [Degree](#degree),
[`AbstractCoordinateDomain`](#abstractcoordinatedomain).

See also: [`TensorProduct`](#tensorproduct), [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum),
[`ExteriorDerivative`](#exteriorderivative).

A concrete subclass of `AbstractExpression` for the wedge product
$\alpha \wedge \beta$ of two differential forms.

| `Product(first, second)` |                                                                                                                |
| ------------------------ | -------------------------------------------------------------------------------------------------------------- |
| Extends                  | `AbstractExpression` |
| Free function            | `multiply(first, second)`, with the overload `*`                                                               |
| Takes                    | `first, second: AbstractExpression` |
| Properties               | `first` and `second` |
| Requires                 | The operands of [Arithmetic](#arithmetic). If both operands have a non-empty degree, each has exactly one slot |
| Static type              | The degrees added slot by slot, where the degree `()` adds 0 to every slot                                     |
| Raises                   | `ValueError` if a requirement fails                                                                            |

For two 0-forms, the wedge product is the ordinary product, and
$f \wedge \omega = f \omega$. With
two slots in an operand, the interior product of [Expanding
components](#expanding-components) removes a slot from one operand only, and
the slots of the operands no longer align. A slot of the result may have a
degree greater than the dimension of the support, and then the result is the
zero form of that degree.

The wedge product is graded commutative,
$\alpha \wedge \beta = (-1)^{kl} \beta \wedge \alpha$ for a $k$-form $\alpha$
and an $l$-form $\beta$. Equality is structural, so `Product` never reorders
the operands, and a later pass may apply the identity. `Product` never adds a
slot, so the product of two 1-forms is a 2-form with the degree `(2,)`, not a
double form with the degree `(1, 1)`.

Examples, with a scalar $f$, the 1-forms $\mathrm{d}x$ and $\mathrm{d}y$ and a
double form $S$ with the degree `(1, 1)`:

| `Product`                                         | Degree              |
| ----------------------------------------------- | ------------------- |
| $f \cdot f$                                     | `()`                |
| $f \, \mathrm{d}x$                              | `(1,)`              |
| $\mathrm{d}x \wedge \mathrm{d}y$                | `(2,)`              |
| $\mathrm{d}x \wedge \mathrm{d}x$, the zero form | `(2,)`              |
| $f \, S$                                        | `(1, 1)`            |
| $\mathrm{d}x \wedge S$                          | raises `ValueError` |

### `Power`

Requires: [Arithmetic](#arithmetic), [Support](#support).

See also: [Applying derivatives](#applying-derivatives).

A concrete subclass of `AbstractExpression` for a scalar raised to a scalar
power.

| `Power(first, second)` |                                                                             |
| ---------------------- | --------------------------------------------------------------------------- |
| Extends                | `AbstractExpression` |
| Free function          | `power(first, second)`, with the overload `**`                              |
| Takes                  | `first, second: AbstractExpression` |
| Properties             | `first` and `second` |
| Requires               | The operands of [Arithmetic](#arithmetic), both of the degree `()`                |
| Raises                 | `ValueError` if a requirement fails                                         |

The exponent may have a support. The derivative of a power whose exponent has
a support needs a logarithm, which the core does not have, so the exterior
derivative of such a power remains after
[`apply_derivatives`](#applying-derivatives), as on a function.

### Negation and subtraction

Requires: [`Sum`](#sum), [`Product`](#product), [Literals](#literals).

Negation and subtraction are free functions, not node classes:

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `negate(f)` | `-f` | `f: AbstractExpression` | `multiply(Integer(-1), f)` |
| `subtract(f, g)` | `f - g` | `f, g: AbstractExpression` | `add(f, negate(g))` |

Both compositions are exact in floating-point arithmetic, so neither
composition changes the arithmetic that the user wrote.

### `Division`

Requires: [Arithmetic](#arithmetic), [`Product`](#product), [`Power`](#power),
[Literals](#literals).

A concrete subclass of `AbstractExpression` for the quotient of a
differential form by a scalar.

| `Division(first, second)` |                                                                                        |
| ------------------------- | -------------------------------------------------------------------------------------- |
| Extends                   | `AbstractExpression` |
| Free function             | `divide(first, second)`, with the overload `/`                                         |
| Takes                     | `first, second: AbstractExpression` |
| Properties                | `first` and `second` |
| Requires                  | The operands of [Arithmetic](#arithmetic), with the degree `()` for the second operand |
| Raises                    | `ValueError` if a requirement fails                                                    |

The quotient has the degree of the first operand, so `dx / Integer(2)` is a
1-form.

### `Absolute`

Requires: [`AbstractExpression`](#abstractexpression), [Scalar
type](#scalar-type).

A concrete subclass of `AbstractExpression` for the absolute value, or the
modulus of a complex scalar.

| `Absolute(f)` |                                           |
| ------------- | ----------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `absolute(f)`, with the overload `abs(f)` |
| Takes         | `f: AbstractExpression` |
| Properties    | `operand` $f$ |
| Requires      | A scalar $f$ of the degree `()`           |
| Static type   | Value shape `()`, real                    |
| Raises        | `ValueError` if a requirement fails       |

### `Conjugate`

Requires: [`AbstractExpression`](#abstractexpression), [Scalar
type](#scalar-type), [Degree](#degree).

See also: [Variational forms](#variational-forms).

A concrete subclass of `AbstractExpression` for the complex conjugate, acting
on the coefficients of a differential form.

| `Conjugate(f)` |                                                                                             |
| -------------- | ------------------------------------------------------------------------------------------- |
| Extends        | `AbstractExpression` |
| Free function  | `conjugate(f)`                                                                              |
| Takes          | `f: AbstractExpression` |
| Properties     | `operand` $f$ |
| Requires       | A scalar differential form $f$ of any degree                                                |
| Static type    | The argument spaces of $f$ with each space conjugated (see [`ConjugateSpace`](#conjugatespace)) |
| Raises         | `ValueError` if a requirement fails                                                         |

A sesquilinear variational form, such as the $L^2$ inner product
$\int u \bar{v}$, needs `Conjugate`.

### `TensorProduct`

Requires: [Arithmetic](#arithmetic), [`Product`](#product), [Free
indices](#free-indices), [Degree](#degree).

See also: [`ExteriorDerivative`](#exteriorderivative),
[`DomainMap`](#domainmap).

A concrete subclass of `AbstractExpression` for the tensor product
$\alpha \otimes \beta$ of two differential forms, which adds a slot.

| `TensorProduct(first, second)` |                                                                                 |
| ------------------------------ | ------------------------------------------------------------------------------- |
| Extends                        | `AbstractExpression` |
| Free function                  | `tensor_product(first, second)`, with no overload                               |
| Takes                          | `first, second: AbstractExpression` |
| Properties                     | `first` and `second` |
| Requires                       | Two operands with value shape `()`, as in [Arithmetic](#arithmetic)             |
| Static type                    | The degree `first.degree + second.degree`, and the free indices as in `Product` |
| Raises                         | `ValueError` if a requirement fails                                             |

An operand with the degree `()` adds no slot, so
$f \otimes \omega = f \omega$.

The tensor product is not commutative, and `TensorProduct` never reorders the
operands. `Product` adds degrees within a slot, and `TensorProduct` adds
slots.

Examples, with the setup of [`ExteriorDerivative`](#exteriorderivative) on
`RD(3)` and `i = index("i")`:

| Expression                                        | Degree              | Meaning                                                    |
| ------------------------------------------------- | ------------------- | ---------------------------------------------------------- |
| `tensor_product(dx, dy)`                          | `(1, 1)`            | $\mathrm{d}x \otimes \mathrm{d}y$                          |
| `index_sum(tensor_product(dX[i], dX[i]), (i,))`   | `(1, 1)`            | The Euclidean metric                                       |
| `tensor_product(dx * dy, dz)`                     | `(2, 1)`            | $(\mathrm{d}x \wedge \mathrm{d}y) \otimes \mathrm{d}z$     |
| `tensor_product(dx, dy) * tensor_product(dy, dz)` | raises `ValueError` | Each operand of `Product` with a slot has exactly one slot |
| `tensor_product(dX, dX)`                          | raises `ValueError` | `dX` has value shape `(3,)`                                |

For a domain map $D$ with the parameter domain $P$ and
`dpsi = exterior_derivative(D.parameterisation, P)`,
`index_sum(tensor_product(dpsi[i], dpsi[i]), (i,))` is the metric that $\psi$
induces on $P$, the pull-back of the Euclidean metric.

### Elementary functions (very incomplete)

Requires: [`AbstractExpression`](#abstractexpression), [`Power`](#power),
[`Division`](#division), [Literals](#literals).

`SquareRoot` is the only elementary function so far, a concrete subclass of
`AbstractExpression` for the square root of a scalar.

| `SquareRoot(f)` |                                     |
| --------------- | ----------------------------------- |
| Extends         | `AbstractExpression` |
| Free function   | `square_root(f)`                    |
| Takes           | `f: AbstractExpression` |
| Properties      | `operand` $f$ |
| Requires        | A scalar $f$ with the degree `()`   |
| Static type     | Value shape `()`                    |
| Raises          | `ValueError` if a requirement fails |

`SquareRoot` is not the composition `power(f, Integer(1) / Integer(2))`,
because IEEE 754 requires a correctly rounded square root and not a correctly
rounded power, so the result of `power` depends on the math library.

### `Tensor`

Requires: [`AbstractExpression`](#abstractexpression), [Degree](#degree).

See also: [`Coordinate`](#coordinate), [`ExteriorDerivative`](#exteriorderivative).

A concrete subclass of `AbstractExpression` that stacks expressions along a
new first axis.

| `Tensor(entries)` |                                                                                                           |
| ----------------- | --------------------------------------------------------------------------------------------------------- |
| Extends           | `AbstractExpression` |
| Free function     | `tensor(entries)`                                                                                         |
| Takes             | `entries: tuple[AbstractExpression, ...]` |
| Properties        | `entries`, a non-empty `tuple` of `AbstractExpression`s with the same value shape `s` and the same degree |
| Static type       | Value shape `(n,) + s`, with $n$ the number of entries                                                    |
| Raises            | `ValueError` if a requirement fails                                                                       |

A matrix is a `Tensor` of `Tensor`s. A `Tensor` of differential forms is an
array-valued differential form, so the new axis is an axis of the value shape
and never a slot.

A value shape may contain the extent 0, for example the value shape `(0,)` of
the coordinate function of `RD(0)`. `Tensor` cannot build an expression with
the extent 0, because an empty `tuple` has no entry to give `s`. An array with
an extent 0 has no entries, so the array is the `Zero` of that value shape.

Examples, with `zero, one, two = Integer(0), Integer(1), Integer(2)`:

- With the 1-forms $\mathrm{d}x$ and $\mathrm{d}y$ on `RD(2)`,
  `tensor((dx, dy))` has value shape `(2,)` and the degree `(1,)`, and is the
  vector-valued 1-form $\mathrm{d}X$.
- With `r0 = Float(0.0)` and `r1 = Float(1.0)`,
  `tensor((tensor((r1, r0)), tensor((r0, r1))))` has `Float` entries, so the
  matrix is not equal to `tensor((tensor((one, zero)), tensor((zero, one))))`.
- `tensor((one, tensor((one, two))))` raises `ValueError`, because the value
  shapes of the entries differ.
- `tensor((one, dx))` raises `ValueError`, because the degrees of the entries
  differ.

### `Indexed`

Requires: [`AbstractExpression`](#abstractexpression), [`Index`](#index),
[`EllipsisIndex`](#ellipsisindex), [Extents](#extents),
[`AbstractIndexValue`](#abstractindexvalue), [Free indices](#free-indices),
[Support](#support), [Degree](#degree).

See also: [`Coordinate`](#coordinate), [`ExteriorDerivative`](#exteriorderivative).

> LLM generated content.

A concrete subclass of `AbstractExpression` that selects one component of a
tensor-valued expression.

| `Indexed(expression, indices)` |                                                                                                                                                                                                                                                 |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                        | `AbstractExpression` |
| Free function                  | `indexed(expression, indices)`, with the overload `expression[indices]`                                                                                                                                                                         |
| Takes                          | `expression: AbstractExpression`, `indices: tuple` |
| Properties                     | `expression`, of rank $r$, and `indices`, a `tuple` of `int`, [`Index`](#index), [`EllipsisIndex`](#ellipsisindex) and [`GatheredIndex`](#gatheredindex) entries, with at most one `EllipsisIndex`                                               |
| Static type                    | Value shape `()`, and the free indices of the expression followed by the `Index` and `EllipsisIndex` entries and the free indices of the `GatheredIndex` entries, with each `Index` taken once                                                  |
| Raises                         | See the table below                                                                                                                                                                                                                             |

Indexing is full indexing only. Each `int`, each `Index` and each
`GatheredIndex` entry covers one axis, and an
`EllipsisIndex` entry covers every axis that the other entries do not cover.
Without an `EllipsisIndex`, the number of entries must equal $r$. With an
`EllipsisIndex`, the number of other entries must be at most $r$. Each `int`
entry $i_j$ must index an axis with an `int` extent $n_j$ and must satisfy
$0 \le i_j < n_j$, where `(n_0, ..., n_{r-1})` is the value shape of the
expression. Each `GatheredIndex` entry must have the extent $n_j$. Each `Index` entry $i_j$ takes the extent $n_j$. An
`EllipsisIndex` takes the `tuple` of the extents of the axes that the
`EllipsisIndex` covers.

An `Index` may occur in two or more entries, or may already be free in the
expression. All the occurrences are the same index, so the indexed axes take
the same value. `A[i, i]` is the diagonal entry of a matrix `A`, and
`index_sum(A[i, i], (i,))` is the trace. All the occurrences must have the same
extent. An `EllipsisIndex` may also already be free in the expression, with the
same rule and the `tuple` of extents of the `EllipsisIndex`. So with an
`EllipsisIndex` `r` that is free with the extents `(2, 3)` in an expression
`B` of value shape `(2, 3)`, `B[r]` is the block diagonal of `B`.

The rank $r$ counts the axes of the value shape only. The indices select axes
of the value shape and never a slot, so the slots of a differential form
cannot be indexed.

`X[0]` and `A[0, 1]` are `indexed(X, (0,))` and `indexed(A, (0, 1))`.
Construction of `Indexed` raises:

| Error        | Cause                                                                                                       |
| ------------ | ----------------------------------------------------------------------------------------------------------- |
| `TypeError`  | An entry that is neither exactly `int`, `Index`, `EllipsisIndex` nor `GatheredIndex`.                       |
| `ValueError` | No entries.                                                                                                 |
| `IndexError` | A number of entries that does not cover exactly the rank. An expression with value shape `()` has rank 0.   |
| `ValueError` | More than one `EllipsisIndex` entry.                                                                        |
| `IndexError` | A negative `int` entry, or an `int` entry $i_j \ge n_j$.                                                    |
| `IndexError` | An `int` entry for an axis whose extent is not an `int`.                                                    |
| `IndexError` | A `GatheredIndex` entry whose extent differs from the extent of the axis.                                   |
| `ValueError` | An index with two extents, in two entries or in an entry and the expression.                                |

Examples:

- For a matrix `A` of value shape `(2, 2)`, `A[1, 0]` has value shape `()`,
  and `A[1]` raises `IndexError`.
- `tensor((Integer(1) - X[0], X[0]))` has value shape `(2,)` and the support
  `RD(2)`.
- For the vector-valued 1-form `dX = tensor((dx, dy))` on `RD(2)`, `dX[0]` is
  the 1-form $\mathrm{d}x$, with value shape `()` and the degree `(1,)`.

### `ComponentTensor` and `IndexSum`

Requires: [`Index`](#index), [`EllipsisIndex`](#ellipsisindex), [Free
indices](#free-indices), [Degree](#degree), [`Indexed`](#indexed).

See also: [`ExteriorDerivative`](#exteriorderivative).

`ComponentTensor` and `IndexSum` are concrete subclasses of
`AbstractExpression` that bind free indices.
`ComponentTensor(expression, indices)` turns the indices into tensor axes.
`IndexSum(expression, indices)` sums over the indices.

| `ComponentTensor(expression, indices)`, `IndexSum(expression, indices)` |                                                                                                                                                                                                                   |
| ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                                                                 | `AbstractExpression` |
| Free function                                                           | `component_tensor(expression, indices)`, `index_sum(expression, indices)`                                                                                                                                         |
| Takes                                                                   | `expression: AbstractExpression`, `indices: tuple[Index, ...]` |
| Properties                                                              | `expression`, a scalar, and `indices`, a non-empty `tuple` of distinct `Index` and `EllipsisIndex` entries free in the expression                                                                                 |
| Static type                                                             | The free indices of the expression without the indices. `ComponentTensor` has the value shape of the extents of the indices in order, with every extent of an `EllipsisIndex`, and `IndexSum` the value shape `()` |
| Raises                                                                  | `ValueError` if a requirement fails                                                                                                                                                                               |

The indices bind axes of the value shape, so a `ComponentTensor` of a
differential form is an array-valued differential form, and an `IndexSum` of
differential forms is a sum of differential forms of the same degree.

Scope is lexical. An occurrence of an `Index` refers to the innermost
`ComponentTensor` or `IndexSum` that encloses the occurrence and binds the
index, and the occurrence is free otherwise. A bound index is not free in the
binding node, so the binding is closed, and a binder may bind an index that a
binder inside the expression of the binder also binds. For example, with
`p = component_tensor(y[gathered_index(cells, (c, w)), k], (c, w, k))`, `p`
has no free index, and `p[c, 0, k]` has the free indices `c` and `k`, which
refer to nothing inside `p`. A pass that substitutes into the expression of a
binder renames the bound indices that would capture a free index of the
substituted expression (see [Passes](#passes)).

One `IndexSum` holds a whole contraction, such as the sum over both indices of
`A[i, j] * B[i, j]`. A pass that chooses the order of summation then finds the
contraction in one node, not in a chain of nested nodes.

Examples, with `i, j = index("i"), index("j")`, two vectors `a` and `b` of value
shape `(3,)`, and a matrix `A` of value shape `(2, 2)`:

- `index_sum(a[i] * b[i], (i,))` is the inner product of `a` and `b`, with no
  free indices.
- `index_sum(A[i, j] * A[i, j], (i, j))` is the sum of the squares of the
  entries of `A`, with no free indices.
- `component_tensor(A[i, j], (j, i))` has value shape `(2, 2)` and is the
  transpose of `A`, written with user indices.
- With the vector-valued 1-form `dX` of value shape `(3,)` and the degree
  `(1,)`, `index_sum(a[i] * dX[i], (i,))` is the 1-form
  $a_0 \mathrm{d}x + a_1 \mathrm{d}y + a_2 \mathrm{d}z$, with the degree
  `(1,)`.
- `index_sum(A[i, i], (i,))` is the trace of `A`, and
  `component_tensor(A[i, i], (i,))` is the diagonal of `A`, for every extent of
  `A`, including an `AbstractExtent`.

### Tensor algebra with `EllipsisIndex`

Requires: [`EllipsisIndex`](#ellipsisindex), [`Indexed`](#indexed), [`ComponentTensor`
and `IndexSum`](#componenttensor-and-indexsum), [`Sum`](#sum), [`Product`](#product),
[`Conjugate`](#conjugate).

An `EllipsisIndex` writes an operation on tensors of every rank as one fixed
composition. The operations in the table below are not part of the core. A
library may provide the operations as free functions (see
[New operations](#new-operations)). With derived ellipsis indices `r` and `s`
and an ordinary index `k`:

| Operation          | Fixed composition                                              |
| ------------------ | -------------------------------------------------------------- |
| `scale(a, T)`      | `component_tensor(a * T[r], (r,))`                             |
| `tensor_sum(A, B)` | `component_tensor(A[r] + B[r], (r,))`                          |
| `outer(u, v)`      | `component_tensor(u[r] * v[s], (r, s))`                        |
| `dot(u, v)`        | `component_tensor(index_sum(u[r, k] * v[k, s], (k,)), (r, s))` |
| `inner(A, B)`      | `index_sum(A[r] * conjugate(B[r]), (r,))`                      |

`A[r] + B[r]` raises `ValueError` if the value shapes of `A` and `B` differ,
because the extents of `r` differ.

### `Inverse`

Requires: [`AbstractExpression`](#abstractexpression), [Degree](#degree),
[`Division`](#division), [Elementary
functions](#elementary-functions-very-incomplete), [Transpose and matrix
product](EXAMPLES.md#transpose-and-matrix-product).

See also: [Expanding pull-backs](#expanding-pull-backs), [Applying
derivatives](#applying-derivatives).

A concrete subclass of `AbstractExpression` for the inverse of a square
matrix, or the pseudo-inverse of a matrix with more rows than columns.

| `Inverse(A)`  |                                                                                         |
| ------------- | --------------------------------------------------------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `inverse(A)`                                                                            |
| Takes         | `A: AbstractExpression` |
| Properties    | `operand` $A$ |
| Requires      | $A$ with value shape `(n, k)`, $k \le n$, and the degree `()`, and $A$ real for $k < n$ |
| Static type   | Value shape `(k, n)`                                                                    |
| Raises        | `ValueError` if a requirement fails                                                     |

`Inverse(A)` is the pseudo-inverse $A^+ = (A^T A)^{-1} A^T$ of a matrix $A$
of rank $k$, the left inverse of $A$ with $A^+ A$ the identity. For $k = n$,
$A^+$ is the inverse matrix $A^{-1}$. If $A$ has rank less than $k$,
`Inverse(A)` is undefined. Construction does not check the rank of $A$.

A pass chooses how to evaluate `Inverse`, in closed form for a small matrix or
by a factorisation. The core writes the determinant as the Leibniz formula,
such as in the density of [`Integral` over a decomposed
domain](#integral-over-a-decomposed-domain), so the core has no node class for
the determinant.

Examples, with `one, two, three = Integer(1), Integer(2), Integer(3)`:

- With the Jacobian matrix
  `J = tensor((tensor((two, one)), tensor((one, three))))` of the [affinely
  mapped triangle](EXAMPLES.md#an-affinely-mapped-triangle), `inverse(J)` has
  value shape `(2, 2)`, and a later pass evaluates `inverse(J)` to
  $\frac{1}{5} \begin{pmatrix} 3 & -1 \\ -1 & 2 \end{pmatrix}$.
- With the Jacobian matrix `Je = tensor((tensor((one,)), tensor((one,))))` of
  value shape `(2, 1)`, of the edge from $(0, 0)$ to $(1, 1)$, `inverse(Je)`
  has value shape `(1, 2)`, and a later pass evaluates `inverse(Je)` to
  $\begin{pmatrix} 1/2 & 1/2 \end{pmatrix}$.
- `inverse(transpose(Je))` with value shape `(1, 2)` raises `ValueError`,
  because $k > n$.
- With the coordinate coframe `dX` of `RD(2)`, `inverse(tensor((dX, dX)))`
  raises `ValueError`, because the degree is `(1,)`.

### Operator overloads

Requires: [`AbstractExpression`](#abstractexpression), [`Sum`](#sum),
[`Product`](#product), [`Division`](#division), [`Power`](#power), [Negation and
subtraction](#negation-and-subtraction), [`Absolute`](#absolute),
[`Indexed`](#indexed), [Transpose and matrix
product](EXAMPLES.md#transpose-and-matrix-product).

See also: [Numbers in expressions](#numbers-in-expressions), [Appendix: Free
functions](#appendix-free-functions).

Every overload in the [table of free functions](#appendix-free-functions) is
a method of `AbstractExpression` (see
[`AbstractExpression`](#abstractexpression)). No other class defines
arithmetic or indexing overloads. The core
has no free function for a matrix product, so `@` is not an overload.

An operator overload calls the free function with the operands unchanged, and
the constructor validates (see [Numbers in expressions](#numbers-in-expressions)).
`AbstractExpression` defines no reflected overloads, such as `__radd__`.

`__getitem__` is the one exception. Python passes `f[i]` to `__getitem__` as the
single entry `i`, and passes `f[i, j]` as the `tuple` `(i, j)`. `__getitem__`
converts a single entry into a `tuple` of one entry and calls `indexed`. The
conversion exists because of Python indexing syntax.

`AbstractExpression` defines `__iter__` to raise `TypeError`, because
expressions are not iterable and Python otherwise iterates through
`__getitem__`.

### Numbers in expressions

Requires: [Literals](#literals), [Operator overloads](#operator-overloads),
[`Tensor`](#tensor), [`Indexed`](#indexed), [`RD`](#rd).

See also: [`Coordinate`](#coordinate).

Every number in an expression is written with literals, such as `Integer(2)`,
`Float(0.5)`, `Integer(1) / Integer(3)` or `ImaginaryUnit()`. A Python `int`, `float`,
`complex`, `Fraction` or NumPy scalar is never converted and is rejected with
`TypeError`. The type of every literal is a choice that the user makes
explicitly.

Examples, with `X = coordinate(RD(2))`:

- `X[0]**Integer(2)`, `Integer(2) * X[0]`, `Integer(1) - X[0]` and
  `tensor((Integer(1), Integer(2)))` are accepted.
- `X[0]**2`, `2 * X[0]`, `1 - X[0]`, `0.5 * X[0]`, `2 + X[0]` and
  `tensor((1, 2))` raise `TypeError`.

Worked examples for this chapter in EXAMPLES.md: [transpose and matrix product](EXAMPLES.md#transpose-and-matrix-product).

## Calculus on a domain

### `Coordinate`

Requires: [`AbstractDomain`](#abstractdomain),
[`AbstractExpression`](#abstractexpression), [Support](#support),
[Degree](#degree).

See also: [`ExteriorDerivative`](#exteriorderivative),
[`DomainMap`](#domainmap).

A concrete subclass of `AbstractExpression` for the coordinate function of a
domain, the terminal that maps each point of the domain to the coordinates of
the point in the ambient coordinate domain.

| `Coordinate(P)` |                                                    |
| --------------- | -------------------------------------------------- |
| Extends         | `AbstractExpression` |
| Free function   | `coordinate(P)`                                    |
| Takes           | `P: AbstractDomain` |
| Properties      | `domain`, an `AbstractDomain` $P$                  |
| Static type     | Value shape `(n,)`, real, the support $P$          |
| Raises          | `TypeError` if $P$ is not an `AbstractDomain`      |

$n$ is the dimension of the ambient coordinate domain of $P$. For a
coordinate domain $P$, the coordinate function is the identity,
$X : P \to R^n$, $X(p) = p$. The coordinate function is a 0-form. The 1-forms
$\mathrm{d}X_i$ come from the exterior derivative, not from `Coordinate`.

The coordinate function of a decomposed domain $D$ is defined through the
pieces of $D$. The pull-back of `coordinate(D)` by a piece is the
parameterisation of the piece (see [`DecomposedDomain`](#decomposeddomain)),
which is built from the coordinate function of the parameter domain of the
piece (see [`DomainMap`](#domainmap)).

Examples:

- `coordinate(RD(3))` is `Coordinate(RD(3))`, with value shape `(3,)` and the
  support `RD(3)`.

### Valid pair

Requires: [`AbstractDomain`](#abstractdomain), [Support](#support).

An expression and a domain $(f, P)$ are a valid pair if $P$ is an
[`AbstractDomain`](#abstractdomain) and $f$ has the support $P$ or no support.
Otherwise construction raises `TypeError` for $P$ or `ValueError` for $f$. An
$f$ with no support is a constant on $P$.

### `ExteriorDerivative`

Requires: [Valid pair](#valid-pair), [`Coordinate`](#coordinate),
[Degree](#degree), [`Product`](#product), [`Indexed`](#indexed).

A concrete subclass of `AbstractExpression` for the exterior derivative
$\mathrm{d}\omega$ of a differential form $\omega$ on a domain $P$.

| `ExteriorDerivative(f, P)` |                                                                     |
| -------------------------- | ------------------------------------------------------------------- |
| Extends                    | `AbstractExpression` |
| Free function              | `exterior_derivative(f, P)`                                         |
| Takes                      | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties                 | `expression` $f$ and `domain` $P$                                   |
| Requires                   | A [valid pair](#valid-pair) $(f, P)$, and $f$ with at most one slot |
| Static type                | The value shape of $f$                                              |
| Raises                     | `TypeError` for $P$ and `ValueError` for $f$ if a requirement fails |
| Evaluated by               | [`apply_derivatives`](#applying-derivatives)                        |

The exterior derivative acts on differential forms with at most one slot. If
$f$ has the degree `()`, the exterior derivative has the degree `(1,)`. If $f$
has the degree `(k,)`, the exterior derivative has the degree `(k + 1,)`. If
$f$ has two or more slots, construction raises `ValueError`, because the
exterior derivative of a double form needs a connection. A degree greater than
the dimension of $P$ gives the zero form, as for [`Product`](#product).

The exterior derivative acts on each component of an array-valued differential
form, so the axes of the value shape pass through unchanged.

The exterior derivative of `Coordinate(P)` is the vector-valued 1-form
$\mathrm{d}X$ and remains after [`apply_derivatives`](#applying-derivatives).

Examples on $R^3$:

```python
P = RD(3)
X = coordinate(P)
x, y, z = X[0], X[1], X[2]

def d(f):
    return exterior_derivative(f, P)

# The coordinate 1-forms, each with value shape () and the degree (1,).
dx, dy, dz = d(x), d(y), d(z)

# The coordinate coframe, with value shape (3,) and the degree (1,).
dX = d(X)
```

After `apply_derivatives`:

| Expression                                   | Degree | After a pass                                                               |
| -------------------------------------------- | ------ | -------------------------------------------------------------------------- |
| `d(x * y * z)`                               | `(1,)` | $yz \, \mathrm{d}x + xz \, \mathrm{d}y + xy \, \mathrm{d}z$, the gradient  |
| `d(-y * dx + x * dy)`                        | `(2,)` | $2 \, \mathrm{d}x \wedge \mathrm{d}y$, the curl $(0, 0, 2)$                |
| `d(x * dy * dz + y * dz * dx + z * dx * dy)` | `(3,)` | $3 \, \mathrm{d}x \wedge \mathrm{d}y \wedge \mathrm{d}z$, the divergence 3 |
| `d(d(x * y * z))`                            | `(2,)` | The zero 2-form, because the curl of a gradient is zero                    |
| `d(d(-y * dx + x * dy))`                     | `(3,)` | The zero 3-form, because the divergence of a curl is zero                  |

- `d(dx * dy * dz)` has the degree `(4,)` and is the zero form, because
  `RD(3)` has dimension 3.
- With a double form $S$ with the degree `(1, 1)`, such as a metric, `d(S)`
  raises `ValueError`.
- `exterior_derivative(X, RD(2))` raises `ValueError`, because `X` has the
  support `RD(3)`.

### `Components`

Requires: [Valid pair](#valid-pair), [`ExteriorDerivative`](#exteriorderivative),
[Degree](#degree), [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum).

See also: [`PullBack`](#pullback), [`DomainMap`](#domainmap).

A concrete subclass of `AbstractExpression` that turns the first slot of a
differential form into an axis of the value shape, with the coordinate frame of
the ambient coordinate domain of a domain $P$, of dimension $n$.

| `Components(f, P)` |                                                                                                        |
| ------------------ | ------------------------------------------------------------------------------------------------------ |
| Extends            | `AbstractExpression` |
| Free function      | `components(f, P)`                                                                                     |
| Takes              | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties         | `expression` $f$ and `domain` $P$, with $n$ the dimension of the ambient coordinate domain of $P$ |
| Requires           | A [valid pair](#valid-pair) $(f, P)$, and $f$ with at least one slot                                   |
| Static type        | Value shape `f.value_shape + (n,)`                                                                     |
| Raises             | `ValueError` if a requirement fails                                                                    |
| Evaluated by       | [`expand_components`](#expanding-components)                                                           |

#### The frame

For a domain $P$ of dimension $k$, the frame is the coordinate frame of the
ambient coordinate domain of $P$, of dimension $n$. A coordinate domain has
$k = n$.

- For $k = n$, that frame is a frame of $P$.
- For $k < n$, such as a decomposed domain that is a surface in $R^3$, the
  coordinate vector fields $\partial_j$ of the ambient coordinate domain are
  not tangent to $P$. $\partial_j$ in the definition below is then the
  orthogonal projection of $\partial_j$ onto the tangent space of $P$, with the
  Euclidean metric of the ambient coordinate domain. The array is then the
  array of the extension of $f$ that is zero on the normal space of $P$, such
  as the surface gradient of a function on $P$.

#### Definition

Let $\partial_j$ be the coordinate vector fields of $P$, dual to the coordinate
1-forms $\mathrm{d}X_j$. The entry $j$ of the new axis is the interior product
of $\partial_j$ with the first slot of $f$:

$$
\mathrm{Components}(f, P)_j(v_2, \dots, v_l) = f(\partial_j, v_2, \dots, v_l).
$$

The degree of the result follows from the degree $l$ of the first slot of $f$:

- If $l > 1$, the first slot of the result has the degree $l - 1$.
- If $l = 1$, the result does not have that slot.

For example, the degree `(1,)` gives `()`, the degree `(2,)` gives `(1,)`, and
the degree `(1, 1)` gives `(1,)`.

The reverse direction is a fixed composition,
`index_sum(v[i] * dX[i], (i,))` for a vector `v` of value shape `(n,)`.

#### Pull-backs of `Components`

The components depend on the coordinates of $P$, so a pull-back does not
commute with `Components`. `PullBack(Components(f, Q), D)` is the array of the
components of $f$ on $Q$, evaluated at $\psi(p)$. [The rule for
`Components`](#the-rule-for-components) of `expand_pull_backs` evaluates the
pull-back.

#### Examples

With `Components`, the arrays of standard differential calculus are fixed
compositions. The examples use the examples of
[`ExteriorDerivative`](#exteriorderivative) on `RD(3)`, and `f = x * y * z`:

| Expression                              | Value shape | Degree | Meaning                                                      |
| --------------------------------------- | ----------- | ------ | ------------------------------------------------------------ |
| `components(d(f), P)`                   | `(3,)`      | `()`   | The gradient $(yz, xz, xy)$, after evaluation                |
| `components(d(components(d(f), P)), P)` | `(3, 3)`    | `()`   | The Hessian of $f$                                           |
| `components(dX, P)`                     | `(3, 3)`    | `()`   | The identity matrix                                          |
| `components(dx * dy, P)`                | `(3,)`      | `(1,)` | $(\mathrm{d}y, -\mathrm{d}x, 0)$                             |
| `components(components(dx * dy, P), P)` | `(3, 3)`    | `()`   | The antisymmetric matrix of $\mathrm{d}x \wedge \mathrm{d}y$ |

- `components(x, P)` raises `ValueError`, because `x` has no slot.
- For a domain map $D$ with the parameter domain $P$,
  `components(exterior_derivative(D.parameterisation, P), P)` is the Jacobian
  matrix of value shape `(m, n)`, with $m$ the dimension of the ambient
  coordinate domain of $D$.
- Let `M` be a decomposed domain of dimension 2 in `RD(3)`, and `u` a function
  on `M`. `components(exterior_derivative(u, M), M)` is the surface gradient
  of `u`, with value shape `(3,)`.
  `components(exterior_derivative(coordinate(M), M), M)` is the orthogonal
  projection onto the tangent plane of `M`, with value shape `(3, 3)`.

### Volume form

Requires: [`Coordinate`](#coordinate), [`ExteriorDerivative`](#exteriorderivative),
[`Product`](#product), [Literals](#literals).

See also: [`Integral`](#integral), [An affinely mapped
triangle](EXAMPLES.md#an-affinely-mapped-triangle).

The free function `volume_form(P)` is the coordinate volume form of an
`AbstractCoordinateDomain` $P$ of dimension $n$,

$$
\mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}.
$$

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `volume_form(P)` | | `P: AbstractCoordinateDomain` | `Product` of the coordinate 1-forms |

With `dX = exterior_derivative(coordinate(P), P)`, `volume_form(P)` is
`dX[0] * ... * dX[n - 1]`, nested to the left as Python evaluates `*`, so
`volume_form(RD(3))` is `(dX[0] * dX[1]) * dX[2]`. The degree is `(n,)`. For $n = 0$,
`volume_form(P)` is `Integer(1)`. The volume form fixes the orientation of $P$ (see [`Integral`](#integral)). The
volume form needs no metric. A Riemannian volume form is a density times
`volume_form(P)`.

For example, with `T = UnitSimplex(2)`, `integral(volume_form(T), T)` is the
area $1/2$ of `T`.

### `Integral`

Requires: [Valid pair](#valid-pair), [Volume form](#volume-form),
[Support](#support), [Degree](#degree).

See also: [An affinely mapped
triangle](EXAMPLES.md#an-affinely-mapped-triangle).

A concrete subclass of `AbstractExpression` for the integral of a differential
form $f$ over a domain $P$ of dimension $k$. A basic integral needs no
topological notion of the domain.

| `Integral(f, P)` |                                                                                                                                                                       |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `integral(f, P)`                                                                                                                                                      |
| Takes            | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties       | `expression` $f$ and `domain` $P$, an `AbstractDomain` of dimension $k$                                                                                               |
| Requires         | A [valid pair](#valid-pair) $(f, P)$, or $f$ with the support `P.ambient_domain`, and $f$ of the degree `(k,)`, or of the degree `()` if $k = 0$ or $P$ is decomposed |
| Static type      | The value shape and the free indices of $f$, with no support                                                                                                          |
| Raises           | `TypeError` for $P$ and `ValueError` for $f$ if a requirement fails                                                                                                   |
| Evaluated by     | [`expand_integrals`](#expanding-integrals), for a decomposed domain $P$                                                                                               |

For a coordinate domain $P$ of dimension $n$, `Integral(f, P)` is the integral
of a top form $f$, with the orientation
$\mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}$ of the coordinates $X$
of $P$. With $f = f_{0 \dots n-1} \, \mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}$,

$$
\int_P f = \int_P f_{0 \dots n-1} \,\mathrm{d}X,
$$

with the Lebesgue measure $\mathrm{d}X$ on $P$. For $n = 0$, the top form is a
0-form, and the integral is the value of $f$ at the single point of $R^0$.

For a decomposed domain $P$ of dimension $k$, `Integral(f, P)` is the integral
of a $k$-form $f$ with the orientation of the pieces, or of a 0-form $f$
against the density that the pieces induce, as the sum of the integrals over
the pieces (see [`Integral` over a decomposed
domain](#integral-over-a-decomposed-domain)). For $k = 0$, a 0-form is a top
form, and the two readings agree.

An array-valued differential form integrates component by component. A later
pass might evaluate an `Integral` over a coordinate domain, for example by
quadrature.

#### Examples

With `T = UnitSimplex(2)` and `X = coordinate(T)`:

- `integral(X[0] * volume_form(T), T)` is $\int_T X_0 \,\mathrm{d}X$. A later
  pass evaluates the integral to $1/6$.
- `integral(X[0], T)` raises `ValueError`, because `X[0]` is a 0-form, not a
  2-form.
- `integral(X[0] * volume_form(T), T) * integral(X[1] * volume_form(T), T)` is a
  product of two integrals, with value shape `()` and no support.

## Maps between domains

### `DomainMap`

Requires: [`AbstractDomain`](#abstractdomain),
[`Coordinate`](#coordinate), [Support](#support), [Valid pair](#valid-pair), [Free indices](#free-indices),
[Scalar type](#scalar-type), [Degree](#degree).

See also: [`PullBack`](#pullback), [`DecomposedDomain`](#decomposeddomain),
[An affinely mapped triangle](EXAMPLES.md#an-affinely-mapped-triangle).

A concrete subclass of `Base` for a map $\psi$ from a parameter domain $P$
into an ambient domain $Q$, with the image $D = \psi(P)$. The sections below
write $D$ for the domain map and for the image of the domain map. A
`DomainMap` is not an `AbstractDomain` and is not the support of an
expression.

| `DomainMap(P, psi, Q)` |                                                                                                                                                                                                             |
| -------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                          | `Base` |
| Free function                    | `domain_map(P, psi, Q)`                                                                                                                                                                           |
| Takes                            | `P: AbstractDomain`, `psi: AbstractExpression`, `Q: AbstractDomain` |
| Properties                       | `parameter_domain` $P$, an [`AbstractDomain`](#abstractdomain) of dimension $k$, `parameterisation` $\psi$, and `ambient_domain` $Q$, an `AbstractDomain`                                                   |
| Requires                         | A [valid pair](#valid-pair) $(\psi, P)$, $\psi$ real with value shape `(n,)` for the dimension $n$ of the ambient coordinate domain of $Q$, the degree `()` and the argument spaces `()`, and $k \le n$ |
| Derived                          | `free_indices`, the free indices of $\psi$ |
| Raises                           | `TypeError` for $P$ or $Q$, and `ValueError` if another requirement fails                                                                                                                                   |

The coordinates of the ambient coordinate domain of $Q$ name the points of
$D$. The 'topological dimension' of $D$ is $k$,
`D.parameter_domain.dimension`, and the 'geometric dimension' is $n$, the
dimension of the ambient coordinate domain of $Q$. Construction does not check
that $\psi(P)$ lies in $Q$.

$P$ is usually a coordinate domain, such as a reference cell. A decomposed
domain as $P$ makes a composite cell, such as a macro element whose parameter
domain is the decomposed domain of the sub-cells of a reference cell. The
pull-back of an expression on $Q$ by $D$ is then an expression on the
decomposed domain $P$, defined through the pieces of $P$ (see
[`DecomposedDomain`](#decomposeddomain)).

A $\psi$ with free indices makes $D$ a family of domains, with one domain for
each value of the free indices. For example, with an array `x` of value shape
`(m, n)` and indices `i` and `k`,
`domain_map(RD(0), component_tensor(x[i, k], (k,)), RD(n))` is a
family of $m$ points in `RD(n)`. A family of domains writes a set of points
once, such as the points of a quadrature rule.

$Q$ is usually `RD(n)`. $Q$ may be smaller, such as the unit triangle for an
edge of the unit triangle, or a decomposed domain for a cell of the decomposed
domain (see [`PullBack`](#pullback)).

Examples:

- The hypotenuse of the unit triangle is the image of a domain map with $k = 1$
  and $n = 2$. $P$ is the unit interval in $R^1$, and
  $\psi(X) = (1 - X_0, X_0)$. $Q$ is `RD(2)` for the hypotenuse in the plane,
  and `UnitSimplex(2)` for the hypotenuse in the unit triangle.
- With `X = coordinate(RD(2))` and `h = X[0]**Integer(2) + X[1]**Integer(2)`,
  `domain_map(RD(2), tensor((X[0], X[1], h)), RD(3))` is the
  paraboloid, with $k = 2$ and $n = 3$.

### `PullBack`

Requires: [`DomainMap`](#domainmap), [Valid pair](#valid-pair),
[Support](#support),
[Free indices](#free-indices), [Degree](#degree),
[`ExteriorDerivative`](#exteriorderivative), [An
affinely mapped triangle](EXAMPLES.md#an-affinely-mapped-triangle).

A concrete subclass of `AbstractExpression` for the pull-back $\psi^* f$ of a
differential form $f$ on $Q$ by the parameterisation $\psi$ of a domain $D$,
with $P$ the parameter domain and $Q$ the ambient domain of $D$.

| `PullBack(f, D)` |                                                                                                            |
| ---------------- | ---------------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `pull_back(f, D)`                                                                                          |
| Takes            | `f: AbstractExpression`, `D: DomainMap \| SideMap` |
| Properties       | `expression` $f$, and `domain` $D$, a `DomainMap` or a [`SideMap`](#sidemap)                 |
| Requires         | A [valid pair](#valid-pair) $(f, Q)$, and an equal extent for an `Index` free in both $f$ and $D$          |
| Static type      | The value shape of $f$, and the free indices of $f$ followed by those of $D$, with each `Index` taken once |
| Raises           | `TypeError` for $D$, and `ValueError` if another requirement fails                                         |
| Evaluated by     | [`expand_pull_backs`](#expanding-pull-backs)                                                               |

The second operand is a domain, not an expression, because a domain states
$Q$.

For a 0-form $f$, the pull-back is the composition $f \circ \psi$, the
expression on $P$ with the value $f(\psi(p))$ at each point $p$ of $P$. For a
differential form with slots, the pull-back acts on each slot,

$$
(\psi^* f)_p(v_1, \dots, v_k) = f_{\psi(p)}(\mathrm{d}\psi \, v_1, \dots, \mathrm{d}\psi \, v_k),
$$

so `PullBack` accepts any degree, including a double form such as a metric.
The axes of the value shape are not transformed. An axis of an expression on
$Q$ is an array in the coordinate frame of $Q$, so the pull-back keeps the
frame of $Q$. A slot with a degree greater
than the dimension of $P$ gives the zero form.

A `PullBack` of a downstream terminal or of a
[`Function`](#function) remains, and the library of the terminal evaluates the
terminal at $\psi(p)$ and transforms the slots. For a function, the library
evaluates the pull-back by a cell domain or a gathered cell domain of the
domain of the space (see [`AbstractFunctionSpace`](#abstractfunctionspace)), and
the pull-back by another domain map into a decomposed domain is rewritten
through those (see [Points and point
evaluation](#points-and-point-evaluation)).

A domain $E$ in the parameter domain of $D$ has the image $\psi(E)$, a subset
of $D$. The image is given by the domain map with the parameterisation
`pull_back(D.parameterisation, E)` on the parameter domain of $E$, in the
ambient domain of $D$, such as `E` in the example below.

Example, with `T`, `X`, `psi`, `one` and `K` of the affinely mapped triangle:

```python
# The unit interval, the parameter domain of the edge, with coordinate t.
S = UnitSimplex(1)
t = coordinate(S)[0]

# The edge of T from (1, 0) to (0, 1), as a domain in T.
e = domain_map(S, tensor((one - t, t)), T)

# PullBack(f, e), with the support UnitSimplex(1), value shape () and the degree ().
f = X[0] * X[1]
fe = pull_back(f, e)

# PullBack(dx, e), with the support UnitSimplex(1) and the degree (1,).
dx = exterior_derivative(X[0], T)
dxe = pull_back(dx, e)

# The edge of K from (3, 1) to (2, 3), as a domain in RD(2).
E = domain_map(S, pull_back(psi, e), RD(2))
```

- A later pass evaluates `fe` to $(1 - t) t$, `dxe` to
  $\mathrm{d}(1 - t) = -\mathrm{d}t$, and the parameterisation of `E` to
  $(3 - t, 1 + 2 t)$.
- `pull_back(dx * exterior_derivative(X[1], T), e)` has the degree `(2,)` and
  is the zero form, because the edge has dimension 1.
- `pull_back(f, domain_map(S, tensor((one - t, t)), RD(2)))`
  raises `ValueError`, because that domain lies in `RD(2)`, not in `T`.

### Concept: orientation of a domain

Requires: [`DomainMap`](#domainmap),
[Volume form](#volume-form), [`PullBack`](#pullback).

An orientation is a choice of direction, such as which way an edge runs or
which side of a surface is up. A coordinate domain $P$ is oriented by the order of
the coordinates of $P$, as fixed by the volume form of $P$ (see
[Volume form](#volume-form)). A domain $D$ = `DomainMap(P, psi, Q)`
takes the orientation of $P$ through $\psi$. For example, an edge runs in the
direction in which $\psi$ moves as the parameter increases. For a map $r$ of
$P$ onto $P$, $\psi$ and $\psi \circ r$ parameterise the same set. The two
orientations are equal if the Jacobian determinant of $r$ is positive, and
opposite if the Jacobian determinant of $r$ is negative. For example,
$r(t) = 1 - t$ reverses an edge, and $r(X) = (X_1, X_0)$ reverses a triangle.
The orientation needs no metric.

If $k = n$, $D$ also has the orientation of $Q$. With $J$ the Jacobian matrix
of $\psi$,

$$
\psi^* (\mathrm{d}Y_0 \wedge \dots \wedge \mathrm{d}Y_{n-1})
= \det J \, \mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1},
$$

with $Y$ the coordinates of $Q$ and $X$ the coordinates of $P$. The two
orientations agree where $\det J > 0$, and differ where $\det J < 0$. The core
has no node class for the determinant, because $\det J$ is the coefficient of
`pull_back(volume_form(Q), D)`.

If $k < n$, such as for a surface in $R^3$, $Q$ gives $D$ no orientation. Each
piece of a decomposed domain of dimension $k < n$ takes the orientation of the
parameter domain of the piece, and the core does not relate the orientations
of two cells. An integral of a $k$-form over such a domain, such as in Stokes'
theorem on a surface, needs the cells to be oriented consistently, as does a
term that turns a $k$-form into a 0-form or a $(k - 1)$-form into a vector
field, such as the scalar curl of a 1-form on a surface. The mesh library
provides a consistent orientation where one is needed, such as an expression
with the sign $+1$ or $-1$ for each cell, and the core needs nothing more (see
[The flux through a surface](EXAMPLES.md#the-flux-through-a-surface)).

An integral of a 0-form uses a density (see [`Integral` over a decomposed
domain](#integral-over-a-decomposed-domain)), and a flux through a facet uses a
side and [`normal(D, E, s)`](#normal), which points out of the cell on the
side. Neither needs an orientation.

### Concept: integral over a domain

Requires: [`PullBack`](#pullback), [`Integral`](#integral),
[`Absolute`](#absolute), [`Division`](#division),
[Concept: orientation of a domain](#concept-orientation-of-a-domain).

The integral of a $k$-form $\omega$ on $Q$ over a domain $D$ =
`DomainMap(P, psi, Q)` of dimension $k$, with the orientation of $D$,
is

$$
\int_D \omega = \int_P \psi^* \omega,
$$

which is `integral(pull_back(omega, D), P)`. Reversing the orientation of $D$
changes the sign of the integral.

An integral that must not depend on the orientation, such as the legacy UFL
integral `f * dx` over a cell whose orientation is arbitrary, multiplies by
the sign $s = \det J / |\det J|$ for $k = n$. Then
$\int_P s \, \psi^* (f \, \mathrm{d}Y_0 \wedge \dots \wedge \mathrm{d}Y_{n-1})
= \int_P (f \circ \psi) \, |\det J| \, \mathrm{d}X$.

### Points and point evaluation

Requires: [`RD`](#rd), [`DomainMap`](#domainmap),
[`PullBack`](#pullback), [`Integral`](#integral), [An affinely mapped
triangle](EXAMPLES.md#an-affinely-mapped-triangle).

A point is the image of a domain map. A point in a domain $Q$ is
`domain_map(RD(0), x, Q)`, with $x$ a constant expression of value shape
`(n,)` for the dimension $n$ of the ambient coordinate domain of $Q$. The core
needs no node class and no free function for a point.

#### Point evaluation

The value of a 0-form $f$ of any value shape at a point $p$ is
`integral(pull_back(f, p), RD(0))` (see [`Integral`](#integral)). The integral
has no support, so the integral combines with any expression.

#### A family of points in a decomposed domain

Notation:

- $D$ is a decomposed domain with a topology $T$.
- $M$ = `domain_map(RD(0), q, D)` is a family of points in $D$, with the free
  index $i$ of the family, of the extent $m$.
- $j$ is the cell type of $T$ of the cells that contain the points, $K_j$ is
  the reference cell of cell type $j$, and $k_j$ is the dimension of $K_j$.

$M$ needs the cell of $D$ that contains each point, and the position of each
point in that cell. These are the side tables `(cells, positions)` of the
topology of the points (see [Sides of a topology](#sides-of-a-topology)). The
topology of the points has one cell of a reference cell type of dimension 0
for each point, and the side `("inside", T, j)`. A family of points in a mesh is
therefore a topology with a side.

The mesh library locates the cells and inverts the map of each cell outside
the language, for example by Newton's method, as the mesh library computes
every side table.

`pull_back(u, M)` for a function `u` on $D$ remains after `expand_pull_backs`.
The reason is that the library of the space of `u` defines `u` on the cells
only (see [`AbstractFunctionSpace`](#abstractfunctionspace)). A pass of the
mesh library rewrites the pull-back as the [restriction](#restriction) to the
topology of the points, with the index $i$ of the family in place of the cell
index:

```python
g = gathered_index(cells, (i,))
R = domain_map(RD(0), component_tensor(positions[i, 0, k], (k,)), K_j)
pull_back(pull_back(u, gathered_cell_domain(D, j, g)), R)
```

Here `cells` has the shape `(m,)` and the extent `T.num_cells[j]`, and
`positions` has the value shape `(m, 1, k_j)`. The result has the free index
$i$. The pull-back of `u` by a gathered cell domain is in the contract of the
space.

The maps below `u` may be a chain. An example of a chain is a family of maps
from a decomposed domain, pulled back by a cell of that domain and then by a
family of points of the reference cell. The pass reads the chain as one family
with all of its free indices (see [A 1D network in a 3D
mesh](EXAMPLES.md#a-1d-network-in-a-3d-mesh)).

A family whose points lie in cells of several types is not yet covered,
because the family would need one piece for each cell type.

#### Examples

The examples use `T` and `X` of the affinely mapped triangle, and
`f = X[0] * X[1]`:

```python
third = Integer(1) / Integer(3)

# DomainMap: the centroid of T, as a point in T.
c = domain_map(RD(0), tensor((third, third)), T)

# Integral: the value of f at the centroid, with no support and value shape ().
# A later pass evaluates the integral to 1/9.
fc = integral(pull_back(f, c), RD(0))
```

- Let `d = domain_map(RD(0), tensor((third, Integer(2) / Integer(3))), T)` be
  a second point. `integral(pull_back(f, c), RD(0)) + integral(pull_back(f, d), RD(0))`
  is the action of a sum of two Dirac measures on $f$.

Worked examples for this chapter in EXAMPLES.md: [an affinely mapped triangle](EXAMPLES.md#an-affinely-mapped-triangle), [orientation of an edge and of a triangle](EXAMPLES.md#orientation-of-an-edge-and-of-a-triangle), [integrating over an edge](EXAMPLES.md#integrating-over-an-edge), [integrating over the boundary of a triangle](EXAMPLES.md#integrating-over-the-boundary-of-a-triangle), [quadrature](EXAMPLES.md#quadrature), [integrating over a cut cell](EXAMPLES.md#integrating-over-a-cut-cell), [the centroid of a triangle](EXAMPLES.md#the-centroid-of-a-triangle).

## Topology

### `AbstractCellType`

Requires: [The base object `Base`](#the-base-object-base).

See also: [Reference cell](#reference-cell),
[`AbstractTopology`](#abstracttopology).

An abstract subclass of the base object `AbstractCellType` for the type of a
cell, with numbered vertices.

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

Requires: [The base object `Base`](#the-base-object-base), [`Index`](#index),
[Extents](#extents), [Free indices](#free-indices),
[`AbstractCellType`](#abstractcelltype).

See also: [Sides of a topology](#sides-of-a-topology),
[Restriction](#restriction),
[`Integral` over a decomposed domain](#integral-over-a-decomposed-domain),
[`DecomposedDomain`](#decomposeddomain).

An abstract subclass of the base object `AbstractTopology` for a set of cells
of given cell types, without coordinates.

| `AbstractTopology` |                                                                                             |
| ------------------ | ------------------------------------------------------------------------------------------- |
| Extends            | `Base` |
| Provides           | `cell_types`, `num_cells` and `sides`, and `side_tables(i, side)` for a topology with sides |
| Derived            | `dimension`, the dimension of the cell types                                                |
| Defined by         | A subclass outside the core, such as the topology of the mesh of a solver                   |

- `cell_types` is a `tuple` of `AbstractCellType`s $C_0, C_1, \dots$ of the
  same dimension, one for each group of cells. A cell of type $j$, a cell of
  group $j$, has the vertices of $C_j$, numbered as in $C_j$.
  The cell types are distinct, except that two groups with different sides may
  have the same cell type (see [Sides of a topology](#sides-of-a-topology)).
- `num_cells` is a `tuple` with the number of cells of each type, an `int`,
  such as for a macro element, or an `AbstractExtent`, such as for the mesh of
  a solver.
- `sides` is `()`, or the sides of [Sides of a topology](#sides-of-a-topology).

Two topologies are the same topology only if the two topologies are equal as
base objects, so a concrete `AbstractTopology` stores the arguments that
identify the topology, such as a label.

The free function `cell_index(T, j)` returns the derived index
`Index(("uflx.cell_index", T, j, 0))` for the cells of type $j$, numbered from
0. The index has the extent `T.num_cells[j]`. The index differs between types
and between topologies, even if the numbers of cells are equal.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_index(T, j)` | | `T: AbstractTopology`, `j: int` | `Index` |

By the rules of [free indices](#free-indices), an expression with the free
index `cell_index(T, j)` has one value for each cell of type $j$.

Examples:

- The Alfeld split of a triangle, a macro element, has the numbers of cells
  `(3,)`.
- A mesh of a solver with triangles and quadrilaterals has two cell types, and
  two `AbstractExtent`s as the numbers of cells.

#### Legacy UFL

In legacy UFL, a `Mesh` holds both the topology and the coordinate element, and the measure `dx` hides the sum over the cells. In UFLx, the topology, the parameterisations of the pieces of a decomposed domain and the `IndexSum` over `cell_index(T, j)` are separate.

### Sides of a topology

Requires: [`AbstractCellType`](#abstractcelltype),
[`AbstractTopology`](#abstracttopology).

A topology $S$ whose cells are sub-entities of cells of other topologies has
sides. A cell of $S$ is seen from one cell of another topology on each side,
and [Restriction](#restriction) picks the value of an expression on a cell of
$S$ by the side name.

#### `sides`

`S.sides` describes the sides of $S$:

- `S.sides` has one entry for each cell type of $S$.
- Entry $i$ is a `tuple` of one or more sides with distinct names.
- Each side is a triple `(name, T, j)` of a side name, a `str`, a topology
  $T$ and a cell type number $j$ of $T$.
- Two entries of the same cell type of $S$ differ. A group of cells of $S$ is
  therefore fixed by the cell type and the sides together.
- A topology whose cells are not seen from other cells has the sides `()`.

Notation, for the group $i$ of $S$ with a side `(name, T, j)`:

- $K_j$ is the reference cell of cell type $j$ of $T$, and $k_j$ is the
  dimension of $K_j$.
- $m$ is the number of vertices of cell type $i$ of $S$.

#### The image of a cell on a side

On each side, a cell of type $i$ of $S$ is the image of cell type $i$ of $S$
under the affine map through the positions of the vertices of the cell of $S$ in $K_j$ (see
[Restriction](#restriction)). The image is one of the following:

- A sub-entity of the cell on the side. The sub-entity is the cell itself, with
  the vertices in the order of the cell, a facet, or a sub-entity of lower
  dimension, such as an edge of a tetrahedron or a vertex.
- A piece of the cell on the side that is not a sub-entity, such as half of a
  facet at a hanging node or a segment of a line through the cell.

The concrete mesh gives the position of each vertex of a cell of $S$ in the
reference cell on each side:

- For a matching facet, the position is that of the same vertex.
- For a periodic boundary, the position is that of the vertex one period away.
- For a piece, the position is that of a point of the cell.

#### Side tables

A concrete topology with sides also provides `side_tables(i, side)`, a pair of
tables `(cells, positions)` for the side named `side` of group $i$, with the
side `(side, T, j)`:

- `cells` is an [`AbstractIndexValue`](#abstractindexvalue) of the shape
  `(S.num_cells[i],)` and the extent `T.num_cells[j]`. Entry $e$ is the cell of
  type $j$ of $T$ on the side of the cell $e$ of $S$.
- `positions` is a constant expression with no support and the value shape
  `(S.num_cells[i], m, k_j)`. Entry $(e, v)$ is the position in $K_j$ of vertex
  $v$ of the cell $e$. For a sub-entity, the position is the coordinates of a
  vertex of $K_j$. For a piece, the position is any point of $K_j$.

[Restriction](#restriction) reads the tables. The core does not check the
sides or the tables.

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
  mesh](EXAMPLES.md#a-1d-network-in-a-3d-mesh)).

The concrete mesh computes which parts of an entity meet which cells, outside
the language, as the concrete mesh locates the cell that contains a point (see
[Points and point evaluation](#points-and-point-evaluation)). The positions of
a piece of a curved cell are an approximation, as a quadrature rule is. The
core does not check that the image of a cell of $S$ lies in the cell on the
side.

#### Examples

The examples use a topology `T` of triangles with the cell types `(K,)`,
`I = Simplex(1)` and `P = Simplex(0)`. The examples show the topology only. A
topology of sub-entities that an expression is restricted to needs reference
cells as cell types, such as `ReferenceSimplex(1)` in place of `I` (see
[Restriction](#restriction)).

- A topology `B` of the boundary facets of `T` has the cell types `(I,)` and
  the sides `((("inside", T, 0),),)`. `cell_index(B, 0)` is the index of the
  facets.
- A topology of the facets shared by two triangles has the cell types `(I,)`
  and the sides `((("+", T, 0), ("-", T, 0)),)`.
- A topology of the corners of the boundary has the cell types `(P,)` and the
  sides `((("+", B, 0), ("-", B, 0)),)`, since each corner is seen from two
  facets of `B`.
- Let `M` be a submesh of the boundary facets of `T`, with the cell types
  `(I,)`. A topology of the facets has the sides
  `((("mesh", T, 0), ("submesh", M, 0)),)`, and each facet is a cell of `M`.
- A topology of the left boundary facets of a periodic `T` has the sides
  `((("left", T, 0), ("right", T, 0)),)`. From `"right"`, the position of each
  vertex is that of the vertex one period to the right.
- Let `M` be a topology of triangles and quadrilaterals with the cell types
  `(Kt, Kq)`. A topology of all the facets shared by two cells has three groups
  of the cell type `I`, with the cell types `(I, I, I)` and the sides
  `((("+", M, 0), ("-", M, 0)), (("+", M, 0), ("-", M, 1)), (("+", M, 1), ("-", M, 1)))`.
  The three groups are the facets between two triangles, between a triangle
  and a quadrilateral, and between two quadrilaterals. An integral over the
  topology has one term for each group.
- For a topology `U3` of prisms with the cell types `(W,)`, a topology of the
  boundary facets has the cell types `(K, Q)`, with the quadrilateral cell type
  `Q`, and the sides `((("inside", U3, 0),), (("inside", U3, 0),))`.

#### Legacy UFL

The legacy UFL measures `ds` and `dS` integrate over a topology of facets with one and two sides, and legacy UFL names the two sides `'+'` and `'-'`. The legacy UFL measure `dx(1)` integrates over a topology of the cells of a subdomain, and the legacy UFL measure `dP` integrates over a topology of vertices.

Worked examples for this chapter in EXAMPLES.md: [the simplex as a cell type](EXAMPLES.md#the-simplex-as-a-cell-type), [a static topology](EXAMPLES.md#a-static-topology).

## Reference cells

### Reference cell

Requires: [`AbstractCoordinateDomain`](#abstractcoordinatedomain), [`RD`](#rd),
[`Tensor`](#tensor), [`DomainMap`](#domainmap),
[`Inverse`](#inverse), [`AbstractCellType`](#abstractcelltype).

See also: [`AbstractTopology`](#abstracttopology), [Integrating over
the boundary of a triangle](EXAMPLES.md#integrating-over-the-boundary-of-a-triangle).

An abstract subclass of `AbstractCoordinateDomain` and of `AbstractCellType`
for a reference cell `AbstractReferenceCell`, a cell type whose vertices have
coordinates.

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `AbstractCoordinateDomain`, `AbstractCellType` |
| Provides                | `vertices`, a constant expression with no support and value shape `(m, n)`, with row $v$ the coordinates of vertex $v$, for dimension $n$ and $m$ vertices |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](EXAMPLES.md#the-unit-simplex-as-a-reference-cell)                                                   |

#### Affine maps

Notation, for `affine_map(F, points)`:

- $F$ is a reference cell of dimension $k$ with $m_F$ vertices, and $V_v$ are
  the vertices of $F$.
- `points` has value shape `(m_F, m)`, and $p_v$ are the rows of `points`.
- $E$ is the matrix with the columns $V_{d+1} - V_0$.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `affine_map(F, points)` | | `F: AbstractReferenceCell`, `points: AbstractExpression` | The affine map of the vertices |

The free function `affine_map(F, points)` returns the affine map on $F$ that
sends vertex $v$ of $F$ to row $v$ of `points` for $v \le k$:

$$
X \mapsto p_0 + \sum_{d < k} (p_{d+1} - p_0) \, \big(E^{-1} (X - V_0)\big)_d.
$$

- The map is an expression on $F$ with value shape `(m,)`, a fixed composition
  for each $k$ with [`Inverse`](#inverse).
- For $k = 0$ the map is the constant $p_0$.
- Only the first $k + 1$ rows enter the map. The other rows are the images of
  the other vertices when the image is affine, such as the facet of a
  reference cell.
- The first $k + 1$ vertices of $F$ must not lie in a hyperplane of $F$.

#### Sub-entities

Let a sub-entity of $K$ have the vertex numbers $(w_0, w_1, \dots)$ and the
reference cell $F$. The sub-entity is
`domain_map(F, affine_map(F, points), K)`. Row $v$ of `points` is the
coordinates `K.vertices[w_v, k]` of the vertex $w_v$ of $K$, stacked with
`component_tensor` over the vertex and the coordinate indices. The core does
not check that the image is a sub-entity of $K$.

#### Examples

- The unit triangle has the vertices $(0, 0)$, $(1, 0)$ and $(0, 1)$. The
  sub-entity with the vertex numbers $(1, 2)$ and the unit interval as $F$ is
  the hypotenuse, with $\psi(t) = (1 - t, t)$, as in the
  [`DomainMap`](#domainmap) examples.
- The unit interval has the vertices $(0)$ and $(1)$. The sub-entity with the
  vertex number $(1)$ and the reference point as $F$ is the end point.

Worked examples for this chapter in EXAMPLES.md: [the unit simplex as a reference cell](EXAMPLES.md#the-unit-simplex-as-a-reference-cell).

## Decomposed domains

### `DecomposedDomain`

Requires: [`AbstractTopology`](#abstracttopology),
[Reference cell](#reference-cell), [`DomainMap`](#domainmap),
[`PullBack`](#pullback), [Free indices](#free-indices).

See also: [Sides of a topology](#sides-of-a-topology),
[`AbstractFunctionSpace`](#abstractfunctionspace).

A concrete subclass of [`AbstractDomain`](#abstractdomain) for a domain that
is the union of the images of pieces, each a family of domain maps.

| `DecomposedDomain(pieces)` |                                                                                                                         |
| -------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Extends                    | `AbstractDomain` |
| Free function              | `decomposed_domain(pieces)`                                                                                             |
| Takes                      | `pieces: tuple[DomainMap, ...]` |
| Properties                 | `pieces`, a non-empty `tuple` of [`DomainMap`](#domainmap)s with the same ambient domain $Q$ and the same dimension $k$ |
| Derived                    | `ambient_domain` $Q$, `dimension` $k$, and `topology`, below                                                            |
| Raises                     | `TypeError` for a piece that is not a `DomainMap`, and `ValueError` if another requirement fails                        |

`topology` is $T$ if the free indices of piece $j$ are `(cell_index(T, j),)`
for every $j$, and otherwise `None`.

Notation, for the piece $j$ of $D$:

- $P_j$ is the parameter domain and $\psi_j$ is the parameterisation.
- For a mesh with the topology $T$, $K_j$ is the reference cell of cell type
  $j$, and $x_j$ are the coordinates of the cells of type $j$.

#### Pieces

A piece with free indices is a family, one domain for each value of the free
indices. $D$ is the union of all the domains of all the pieces. $D$ may have a
lower dimension than the ambient coordinate domain, such as a surface in $R^3$.

- A mesh has one piece for each cell type $j$, `domain_map(K_j, x_j, Q)`, with
  the free index `cell_index(T, j)`. Then the topology of $D$ is $T$.
- A piece with an ordinary `Index`, such as a cloud of points, gives no
  topology.
- A single `DomainMap` with no free indices is a decomposed domain of one
  piece.

The core does not check that the domains meet only in sets of measure zero.

#### Nesting

- The parameter domain of a piece may itself be a decomposed domain, so
  decomposed domains nest. An expression on $D$ pulled back by such a piece is
  an expression on the parameter domain of the piece, defined through the
  pieces of the parameter domain in turn.
- The ambient domain $Q$ may also be a decomposed domain, such as a mesh in
  which a family of points lies.

#### Cell domains

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_domain(D, j)` | | `D: DecomposedDomain`, `j: int` | `DomainMap` of the cells of type `j` |

The free function `cell_domain(D, j)` returns piece $j$ with the ambient domain
$D$, `domain_map(P_j, psi_j, D)`, with the free indices of piece $j$.

An expression $f$ with the support $D$ is defined through the pieces.
`pull_back(f, cell_domain(D, j))` is an expression on $P_j$, the value of $f$
on each domain of piece $j$. The pull-back of `coordinate(D)` is $\psi_j$.

#### Gathered cell domains

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `gathered_cell_domain(D, j, g)` | | `D: DecomposedDomain`, `j: int`, `g: GatheredIndex` | `DomainMap` of the cells of type `j` selected by `g` |

The free function `gathered_cell_domain(D, j, g)` applies to a decomposed
domain $D$ with a topology $T$. The function returns the domains of piece $j$
selected by a [`GatheredIndex`](#gatheredindex) $g$ of the extent
`T.num_cells[j]`. With $c$ = `cell_index(T, j)`, the one free index of piece
$j$, the result is `domain_map(P_j, chi, D)` with

```python
chi = component_tensor(component_tensor(psi_j[k], (c, k))[g, k], (k,))
```

and $k$ the derived index `Index(("uflx.gathered_cell_domain", D, j, g, 0))`.

The library of a function in a space on $D$ evaluates the pull-back by
`gathered_cell_domain(D, j, g)` as the pull-back by `cell_domain(D, j)` on the
cell $g$ (see [`AbstractFunctionSpace`](#abstractfunctionspace)).

#### Examples

The examples use the following objects:

- `K = ReferenceSimplex(2)`.
- The topology
  `square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())` of the
  [static topology example](EXAMPLES.md#a-static-topology), on
  reference cells.
- An array `y` of value shape `(4, 2)` with the positions $(0, 0)$, $(1, 0)$,
  $(0, 1)$ and $(1, 1)$ of the vertices.
- A table `cell_vertices` of the shape `(2, 3)`, the extent 4 and the values
  `((0, 1, 2), (1, 3, 2))`, the vertices of each triangle.

```python
X = coordinate(K)
c = cell_index(square, 0)
k, v = index("k"), index("v")

# The position of each local vertex of each triangle, of value shape (2, 3, 2),
# gathered from the vertex positions.
p = component_tensor(y[gathered_index(cell_vertices, (c, v)), k], (c, v, k))

# The affine map of each triangle.
x = component_tensor(
    p[c, 0, k] + (p[c, 1, k] - p[c, 0, k]) * X[0] + (p[c, 2, k] - p[c, 0, k]) * X[1],
    (k,),
)
D = decomposed_domain((domain_map(K, x, RD(2)),))
```

- `D.topology` is `square`, because `x` has the free index `c`.
  `cell_domain(D, 0)` is the family of the two triangles, with the free index
  `c`.
- `pull_back(coordinate(D), cell_domain(D, 0))` is `x`. A later pass evaluates
  `x` on triangle 1 to $(1 - X_1, X_0 + X_1)$.

For a topology `M` of triangles and quadrilaterals, with the cell types
`(Kt, Kq)`, $D$ has one piece for each cell type. A function `u` on `D` is one
expression, pulled back through each cell type:

```python
D = decomposed_domain(
    (domain_map(Kt, x_tri, RD(2)), domain_map(Kq, x_quad, RD(2)))
)

# On Kt with the free index cell_index(M, 0), and on Kq with cell_index(M, 1).
u_tri = pull_back(u, cell_domain(D, 0))
u_quad = pull_back(u, cell_domain(D, 1))
```

#### Legacy UFL

In legacy UFL, a `Mesh` is a decomposed domain, and `SpatialCoordinate(mesh)` is `coordinate(D)`.

### Restriction

Requires: [`DecomposedDomain`](#decomposeddomain),
[Sides of a topology](#sides-of-a-topology), [`GatheredIndex`](#gatheredindex),
[`PullBack`](#pullback), [Reference cell](#reference-cell).

See also: [`Integral` over a decomposed domain](#integral-over-a-decomposed-domain).

The restriction of an expression $f$ on a decomposed domain $D$ to the cells
of a topology $S$ with sides, from the side $s$, is the pull-back of $f$ along
the inclusion of each cell of $S$ into the cell on the side $s$. The side
tables of $S$ give the cell on the side, through a
[`GatheredIndex`](#gatheredindex) entry, and the positions of the vertices of
the inclusion.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `restriction(f, E, s)` | | `f: AbstractExpression`, `E: DecomposedDomain`, `s: str` | `pull_back(f, side_map(D, E, s))`, with $D$ the support of $f$, or $f$ if $f$ has no support |
| `side_domain(D, S)` | | `D: DecomposedDomain`, `S: AbstractTopology` | The decomposed domain of the cells of $S$, with the geometry of the first side of each group |

Notation, for the group $i$ of $S$ with the side `(s, T, j)`:

- $T$ is the topology of $D$, and $K_j$ is the reference cell of cell type $j$
  of $T$.
- `(cells, positions)` are the side tables of the side $s$ of group $i$.
- $F_i$ is the cell type of group $i$.
- `e = cell_index(S, i)` is the free index of the cells of group $i$.
- $Q$ is the ambient domain of $D$.

#### The inclusion

Two families of domain maps with the free index `e` make the inclusion of the
cells of group $i$:

- $B_i$, the cell on the side, is
  `gathered_cell_domain(D, j, gathered_index(cells, (e,)))`. $B_i$ is the cell
  of type $j$ of $T$ on the side of each cell $e$ of $S$ (see
  [`DecomposedDomain`](#decomposeddomain)).
- $A_i$, the inclusion into the cell on the side, is
  `domain_map(F_i, affine_map(F_i, points), K_j)`. Here
  `points = component_tensor(positions[e, v, k], (v, k))`, with the derived
  indices $v$ and $k$, is the position in $K_j$ of each vertex $v$ of the cell
  $e$. $F_i$ must be a reference cell. $A_i$ sends each vertex of the cell $e$
  to the position of that vertex in the cell on the side:
  - For a vertex, $A_i$ is the constant point.
  - For the cell itself with the vertices in order, $A_i$ is the identity map.
  - For a piece of the cell, $A_i$ is the affine map onto the piece.

The restriction of $f$ to group $i$ from the side $s$ is
`pull_back(pull_back(f, B_i), A_i)`, an expression on $F_i$ with the free index
`e`.

Each side uses the vertex numbering of the cell $e$. The restrictions from two
sides are therefore expressions on the same $F_i$, with the same coordinates,
and the difference of the two restrictions is the jump. A cell of $S$ with more
than two sides, such as an edge where three triangles meet, has no jump, and a
form uses the restriction from each side.

#### `SideMap`

`SideMap` is a concrete subclass of `Base`. A `SideMap` is the lazy map of the
cells of $E$ into $D$ through the side $s$, where $E$ is a decomposed domain of
the cells of $S$.

| `SideMap(D, E, s)` |                                                                                                         |
| ------------------ | ------------------------------------------------------------------------------------------------------- |
| Extends            | `Base` |
| Free function      | `side_map(D, E, s)`                                                                                     |
| Takes              | `D, E: DecomposedDomain`, `s: str` |
| Properties         | `ambient_domain` $D$, `parameter_domain` $E$, a decomposed domain with the topology $S$, and `side` $s$ |
| Derived            | `free_indices`, always `()` |
| Requires           | Conditions 1 to 3 below                                                                                 |
| Raises             | `ValueError` if a requirement fails                                                                     |

1. $D$ has a topology.
2. Every cell type of $S$ is an `AbstractReferenceCell`.
3. Every group of $S$ has a side named $s$ with the topology of $D$.

`pull_back(f, SideMap(D, E, s))` is an expression with the support $E$, the
restriction of $f$ to all the groups at once. The pull-back of that expression
by `cell_domain(E, i)` is `pull_back(pull_back(f, B_i), A_i)` (see
[Expanding pull-backs](#expanding-pull-backs)).

#### Side domains

Piece $i$ of `side_domain(D, S)` is
`domain_map(F_i, pull_back(pull_back(coordinate(D), B_i), A_i), Q)`, with
$B_i$ and $A_i$ of the first side of group $i$. $F_i$ must be an
`AbstractReferenceCell`. The topology of the result is $S$. The core does not
check that the other sides have the same geometry.

#### Traces and normals

The pull-back gives the trace of each degree:

- The restriction of a 0-form is the value on the cell of $S$.
- The restriction of the 1-form $\mathrm{d}u$ keeps only the derivatives along
  the cell of $S$. To keep the whole array of $\mathrm{d}u$, as for a normal
  derivative, restrict `components(du, D)`. The axes of `components(du, D)`
  are in the frame of the ambient coordinate domain from every side.

The outward normal of the cell on a side is [`normal(D, E, s)`](#normal).

#### Example

The example uses `D` and `square` of the
[`DecomposedDomain`](#decomposeddomain) example, `I = ReferenceSimplex(1)`, a
function `u` on `D`, and the diagonal of `square`:

```python
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
E = side_domain(D, diagonal)
jump = restriction(u, E, "+") - restriction(u, E, "-")
```

- From `"+"`, the side tables are the cells `(0,)` and the positions
  `(((1, 0), (0, 1)),)`, those of the local vertices 1 and 2 of `K`. From
  `"-"`, the side tables are the cells `(1,)` and the positions
  `(((0, 0), (0, 1)),)`, those of the local vertices 0 and 2. So
  $A_+(t) = (1 - t, t)$ and $A_-(t) = (0, t)$. The coordinates of `E` are
  $(1 - t, t)$.
- Let `pull_back(u, cell_domain(D, 0))` be
  `a[c, 0] + a[c, 1] * X[0] + a[c, 2] * X[1]`. The two restrictions pulled
  back to the diagonal are $a_{0,0} + a_{0,1}(1 - t) + a_{0,2} t$ and
  $a_{1,0} + a_{1,2} t$. The two are equal for all $t$ if `u` is continuous at
  the vertices 1 and 2.
- With the diagonal listed as `(2, 1)`, the positions from `"+"` are
  `(((0, 1), (1, 0)),)` and $A_+(t) = (t, 1 - t)$. The restriction of
  $\mathrm{d}x_0$ from `"+"` then changes sign, from $-\mathrm{d}t$ to
  $\mathrm{d}t$.
- For the facets `S` of a triangle and quadrilateral mesh, with the three
  groups of [Sides of a topology](#sides-of-a-topology), the same two lines
  with `E = side_domain(D, S)` give one jump for all three groups.

#### Legacy UFL

The legacy UFL restriction `f('+')` is `restriction(f, E, "+")` on the decomposed domain $E$ of the interior facets, whose sides are named `'+'` and `'-'`.

### `Integral` over a decomposed domain

Requires: [`Integral`](#integral),
[Concept: orientation of a domain](#concept-orientation-of-a-domain),
[Concept: integral over a domain](#concept-integral-over-a-domain),
[`DecomposedDomain`](#decomposeddomain).

See also: [Expanding integrals](#expanding-integrals).

`Integral(f, D)` for a decomposed domain $D$ of dimension $k$ integrates a
0-form against the density that the pieces induce, and integrates a $k$-form
with the orientation of the pieces. $f$ has the support $D$, the support
`D.ambient_domain` or no support.

Notation, for the piece $j$ of $D$:

- $P_j$ is the parameter domain and $\psi_j$ is the parameterisation.
- `indices_j` are the free indices of piece $j$.
- $n$ is the dimension of the ambient coordinate domain.

#### Expansion into pieces

The integral is the sum over the pieces $j$ of the integrals over the domains
of piece $j$. [`expand_integrals`](#expanding-integrals) rewrites
`Integral(f, D)` as that sum:

```python
# A 0-form f, with the density rho_j of piece j.
index_sum(integral(pull_back(f, cell_domain(D, j)) * rho_j * vol_j, P_j), indices_j)

# A k-form f, with the orientation of piece j.
index_sum(integral(pull_back(f, cell_domain(D, j)), P_j), indices_j)
```

- A piece with no free indices has no `IndexSum`.
- For a coordinate domain $P_j$, `vol_j` is `volume_form(P_j)`.
- For a decomposed domain $P_j$, `vol_j` is absent. The inner integral is
  again the integral of a 0-form against a density, and `expand_integrals`
  rewrites the inner integral in turn.

#### The density

The Jacobian matrix of piece $j$ is
$J_j$ = `components(exterior_derivative(psi_j, P_j), P_j)`, of value shape
`(n, k)`. The density is $\rho_j = \sqrt{\det(J_j^T J_j)}$, the length, the
area or the volume of the image of a small region of $P_j$.

- The density is written as `square_root` of the Leibniz formula for the
  determinant of the $k \times k$ matrix $J_j^T J_j$, a fixed composition for
  each $k$.
- For $k = 0$, the empty determinant is `Integer(1)`, so the integral over a
  set of points is the sum of the values.

The density needs no orientation, so an integral of a positive 0-form is
positive.

#### Integrands on the ambient domain

An $f$ with the support `D.ambient_domain` is pulled back by piece $j$ itself
in place of `cell_domain(D, j)`. $D$ may therefore be a set of points or of
curves in a mesh. A free index of $f$ that is also a free index of piece $j$
takes the value of the member. So `d[i]` pairs a datum with each point of a
family with the index `i`.

#### Examples

The examples use the following objects:

- `D` of the [`DecomposedDomain`](#decomposeddomain) example.
- A function `u` on `D`.
- The decomposed domain `B` of the boundary facets of `D`, with the side
  `"inside"`.

Examples:

- `integral(Integer(1), D)` is the area 1 of the square, the sum of two
  integrals over `K`.
- `integral(Integer(1), B)` is the perimeter 4.
- `integral(restriction(u, B, "inside"), B)` integrates `u` over the boundary.
- `integral(restriction(du, B, "inside"), B)` integrates the 1-form
  $\mathrm{d}u$ along the boundary facets, with the orientation of each facet.

#### Legacy UFL

`f * dx`, `f * ds` and `f * dS` are `integral(f, D)` over the decomposed domain of the cells, of the boundary facets and of the interior facets, with the fields restricted to the facets in `f`. The three are one integral over three domains.

### `Normal`

Requires: [Restriction](#restriction), [Reference cell](#reference-cell),
[`Inverse`](#inverse).

A concrete subclass of `AbstractExpression` for the unit outward normal of the
cell on a side, on the cells of a topology with sides.

| `Normal(D, E, s)` |                                                                                                       |
| ----------------- | ----------------------------------------------------------------------------------------------------- |
| Extends           | `AbstractExpression` |
| Free function     | `normal(D, E, s)`                                                                                     |
| Takes             | `D, E: DecomposedDomain`, `s: str` |
| Properties        | `ambient_domain` $D$, `parameter_domain` $E$ and `side` $s$, as for [`SideMap`](#sidemap)         |
| Requires          | The requirements of `SideMap(D, E, s)`, and $E$ of the dimension $k - 1$ for the dimension $k$ of $D$ |
| Static type       | Value shape `(n,)` for the dimension $n$ of the ambient coordinate domain, the support $E$, real      |
| Raises            | `ValueError` if a requirement fails                                                                   |
| Evaluated by      | [`expand_pull_backs`](#expanding-pull-backs)                                                          |

Notation, for the piece $i$ of $E$:

- $B_i$ and $A_i$ are the maps of [Restriction](#restriction) for the side $s$.
- $K_j$ is the reference cell of the cell on the side.
- $k$ is the dimension of $D$, and $n$ is the dimension of the ambient
  coordinate domain.

#### The pull-back of `Normal`

The pull-back of `Normal(D, E, s)` by `cell_domain(E, i)` is

$$
\nu = \frac{K^T \hat m}{\lvert K^T \hat m \rvert},
$$

where:

- $K$ is the [inverse](#inverse) of the Jacobian matrix of $B_i$ at the points
  of $A_i$. For $k < n$, $K$ is the pseudo-inverse.
- $\hat m$ is a normal of the cell of $E$ in $K_j$. $\hat m$ is orthogonal to
  the columns of the Jacobian matrix of $A_i$, and $\hat m$ points away from
  the centroid of the vertices of $K_j$.

$K$ and $\hat m$ are fixed compositions for each $k$ and $n$, from the side
tables and the vertices of $K_j$, so the core needs no numbering of the
facets.

#### The direction of the normal

- For a cell of $E$ on the boundary of $K_j$, such as a facet or a piece of a
  facet, $\hat m$ points out of $K_j$. The reason is that the centroid lies
  inside a convex reference cell.
- For a cell of $E$ inside $K_j$, such as a cut through the cell, $\hat m$
  points away from the centroid, by convention.
- For $k = n$, the normals from two sides of one cell of $E$ are opposite.
- For $k < n$, such as on the edges of a surface in $R^3$, $\nu$ is the unit
  conormal. The conormal is tangent to the cell on the side and orthogonal to
  the cell of $E$. The conormals from two sides are opposite only where the
  cells on the two sides have the same tangent space. On a surface, a jump
  therefore weights each side by the normal of that side, as in
  `u_p * n_p + u_m * n_m`.

`Normal` is the normal of a cell on a side, not the unit normal field of a
domain of dimension $n - 1$, such as of a surface in $R^3$. The unit normal
field needs an orientation of the domain (see [Concept: orientation of a
domain](#concept-orientation-of-a-domain)) and belongs to a downstream library.

#### Examples

The examples use the following objects:

- `B` of the [`Integral` over a decomposed
  domain](#integral-over-a-decomposed-domain) examples.
- The decomposed domain `E` of the interior facets of `D`, with the sides
  `"+"` and `"-"`.
- `grad_u = components(exterior_derivative(u, D), D)` and `i = index("i")`.

```python
n_B = normal(D, B, "inside")

# The flux of grad u out of the square, the legacy UFL dot(grad(u), n) * ds.
flux = integral(index_sum(restriction(grad_u, B, "inside")[i] * n_B[i], (i,)), B)

# The jump of the normal derivative across the diagonal, from the side "+".
n_E = normal(D, E, "+")
jump = restriction(grad_u, E, "+") - restriction(grad_u, E, "-")
penalty = integral(index_sum(jump[i] * n_E[i], (i,)) ** Integer(2), E)
```

Worked examples for this chapter in EXAMPLES.md: [an interior facet integral](EXAMPLES.md#an-interior-facet-integral), [a Bank-Weiser error estimator](EXAMPLES.md#a-bank-weiser-error-estimator), [a fluid and a solid](EXAMPLES.md#a-fluid-and-a-solid), [a 1D network in a 3D mesh](EXAMPLES.md#a-1d-network-in-a-3d-mesh).

## Function spaces

A function is a member of a function space, such as a coefficient or an
argument of a variational form.

### `AbstractFunctionSpace`

Requires: [`AbstractDomain`](#abstractdomain), [Degree](#degree),
[Scalar type](#scalar-type).

An abstract subclass of the base object `AbstractFunctionSpace` for a space of
functions on a domain.

| `AbstractFunctionSpace` |                                                                                                                                                       |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `Base` |
| Provides                | `domain`, an [`AbstractDomain`](#abstractdomain) $D$ or `None`, `degree`, `value_shape`, `scalar_type` and `argument_spaces`, `()` for a space of fields |
| Defined by              | A subclass outside the core, such as a subclass for the finite element spaces of Basix, and [`DualSpace`](#dualspace) in the core                           |

A space with the domain `None` has functions that do not depend on a domain.

#### Functions in a space

A function in $V$ is an expression with the support $D$ and with the degree,
the value shape and the scalar type of $V$.

For a decomposed domain $D$, the function is defined through the cells of $D$
(see [`DecomposedDomain`](#decomposeddomain)). The concrete subclass does the
following:

- The concrete subclass defines the pull-back of a function in $V$ by
  `cell_domain(D, j)` for every cell type $j$ of $D$, such as with one finite
  element for each cell type.
- The concrete subclass evaluates the pull-back by
  `gathered_cell_domain(D, j, g)` as the same function on the cell $g$ (see
  [`DecomposedDomain`](#decomposeddomain)).

The pull-back of a function by any other domain map into $D$, such as a family
of points, is rewritten through these two (see [Points and point
evaluation](#points-and-point-evaluation)).

For a coordinate domain $D$, such as a reference cell, a function in $V$ is a
terminal on $D$ that the concrete subclass evaluates.

The core does not represent the polynomials, the degrees of freedom or the
basis of $V$, which belong to the concrete subclass. Two spaces are the same
space only if the two spaces are equal as base objects. A concrete
`AbstractFunctionSpace` therefore stores the arguments that identify the
space, such as the finite element.

#### Spaces of constants

A space with the domain `None` and the argument spaces `()` is a space of
constants. A function in the space is an expression with no support and the
degree `()`, a single array with the value shape and the scalar type of $V$.

The integrals of a form decide the cells on which a constant acts.

A space of constants is needed for a constant unknown, such as a Lagrange
multiplier, because an argument comes from a space. Constant data needs no
space. An expression with no support, such as a literal or an array that a
kernel receives, is already a constant.

#### Representation on the reference cell

The pull-back by a cell fixes what a function in $V$ is on the cell. The
pull-back transforms each slot and keeps the frame of each axis of the value
shape (see [`PullBack`](#pullback)). The library of $V$ fixes how the function
is represented on the reference cell.

- For a space of differential forms, the natural representation is the
  pull-back itself, and the legacy UFL map is the pull-back of the slots.
- For a space of vector fields, the library chooses the representation, and
  the legacy UFL map is the choice of the library.

In the table, $k$ is the dimension of $D$ and $n$ is the dimension of the
ambient coordinate domain:

| Space                             | Degree     | Value shape | Legacy UFL map                   |
| --------------------------------- | ---------- | ----------- | -------------------------------- |
| Lagrange                          | `()`       | `()`        | identity                         |
| Vector Lagrange                   | `()`       | `(n,)`      | identity                         |
| Nédélec, $H(\mathrm{curl})$       | `(1,)`     | `()`        | covariant Piola                  |
| Raviart–Thomas, $H(\mathrm{div})$ | `(k - 1,)` | `()`        | contravariant Piola              |
| Raviart–Thomas, $H(\mathrm{div})$ | `()`       | `(n,)`      | contravariant Piola of a library |
| Discontinuous densities           | `(k,)`     | `()`        | $L^2$ Piola                      |
| Regge                             | `(1, 1)`   | `()`        | double covariant                 |

The two rows of Raviart–Thomas differ on a domain of dimension $k < n$, such as
a surface in $R^3$:

- The space of $(k - 1)$-forms has the divergence $\mathrm{d}u$, a $k$-form.
  An integral of $\mathrm{d}u$ needs an orientation of the domain (see
  [Concept: orientation of a domain](#concept-orientation-of-a-domain)).
- The space of vector fields has the divergence
  `index_sum(components(exterior_derivative(u, D), D)[i, i], (i,))`, a 0-form,
  and needs no orientation. The library represents a function in the space of
  vector fields on a cell as $J \hat u / \rho$. Here $\hat u$ is a function on
  the reference cell, $J$ is the Jacobian matrix of the cell, and
  $\rho = \sqrt{\det(J^T J)}$ is the density of [`Integral` over a decomposed
  domain](#integral-over-a-decomposed-domain).

#### Mixed and Sobolev spaces

A mixed problem needs no mixed space, because an argument of a [form](#form)
may be a `tuple` of functions in several spaces.

The Sobolev space of $V$ is not yet covered. The Sobolev space is expected to
be stated through the [restrictions](#restriction) of the functions of $V$ to
the facets of the topology of $D$, from the two sides of each facet.

#### Examples

The examples use a decomposed domain `D` of `ReferenceSimplex(2)` cells in
`RD(2)`, so $k = n = 2$.

- The Lagrange space of polynomial degree 1 on `D` has the degree `()` and the
  value shape `()`.
- The Nédélec and the Raviart–Thomas spaces on `D` both have the degree `(1,)`,
  because $k - 1 = 1$. A Raviart–Thomas function is the vector field of legacy
  UFL contracted with the volume form.
- The discontinuous Lagrange space of polynomial degree 0 on `D`, as densities,
  has the degree `(2,)`.

#### Legacy UFL

The legacy UFL `Real` space is a space of constants, and the legacy UFL `Constant` is constant data.

In legacy UFL, a `FunctionSpace` holds a mesh and an element, and the element states the value shape, the map to the reference cell and the Sobolev space. In UFLx, the pull-back states the meaning of a function on a cell, the library states the representation, and the space holds a decomposed domain, the legacy mesh.

### `Function`

Requires: [`AbstractFunctionSpace`](#abstractfunctionspace),
[`AbstractExpression`](#abstractexpression), [Support](#support),
[Degree](#degree).

See also: [`PullBack`](#pullback), [Restriction](#restriction),
[Variational forms](#variational-forms).

A concrete subclass of `AbstractExpression` for a function, a member of a
function space.

| `Function(V, label)` |                                                                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| Extends              | `AbstractExpression` |
| Free function        | `function(V, label)`                                                                                                         |
| Takes                | `V: AbstractFunctionSpace`, `label: str \| tuple` |
| Properties           | `space` $V$, an `AbstractFunctionSpace`, and `label`, a `str` or a `tuple` with the rules of the label of an [`Index`](#index) |
| Static type          | The value shape of $V$, and no free indices                                                                                  |
| Raises               | `TypeError` if a requirement fails                                                                                           |

The label identifies the function, so two functions are the same only if the
spaces and the labels are equal. A function in a [dual space](#dualspace) is an
antilinear form, such as a known load.

The exterior derivative of a function
remains after [applying derivatives](#applying-derivatives).

A function is not an argument or a coefficient by itself. Whether a function is
an argument of a form is stated by the form (see [`Form`](#form)).

Examples, with a decomposed domain `D` and a space `V` on `D` of the degree
`()` and the value shape `()`:

- `u = function(V, "u")` has the support `D`, the degree `()` and the value
  shape `()`, and `function(V, "u") == u` is `True`.

#### Legacy UFL

The legacy UFL `Coefficient`, `TestFunction` and `TrialFunction` are all functions. A test function and a trial function are the functions that a form states as the arguments at the positions 0 and 1.

## Variational forms

A variational form is an expression that is abstracted over some functions, the
arguments, as a multilinear map is over the arguments of the map. A form binds
the arguments, as a [`ComponentTensor`](#componenttensor-and-indexsum) binds
indices, and an action gives a function to an argument of a form, as an
[`Indexed`](#indexed) gives an index to an axis.

### Argument spaces

Requires: [`AbstractExpression`](#abstractexpression),
[`AbstractFunctionSpace`](#abstractfunctionspace).

Every `AbstractExpression` derives the argument spaces, a `tuple` with the
space of each argument of a form, in order of position, and `()` for an
expression that is not a form. The space of an argument is an
`AbstractFunctionSpace` for an argument of one function, or a `tuple` of at
least two spaces for an argument of a `tuple` of functions, the product of the
spaces. The arity is the number of
argument spaces. A functional has the arity 0, a linear form the arity 1 and a
bilinear form the arity 2, whatever the number of functions in each argument.

A node passes on the argument spaces of an operand in which the node is linear
or antilinear, and every other operand must have the argument spaces `()`, or
construction raises `ValueError`. The [static type
table](#appendix-static-type-table) gives the rule of each node class. A leaf
class states the argument spaces as part of the contract of the class, as the
degree.

A form is linear in an argument with the space $V$, and antilinear in an
argument with the space `conjugate_space(V)` (see
[`ConjugateSpace`](#conjugatespace)). The form states the linearity when the
form binds the argument (see [`Form`](#form)).

The argument spaces are static, as the value shape and the degree are, so a
form with a wrong arity is rejected during construction, and a load that is
antilinear in `v` and a form that is linear in `v` do not add. The arity and
the degree are independent. The degree counts the slots of a differential
form, and the arity counts the arguments that a form binds.

### `ConjugateSpace`

Requires: [Argument spaces](#argument-spaces),
[`AbstractFunctionSpace`](#abstractfunctionspace), [Scalar type](#scalar-type),
[`Conjugate`](#conjugate).

See also: [`DualSpace`](#dualspace), [`Form`](#form).

A concrete subclass of `AbstractFunctionSpace` for the complex conjugate space
$\bar V$ of a complex space $V$.

| `ConjugateSpace(V)` |                                                                                                   |
| ------------------- | ------------------------------------------------------------------------------------------------- |
| Extends             | `AbstractFunctionSpace` |
| Free function       | `conjugate_space(V)` |
| Takes               | `V: AbstractFunctionSpace` |
| Returns             | `V` for a real `V`, and `ConjugateSpace(V)` for a complex `V` |
| Properties          | `space` $V$, an `AbstractFunctionSpace` with the complex scalar type, not a `ConjugateSpace`       |
| Derived             | `domain`, `degree`, `value_shape`, `scalar_type` and `argument_spaces`, those of $V$ |
| Raises              | `TypeError` for a $V$ that is not a space, and `ValueError` for a real $V$ or a `ConjugateSpace`  |

$\bar V$ has the functions of $V$ with the conjugate scalar multiplication, so
a function in $\bar V$ is the conjugate of a function in $V$, and a linear
form on $\bar V$ is an antilinear form on $V$. The core uses $\bar V$ only as
an argument space, to record that a form is antilinear in an argument. The conjugate
of an argument space is `conjugate_space(V)` for a space $V$ and $V$ for
`ConjugateSpace(V)`, so conjugating twice gives the space back (see
[`Conjugate`](#conjugate)).

### `DualSpace`

Requires: [Argument spaces](#argument-spaces),
[`AbstractFunctionSpace`](#abstractfunctionspace), [`Function`](#function).

See also: [`Action`](#action).

A concrete subclass of `AbstractFunctionSpace` for the space of antilinear
forms on a space or on the product of spaces.

| `DualSpace(S)` |                                                                                                                            |
| -------------- | -------------------------------------------------------------------------------------------------------------------------- |
| Extends        | `AbstractFunctionSpace` |
| Free function  | `dual(S)`                                                                                                                  |
| Takes          | `S: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]` |
| Properties     | `space` $S$, an `AbstractFunctionSpace`, or a `tuple` of at least two spaces                                               |
| Derived        | `domain`, `degree`, `value_shape`, `scalar_type` and `argument_spaces`, below |
| Raises         | `TypeError` for an $S$ that is not a space or a `tuple` of spaces, and `ValueError` for a `tuple` of fewer than two spaces |

`dual(S)` provides:

- The domain `None`, the degree `()` and the value shape `()`.
- The argument spaces `(S',)`, with $S'$ the conjugate of $S$ as in
  [`Conjugate`](#conjugate), member by member for a `tuple` $S$.
- The scalar type of $S$, complex if any space of a `tuple` is complex.

#### Antilinear forms

A function in `dual(S)` is an antilinear form on $S$, a linear form on
`conjugate_space(S)` (see [`ConjugateSpace`](#conjugatespace)). A form is
antilinear in the same way in an argument bound as `conjugate(v)`, such as the
argument $v$ of $\int u \bar v$. A residual bound as `(conjugate(v),)` and an
assembled load therefore have the same argument spaces, and the difference of
the two is a form.

The space of linear forms on a complex $S$ is `dual(conjugate_space(S))`,
with the argument spaces `(S,)`. The core therefore fixes no convention on the
argument in which a form is antilinear.

The pairing of a function `l` in `dual(S)` with a function `w` in $S$ is
`action(l, 0, w)`, antilinear in `w`. The library that provides the values of
`l` evaluates the pairing.

#### The dual of a product

For a `tuple` $S$, a function in `dual(S)` is one antilinear form on the
product of the spaces, such as the load of a mixed problem. No single function
is in the product.

#### The double dual

`dual(dual(V))` is the space of antilinear forms on `dual(V)`, with the
argument spaces `(conjugate_space(dual(V)),)`. `dual(dual(V))` is not `V`,
because `V` has the argument spaces `()`.

The canonical map takes a function `w` in `V` to
`form(conjugate(action(l, 0, w)), (conjugate(l),))`, with a function `l` in
`dual(V)`. The image is linear in `w` and antilinear in `l`. A library may
store `w` and the image with the same values.

#### Examples

The examples use the spaces `V` and `Q` and the functions `v` in `V` and `q`
in `Q` of the [`Form`](#form) examples.

- `l = function(dual(V), "l")` has the argument spaces
  `(conjugate_space(V),)`. `action(l, 0, v)` is the pairing
  $\langle l, v \rangle$.
- `m = function(dual((V, Q)), "m")` has the argument spaces
  `((conjugate_space(V), conjugate_space(Q)),)`. `action(m, 0, (v, q))` is the
  pairing of `m` with `(v, q)`.
- Let `w` be a function in `V` and `l` a function in `dual(V)`.
  `form(conjugate(action(l, 0, w)), (conjugate(l),))` has the argument spaces
  `(conjugate_space(dual(V)),)`, the static type of a function in
  `dual(dual(V))`. The form evaluates every antilinear form on `V` at `w`.

#### Legacy UFL

The legacy UFL `Cofunction` and `Coargument` are a fixed and a bound function in `dual(V)`. In legacy UFL, the dual of the dual space is $V$. In UFLx, the canonical map is a form, so the spaces keep distinct argument spaces.

### `Form`

Requires: [Argument spaces](#argument-spaces), [`Function`](#function).

A concrete subclass of `AbstractExpression` for a variational form, an
expression with some of its functions bound as arguments.

| `Form(expression, arguments)` |                                                                                                                                                                                                    |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                       | `AbstractExpression` |
| Free function                 | `form(expression, arguments)`                                                                                                                                                                      |
| Takes                         | `expression: AbstractExpression`, `arguments: tuple[Function \| tuple[Function, ...], ...]` |
| Properties                    | `expression`, an expression with the argument spaces `()`, and `arguments`, a `tuple` of entries, each a [function](#function) `v`, the [conjugate](#conjugate) `conjugate(v)` of a function, or a `tuple` of at least two of these |
| Requires                      | Every function in `arguments` distinct                                                                                                                                                             |
| Static type                   | That of the expression, with the argument spaces the space $V$ of each entry `v`, `conjugate_space(V)` of each entry `conjugate(v)`, or the `tuple` of the spaces of the members of a `tuple` entry, in order |
| Raises                        | `TypeError` for an entry that is not a `Function`, a `Conjugate` of a `Function` or a `tuple` of these, and `ValueError` for a repeated function, a `tuple` of fewer than two members or an expression with argument spaces |

The position
of an argument is the index of the entry in `arguments`. The entry states the
linearity: for an entry `v`, the expression is linear in `v`, and for an entry
`conjugate(v)`, the expression is linear in $\bar v$, so antilinear in `v`,
as in $\int u \bar v$ (see [`ConjugateSpace`](#conjugatespace)). For a real
space, the two entries give the same argument space.
For a `tuple` entry, the expression is linear in the members jointly, so a term
may depend on one member only. Construction does not check the linearity, and
a pass may check the linearity (see [Checking linearity](#checking-linearity)).
Every other function in the expression is fixed, so one function may be an argument
of one form and fixed in another. The scope of a bound function is lexical, as
the scope of an index is (see [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum)): an occurrence of a function refers
to the innermost enclosing `Form` that binds the function, and is fixed
otherwise.

An argument that is a `tuple` of functions gives the forms of mixed problems
with no mixed space. Each member is a function in an ordinary space, so the
members may have different degrees, such as a Raviart–Thomas and a
discontinuous function, and different domains, such as a Lagrange multiplier
on the boundary. A library may assemble the arguments as one system or in
blocks. The block of a member of each argument is the form of the expression
with the other members replaced by `Zero`.

The expression may have any value shape and support. A form whose expression
is a field is a linear operator. The support is that of the expression, so
a form of integrals has no support, although the arguments have one. A
substitution renames a bound argument to avoid capture, as for indices (see
[Passes](#passes)).

Examples, with `u`, `v` and a fixed `w` in a space `V` on `D`, and an
expression `a` linear in `u` and `v`:

- `form(a, (v, u))` is a bilinear form with the argument spaces `(V, V)`.
- `form(a, (v,))` is a linear form in `v`.
- `form(a, (v, u)) + form(a, (v,))` raises `ValueError`.
- With `m = integral(u * conjugate(v), D)`, `form(m, (conjugate(v), u))` is
  the sesquilinear mass form, with the argument spaces
  `(conjugate_space(V), V)`.
- With `r = integral(f * conjugate(v), D)` and `l = function(dual(V), "l")`,
  `form(r, (conjugate(v),)) - l` is a residual, with the argument spaces
  `(conjugate_space(V),)`. For a complex `V`, `form(integral(f * v, D), (v,)) - l`
  raises `ValueError`, because the argument spaces `(V,)` and
  `(conjugate_space(V),)` differ.
- With a velocity space `V` and a pressure space `Q`, functions `u`, `v` in `V`
  and `p`, `q` in `Q`, and the Stokes expression `e` of the terms
  $\nabla u : \nabla v$, $-p \, \mathrm{div}\, v$ and $-q \, \mathrm{div}\, u$,
  `form(e, ((v, q), (u, p)))` is a bilinear form with the argument spaces
  `((V, Q), (V, Q))`. `form(e, (v, q, u, p))` states a wrong form of arity
  4, because no term of `e` is linear in all four functions. Construction
  accepts the form, and a linearity check rejects the form.

#### Checking linearity

A linearity check is a pass that asks whether an expression is linear in a
function `v`, or jointly in a `tuple` of functions (see [New
questions](#new-questions)). The pass gives every node one of five classes:
zero, constant, linear, antilinear or nonlinear in `v`. A `Form` is checked
one entry at a time. The expression of the `Form` must be linear or zero in
`v` for an entry `v`, and antilinear or zero in `v` for an entry
`conjugate(v)`.

#### Legacy UFL

In legacy UFL, the arguments of a form are the test functions and the trial functions in the integrand, with the argument numbers of the functions, and a form is not an expression. In UFLx, a function has no argument number, the form states the arguments, and a functional is any expression with the argument spaces `()`. The legacy UFL `MixedElement` and `split` are an argument that is a `tuple` of functions.

### `Action`

Requires: [`Form`](#form).

A concrete subclass of `AbstractExpression` for the action of a form, an
expression given to one argument of the form.

| `Action(f, position, g)` |                                                                                                                                                        |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Extends                  | `AbstractExpression` |
| Free function            | `action(f, position, g)`                                                                                                                               |
| Takes                    | `f: AbstractExpression`, `position: int`, `g: AbstractExpression \| tuple[AbstractExpression, ...]` |
| Properties               | `form`, an expression with more than `position` argument spaces, `position`, an `int`, and `argument`, the expression $g$, or a `tuple` of expressions |
| Requires                 | The conditions below                                                                                                                                   |
| Static type              | That of $f$, with the argument spaces of $f$ without the entry at `position`.                                                                          |
| Raises                   | `ValueError` if a requirement fails                                                                                                                    |
| Evaluated by             | [`apply_actions`](#applying-actions)                                                                                                                   |

Let the argument space at `position` be $V$ or `ConjugateSpace(V)`.

- For a space $V$, $g$ has the static type of a function in $V$, with the
  support of that function or no support.
- For a `tuple` of spaces, $g$ is a `tuple` with one such member per space.

`Action` is lazy. An action on a leaf, such as a function in a dual space,
remains for the library of the leaf.

#### Examples

The examples use `a`, `u`, `v` and `w` of the [`Form`](#form) examples, and
`l = function(dual(V), "l")`.

- `action(form(a, (v, u)), 1, w)` is a linear form in `v`. A pass evaluates
  the action to `form(b, (v,))`, with `b` the expression `a` with `u` replaced
  by `w`.
- `action(action(form(a, (v, u)), 1, w), 0, w)` is a functional.
- `action(l, 0, w)` is the pairing of `l` with `w`, with no support.
- `action(form(a, (v, u)) + 2 * form(a, (v, u)), 1, w)` moves into both
  operands of the sum, and into the form of the product.
- Let `w` in `V` and `r` in `Q` be fixed functions, with the Stokes form of the
  [`Form`](#form) examples. `action(form(e, ((v, q), (u, p))), 1, (w, r))` is
  a linear form in `(v, q)`, the residual at `(w, r)`.

#### Legacy UFL

Legacy UFL operations are actions and rebindings. The action on the last argument of a bilinear form is `action(f, 1, w)`. The adjoint of `form(e, (conjugate(v), u))` is `form(conjugate(e), (conjugate(u), v))`, because `conjugate(e)` is linear in `v` and antilinear in `u`. The Gateaux derivative with respect to a function adds an argument (see [Gateaux derivative](#gateaux-derivative)).

Worked examples for this chapter in EXAMPLES.md: [the Laplace-Beltrami operator](EXAMPLES.md#the-laplace-beltrami-operator), [the flux through a surface](EXAMPLES.md#the-flux-through-a-surface), [a classic FEniCS layer for the Poisson problem](EXAMPLES.md#a-classic-fenics-layer-for-the-poisson-problem).

## Core transformations

Requires: [Passes](#passes), [`Integral`](#integral), [`PullBack`](#pullback),
[`ExteriorDerivative`](#exteriorderivative), [`Components`](#components),
[`Action`](#action).

The core has six transformations, each a pass (see [Passes](#passes)) that
evaluates one node class. `expand_integrals` dispatches on the domain of the
`Integral` and not on the operand, and `reduce_component_tensors` has one rule,
for an `Indexed` whose operand is a `ComponentTensor`. Each of the other four
has a default for the node classes of the operand, and a table of the operand
classes with a rule of their own. The [pass table](#appendix-pass-table) gives
every node class for each of the four.

| Core transformation                                         | Evaluates            | Default                                                     |
| ----------------------------------------------------------- | -------------------- | ----------------------------------------------------------- |
| [`expand_integrals(f)`](#expanding-integrals)               | `Integral`           | One rule, by the domain of the `Integral`                   |
| [`expand_pull_backs(f)`](#expanding-pull-backs)             | `PullBack`           | Commutes with every node class                              |
| [`apply_derivatives(f)`](#applying-derivatives)             | `ExteriorDerivative` | Linear for the linear node classes, the chain rule for `()` |
| [`expand_components(f)`](#expanding-components)             | `Components`         | Linear for the linear node classes                          |
| [`apply_actions(f)`](#applying-actions)                     | `Action`             | Moves into the operand that carries the argument            |
| [`reduce_component_tensors(f)`](#reducing-component-tensors) | `Indexed`            | One rule, for a `ComponentTensor` operand, and otherwise remains |

The linear node classes are the following:

- `Sum`, `Tensor`, `Indexed`, `ComponentTensor` and `IndexSum`, linear in
  every operand.
- `Division`, linear in the first operand.

An evaluated node that is linear moves into the operands of a linear node
class in which the node is linear.

A rule obeys the following:

- A rule rewrites a node of the evaluated class to a mathematically equal
  expression of node classes of the core, and never simplifies.
- The rewritten expression has the free indices of the node, with the same
  extents, possibly in another order, and the argument spaces of the node.
- The rule dispatches on the class of the operand and moves the evaluated node
  towards the leaves.

The pass rewrites inner nodes first, until no rule applies. A node without a
rule is rebuilt on the rewritten edges, and the evaluated node remains on the
rebuilt node.

In a `Form`, an evaluated node moves only into the expression. In an `Action`,
an evaluated node moves only into the form $f$. The reason is that the
arguments and $g$ are bound, not on the support. For example,
$\psi^* \mathrm{action}(f, s, g) = \mathrm{action}(\psi^* f, s, g)$, because
the action replaces $\psi^* u$ in $\psi^* f$ by $\psi^* g$.

Each core transformation rewrites only the evaluated class, so the six run
separately or together. Evaluation runs the six together until no rule
applies. Evaluation terminates for the following reasons:

- An `Integral` over a decomposed domain is replaced by `Integral`s over the
  parameter domains of the pieces, which nest less deeply.
- A `PullBack` moves inside an `ExteriorDerivative`, and never the reverse.
- An action removes an `Action`.
- A reduction removes an `Indexed` and a `ComponentTensor`.

A downstream lazy node class may join evaluation if the rules of the class
satisfy the following:

- The rules rewrite only that class.
- The rules move the node towards the leaves.
- The rules never reverse a core rule, as moving an `ExteriorDerivative`
  inside a `PullBack` would.

The terminals are `Coordinate`, `Function`, `Normal` and the leaf expressions
with a support of downstream libraries (see [Leaves and
subclassing](#leaves-and-subclassing)).

After evaluation, every `PullBack`, `ExteriorDerivative`, `Components` and
`Action` applies to an expression that contains a downstream terminal, with
two exceptions:

- the coordinate frame $\mathrm{d}X$,
- the components of $\mathrm{d}X$ on a decomposed domain of lower dimension
  than the ambient coordinate domain.

An expression without a downstream terminal evaluates to a sum of products and
tensor products of 0-forms and the coordinate 1-forms $\mathrm{d}X_i$, inside
any `Form`.

### Expanding integrals

Requires: [Core transformations](#core-transformations), [`Integral`](#integral),
[`Integral` over a decomposed domain](#integral-over-a-decomposed-domain).

`expand_integrals` rewrites every `Integral` whose domain is a decomposed
domain as the sum over the pieces of [`Integral` over a decomposed
domain](#integral-over-a-decomposed-domain), with these rules:

| Domain of the `Integral`                    | Rewrite of `Integral(f, D)`                                                                          |
| ------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| A coordinate domain                         | Remains                                                                                              |
| A decomposed domain, $f$ of the degree `()` | The sum over the pieces of the integrals of the pull-back of $f$ times the density $\rho_j$ `vol_j` |
| A decomposed domain, $f$ of the degree `(k,)`, $k \ge 1$ | The sum over the pieces of the integrals of the pull-back of $f$                             |

The pull-back is by `cell_domain(D, j)` for $f$ with the support $D$, and by
piece $j$ itself for $f$ with the support `D.ambient_domain`. An $f$ with no
support takes either, with the same result. The pass rewrites inner integrals
first, so an `Integral` over a nested decomposed domain is rewritten from the
inside out.

### Expanding pull-backs

Requires: [Core transformations](#core-transformations), [`PullBack`](#pullback),
[`Coordinate`](#coordinate), [`ExteriorDerivative`](#exteriorderivative),
[`Components`](#components), [`Inverse`](#inverse).

For a domain $D$ with the parameterisation $\psi$, the parameter domain $P$
and the ambient domain $Q$, $\psi^*$ commutes with every node class of the
core, with these rules:

| Operand                                                                   | Rewrite of $\psi^*$ applied to the operand                                           |
| ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| An expression $c$ with no support                                         | $c$, or $c$ `+ Zero(...)` (see [Expressions with no support](#expressions-with-no-support)) |
| `Coordinate(Q)`                                                           | $\psi$                                                                               |
| `ExteriorDerivative(f, Q)`                                                | $\mathrm{d}(\psi^* f)$, with the exterior derivative on $P$                          |
| `PullBack(g, D2)`, for a `DomainMap` `D2`                                 | Remains, so pull-backs nest                                                          |
| `PullBack(f, SideMap(D2, E, s))`, for `cell_domain(E, i)`                 | `pull_back(pull_back(f, B_i), A_i)` (see [Restriction](#restriction))                |
| `Normal(D2, E, s)`, for `cell_domain(E, i)`                               | $K^T \hat m / \lvert K^T \hat m \rvert$ (see [`Normal`](#normal))                      |
| `PullBack(f, SideMap(D2, E, s))` or `Normal(D2, E, s)`, for any other $D$ | Remains                                                                              |
| `Components(f, Q)`, $P$ of the dimension of $Q$                           | $\sum_i K_{ij} \, \mathrm{Components}(\psi^* f, P)_i$ for entry $j$ of the last axis |
| `Components(f, Q)`, $P$ of lower dimension than $Q$                       | Remains                                                                              |
| `Function` or a downstream terminal                                       | Remains, for the library of the terminal                                             |

#### Expressions with no support

The rule for an expression with no support comes before every other rule and
before commuting, whatever the class of the expression. A pull-back therefore
never moves inside an `Integral`. The integral has no support, and the
integrand stays on the domain of the integral.

`PullBack(c, D)` has the free indices of the family $D$, and $c$ may not.

- If every free index of $D$ is free in $c$, the rewrite is $c$.
- Otherwise the rewrite is $c$ `+ Zero(...)`. The `Zero` has the value shape,
  the degree and the argument spaces of $c$, and the missing free indices with
  the extents of the family (see [`Zero`](#zero)).

For example, take the family of three quadrature points `xs` of the [quadrature example](EXAMPLES.md#quadrature).
`pull_back(Integer(1), xs)` is rewritten as
`Integer(1) + Zero((), ((i, 3),), (), ())`, and
`index_sum(integral(pull_back(Integer(1), xs), RD(0)), (i,))` counts the
points.

#### Nested pull-backs

The pass rewrites inner nodes first. A `PullBack(g, D2)` that the outer
pull-back meets is therefore one that remains, with $g$ one of the following
(see the rows above):

- a terminal,
- a `Components` of lower dimension,
- a restriction by a `SideMap`.

The outer pull-back then remains outside the inner one, so pull-backs nest.
The library of a terminal evaluates the inner pull-back, such as a function on
the cell on the side of a facet. The outer pull-back applies to the expression
the library gives, such as the pull-back to the facet of a basis function on
the reference cell. That outer pull-back remains in turn, for the library to
evaluate at the points of the facet (see [Restriction](#restriction) and
[An interior facet
integral](EXAMPLES.md#an-interior-facet-integral)).

#### The rule for `Components`

In the rule for `Components`:

- $J$ is `components(exterior_derivative(D.parameterisation, P), P)`, of value
  shape `(n, k)`, with $k$ the dimension of $P$ and $n$ the dimension of the
  frame of `Components(f, Q)`.
- $K$ is `inverse(J)`, the pseudo-inverse for $k < n$ (see
  [`Inverse`](#inverse)).

The vector fields $\sum_i K_{ij} \, \partial_i$ on $P$ push forward to the
orthogonal projections of the $\partial_j$ of the frame onto the tangent space
of $D$. For $k = n$, the projections are the $\partial_j$ themselves. For a
decomposed $Q$ of dimension $k < n$, `Components(f, Q)` is defined with those
projections (see [`Components`](#components)), so the rule holds for $f$ on
$Q$.

If $P$ has lower dimension than $Q$, the same composition gives the components
of the part of $f$ tangent to $D$, not the components of $f$. The pull-back
therefore remains.

### Applying derivatives

Requires: [Core transformations](#core-transformations),
[`ExteriorDerivative`](#exteriorderivative), [`Product`](#product),
[`TensorProduct`](#tensorproduct), [`Zero`](#zero), [`Conjugate`](#conjugate).

The exterior derivative is linear for the linear node classes and follows the
chain rule for the node classes of the degree `()`, with these rules:

| Operand                                                                     | Rewrite of $\mathrm{d}$ applied to the operand                                                                     |
| --------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `Product` $\alpha \wedge \beta$, with $\alpha$ a $k$-form                   | $\mathrm{d}\alpha \wedge \beta + (-1)^k \alpha \wedge \mathrm{d}\beta$                                             |
| `TensorProduct` $\alpha \otimes \beta$                                      | The rule of `Product`, written with `Product`                                                                      |
| `Division` $a / b$                                                          | $\mathrm{d}a / b - (\mathrm{d}b \wedge a) / (b b)$                                                                 |
| `Inverse` $A^+$                                                             | $-A^+ \, \mathrm{d}A \, A^+ + A^+ (A^+)^T \, \mathrm{d}A^T \, (I - A A^+)$                                         |
| `Absolute` $\lvert a \rvert$                                                | $(\bar{a} \, \mathrm{d}a + a \, \overline{\mathrm{d}a}) / (2 \lvert a \rvert)$                                     |
| `Conjugate` $\bar{\omega}$                                                  | $\overline{\mathrm{d}\omega}$                                                                                      |
| `Power` $a^b$ with an exponent $b$ with no support                          | $b \, a^{b - 1} \, \mathrm{d}a$, the chain rule, with `b - Integer(1)` the exponent                                |
| `Power` $a^b$ with an exponent $b$ with a support                           | Remains                                                                                                            |
| `ExteriorDerivative`, an expression with no support, or a `PullBack` of one | `Zero`, with the value shape, the free indices and the argument spaces of the operand and the degree raised by one |
| `Coordinate(P)`, `Normal`, `Function` or a downstream terminal              | Remains                                                                                                            |
| `Components`, or a `PullBack` of an expression with a support               | Remains                                                                                                            |

`Components` remains because the exterior derivative and the interior product
do not commute. A `PullBack` with a support remains because moving the
exterior derivative inside the pull-back would reverse a rule of
[`expand_pull_backs`](#expanding-pull-backs).

A `TensorProduct` under an exterior derivative has at most one slot, so
$\alpha \otimes \beta = \alpha \wedge \beta$. Written with `TensorProduct`,
the two terms would have different numbers of slots. In the rule for
`Inverse`, the matrix products are in the value shape, with the wedge product
of the 0-forms $A^+$ and $(A^+)^T$ with the 1-forms $\mathrm{d}A$ and
$\mathrm{d}A^T$, and $I$ is the identity matrix of value shape `(n, n)`, the
`Tensor` of `Integer(1)` and `Integer(0)` entries. For $k = n$,
$I - A A^+$ is zero and the rule is $-A^{-1} \, \mathrm{d}A \, A^{-1}$, but
the rule does not drop the term, because the core never simplifies.
`Absolute` is not complex differentiable, and for a real $a$ the rule is the
chain rule.

### Expanding components

Requires: [Core transformations](#core-transformations),
[`Components`](#components), [`Product`](#product), [`TensorProduct`](#tensorproduct),
[`Zero`](#zero), [`Conjugate`](#conjugate).

The interior product $\iota_j$ with $\partial_j$ is linear for the linear node
classes, and the new axis stays the last axis of the value shape, with these
rules:

| Operand                                                                                                                              | Rewrite of $\iota_j$ applied to the operand                                                                |
| ------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------- |
| `Product` $\alpha \wedge \beta$, with the first slot of $\alpha$ of the degree $k$                                                   | $\iota_j\alpha \wedge \beta + (-1)^k \alpha \wedge \iota_j\beta$ when both operands have a slot, $f \, \iota_j\beta$ for a 0-form $f = \alpha$, and $\iota_j\alpha \, f$ for a 0-form $f = \beta$ |
| `TensorProduct` $\alpha \otimes \beta$                                                                                               | $(\iota_j\alpha) \otimes \beta$ if $\alpha$ has a slot, and $f \, \iota_j\beta$ for a 0-form $f$           |
| $\mathrm{d}X_i$, on a support of the dimension of the frame                                                                          | $\delta_{ij}$                                                                                              |
| `Zero`                                                                                                                               | `Zero` with the new value shape and degree                                                                 |
| `Conjugate` $\bar{\alpha}$                                                                                                           | $\overline{\iota_j \alpha}$                                                                                |
| `Function` or a downstream terminal                                                                                                  | Remains                                                                                                    |
| $\mathrm{d}X_i$ on a decomposed domain of lower dimension than the frame, any other `ExteriorDerivative`, `Components` or `PullBack` | Remains                                                                                                    |

The rule for $\mathrm{d}X_i$ applies to `ExteriorDerivative(Coordinate(P), P)`,
reached through `Indexed`. On a decomposed domain of lower dimension than the
frame, the components of $\mathrm{d}X$ are the orthogonal projection onto the
tangent space, which `expand_pull_backs` evaluates on each piece as $J K$.
The components of the exterior derivative of any other expression remain, such
as the gradient `components(exterior_derivative(u, P), P)` of a function `u`.

The components of `dX` are the identity matrix, the `Tensor` of `Integer(1)`
and `Integer(0)` entries, or `Zero((0, 0), (), (), ())` for $n = 0$.

### Applying actions

Requires: [Core transformations](#core-transformations), [`Action`](#action),
[`Form`](#form), [Passes](#passes).

An action moves into the operand that carries the argument, with these rules:

| Operand                                                  | Rewrite of `Action(f, position, g)` with $f$ the operand                                                                               |
| -------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `Form(e, arguments)`                                     | $e$ with the entry at `position` replaced by $g$, member by member for a `tuple`, in a `Form` of the other entries if any remain        |
| `Sum`, `Tensor`                                          | The action of every operand                                                                                                            |
| `Product`, `TensorProduct`                               | The action of the operand that carries the argument at `position`, with the position counted within that operand                       |
| `Zero`                                                   | `Zero` with the argument spaces without the entry at `position`                                                                        |
| `Function` or a downstream leaf, such as in a dual space | Remains, for the library of the leaf                                                                                           |
| `Action`                                                 | Remains                                                                                                                        |

The substitution into a `Form` avoids capture (see [Passes](#passes)). The
pass rewrites inner actions first, so an outer action meets an `Action` only on
a leaf, and the outer action remains.

### Reducing component tensors

Requires: [Core transformations](#core-transformations), [`Indexed`](#indexed),
[`ComponentTensor` and `IndexSum`](#componenttensor-and-indexsum),
[Passes](#passes).

`reduce_component_tensors` evaluates an `Indexed` whose operand is a
`ComponentTensor`, with these rules:

| Operand of the `Indexed`      | Rewrite of `Indexed(operand, entries)`                                                            |
| ----------------------------- | ------------------------------------------------------------------------------------------------- |
| `ComponentTensor(e, indices)` | $e$ with each index of `indices` replaced by the entries of the axes of the index, when possible |
| Any other node class          | Remains                                                                                           |

The entries are matched to the indices axis by axis. An `Index` of `indices`
covers one axis and is replaced by the entry of that axis, an `int`, an
`Index`, a `GatheredIndex` or an `EllipsisIndex` that covers one axis. An
`EllipsisIndex` of `indices` covers the axes of the extents of the
`EllipsisIndex`, and is replaced by the entries of those axes, or by one
`EllipsisIndex` entry that covers exactly those axes. An occurrence of the
replaced index in the entries of an `Indexed` or of a `GatheredIndex` takes
the replacing entries in place, and an occurrence in the free indices of a
`Zero` is removed or replaced in the same way. An `EllipsisIndex` cannot occur
in a `GatheredIndex` and cannot replace a part of the axes of an
`EllipsisIndex` of `indices`. When a replacement is not possible, the
`Indexed` remains. The substitution
avoids capture (see [Passes](#passes)).

The rewrite is the identity `component_tensor(e, (i,))[j]` $= e[i := j]$, so
the rewrite never changes the value, and the result has the free indices of
the `Indexed`. For example, `component_tensor(points[i, k], (k,))[0]` is
rewritten as `points[i, 0]`, and with `p` of the
[`DecomposedDomain`](#decomposeddomain) examples, `p[c, 0, k]` is rewritten as
`y[gathered_index(cell_vertices, (c, 0)), k]`.

### Gateaux derivative

Requires: [Core transformations](#core-transformations), [`Form`](#form),
[`Function`](#function).

> Not yet covered. This section states only the shape of the transformation.

The Gateaux derivative of an expression with respect to a function $u$ in a
space $V$ is expected to be a core transformation that adds an argument. The
result is linear in a new function in $V$, the direction, bound by a `Form` at
a new position, and the rule for the class of $u$ replaces $u$ by the direction
in the derivative of each occurrence of $u$. The rules for the node classes of
the core, the derivative with respect to a function bound by an enclosing
`Form`, and the derivative of a downstream terminal through the pull-back by a
cell are not yet decided.

## Extending the core

Requires: [Passes](#passes), [Scalar type](#scalar-type), [`Absolute`](#absolute).

The core has a small, fixed set of node classes and an open set of passes. A
downstream library extends UFLx mainly by adding passes, not by adding
properties to the node classes of the core. A new pass needs one rule for each
node class of the core and does not change the core. Only a language extension
adds node classes.

The core supports five kinds of extension:

- [New questions](#new-questions), such as physical units, polynomial degree
  or linearity checks.
- [New subclasses](#new-subclasses), such as the elements and functions of
  Basix, or new coordinate domains.
- [New operations](#new-operations), such as the trace or the divergence.
- [New transformations](#new-transformations), such as the Gateaux
  derivative.
- [Language extensions](#language-extensions), such as external operators.

A property that the rules of the core depend on is part of the core, not an
extension. For example, the scalar type is part of the core, because
`Absolute` is real for a complex operand. Adding such a property changes the
core.

#### Legacy UFL

In legacy UFL, properties are stored on the node classes, and the set of properties is fixed. Every new property touches every node class, so features such as complex numbers and dual spaces were patched in across the node classes.

### Language extensions

Requires: [Leaves and subclassing](#leaves-and-subclassing), [Passes](#passes),
[Core transformations](#core-transformations), [Transpose and matrix
product](EXAMPLES.md#transpose-and-matrix-product).

A language extension adds node classes that are not leaves, in a package
outside the core. A node class of a language extension subclasses
`AbstractExpression`, never a concrete class of the core (see [Leaves and
subclassing](#leaves-and-subclassing)). Each node class of a language
extension provides, as part of the class:

- a rule for each [core transformation](#core-transformations) that keeps
  the node, which produces only node classes of the core or of the language
  extension, and
- a lowering to an equal expression of node classes of the core, if such an
  expression exists.

A pass without a rule for a node class of a language extension lowers the node
and continues, so a pass such as `units` needs no rule for the node class. For
a node class without a lowering, every pass without a rule raises
`TypeError`. A rule must agree with the lowering followed by the rules for the
lowered expression, up to mathematical equality.

For example, a language extension for local external operators adds a node
class for an operation known only by name, such as `exp`, a material law or a
neural network. The rule for `apply_derivatives` differentiates the node class.
The rule for `expand_pull_backs` pulls back the operands, because the
operation is pointwise. The node class has no lowering.

A named composition that must survive a pass, or print as written, is a node
class of a language extension whose lowering is the composition. For example,
a `Transpose` node class lowers to `component_tensor(A[k, l], (l, k))`. The
rule for `expand_pull_backs` moves the pull-back inside `Transpose`, because
the transpose is pointwise, so `Transpose` survives the pull-back.

### New questions

Requires: [Passes](#passes), [Language extensions](#language-extensions),
[`Coordinate`](#coordinate), [`ExteriorDerivative`](#exteriorderivative),
[`Components`](#components).

A downstream library asks a new question of an expression with a pass. The pass
has one rule for each node class of the core and one rule for each leaf class
of the library. The pass needs no rule for a node class of a language extension
that has a lowering (see [Language extensions](#language-extensions)).

For example, a library that checks physical units, such as the work of Habera
et al.:

1. Assigns a unit to leaves in a table of the library. For example, with
   `X = coordinate(RD(2))`, the table `{X: "m"}` states that the coordinates
   are in metres. A leaf that is not in the table, such as a literal, is
   dimensionless.
2. Writes passes, such as `units(expression, table)`, which returns the unit of
   any expression. `units` derives the unit of an expression from the units of
   the leaves, in the same way as the core derives the value shape.
   `X[0] * X[1]` is in square metres. The gradient
   `components(exterior_derivative(X[0], RD(2)), RD(2))` is in metres divided
   by metres, so the gradient is dimensionless.
   `X[0] + X[0] * X[1]` is an error of the units library.
3. Calls `units` on an expression whenever the library needs the check.

### New operations

Requires: [Node class or free function](#node-class-or-free-function), [`Tensor`
algebra with `EllipsisIndex`](#tensor-algebra-with-ellipsisindex), [Language
extensions](#language-extensions).

A downstream library builds a new operation from existing node classes, as a
free function. As in the core, the free function may read the static type of
the operands, such as a value shape, to choose the composition (see
[Node class or free function](#node-class-or-free-function)). For example,
`determinant(A)` reads the extent $n$ of `A` and returns the Leibniz formula
for an $n \times n$ matrix. So `determinant(A)` requires an `int` extent $n$
(see [Extents](#extents)). `trace(A)` reads no extent, because
`index_sum(A[k, k], (k,))`, with a derived index `k`, is the trace for every
extent. An operation that an `EllipsisIndex` expresses, such as
`scale(a, T)`, needs no value shape (see
[`EllipsisIndex`](#tensor-algebra-with-ellipsisindex)).

### New subclasses

Requires: [Leaves and subclassing](#leaves-and-subclassing), [New
operations](#new-operations), [Language extensions](#language-extensions),
[`Coordinate`](#coordinate), [An affinely mapped
triangle](EXAMPLES.md#an-affinely-mapped-triangle).

Every abstract class is a contract. The contract states what a subclass
provides, so that every core transformation handles the subclass through the
abstract class alone. A downstream library may subclass an abstract class by
fulfilling the contract, without changing the core. A contract is part of the
class, never a global registry. For example, a downstream subclass of
`AbstractCoordinateDomain`, such as
[`UnitSimplex`](EXAMPLES.md#the-unit-simplex), provides the dimension, and every
core transformation handles the coordinate function of the subclass through
`Coordinate`.

Information that a downstream library attaches to the graph goes in one of two
places:

- Information that defines the variational problem goes in the construction
  arguments of a leaf. For example, the space that a function belongs to
  changes the meaning of an expression, so the space is part of the leaf.
- Information that only checks or analyses the variational problem may stay in
  a table of the library, keyed by leaves. For example, units check an
  expression but do not change the meaning of the expression.

### New transformations

Requires: [Language extensions](#language-extensions).

A new transformation is a pass with one rule for each node class of the core
and one rule for each leaf class that the transformation rewrites. For example,
the [Gateaux derivative](#gateaux-derivative) rewrites the functions that the
derivative is taken with respect to, so the Gateaux derivative needs a rule for
the leaf classes of those functions.

## Appendix: Class hierarchy

Every object in the graph extends `Base`. Each row gives the classes of a
chapter that extend the classes in the Extends column.

<!-- generated: class-hierarchy, run design/appendices.py to update -->
| Group                       | Extends                                                                                          | Classes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| --------------------------- | ------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Domains                     | [`Base`](#the-base-object-base)                                                                  | [`AbstractDomain`](#abstractdomain)                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|                             | [`AbstractDomain`](#abstractdomain)                                                              | [`AbstractCoordinateDomain`](#abstractcoordinatedomain)                                                                                                                                                                                                                                                                                                                                                                                                                                         |
|                             | [`AbstractCoordinateDomain`](#abstractcoordinatedomain)                                          | [`RD`](#rd)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| Expressions and static type | [`Base`](#the-base-object-base)                                                                  | [`AbstractExpression`](#abstractexpression), [`Index`](#index), [`EllipsisIndex`](#ellipsisindex), [`AbstractExtent`](#extents), [`AbstractIndexValue`](#abstractindexvalue), [`GatheredIndex`](#gatheredindex)                                                                                                                                                                                                                                                                                 |
| Algebra                     | [`AbstractExpression`](#abstractexpression)                                                      | [`Integer`](#literals), [`Float`](#literals), [`ImaginaryUnit`](#literals), [`Zero`](#zero), [`Sum`](#sum), [`Product`](#product), [`Power`](#power), [`Division`](#division), [`Absolute`](#absolute), [`Conjugate`](#conjugate), [`TensorProduct`](#tensorproduct), [`SquareRoot`](#elementary-functions-very-incomplete), [`Tensor`](#tensor), [`Indexed`](#indexed), [`ComponentTensor`](#componenttensor-and-indexsum), [`IndexSum`](#componenttensor-and-indexsum), [`Inverse`](#inverse) |
| Calculus on a domain        | [`AbstractExpression`](#abstractexpression)                                                      | [`Coordinate`](#coordinate), [`ExteriorDerivative`](#exteriorderivative), [`Components`](#components), [`Integral`](#integral)                                                                                                                                                                                                                                                                                                                                                                  |
| Maps between domains        | [`Base`](#the-base-object-base)                                                                  | [`DomainMap`](#domainmap)                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|                             | [`AbstractExpression`](#abstractexpression)                                                      | [`PullBack`](#pullback)                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Topology                    | [`Base`](#the-base-object-base)                                                                  | [`AbstractCellType`](#abstractcelltype), [`AbstractTopology`](#abstracttopology)                                                                                                                                                                                                                                                                                                                                                                                                                |
| Reference cells             | [`AbstractCoordinateDomain`](#abstractcoordinatedomain), [`AbstractCellType`](#abstractcelltype) | [`AbstractReferenceCell`](#reference-cell)                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Decomposed domains          | [`AbstractDomain`](#abstractdomain)                                                              | [`DecomposedDomain`](#decomposeddomain)                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
|                             | [`Base`](#the-base-object-base)                                                                  | [`SideMap`](#sidemap)                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|                             | [`AbstractExpression`](#abstractexpression)                                                      | [`Normal`](#normal)                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Function spaces             | [`Base`](#the-base-object-base)                                                                  | [`AbstractFunctionSpace`](#abstractfunctionspace)                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|                             | [`AbstractExpression`](#abstractexpression)                                                      | [`Function`](#function)                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Variational forms           | [`AbstractFunctionSpace`](#abstractfunctionspace)                                                | [`ConjugateSpace`](#conjugatespace), [`DualSpace`](#dualspace)                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|                             | [`AbstractExpression`](#abstractexpression)                                                      | [`Form`](#form), [`Action`](#action)                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
<!-- end generated: class-hierarchy -->

## Appendix: Free functions

<!-- generated: free-functions, run design/appendices.py to update -->
| Group                       | Free function                                                            | `AbstractExpression` overload | Takes                                                                                               | Returns                                                                                      |
| --------------------------- | ------------------------------------------------------------------------ | ----------------------------- | --------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| Expressions and static type | [`index(label)`](#index)                                                 |                               | `label: str \| tuple`                                                                               | `Index`                                                                                      |
|                             | [`ellipsis_index(label)`](#ellipsisindex)                                |                               | `label: str \| tuple`                                                                               | `EllipsisIndex`                                                                              |
|                             | [`gathered_index(table, indices)`](#gatheredindex)                       |                               | `table: AbstractIndexValue`, `indices: tuple`                                                       | `GatheredIndex`                                                                              |
| Algebra                     | [`add(first, second)`](#sum)                                             | `+`                           | `first, second: AbstractExpression`                                                                 | `Sum`                                                                                        |
|                             | [`multiply(first, second)`](#product)                                    | `*`                           | `first, second: AbstractExpression`                                                                 | `Product`                                                                                    |
|                             | [`power(first, second)`](#power)                                         | `**`                          | `first, second: AbstractExpression`                                                                 | `Power`                                                                                      |
|                             | [`negate(f)`](#negation-and-subtraction)                                 | `-f`                          | `f: AbstractExpression`                                                                             | `multiply(Integer(-1), f)`                                                                   |
|                             | [`subtract(f, g)`](#negation-and-subtraction)                            | `f - g`                       | `f, g: AbstractExpression`                                                                          | `add(f, negate(g))`                                                                          |
|                             | [`divide(first, second)`](#division)                                     | `/`                           | `first, second: AbstractExpression`                                                                 | `Division`                                                                                   |
|                             | [`absolute(f)`](#absolute)                                               | `abs(f)`                      | `f: AbstractExpression`                                                                             | `Absolute`                                                                                   |
|                             | [`conjugate(f)`](#conjugate)                                             |                               | `f: AbstractExpression`                                                                             | `Conjugate`                                                                                  |
|                             | [`tensor_product(first, second)`](#tensorproduct)                        |                               | `first, second: AbstractExpression`                                                                 | `TensorProduct`                                                                              |
|                             | [`square_root(f)`](#elementary-functions-very-incomplete)                |                               | `f: AbstractExpression`                                                                             | `SquareRoot`                                                                                 |
|                             | [`tensor(entries)`](#tensor)                                             |                               | `entries: tuple[AbstractExpression, ...]`                                                           | `Tensor`                                                                                     |
|                             | [`indexed(expression, indices)`](#indexed)                               | `expression[indices]`         | `expression: AbstractExpression`, `indices: tuple`                                                  | `Indexed`                                                                                    |
|                             | [`component_tensor(expression, indices)`](#componenttensor-and-indexsum) |                               | `expression: AbstractExpression`, `indices: tuple[Index, ...]`                                      | `ComponentTensor`                                                                            |
|                             | [`index_sum(expression, indices)`](#componenttensor-and-indexsum)        |                               | `expression: AbstractExpression`, `indices: tuple[Index, ...]`                                      | `IndexSum`                                                                                   |
|                             | [`inverse(A)`](#inverse)                                                 |                               | `A: AbstractExpression`                                                                             | `Inverse`                                                                                    |
| Calculus on a domain        | [`coordinate(P)`](#coordinate)                                           |                               | `P: AbstractDomain`                                                                                 | `Coordinate`                                                                                 |
|                             | [`exterior_derivative(f, P)`](#exteriorderivative)                       |                               | `f: AbstractExpression`, `P: AbstractDomain`                                                        | `ExteriorDerivative`                                                                         |
|                             | [`components(f, P)`](#components)                                        |                               | `f: AbstractExpression`, `P: AbstractDomain`                                                        | `Components`                                                                                 |
|                             | [`volume_form(P)`](#volume-form)                                         |                               | `P: AbstractCoordinateDomain`                                                                       | `Product` of the coordinate 1-forms                                                          |
|                             | [`integral(f, P)`](#integral)                                            |                               | `f: AbstractExpression`, `P: AbstractDomain`                                                        | `Integral`                                                                                   |
| Maps between domains        | [`domain_map(P, psi, Q)`](#domainmap)                                    |                               | `P: AbstractDomain`, `psi: AbstractExpression`, `Q: AbstractDomain`                                 | `DomainMap`                                                                                  |
|                             | [`pull_back(f, D)`](#pullback)                                           |                               | `f: AbstractExpression`, `D: DomainMap \| SideMap`                                                  | `PullBack`                                                                                   |
| Topology                    | [`cell_index(T, j)`](#abstracttopology)                                  |                               | `T: AbstractTopology`, `j: int`                                                                     | `Index`                                                                                      |
| Reference cells             | [`affine_map(F, points)`](#affine-maps)                                  |                               | `F: AbstractReferenceCell`, `points: AbstractExpression`                                            | The affine map of the vertices                                                               |
| Decomposed domains          | [`decomposed_domain(pieces)`](#decomposeddomain)                         |                               | `pieces: tuple[DomainMap, ...]`                                                                     | `DecomposedDomain`                                                                           |
|                             | [`cell_domain(D, j)`](#cell-domains)                                     |                               | `D: DecomposedDomain`, `j: int`                                                                     | `DomainMap` of the cells of type `j`                                                         |
|                             | [`gathered_cell_domain(D, j, g)`](#gathered-cell-domains)                |                               | `D: DecomposedDomain`, `j: int`, `g: GatheredIndex`                                                 | `DomainMap` of the cells of type `j` selected by `g`                                         |
|                             | [`restriction(f, E, s)`](#restriction)                                   |                               | `f: AbstractExpression`, `E: DecomposedDomain`, `s: str`                                            | `pull_back(f, side_map(D, E, s))`, with $D$ the support of $f$, or $f$ if $f$ has no support |
|                             | [`side_domain(D, S)`](#restriction)                                      |                               | `D: DecomposedDomain`, `S: AbstractTopology`                                                        | The decomposed domain of the cells of $S$, with the geometry of the first side of each group |
|                             | [`side_map(D, E, s)`](#sidemap)                                          |                               | `D, E: DecomposedDomain`, `s: str`                                                                  | `SideMap`                                                                                    |
|                             | [`normal(D, E, s)`](#normal)                                             |                               | `D, E: DecomposedDomain`, `s: str`                                                                  | `Normal`                                                                                     |
| Function spaces             | [`function(V, label)`](#function)                                        |                               | `V: AbstractFunctionSpace`, `label: str \| tuple`                                                   | `Function`                                                                                   |
| Variational forms           | [`conjugate_space(V)`](#conjugatespace)                                  |                               | `V: AbstractFunctionSpace`                                                                          | `V` for a real `V`, and `ConjugateSpace(V)` for a complex `V`                                |
|                             | [`dual(S)`](#dualspace)                                                  |                               | `S: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]`                                     | `DualSpace`                                                                                  |
|                             | [`form(expression, arguments)`](#form)                                   |                               | `expression: AbstractExpression`, `arguments: tuple[Function \| tuple[Function, ...], ...]`         | `Form`                                                                                       |
|                             | [`action(f, position, g)`](#action)                                      |                               | `f: AbstractExpression`, `position: int`, `g: AbstractExpression \| tuple[AbstractExpression, ...]` | `Action`                                                                                     |
<!-- end generated: free-functions -->

A free function and an operator overload accept only the
types in the "Takes" column.

A free function whose "Returns" column names only a node class `C` returns
`C(...)` of the arguments. Any other free function returns the composition in
the "Returns" column.

## Appendix: Static type table

The degree and the argument spaces of each node class of the core. An operand
not named in the argument spaces column must have the argument spaces `()`.

| Class                                    | Edges that are expressions must have                      | Degree of the node                            | Argument spaces of the node                |
| ---------------------------------------- | --------------------------------------------------------- | --------------------------------------------- | ------------------------------------------ |
| `Integer`, `Float`, `ImaginaryUnit`      | No such edges                                             | `()`                                          | `()`                                       |
| `Zero`                                   | No such edges                                             | The stored degree                             | The stored argument spaces                 |
| `Sum`                                    | Equal degrees                                             | The degree of the operands                    | Those of the operands, which must be equal |
| `Product`                                | The degree `()`, or exactly one slot each                 | The slotwise sum (see [`Product`](#product))    | Those of the operands, concatenated        |
| `Power`, `Absolute`                      | The degree `()`                                           | `()`                                          | `()`                                       |
| `Division`                               | The degree `()` for the second operand                    | The degree of the first operand               | Those of the first operand                 |
| `Conjugate`                              | Any degree                                                | The degree of the operand                     | Those of the operand, each conjugated      |
| `TensorProduct`                          | Any degree                                                | The concatenated degrees                      | Those of the operands, concatenated        |
| `SquareRoot`, `Inverse`                  | The degree `()`                                           | `()`                                          | `()`                                       |
| `Tensor`                                 | Equal degrees                                             | The degree of the entries                     | Those of the entries, which must be equal  |
| `Indexed`, `ComponentTensor`, `IndexSum` | Any degree                                                | The degree of the expression                  | Those of the expression                    |
| `Coordinate`                             | No such edges                                             | `()`                                          | `()`                                       |
| `Normal`                                 | No such edges                                             | `()`                                          | `()`                                       |
| `Function`                               | No such edges                                             | The degree of the space                       | Those of the space                         |
| `Form`                                   | Any degree                                                | The degree of the expression                  | The space of each entry, conjugated for a `conjugate(v)` entry |
| `Action`                                 | Any degree                                                | The degree of the form                        | Those of the form without the slot         |
| `ExteriorDerivative`                     | The degree `()` or `(k,)`                                 | `(1,)` for `()`, and `(k + 1,)` for `(k,)`    | Those of the operand                       |
| `Components`                             | At least one slot                                         | The first slot lowered by 1, and removed at 0 | Those of the operand                       |
| `Integral`                               | `(n,)` with $n$ the dimension of $P$, or `()` for $n = 0$ | `()`                                          | Those of the integrand                     |
| `PullBack`                               | Any degree                                                | The degree of the expression                  | Those of the expression                    |

The value shape and the free indices of each node class are in the section
of the class. The support and the scalar type of the node classes that differ
from the default rules of [Support](#support) and [Scalar type](#scalar-type):

| Class                                            | Support of the node                | Scalar type of the node           |
| ------------------------------------------------ | ---------------------------------- | --------------------------------- |
| `ImaginaryUnit`                                  | None                               | Complex                           |
| `Zero`                                           | None                               | Real                              |
| `Absolute`                                       | The default                        | Real                              |
| `Coordinate`, `ExteriorDerivative`, `Components` | The domain $P$                     | The default                       |
| `Integral`                                       | None                               | The default                       |
| `PullBack`                                       | The parameter domain of the domain | The default                       |
| `Normal`                                         | $E$                                | Real                              |
| `Function`                                       | The domain of the space, or none   | The scalar type of the space      |
| `Form`                                           | The support of the expression      | The scalar type of the expression |
| `Action`                                         | The support of the form            | The scalar type of the form       |

## Appendix: Pass table

What each core transformation does with each node class as the operand of the
evaluated node. `expand_integrals` and `reduce_component_tensors` have no column, because
neither has a rule for each operand class (see [Core
transformations](#core-transformations)). "In" moves the evaluated node into the operands, by linearity
or commuting. "Chain" is the chain rule. "Rule" is a rule of the
transformation. "Bound" moves the evaluated node into the expression or the
form only, past the bound arguments. "Remains" keeps the
evaluated node on the operand. "Rule, Remains" is a rule for some operands of
the class and "Remains" for the others, as the table of the transformation
states. "—" cannot occur, because the static type rules out the pair.

For `expand_pull_backs`, the entry is for an operand with a support. An
operand with no support takes the rule for an expression with no support,
whatever the class of the operand.

| Class                                    | `expand_pull_backs` | `apply_derivatives` | `expand_components` | `apply_actions` |
| ---------------------------------------- | ------------------- | ------------------- | ------------------- | --------------- |
| `Integer`, `Float`, `ImaginaryUnit`      | Rule                | Rule                | —                   | —               |
| `Zero`                                   | Rule                | Rule                | Rule                | Rule            |
| `Sum`, `Tensor`                          | In                  | In                  | In                  | Rule            |
| `Product`, `TensorProduct`               | In                  | Rule                | Rule                | Rule            |
| `Division`                               | In                  | Rule                | In                  | In              |
| `Power`                                  | In                  | Chain, Remains      | —                   | —               |
| `SquareRoot`                             | In                  | Chain               | —                   | —               |
| `Absolute`, `Inverse`                    | In                  | Rule                | —                   | —               |
| `Conjugate`                              | In                  | Rule                | Rule                | In              |
| `Indexed`, `ComponentTensor`, `IndexSum` | In                  | In                  | In                  | In              |
| `Coordinate`                             | Rule                | Remains             | —                   | —               |
| `Normal`                                 | Rule, Remains       | Remains             | —                   | —               |
| `Function`                               | Remains             | Remains             | Remains             | Remains         |
| `Form`                                   | Bound               | Bound               | Bound               | Rule            |
| `Action`                                 | Bound               | Bound               | Bound               | Remains         |
| `ExteriorDerivative`                     | Rule                | Rule                | Rule, Remains       | In              |
| `Components`                             | Rule, Remains       | Remains             | Remains             | In              |
| `Integral`                               | Rule                | Rule                | —                   | In              |
| `PullBack`                               | Rule, Remains       | Rule, Remains       | Remains             | In              |

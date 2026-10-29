# Core UFLx language definition

This note describes the smallest core language behind UFLx. The note leaves
out implementation details, compiler algorithms and broad redesign ideas.
The library classes of the examples and further worked examples are in
[EXAMPLES.md](EXAMPLES.md). The note assumes knowledge of differential
geometry and functional analysis, and explains some advanced concepts of
differential geometry in more detail.

This document is written for users and developers who know domain-specific
languages for variational forms, such as UFL and ngsolve.

## Conventions

The note uses "must" for a rule that an implementation or a library has to
satisfy, "should" for a rule that an implementation or a library satisfies
unless the note states a reason not to, and "may" for a choice that the note
leaves open. A statement in the present tense, such as "the value shape is a
`tuple`", is a rule of the language. "Not yet covered" marks a decided gap, a
part of the language that the authors intend to specify later. "Open" marks a
question that the authors have not decided. A "Legacy UFL" dropdown relates a
section to legacy UFL and is not a rule.

## Scope

### What UFLx is and isn't

In one sentence:

> UFLx is a typed symbolic language for writing variational forms.

The UFLx core language, henceforth the core, should not use finite element or
numerical analysis terminology. Differential geometry and functional analysis
should be all that the core requires. This draft requires vector calculus on
domains in $R^n$ and on manifolds embedded in $R^n$, in the coordinates of
$R^n$. One gradient serves every dimension, so the statement of each concept
does not depend on the dimension.

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

```{dropdown} Legacy UFL
TSFC lowers a UFL form to GEM, which repeats much of the algebra of UFL, partly because UFL simplifies during construction. The UFLx core does not simplify during construction, so a form compiler may be able to use a lowered UFLx form directly, adding only leaves for tables and kernel inputs. This has not been confirmed with the TSFC or FFCx developers.
```

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
   written at construction.
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
[Transpose](EXAMPLES.md#transpose).

A property reads and a free function builds. Data that an object stores, or
derives from the class and the construction arguments of the object, is a
property of the object, such as `value_shape` of an expression or
`parameter_domain` of a domain. A property never constructs an object. Any
operation that constructs objects is a pure free function, such as
`coordinate(P)` or `density(D)`. A class method, where
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
`IndexSum`](#componenttensor-and-indexsum), [Derived indices](#derived-indices).

A pass is a pure free function that walks a graph after construction of the
graph. The result of a pass depends only on the arguments of the pass. A pass
returns a value, such as a unit, or a new graph. A pass builds a new graph by
constructing an object for each node that the pass replaces, and by
reconstructing every object that depends on a replaced node, up to the root
(see [The base object](#the-base-object-base)). An object depends on a
replaced node through any edge, so a pass that replaces an object also
reconstructs an [`Index`](#index) whose label contains the object, other than
a renamed name (see [Derived indices](#derived-indices)). Hence a
replacement commutes with the free functions: replacing `A` by `B` in
`transpose(A)` gives `transpose(B)`, with the derived indices of `B`. A pass
may cache results for each object, because objects are immutable.

Substitution respects scope. A binding node binds names in some of its
operands:

| Binding node | Bound names | Bound in |
| ------------ | ----------- | -------- |
| [`ComponentTensor(e, indices)`, `IndexSum(e, indices)`](#componenttensor-and-indexsum) | The indices of `indices` | `e` and `indices` |
| [`Form(e, arguments)`](#form) | The functions of `arguments` | `e` and `arguments` |
| [`DirectionalDerivative(f, x, a, h)`](#directionalderivative) | $x$ | $f$ and $x$ |
| [`GateauxDerivative(F, x, a, w)`](#gateaux-derivative) | The functions of $x$ and of $w$ | $F$, $x$ and $w$ |

A binding node is made fresh for expressions $s_1, \dots, s_n$ by renaming,
in the operands in which the names are bound, each bound name that is free in
some $s_k$, so that the binding node does not capture the name.

A substitution replaces expressions $r_1, \dots, r_n$ by $s_1, \dots, s_n$
simultaneously, so a substituted $s_k$ is not searched for the other $r_l$.
Inside a binding node, the substitution goes into the other operands as
anywhere else. For the operands in which the names are bound, the substitution
drops each pair whose $r_k$ has a bound name free in it, because every
occurrence of that $r_k$ there refers to the bound name, makes the binding node
fresh for the remaining $s_k$, and substitutes the remaining pairs.

When a binding node is made fresh for $s_1, \dots, s_n$, the new name of a
bound index `i` is the derived index
`Index(("uflx.rename", i, body, (s_1, ..., s_n)))`, or the derived
`EllipsisIndex` with the same label for a bound `EllipsisIndex`, and the new
name of a bound function `v` in a space $V$ is
`Function(V, ("uflx.rename", v, body, (s_1, ..., s_n)))`, with `body` the
`tuple` of the operands in which the names are bound. The label refers to
`body` and to each $s_k$, so the new name occurs in neither (see [Derived
indices](#derived-indices)). The bound names of a binding node are distinct,
so the new names are distinct. For example, replacing
`c` by `b[i]` in `index_sum(a[i] * c, (i,))` gives
`index_sum(a[j] * b[i], (j,))`, with `j` the renamed index, and the result has
the free index `i`. Without the renaming, the result would be
`index_sum(a[i] * b[i], (i,))`, with no free index. Replacing `v` by `w` in
`form(a, (v,))` gives `form(a, (v,))`, because the form binds `v`. Replacing
`u` by `u0` in `DirectionalDerivative(f, u, u, h)` gives
`DirectionalDerivative(f, u, u0, h)`, the derivative at `u0`, because the
point is not bound.

The pass `replace(e, mapping)` of the core substitutes the value of each key of
the `dict` `mapping` for the key in $e$, simultaneously and with the scope
above.

A pass does one thing, such as differentiating or lowering. Passes compose,
and the user of the core chooses which passes run and in what order.

### Kinds of code that produce objects

Requires: [The base object `Base`](#the-base-object-base), [Properties and free
functions](#properties-and-free-functions), [Passes](#passes).

See also: [Transpose](EXAMPLES.md#transpose),
[`Indexed`](#indexed), [New questions](#new-questions).

Three kinds of code produce objects:

| Kind                                               | When the code runs              | Inspects arguments to choose what to build                                | Validates                     |
| -------------------------------------------------- | ------------------------------- | ------------------------------------------------------------------------- | ----------------------------- |
| Constructor                                        | Whenever an object is made      | No, a constructor inspects arguments only to validate                     | Yes, and only constructors do |
| Construction free function, such as `determinant`  | While the user builds the graph | Only the static type, to build one fixed composition for each static type | No                            |
| Pass, such as differentiation, lowering or `units` | After construction              | Yes, and a pass may rewrite the graph                                     | No                            |

A construction free function fails when a constructor rejects the
arguments, for example `side_map(D, E, s)` raises `ValueError` from the
constructor of `SideMap` when no group of the topology of `E` has a side named
`s`. A pass fails with `TypeError` for a node class that
the pass has no rule for, or with an error of the pass, such as inconsistent
units.

## Domains

The function spaces and operators of UFLx are defined on a domain and need the
topological structure of the domain, which is the dimension, boundary and
decomposition into cells of the domain. The core needs no metric beyond the
Euclidean metric of the ambient coordinate domain, which gives the gradient on
a domain of lower dimension (see [`Gradient`](#gradient)), the density of an
integral (see [Jacobian and density](#jacobian-and-density)) and the unit
normal (see [`Normal`](#normal)).

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
takes a domain accepts either kind, and the section of the node class gives the
meaning of the node on a decomposed domain through the pieces of the decomposed
domain.

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
(see [`Coordinate`](#coordinate)), and the integral is the Lebesgue integral in
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
| Provides             | `value_shape`, `free_indices`, `scalar_type`, `support`, `argument_spaces` |

Every `AbstractExpression` derives the value shape, the free indices, the
scalar type, the support and the [argument spaces](#argument-spaces) of the
expression. The value shape is
a `tuple` of [extents](#extents), one for each axis. An expression with value
shape `()` is a scalar.

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

```{dropdown} Open: ragged indexing
An axis whose length depends on another index, such as the quadrature points
of each cut cell, cannot be written. A dependent extent, an `AbstractExtent`
whose label contains the other index, would cover the case. An index value
read from a table is a [`GatheredIndex`](#gatheredindex).
```

### `Index`

Requires: [The base object `Base`](#the-base-object-base), [Properties and free
functions](#properties-and-free-functions), [Extents](#extents).

See also: [Free indices](#free-indices), [Index notation](#index-notation),
[`Indexed`](#indexed), [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum), [Derived
indices](#derived-indices).

> LLM generated content.

A concrete subclass of the base object for a free index. An `Index` is not an
expression. An `Index` names an axis of a tensor, as `i` and `j` in `A[i, j]`
(see [Index notation](#index-notation)).

| `Index(label)` |                                                                                                 |
| -------------- | ----------------------------------------------------------------------------------------------- |
| Extends        | `Base` |
| Free function  | `index(label)`                                                                                  |
| Takes          | `label: str \| tuple` |
| Properties     | `label`, a `str` for a user index, or a `tuple` that contains graph objects for a derived index |
| Raises         | `TypeError` if the label is neither a `str` nor a `tuple`                                       |

An `Index` has no extent of its own. An `Index` takes the extent of the axis
that the `Index` indexes.

A user index has a `str` label, which the user chooses. Two free occurrences
of one `Index` in one scope are the same index, so an index that the user
means to be new needs a new name. The scope of an index is given in
[`ComponentTensor` and `IndexSum`](#componenttensor-and-indexsum).

A derived index has a `tuple` label, which a free function or a pass builds
from the operands, so that the index is fresh without a global counter (see
[Derived indices](#derived-indices)).

### Free indices

Requires: [`AbstractExpression`](#abstractexpression), [`Index`](#index),
[Extents](#extents).

See also: [`EllipsisIndex`](#ellipsisindex), [`Sum`](#sum),
[`Tensor`](#tensor), [`Indexed`](#indexed), [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum).

Every `AbstractExpression` derives the free indices, a `tuple` of pairs. Each
pair is an `Index` and an extent, or an [`EllipsisIndex`](#ellipsisindex) and
a `tuple` of extents. The order is the order of first occurrence, never the
order of a `set`.

Unless a node class states otherwise, the free indices of a node are derived
from the edges of the node that are expressions:

- The free indices of the edges are taken in edge order, with each `Index`
  taken once, at the first occurrence.
- An `Index` that occurs in two edges must have the same extent in both edges,
  or construction raises `ValueError`.

Summation is never implicit. `a[i] * b[i]` has the free index `i` and denotes
the component-wise product. The sum over `i` is `index_sum(a[i] * b[i], (i,))`.
The operands of a node may have different free indices, so `a[i] + b[j]` has
the free indices `i` and `j` and denotes the outer sum.

### Scalar type

Requires: [`AbstractExpression`](#abstractexpression).

See also: [Literals](#literals), [`Absolute`](#absolute), [`RealPart` and
`ImaginaryPart`](#realpart-and-imaginarypart).

The scalar type is real or complex. Unless the [static type
table](#appendix-static-type-table) states otherwise, an expression is complex
if and only if at least one edge is a complex expression.

The scalar type is the arithmetic in which the expression is evaluated, not a
statement about the value. A real `SquareRoot`, `Logarithm` or `Power` whose
value is not a real number, such as `square_root(Integer(-1))`, is undefined.

### Support

Requires: [`AbstractExpression`](#abstractexpression),
[`AbstractDomain`](#abstractdomain).

The support is an [`AbstractDomain`](#abstractdomain). The support is absent
for an expression with no support, such as a literal. Unless the [static type
table](#appendix-static-type-table) states otherwise, the support
of a node is the support of the edges of the node that are expressions with a support.
All such edges must have the same support, or construction raises `ValueError`. A
node with no such edge has no support.

### Expressions as functions

Requires: [Support](#support).

See also: [`Zero`](#zero).

An expression denotes an array-valued function on the support of the
expression. An expression with no support denotes a single array. An
expression with argument spaces other than `()` denotes a multilinear map from
the argument spaces to such functions (see [Variational
forms](#variational-forms)).

### Node class or free function

Requires: [Properties and free functions](#properties-and-free-functions),
[`AbstractExpression`](#abstractexpression), [Free indices](#free-indices),
[Scalar type](#scalar-type), [Support](#support),
[`AbstractDomain`](#abstractdomain).

See also: [Transpose](EXAMPLES.md#transpose),
[Determinant](#determinant), [`Gradient`](#gradient),
[`Division`](#division), [Elementary functions (very
incomplete)](#elementary-functions-very-incomplete),
[Tensor algebra with
`EllipsisIndex`](EXAMPLES.md#tensor-algebra-with-ellipsisindex),
[Language extensions](#language-extensions).

A concept must be a free function, not a node class, if the free function can
build the concept as one fixed composition of existing node classes for each
static type of the arguments. The static type of an expression is the value
shape, the free indices, the scalar type, the support and the argument
spaces. The static type of a domain is the stored and the derived
properties of the domain, such as the dimension of a coordinate domain or the pieces of a
decomposed domain. A free function may read the static type to choose the
composition, and never reads the structure of an argument, such as the class,
the edges or the values. Otherwise the concept gets its own node class. A
concept that is evaluated by a numerical algorithm that the user supplies also
gets its own node class, even if a fixed composition exists, because the
composition would fix the algorithm. A concept that a pass must
meet as the user wrote the concept also gets its own node class, even if a
fixed composition exists, and the section of the node class states the
composition.
`Integral` over a decomposed domain is such a concept (see
[`Integral`](#integral)).

For example, the transpose is always `component_tensor(A[k, l], (l, k))`, and
`determinant(A)` is the Leibniz formula for each extent $n$ of a square matrix
`A`, so both are free functions.

A free function is part of the core only if a definition or a rule of the
core uses the free function, such as `density` in [`Integral` over a
decomposed domain](#integral-over-a-decomposed-domain) and `cell_domain` in
[Restriction](#restriction), or if a Python operator calls the free function
as an overload, such as `matrix_product` for `@` (see [Operator
overloads](#operator-overloads)). `transpose` is in the core because the rule
for `Inverse` of [Applying derivatives](#applying-derivatives) uses it. Every
other free function belongs to a library, such as the operations of [Tensor
algebra with `EllipsisIndex`](EXAMPLES.md#tensor-algebra-with-ellipsisindex).

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
| `Zero`                                           | No composition gives a zero with free indices without edges to other expressions                              |
| `Sum`, `Product`, `Power`, `Absolute`             | Operations of the scalar field, which nothing else composes                                                 |
| `Conjugate`                                      | An operation of the scalar field, and `Form` meets `conjugate(v)` of a function of any value shape as written |
| `RealPart`                                       | `(f + conjugate(f)) / Integer(2)` has the value of the real part but is complex, and the real part is real    |
| `ImaginaryPart`                                  | `real_part(ImaginaryUnit() * conjugate(f))` multiplies the real part by zero, which is NaN for an infinite real part |
| `Division`                                       | `multiply(f, power(g, Integer(-1)))` rounds twice, and the quotient the user wrote rounds once                |
| `SquareRoot`                                     | IEEE 754 rounds the square root correctly and the power not, so `power(f, Integer(1) / Integer(2))` differs   |
| `Logarithm`                                      | An elementary function, which nothing else composes                                                           |
| `Inverse`                                        | No one composition covers every $k \le n$, and a pass chooses the algorithm, as for `Division`                |
| `Tensor`, `Indexed`, `ComponentTensor`, `IndexSum` | The primitives of index notation                                                                            |
| `Coordinate`, `Function`, `Normal`               | Terminals: the coordinate function of a domain, a member of a space, the outward normal of a cell on a side   |
| `Gradient`                                       | Nothing composes to differentiate                                                                             |
| `DirectionalDerivative`                          | Nothing composes to differentiate with respect to an expression                                               |
| `GateauxDerivative`                              | The value depends on the forms in the expression, which a free function cannot read                           |
| `Integral`                                       | Nothing composes to integrate over a coordinate domain, and over a decomposed domain a pass must meet the integral as written |
| `PullBack`                                       | Nothing composes to move an expression between domains                                                        |
| `Form`, `Action`                                 | Binding arguments and giving an expression to an argument are not compositions                                |

## Algebra

### Literals

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices), [Scalar type](#scalar-type), [Support](#support).

See also: [`Division`](#division).

A literal is a number written down in an expression, such as
`Integer(2) * u`. Every literal has value shape
`()`, no free indices and no support.

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
indices](#free-indices), [Support](#support).

See also: [`Indexed`](#indexed).

A concrete subclass of `AbstractExpression` for the zero of any static type.

| `Zero(value_shape, free_indices, argument_spaces)` |                                                                                                                                |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| Extends                                            | `AbstractExpression` |
| Properties                                         | `value_shape`, `free_indices` and `argument_spaces`, as in [Extents](#extents), [Free indices](#free-indices) and [Argument spaces](#argument-spaces) |
| Static type                                        | The stored properties, with no support                                                                                         |
| Raises                                             | `TypeError` for an argument of another type, and `ValueError` for a negative extent or an index that occurs twice             |

A `Zero` is real and has no support, and is also the zero of a complex
expression and of an expression with a support (see
[Evaluation](#evaluation)). The derivative of an expression with no support is
a `Zero`. For example, `Zero((3,), (), ())` is the
zero array of value shape `(3,)`, and `Zero((), (), (V,))` is the zero linear
form on a space `V`.

### Arithmetic

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices).

See also: [Tensor arithmetic](#tensor-arithmetic).

`Sum`, `Product`, `Division` and `Power` are concrete subclasses of
`AbstractExpression` that each store two expressions. The node classes are
scalar. Both operands must have value shape `()`, or construction raises
`ValueError`, and the result has value shape `()`. Free indices pass through
arithmetic (see [Free indices](#free-indices)). On scalars, the free functions
`add`, `multiply` and `divide` construct `Sum`, `Product` and `Division`. On
tensors, the free functions build compositions of the node classes (see
[Tensor arithmetic](#tensor-arithmetic)), so `A + B` and `a * T` are tensor
arithmetic.

| Class      | Free function | Overload |
| ---------- | ------------- | -------- |
| `Sum`      | `add`         | `+`      |
| `Product`  | `multiply`    | `*`      |
| `Division` | `divide`      | `/`      |
| `Power`    | `power`       | `**`     |

### `Sum`

Requires: [Arithmetic](#arithmetic).

A concrete subclass of `AbstractExpression` for the sum of two scalars.

| `Sum(first, second)` |                                             |
| -------------------- | ------------------------------------------- |
| Extends              | `AbstractExpression` |
| Free function        | [`add`](#tensor-arithmetic), with the overload `+` |
| Takes                | `first, second: AbstractExpression` |
| Properties           | `first` and `second` |
| Requires             | The operands of [Arithmetic](#arithmetic)   |
| Raises               | `ValueError` if a requirement fails         |

### `Product`

Requires: [Arithmetic](#arithmetic), [Free indices](#free-indices),
[Support](#support).

See also: [`ComponentTensor` and `IndexSum`](#componenttensor-and-indexsum).

A concrete subclass of `AbstractExpression` for the product of two scalars.

| `Product(first, second)` |                                                  |
| ------------------------ | ------------------------------------------------ |
| Extends                  | `AbstractExpression` |
| Free function            | [`multiply`](#tensor-arithmetic), with the overload `*` |
| Takes                    | `first, second: AbstractExpression` |
| Properties               | `first` and `second` |
| Requires                 | The operands of [Arithmetic](#arithmetic)        |
| Raises                   | `ValueError` if a requirement fails              |

Equality is structural, so `Product` never reorders the operands, and a later
pass may apply commutativity.

### `Power`

Requires: [Arithmetic](#arithmetic), [Support](#support).

A concrete subclass of `AbstractExpression` for a scalar raised to a scalar
power.

| `Power(first, second)` |                                                                             |
| ---------------------- | --------------------------------------------------------------------------- |
| Extends                | `AbstractExpression` |
| Free function          | `power(first, second)`, with the overload `**`                              |
| Takes                  | `first, second: AbstractExpression` |
| Properties             | `first` and `second` |
| Requires               | The operands of [Arithmetic](#arithmetic)                                   |
| Raises                 | `ValueError` if a requirement fails                                         |

The exponent may have a support. The derivative of a power $a^b$ whose
exponent has a support is $a^b (\log(a) \, \mathrm{d}b + b \, \mathrm{d}a / a)$,
with `Logarithm` (see [Elementary
functions](#elementary-functions-very-incomplete)).

### Negation and subtraction

Requires: [`Sum`](#sum), [`Product`](#product), [Literals](#literals).

Negation and subtraction are free functions, not node classes:

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `negate(f)` | `-f` | `f: AbstractExpression` | `multiply(Integer(-1), f)` |
| `subtract(f, g)` | `f - g` | `f, g: AbstractExpression` | `add(f, negate(g))` |

Both compositions are exact in floating-point arithmetic, so neither
composition changes the arithmetic that the user wrote. On tensors, `negate`
and `subtract` follow `multiply` and `add` (see [Tensor
arithmetic](#tensor-arithmetic)).

### `Division`

Requires: [Arithmetic](#arithmetic), [`Product`](#product), [`Power`](#power),
[Literals](#literals).

A concrete subclass of `AbstractExpression` for the quotient of two
scalars.

| `Division(first, second)` |                                                                                        |
| ------------------------- | -------------------------------------------------------------------------------------- |
| Extends                   | `AbstractExpression` |
| Free function             | [`divide`](#tensor-arithmetic), with the overload `/` |
| Takes                     | `first, second: AbstractExpression` |
| Properties                | `first` and `second` |
| Requires                  | The operands of [Arithmetic](#arithmetic)                                              |
| Raises                    | `ValueError` if a requirement fails                                                    |

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
| Requires      | A scalar $f$                              |
| Static type   | Value shape `()`, real                    |
| Raises        | `ValueError` if a requirement fails       |

### `Conjugate`

Requires: [`AbstractExpression`](#abstractexpression), [Scalar
type](#scalar-type).

See also: [Variational forms](#variational-forms).

A concrete subclass of `AbstractExpression` for the complex conjugate of an
expression of any value shape, entry by entry.

| `Conjugate(f)` |                                                                                             |
| -------------- | ------------------------------------------------------------------------------------------- |
| Extends        | `AbstractExpression` |
| Free function  | `conjugate(f)`                                                                              |
| Takes          | `f: AbstractExpression` |
| Properties     | `operand` $f$ |
| Static type    | The value shape and the free indices of $f$, and the argument spaces of $f$ with each space conjugated (see [`ConjugateSpace`](#conjugatespace)) |

A sesquilinear variational form, such as the $L^2$ inner product
$\int u \bar{v}$, needs `Conjugate`. `Conjugate` is the one operation of the
scalar field that takes an operand of any value shape. A form binds an
argument as `conjugate(v)` (see [`Form`](#form)), and `v` may be a vector,
such as the argument of $\int \operatorname{curl} u \cdot \overline{\operatorname{curl} v}$.
`conjugate(f)[e]` has the value of `conjugate(f[e])`, so the entries of a
tensor are conjugated as scalars.



### `RealPart` and `ImaginaryPart`

Requires: [`AbstractExpression`](#abstractexpression), [Scalar
type](#scalar-type), [`Conjugate`](#conjugate).

`RealPart` and `ImaginaryPart` are concrete subclasses of `AbstractExpression`
for the real part and the imaginary part of a scalar.

| `RealPart(f)`, `ImaginaryPart(f)` |                                       |
| --------------------------------- | ------------------------------------- |
| Extends                           | `AbstractExpression` |
| Free function                     | [`real_part`](#tensor-arithmetic), [`imaginary_part`](#tensor-arithmetic) |
| Takes                             | `f: AbstractExpression` |
| Properties                        | `operand` $f$ |
| Requires                          | A scalar $f$                          |
| Static type                       | Value shape `()`, real                |
| Raises                            | `ValueError` if a requirement fails   |

Apart from `Absolute`, `RealPart` and `ImaginaryPart` are the only node
classes that take a complex operand to a real expression. A real quantity of
a complex field, such as $\operatorname{Re} \int u \bar{v}$, needs
`RealPart`. `(f + conjugate(f)) / Integer(2)` has the value of the real part,
but is complex by the rule of [Scalar type](#scalar-type).

The free functions `real_part` and `imaginary_part` take a tensor and apply
the node class entry by entry (see [Tensor arithmetic](#tensor-arithmetic)).

The operand may be real. Construction does not simplify, so `real_part(f)`
of a real $f$ is not $f$, and `imaginary_part(f)` of a real $f$ is not
`Zero`.

Examples, with `one, two = Integer(1), Integer(2)`:

- `real_part(one + ImaginaryUnit() * two)` is real, and a later pass
  evaluates the expression to 1.
- `imaginary_part(one + ImaginaryUnit() * two)` is real, and a later pass
  evaluates the expression to 2.

### Elementary functions (very incomplete)

Requires: [`AbstractExpression`](#abstractexpression), [`Power`](#power),
[`Division`](#division), [Literals](#literals).

The elementary functions so far are `SquareRoot` and `Logarithm`, concrete
subclasses of `AbstractExpression` for the square root and the natural
logarithm of a scalar.

| `SquareRoot(f)`, `Logarithm(f)` |                                     |
| ------------------------------- | ----------------------------------- |
| Extends         | `AbstractExpression` |
| Free function   | `square_root(f)` and `logarithm(f)` |
| Takes           | `f: AbstractExpression` |
| Properties      | `operand` $f$ |
| Requires        | A scalar $f$                        |
| Static type     | Value shape `()`                    |
| Raises          | `ValueError` if a requirement fails |

For a complex operand, `Logarithm` is the principal branch, with the imaginary
part in $(-\pi, \pi]$, and `SquareRoot` is the principal square root.

~~~{dropdown} Why `SquareRoot` is a node class
`SquareRoot` is not the composition `power(f, Integer(1) / Integer(2))`,
because IEEE 754 requires a correctly rounded square root and not a correctly
rounded power, so the result of `power` depends on the math library.
~~~

### `Tensor`

Requires: [`AbstractExpression`](#abstractexpression).

See also: [`Coordinate`](#coordinate).

A concrete subclass of `AbstractExpression` that stacks expressions along a
new first axis.

| `Tensor(entries)` |                                                                                                           |
| ----------------- | --------------------------------------------------------------------------------------------------------- |
| Extends           | `AbstractExpression` |
| Free function     | `tensor(entries)`                                                                                         |
| Takes             | `entries: tuple[AbstractExpression, ...]` |
| Properties        | `entries`, a non-empty `tuple` of `AbstractExpression`s with the same value shape `s`                     |
| Static type       | Value shape `(n,) + s`, with $n$ the number of entries                                                    |
| Raises            | `ValueError` if a requirement fails                                                                       |

A matrix is a `Tensor` of `Tensor`s.

A value shape may contain the extent 0, for example the value shape `(0,)` of
the coordinate function of `RD(0)`. `Tensor` cannot build an expression with
the extent 0, because an empty `tuple` has no entry to give `s`. An array with
an extent 0 has no entries, so the array is the `Zero` of that value shape.

Examples, with `zero, one, two = Integer(0), Integer(1), Integer(2)`:

- With `r0 = Float(0.0)` and `r1 = Float(1.0)`,
  `tensor((tensor((r1, r0)), tensor((r0, r1))))` has `Float` entries, so the
  matrix is not equal to `tensor((tensor((one, zero)), tensor((zero, one))))`.
- `tensor((one, tensor((one, two))))` raises `ValueError`, because the value
  shapes of the entries differ.

### Index notation

Requires: [`Tensor`](#tensor), [`Index`](#index), [Free
indices](#free-indices), [Node class or free
function](#node-class-or-free-function).

See also: [Transpose](EXAMPLES.md#transpose).

> LLM generated content.

Index notation writes tensor algebra with four node classes. `Tensor` builds a
tensor from scalars, `Indexed` selects one component, and `ComponentTensor`
and `IndexSum` bind the free indices of a scalar, as axes or as a sum. Every
other tensor operation is a fixed composition of the four, so `transpose`,
`dot` and `inner` are free functions, not node classes (see [Node class or
free function](#node-class-or-free-function)).

With `i, j = index("i"), index("j")` and a matrix `A` of value shape `(2, 3)`:

- `A[i, j]` is a scalar with the free indices `i`, of extent 2, and `j`, of
  extent 3. The scalar denotes one value for each value of `i` and `j`.
- `component_tensor(A[i, j], (j, i))` binds `j` and `i` as axes, in that
  order, and is the transpose of `A`, of value shape `(3, 2)`.
- `index_sum(A[i, j] * A[i, j], (i, j))` sums over `i` and `j`, and is the sum
  of the squares of the entries of `A`, with no free index.

The sections below give the rules. An [`EllipsisIndex`](#ellipsisindex)
covers any number of axes, so one composition serves tensors of every rank. A
[`GatheredIndex`](#gatheredindex) reads an index value from a table that the
solver gives, such as the cell of each facet. A free function or a pass that
needs a fresh index derives the index from the operands (see [Derived
indices](#derived-indices)).

### `Indexed`

Requires: [`AbstractExpression`](#abstractexpression), [`Index`](#index),
[Extents](#extents), [Free indices](#free-indices), [Support](#support),
[Index notation](#index-notation).

See also: [`EllipsisIndex`](#ellipsisindex),
[`GatheredIndex`](#gatheredindex), [`Coordinate`](#coordinate),
[`Gradient`](#gradient).

> LLM generated content.

A concrete subclass of `AbstractExpression` that selects one component of a
tensor-valued expression. `X[0]` and `A[0, 1]` are `indexed(X, (0,))` and
`indexed(A, (0, 1))`, through the free function `component` below.

| `Indexed(expression, indices)` |                                                                                   |
| ------------------------------ | --------------------------------------------------------------------------------- |
| Extends                        | `AbstractExpression` |
| Free function                  | `indexed(expression, indices)`                                                    |
| Takes                          | `expression: AbstractExpression`, `indices: tuple` |
| Properties                     | `expression`, and `indices`, a non-empty `tuple` of entries of the kinds below    |
| Static type                    | Value shape `()`, and the free indices below                                      |
| Raises                         | `TypeError`, `ValueError` and `IndexError`, below                                 |

Indexing is full indexing. Each entry covers one or more axes of the
expression, and the entries together cover every axis exactly once. With
`(n_0, ..., n_{r-1})` the value shape of the expression, an entry is one of
the following kinds:

| Entry                                | Covers                                | Rule                                                                                               | Raises       |
| ------------------------------------ | ------------------------------------- | -------------------------------------------------------------------------------------------------- | ------------ |
| An `int` $k$                         | One axis, of extent $n_j$             | $n_j$ is an `int`, and $0 \le k < n_j$                                                             | `IndexError` |
| An [`Index`](#index)                 | One axis, of extent $n_j$             | The `Index` takes the extent $n_j$                                                                 |              |
| A [`GatheredIndex`](#gatheredindex)  | One axis, of extent $n_j$             | The extent of the `GatheredIndex` equals $n_j$                                                     | `IndexError` |
| An [`EllipsisIndex`](#ellipsisindex) | Every axis that no other entry covers | At most one per `Indexed`. The `EllipsisIndex` takes the `tuple` of the extents of the covered axes, which may be empty | `ValueError` |

The free indices of the `Indexed` are the free indices of the expression,
followed by each `Index` and `EllipsisIndex` entry with the extent that the
entry takes, followed by the free indices of each `GatheredIndex` entry, with
each index taken once.

An `Index` may occur in two or more entries, or may already be free in the
expression. All the occurrences are the same index and must have the same
extent, so the indexed axes take the same value. `A[i, i]` is the diagonal
entry of a matrix `A`, and `index_sum(A[i, i], (i,))` is the trace.

Construction also raises `TypeError` for an entry of another kind,
`ValueError` for no entries or for an index with two extents, in two entries
or in an entry and the expression, and `IndexError` if the entries do not
cover exactly the rank. An expression with value shape `()` has rank 0.

Examples:

- For a matrix `A` of value shape `(2, 2)`, `A[1, 0]` has value shape `()`,
  and `indexed(A, (1,))` raises `IndexError`.
- `tensor((Integer(1) - X[0], X[0]))` has value shape `(2,)` and the support
  `RD(2)`.

#### Partial indexing

The overload `[]` calls the free function `component`, which also takes fewer
entries than the rank and indexes the leading axes:

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `component(f, entries)` | `[]` | `f: AbstractExpression`, `entries: tuple` | `component_tensor(indexed(f, entries + K), K)` if the entries are fewer than the rank of $f$ and contain no `EllipsisIndex`, otherwise `indexed(f, entries)` |

`K` is the `tuple` of the derived indices
`Index(("uflx.component", f, entries, m))`, one for each axis that the entries
do not cover (see [Derived indices](#derived-indices)). `component` reads only
the rank of $f$ and the number of entries. `Indexed` validates the entries, so
too many entries raise `IndexError` from `Indexed`.

A free function of the core calls `indexed`, never `component`, so a
composition in this document that writes `f[K]` is `indexed(f, K)`, and an
operand of the wrong rank raises from `Indexed`.

For the matrix `A` of the examples, `A[1]` is row 1 of `A`, the composition
`component_tensor(indexed(A, (1, k)), (k,))` of value shape `(2,)`, with `k`
the derived index.

### `ComponentTensor` and `IndexSum`

Requires: [`Index`](#index), [Free indices](#free-indices),
[`Indexed`](#indexed).

See also: [`EllipsisIndex`](#ellipsisindex), [Passes](#passes).

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

Scope is lexical. An occurrence of an `Index` refers to the innermost
`ComponentTensor` or `IndexSum` that encloses the occurrence and binds the
index. An occurrence that no such binder encloses is free. A bound index is
not free in the binding node, so the binding is closed, and a binder may bind
an index that a binder inside the expression of the binder also binds. For
example, with
`p = component_tensor(y[gathered_index(cells, (c, w)), k], (c, w, k))`, `p`
has no free index, and `p[c, 0, k]` has the free indices `c` and `k`, which
refer to nothing inside `p`.

A pass that substitutes into the expression of a binder renames the bound
indices that would capture a free index of the substituted expression (see
[Passes](#passes)).

Examples, with `i, j = index("i"), index("j")`, two vectors `a` and `b` of value
shape `(3,)`, and a matrix `A` of value shape `(2, 2)`:

- `index_sum(a[i] * b[i], (i,))` is the inner product of `a` and `b`, with no
  free indices.
- `index_sum(A[i, i], (i,))` is the trace of `A`, and
  `component_tensor(A[i, i], (i,))` is the diagonal of `A`, for every extent of
  `A`, including an `AbstractExtent`.

### `EllipsisIndex`

Requires: [`Index`](#index), [`Indexed`](#indexed), [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum).

See also: [Tensor algebra with
`EllipsisIndex`](EXAMPLES.md#tensor-algebra-with-ellipsisindex).

> LLM generated content.

A concrete subclass of the base object for an index that covers any number of
axes. An `EllipsisIndex` is not an expression. With an `EllipsisIndex`, one
composition serves tensors of every rank. For example,
`component_tensor(a * T[r], (r,))` scales a tensor `T` of any rank by a scalar
`a`, with the value of `a * T` (see [Tensor arithmetic](#tensor-arithmetic)).
An `EllipsisIndex` entry does not reduce against a binder with one index for
each axis (see [Derived indices](#derived-indices)).

| `EllipsisIndex(label)` |                                                           |
| ---------------------- | --------------------------------------------------------- |
| Extends                | `Base` |
| Free function          | `ellipsis_index(label)`                                   |
| Takes                  | `label: str \| tuple` |
| Properties             | `label`, with the rules of the label of an `Index`        |
| Raises                 | `TypeError` if the label is neither a `str` nor a `tuple` |

As an entry of an `Indexed`, an `EllipsisIndex` covers every axis that the
other entries do not cover, so an `EllipsisIndex` may cover no axis. The
extent of an `EllipsisIndex` is the `tuple` of the extents of the covered
axes. A `ComponentTensor` that binds an `EllipsisIndex` gets one axis for
each extent of the `EllipsisIndex`. In every other rule of index notation, an
`EllipsisIndex` behaves as an `Index`.

### `AbstractIndexValue`

Requires: [The base object `Base`](#the-base-object-base), [Extents](#extents).

See also: [`GatheredIndex`](#gatheredindex), [Sides of a
topology](#sides-of-a-topology).

An abstract subclass of the base object for a table of index values that the
solver gives, such as the cell on one side of each facet, or the local number
of the facet of a cell in a facet integral. The values are fixed but not known
during construction. An `AbstractIndexValue` is not an expression and not a
free index.

| `AbstractIndexValue` |                                                        |
| -------------------- | ------------------------------------------------------ |
| Extends              | `Base` |
| Provides             | `extent`, an extent, and `shape`, a `tuple` of extents |
| Defined by           | A subclass outside the core                            |

The table has the shape `shape`, and each value is an `int` $v$ with
$0 \le v <$ `extent`. An `AbstractIndexValue` with the shape `()` is one
value. An expression reads the table only through a
[`GatheredIndex`](#gatheredindex), so one value is read as
`gathered_index(v, ())`.

### `GatheredIndex`

Requires: [`AbstractIndexValue`](#abstractindexvalue), [`Index`](#index),
[Extents](#extents), [`Indexed`](#indexed).

See also: [Restriction](#restriction), [Gathered cell
domains](#gathered-cell-domains).

A concrete subclass of the base object for an index whose value is read from
a table. A `GatheredIndex` is not an expression and not a free index. With
`cells` the table of the cell of each facet, `a[gathered_index(cells, (e,))]`
is the value of `a` on the cell of facet `e`, with the free index `e`.

| `GatheredIndex(table, indices)` |                                                                                                                                    |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| Extends                         | `Base` |
| Free function                   | `gathered_index(table, indices)`                                                                                                   |
| Takes                           | `table: AbstractIndexValue`, `indices: tuple` |
| Properties                      | `table`, an `AbstractIndexValue`, and `indices`, an `int`, an `Index` or a `GatheredIndex` for each axis of the shape of the table |
| Derived                         | `extent` and `free_indices`, below |
| Raises                          | The errors of the entries of an [`Indexed`](#indexed)                                                                              |

The value of a `GatheredIndex` is the entry of the table at the indices. The
indices follow the rules and raise the errors of the entries of an
[`Indexed`](#indexed), without an `EllipsisIndex` entry. A `GatheredIndex`
entry reads through another table, so gathers nest. The extent of the
`GatheredIndex` is the extent of the table. The free indices are the `Index`
entries, each with the extent of the axis, followed by the free indices of
the `GatheredIndex` entries, with each `Index` taken once.

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

### Derived indices

Requires: [`Index`](#index), [`EllipsisIndex`](#ellipsisindex),
[Passes](#passes), [Node class or free function](#node-class-or-free-function).

See also: [Transpose](EXAMPLES.md#transpose), [New operations](#new-operations).

> LLM generated content.

A free function that needs fresh indices derives the indices from the
operands. The label of a derived index is a `tuple` of the identifier of the
free function, the arguments of the free function and a number. The
identifier is the public name of the free function, qualified by the package
name, such as `"uflx.transpose"` for the core or `"library.outer"` for a
library. For example, `transpose(A)` (see [Tensor
arithmetic](#tensor-arithmetic)) uses `Index(("uflx.transpose", A, 0))` and
`Index(("uflx.transpose", A, 1))`.
A derived index is fresh with respect to the operands, so the free function
does not inspect the operands to find a fresh index. A global counter is not
allowed, because a counter makes construction depend on global state.

A pass that needs fresh indices derives the indices in the same way, with the
identifier of the pass, the arguments of the rewritten node and a number. The
renamed index of [Passes](#passes) has the label
`("uflx.rename", i, body, (s_1, ..., s_n))`.

A label refers only to objects that the index is fresh for, and never to an
object that contains the index. A replacement rebuilds a derived index whose
label contains a replaced object, so that the replacement commutes with the
free functions. A renamed name keeps its label. The substitution that renames
a name continues into the operands that contain the new name, and rebuilding
the label there would substitute into $s_1, \dots, s_n$ again, which may
contain the binding node, so the substitution would not terminate.

The edges of an index to the objects in the label of the index identify the
index and are not operands.

A free function or a pass of the core derives one `Index` for each axis, and
never an `EllipsisIndex`. `component_tensor(e, K)[J]` is $e$ with each index
of `K` replaced by the entry of `J` for the axis of the index. An
`EllipsisIndex` of `K` may be replaced by the entries of the axes that the
`EllipsisIndex` covers, but an `EllipsisIndex` entry of `J` cannot be split
over the indices of `K`. So a composition with one `Index` for each axis
reduces against every binder.

### Tensor arithmetic

Requires: [Arithmetic](#arithmetic), [Negation and
subtraction](#negation-and-subtraction), [`Conjugate`](#conjugate),
[`RealPart` and `ImaginaryPart`](#realpart-and-imaginarypart), [`Indexed`](#indexed),
[`ComponentTensor` and `IndexSum`](#componenttensor-and-indexsum), [Derived
indices](#derived-indices), [Node class or free
function](#node-class-or-free-function).

See also: [Operator overloads](#operator-overloads), [Tensor algebra with
`EllipsisIndex`](EXAMPLES.md#tensor-algebra-with-ellipsisindex).

The node classes of [Arithmetic](#arithmetic) are scalar. The free functions
`add`, `multiply`, `divide`, `real_part` and `imaginary_part` read the ranks
of the operands, and on tensors build the fixed compositions below. The free function `matrix_product`
contracts one axis, and `transpose` swaps the two axes of a matrix. `add`, `multiply`, `divide` and `matrix_product` are the
overloads `+`, `*`, `/` and `@` (see [Operator overloads](#operator-overloads)),
so tensor arithmetic is written as `A + B`, `a * T`, `T / a` and `A @ x`.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `add(first, second)` | `+` | `first, second: AbstractExpression` | `component_tensor(first[K] + second[K], K)` if both operands have rank 1 or more, otherwise `Sum(first, second)` |
| `multiply(first, second)` | `*` | `first, second: AbstractExpression` | `component_tensor(first * second[K], K)` if only `first` is a scalar, `component_tensor(first[K] * second, K)` if only `second` is a scalar, otherwise `Product(first, second)` |
| `divide(first, second)` | `/` | `first, second: AbstractExpression` | `component_tensor(first[K] / second, K)` if `first` has rank 1 or more and `second` is a scalar, otherwise `Division(first, second)` |
| `real_part(f)` | | `f: AbstractExpression` | `component_tensor(RealPart(f[K]), K)` if $f$ has rank 1 or more, otherwise `RealPart(f)` |
| `imaginary_part(f)` | | `f: AbstractExpression` | `component_tensor(ImaginaryPart(f[K]), K)` if $f$ has rank 1 or more, otherwise `ImaginaryPart(f)` |
| `matrix_product(first, second)` | `@` | `first, second: AbstractExpression` | The contraction of the last axis of `first` with the first axis of `second`, below |
| `transpose(A)` | | `A: AbstractExpression` | `component_tensor(A[k, l], (l, k))`, with `(k, l)` = `K`, for $A$ of rank 2 |

In each composition, `K` is the `tuple` of the derived indices
`Index((name, first, second, m))`, or `Index((name, f, m))`, for
$m = 0, \dots, n - 1$, with `name` the identifier of the free function, such
as `"uflx.add"`, and $n$ the rank of the first operand of rank 1 or more (see
[Derived indices](#derived-indices)). The operands of `+`, `*` and `/` inside
the composition are scalars, so the composition constructs `Sum`, `Product`
and `Division`. The indices of `K` are fresh with respect to the operands, so
the free indices of the operands are free in the result.

The free functions read only the ranks, and never validate (see [The base
object `Base`](#the-base-object-base)). A constructor in the composition
raises for operands that do not fit. There is no broadcasting, and no product
of two tensors other than `@`:

- `A + B` for ranks that differ raises `IndexError` from `Indexed`, and for
  extents that differ raises `ValueError` from `Sum`, because an index of `K`
  has two extents.
- `a + T` for a scalar `a` and a tensor `T` raises `ValueError` from `Sum`.
- `A * B` for two tensors raises `ValueError` from `Product`, and `a / T`
  raises `ValueError` from `Division`.

`negate` and `subtract` compose `multiply` and `add` (see [Negation and
subtraction](#negation-and-subtraction)), so `-T` is
`component_tensor(Integer(-1) * T[K], K)` and `A - B` is
`add(A, negate(B))` on tensors.

A free function of a scalar operation takes a tensor only if the operation
has no other meaning on a tensor. `conjugate`, `real_part` and
`imaginary_part` have one meaning on a tensor, entry by entry. `Conjugate` itself takes
any value shape (see [`Conjugate`](#conjugate)). The following stay scalar,
and raise `ValueError` for a tensor, because the tensor has a second meaning:

| Free function | Second meaning on a tensor |
| ------------- | -------------------------- |
| `power`, `**` | The matrix power |
| `absolute`, `abs` | The norm of a vector |
| `square_root` | The matrix square root |
| `logarithm` | The matrix logarithm |

An elementary function added later, such as the exponential, stays scalar for
the same reason, because the matrix exponential is a second meaning. The
entry-by-entry application of a scalar operation `op` is written as
`component_tensor(op(T[K]), K)`, with one index in `K` for each axis of `T`.

The graph of `A + B` is the composition, with no node for the sum of tensors.
A language extension adds such a node if a pass must meet `A + B` as written
(see [Language extensions](#language-extensions)).

`matrix_product(A, B)` takes operands `A` and `B` each of rank 1 or 2. With
the derived indices $i$, $k$ and $j$ =
`Index(("uflx.matrix_product", A, B, m))` for $m = 0, 1, 2$, the free
function reads the ranks and builds:

| Value shapes of `A` and `B` | `matrix_product(A, B)` |
| --------------------------- | ---------------------- |
| `(n,)` and `(n,)` | `index_sum(A[k] * B[k], (k,))` |
| `(m, n)` and `(n,)` | `component_tensor(index_sum(A[i, k] * B[k], (k,)), (i,))` |
| `(n,)` and `(n, p)` | `component_tensor(index_sum(A[k] * B[k, j], (k,)), (j,))` |
| `(m, n)` and `(n, p)` | `component_tensor(index_sum(A[i, k] * B[k, j], (k,)), (i, j))` |

`matrix_product` indexes an operand of rank 2 with two indices and every other
operand with one index, so `Indexed` raises `IndexError` for an operand of
rank 0 or of rank 3 or more. Unequal inner extents raise `ValueError` from
`Product`, because $k$ has two extents. Neither operand is conjugated, so
`x @ y` of two complex vectors is not the inner product.

For ranks 1 and 2, `matrix_product` agrees with `@` of NumPy. NumPy treats an
operand of rank 3 or more as a stack of matrices, so the core rejects such an
operand and does not give `@` a second meaning. A contraction of higher rank
is written in index notation, or with a free function of a library, such as
`dot` (see [Tensor algebra with
`EllipsisIndex`](EXAMPLES.md#tensor-algebra-with-ellipsisindex)).

Examples, with `i = index("i")`, a scalar `a`, vectors `x` and `y` of value
shape `(3,)`, a matrix `A` of value shape `(2, 3)` and a matrix `B` of value
shape `(3, 3)`:

- `A + A` has value shape `(2, 3)`, and `A + B` raises `ValueError`.
- `a * x - x / a` has value shape `(3,)`.
- `x[i] * y` has value shape `(3,)` and the free index `i`.
- `A @ x` has value shape `(2,)`, `x @ y` has value shape `()`, and `A @ B`
  has value shape `(2, 3)`.
- `B @ A` raises `ValueError`, because the inner extents 3 and 2 differ, and
  `A * B` raises `ValueError`.

### `Inverse`

Requires: [`AbstractExpression`](#abstractexpression),
[`Division`](#division), [Elementary
functions](#elementary-functions-very-incomplete), [Tensor
arithmetic](#tensor-arithmetic).

A concrete subclass of `AbstractExpression` for the inverse of a square
matrix, or the pseudo-inverse of a matrix with more rows than columns.

| `Inverse(A)`  |                                                                                         |
| ------------- | --------------------------------------------------------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `inverse(A)`                                                                            |
| Takes         | `A: AbstractExpression` |
| Properties    | `operand` $A$ |
| Requires      | $A$ with value shape `(n, k)`, $k \le n$, and $A$ real for $k < n$                      |
| Static type   | Value shape `(k, n)`                                                                    |
| Raises        | `ValueError` if a requirement fails                                                     |

`Inverse(A)` is the pseudo-inverse $A^+ = (A^T A)^{-1} A^T$ of a matrix $A$
of rank $k$, the left inverse of $A$ with $A^+ A$ the identity. For $k = n$,
$A^+$ is the inverse matrix $A^{-1}$. If $A$ has rank less than $k$,
`Inverse(A)` is undefined. Construction does not check the rank of $A$.

A pass chooses how to evaluate `Inverse`, in closed form for a small matrix or
by a factorisation. The definition fixes the value and not the algorithm: the
normal equations square the condition number of $A$, and a factorisation such
as QR does not. The core writes the determinant as the free function
[`determinant`](#determinant), so the core has no node class for the
determinant.

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

### Determinant

Requires: [`Indexed`](#indexed), [`Sum`](#sum), [`Product`](#product),
[Negation and subtraction](#negation-and-subtraction), [Extents](#extents).

See also: [`Inverse`](#inverse), [Jacobian and density](#jacobian-and-density).

The free function `determinant(A)` is the Leibniz formula for the determinant
of a square matrix $A$.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `determinant(A)` | | `A: AbstractExpression` | The Leibniz formula for $\det A$ |

`determinant(A)` requires $A$ with value shape `(n, n)` and an `int` extent
$n$ (see [Extents](#extents)), and reads $n$ to choose the composition. A free
function does not validate, so construction does not check the requirement,
and `determinant(A)` of another $A$ is undefined, as `Inverse` of a matrix of
lower rank is. The Leibniz formula is

$$
\det A = \sum_\sigma \operatorname{sgn}(\sigma) \prod_{i=0}^{n-1} A_{i, \sigma(i)},
$$

with $\sigma$ over the permutations of $(0, \dots, n - 1)$. The product for
$\sigma$ is `A[0, σ(0)] * ... * A[n - 1, σ(n - 1)]`, nested to the left. The
products are joined in the lexicographic order of $\sigma$, nested to the
left, with `+` for an even $\sigma$ and `-` for an odd $\sigma$. For example,
for $n = 2$, `determinant(A)` is `A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]`. For
$n = 1$, `determinant(A)` is `A[0, 0]`, and for $n = 0$, `determinant(A)` is
`Integer(1)`.

The composition is fixed, so two calls of `determinant` on equal matrices
build equal expressions, and a pass that computes $\det A$ once for each
$A$ finds every use of $\det A$.

### Operator overloads

Requires: [`AbstractExpression`](#abstractexpression), [`Sum`](#sum),
[`Product`](#product), [`Division`](#division), [`Power`](#power), [Negation and
subtraction](#negation-and-subtraction), [`Absolute`](#absolute),
[`Indexed`](#indexed), [Tensor arithmetic](#tensor-arithmetic).

See also: [Numbers in expressions](#numbers-in-expressions), [Appendix: Free
functions](#appendix-free-functions).

Every overload in the [table of free functions](#appendix-free-functions) is
a method of `AbstractExpression` (see
[`AbstractExpression`](#abstractexpression)). No other class defines
arithmetic or indexing overloads. On
tensors, `+`, `-`, `*`, `/` and `@` are the free functions of [Tensor
arithmetic](#tensor-arithmetic).

An operator overload calls the free function with the operands unchanged, and
the constructor validates (see [Numbers in expressions](#numbers-in-expressions)).
`AbstractExpression` defines no reflected overloads, such as `__radd__`.

An overload applies one free function to the operands, so the graph follows
the grouping and the precedence of Python, and construction never regroups
(see [No simplification during
construction](#no-simplification-during-construction)). For example,
`a + b + c` is `Sum(Sum(a, b), c)` and `a + (b + c)` is `Sum(a, Sum(b, c))`,
and the two are not equal. `a ** b ** c` is `power(a, power(b, c))`, and
`-a ** b` is `negate(power(a, b))`.

`__getitem__` is the one exception. Python passes `f[i]` to `__getitem__` as the
single entry `i`, and passes `f[i, j]` as the `tuple` `(i, j)`. `__getitem__`
converts a single entry into a `tuple` of one entry and calls `component` (see
[Partial indexing](#partial-indexing)). The
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

Worked examples for this chapter in EXAMPLES.md: [transpose](EXAMPLES.md#transpose), [tensor algebra with `EllipsisIndex`](EXAMPLES.md#tensor-algebra-with-ellipsisindex).

## Calculus on a domain

### `Coordinate`

Requires: [`AbstractDomain`](#abstractdomain),
[`AbstractExpression`](#abstractexpression), [Support](#support).

See also: [`Gradient`](#gradient), [`DomainMap`](#domainmap).

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
$X : P \to R^n$, $X(p) = p$, and the gradient of $X$ is the identity matrix
(see [`Gradient`](#gradient)).

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

### `Gradient`

Requires: [Valid pair](#valid-pair), [`Coordinate`](#coordinate),
[`Indexed`](#indexed), [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum).

See also: [`PullBack`](#pullback), [`DomainMap`](#domainmap), [New
operations](#new-operations).

A concrete subclass of `AbstractExpression` for the gradient of an expression
on a domain $P$. The gradient adds a last axis in the coordinate frame of the
ambient coordinate domain of $P$, of dimension $n$.

| `Gradient(f, P)` |                                                                                                   |
| ---------------- | ------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `gradient(f, P)`                                                                                  |
| Takes            | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties       | `expression` $f$ and `domain` $P$, with $n$ the dimension of the ambient coordinate domain of $P$ |
| Requires         | A [valid pair](#valid-pair) $(f, P)$                                                              |
| Static type      | Value shape `f.value_shape + (n,)`                                                                |
| Raises           | `TypeError` for $P$ and `ValueError` for $f$ if a requirement fails                               |

The gradient acts on each component of an array-valued $f$. Entry `j` of the
last axis of `Gradient(f, P)` is the partial derivative $\partial_j f$ in the
coordinate $X_j$ of the frame.

#### The frame

For a domain $P$ of dimension $k$, the frame is the coordinate frame of the
ambient coordinate domain of $P$, of dimension $n$. A coordinate domain has
$k = n$.

- For $k = n$, $\partial_j$ is the partial derivative in the coordinate
  $X_j$ of $P$.
- For $k < n$, such as a decomposed domain that is a surface in $R^3$, the
  coordinate directions of the ambient coordinate domain are not tangent to
  $P$. `Gradient(f, P)` is then the tangential gradient, the gradient of an
  extension of $f$ off $P$, projected orthogonally onto the tangent space of
  $P$ with the Euclidean metric of the ambient coordinate domain. The
  tangential gradient does not depend on the extension.

The divergence, the curl and other operators of vector calculus are fixed
compositions of `Gradient`, which a library provides as free functions (see
[New operations](#new-operations)).

#### Pull-backs of `Gradient`

The gradient depends on the coordinates of $P$, so a pull-back does not
commute with `Gradient`. `PullBack(Gradient(f, Q), D)` is the gradient of $f$
on $Q$, evaluated at $\psi(p)$.

#### Examples

The examples use `P = RD(3)`, `X = coordinate(P)`, `x, y, z = X[0], X[1], X[2]`,
`f = x * y * z`, a function `u` on `P` of value shape `(3,)`, and
`i = index("i")`:

| Expression                              | Value shape | Meaning                                       |
| --------------------------------------- | ----------- | --------------------------------------------- |
| `gradient(f, P)`                        | `(3,)`      | The gradient $(yz, xz, xy)$, after evaluation |
| `gradient(gradient(f, P), P)`           | `(3, 3)`    | The Hessian of $f$                            |
| `gradient(X, P)`                        | `(3, 3)`    | The identity matrix, after evaluation         |
| `index_sum(gradient(u, P)[i, i], (i,))` | `()`        | The divergence of `u`                         |

- For a domain map $D$ with the parameter domain $P$,
  `gradient(D.parameterisation, P)` is the Jacobian matrix of value shape
  `(n, m)`, with $n$ the dimension of the ambient coordinate domain of
  `D.ambient_domain` and $m$ the dimension of the ambient coordinate domain of
  $P$ (see [Jacobian and density](#jacobian-and-density)).
- Let `M` be a decomposed domain of dimension 2 in `RD(3)`, and `w` a function
  on `M`. `gradient(w, M)` is the surface gradient of `w`, with value shape
  `(3,)`. `gradient(coordinate(M), M)` is the orthogonal projection onto the
  tangent plane of `M`, with value shape `(3, 3)`.
- `gradient(X, RD(2))` raises `ValueError`, because `X` has the support
  `RD(3)`.

### `DirectionalDerivative`

Requires: [`Gradient`](#gradient), [Free indices](#free-indices),
[Support](#support), [Passes](#passes).

See also: [Gateaux
derivative](#gateaux-derivative), [Argument spaces](#argument-spaces),
[Hyperelasticity and the stress](EXAMPLES.md#hyperelasticity-and-the-stress),
[Symmetric tensors](EXAMPLES.md#symmetric-tensors).

A concrete subclass of `AbstractExpression` for the derivative of an
expression $f$ with respect to a variable $x$ at a point $a$ in a direction
$h$,

$$
\mathrm{D}(f, x, a, h) = \frac{d}{d\varepsilon} f[x := a + \varepsilon h] \Big|_{\varepsilon = 0},
$$

with $\varepsilon$ real and $f[x := s]$ the substitution of $s$ for $x$ in $f$
that does not enter domains, domain maps or the spaces of functions. The
variable $x$ is a function or a `tuple` of functions. For a `tuple`, $a$ and
$h$ are `tuple`s with one member for each member of $x$, and the substitution
replaces every member $x_k$ by $a_k + \varepsilon h_k$ at once. The functions
of $x$ are bound in $f$ (see [Passes](#passes)).

| `DirectionalDerivative(f, x, a, h)` |                                                                                                       |
| -------------------------------- | ----------------------------------------------------------------------------------------------------- |
| Extends                          | `AbstractExpression` |
| Takes                            | `f: AbstractExpression`, `x: Function \| tuple[Function, ...]`, `a, h: AbstractExpression \| tuple[AbstractExpression, ...]` |
| Properties                       | `expression` $f$, `variable` $x$, `point` $a$ and `direction` $h$ |
| Requires                         | $x$ a `Function`, or a `tuple` of at least two distinct functions, each with the argument spaces `()`, and member by member: $a$ and $h$ with the value shape of the member of $x$ and the argument spaces `()`, $h$ with the support of the member of $x$ or no support, and $a$ and $h$ real for a real member of $x$ |
| Static type                      | The value shape of $f$, the free indices of $f$ followed by those of the members of $a$ and $h$, each `Index` taken once, and the argument spaces and the scalar type of $f$ |
| Raises                           | `TypeError` for an $x$ that is not a `Function` or a `tuple` of functions, and `ValueError` if another requirement fails |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `directional_derivative(f, x, a, h)` | | `f: AbstractExpression`, `x: Function \| tuple`, `a, h: AbstractExpression \| tuple` | `DirectionalDerivative(f, x, a, h)` |

The variable is a function, and the point is any expression of the value shape
of the function, so the derivative of an expression of a symbol $x$ at a point
$a$ is the derivative of a function of $x$ at $a$. The substitution replaces only the occurrences of $x$, so
`DirectionalDerivative` is a partial derivative. A `tuple` variable gives the
derivative in all of its members at once, at the point $a$ of all the
members, such as for the unknowns of a mixed problem.

The substitution does not enter domains, so the derivative holds the geometry
fixed. A shape derivative, with respect to the coordinates of a mesh, is a
language extension (see [Language extensions](#language-extensions)).

The derivative is real-linear in $h$. The gradient of a direction with no
support is zero, so the derivative with respect to a function `u` in a constant
direction holds `gradient(u, P)` fixed. A direction that is a function in the
space of `u` gives the Gateaux derivative (see [Gateaux
derivative](#gateaux-derivative)).

The derivative with respect to a tensor is a composition over unit
directions. The direction states the symmetry, and the core has none. With a
Python function `strain_energy` that builds the strain energy from a tensor,
the right Cauchy–Green tensor `C`, a function `C_` in a real space of
constants of value shape `(3, 3)`, the identity matrix `I` of value shape
`(3, 3)` and indices `i`, `j`, `k` and `l`:

```python
E = component_tensor((I[i, k] * I[j, l] + I[i, l] * I[j, k]) / Integer(2), (i, j))
S = Integer(2) * component_tensor(directional_derivative(strain_energy(C_), C_, C, E), (k, l))
```

`E` has the free indices `k` and `l`, and is symmetric for each value of `k`
and `l`, so `S` is $2 \operatorname{sym}(\partial \psi / \partial C)$ at `C`,
the second Piola–Kirchhoff stress. For a strain energy `psi` already built
from `C`, `replace(psi, {C: C_})` gives the strain energy of the symbol (see
[Passes](#passes)). The direction `component_tensor(I[i, k] * I[j, l], (i, j))`
gives the full derivative, without `sym`.

### `Integral`

Requires: [Valid pair](#valid-pair), [Support](#support).

See also: [Concept: integral over a domain](#concept-integral-over-a-domain),
[An affinely mapped triangle](EXAMPLES.md#an-affinely-mapped-triangle).

A concrete subclass of `AbstractExpression` for the integral of an expression
$f$ over a domain $P$. A basic integral needs no topological notion of the
domain.

| `Integral(f, P)` |                                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `integral(f, P)`                                                                            |
| Takes            | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties       | `expression` $f$ and `domain` $P$                                                           |
| Requires         | A [valid pair](#valid-pair) $(f, P)$, or $f$ with the support `P.ambient_domain`, and an `Index` free in both $f$ and a piece of $P$ with the same extent in both |
| Static type      | The value shape of $f$, and the free indices of $f$ other than the free indices of the pieces of $P$, with no support |
| Raises           | `TypeError` for $P$ and `ValueError` for $f$ if a requirement fails                         |

For a coordinate domain $P$ of dimension $n$, `Integral(f, P)` is the
Lebesgue integral in the coordinates $X$ of $P$,

$$
\int_P f \,\mathrm{d}X.
$$

For $n = 0$, the integral is the value of $f$ at the single point of $R^0$.

For a decomposed domain $P$, `Integral(f, P)` is the integral of $f$ against
the density that the pieces induce, as the sum of the integrals over the
pieces (see [`Integral` over a decomposed
domain](#integral-over-a-decomposed-domain)).

An array-valued $f$ integrates component by component. A later pass might
evaluate an `Integral` over a coordinate domain, for example by quadrature.

#### Examples

With `T = UnitSimplex(2)` and `X = coordinate(T)`:

- `integral(X[0], T)` is $\int_T X_0 \,\mathrm{d}X$. A later pass evaluates
  the integral to $1/6$.
- `integral(Integer(1), T)` is the area $1/2$ of `T`.
- `integral(X[0], T) * integral(X[1], T)` is a product of two integrals, with
  value shape `()` and no support.

## Maps between domains

### `DomainMap`

Requires: [`AbstractDomain`](#abstractdomain),
[`Coordinate`](#coordinate), [Support](#support), [Valid pair](#valid-pair), [Free indices](#free-indices),
[Scalar type](#scalar-type).

See also: [`PullBack`](#pullback), [`DecomposedDomain`](#decomposeddomain),
[Jacobian and density](#jacobian-and-density), [An affinely mapped
triangle](EXAMPLES.md#an-affinely-mapped-triangle).

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
| Requires                         | A [valid pair](#valid-pair) $(\psi, P)$, $\psi$ real with value shape `(n,)` for the dimension $n$ of the ambient coordinate domain of $Q$ and the argument spaces `()`, and $k \le n$ |
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
[`DecomposedDomain`](#decomposeddomain)). An integral over such a cell expands
into the sub-cells, which is where the composite quadrature rule of a macro
element comes from.

A $\psi$ with free indices makes $D$ a family of domains, with one domain for
each value of the free indices. For example, with an array `x` of value shape
`(m, n)` and an index `i`, `domain_map(RD(0), x[i], RD(n))` is a family of $m$
points in `RD(n)`. A family of domains writes a set of points
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
[Free indices](#free-indices), [An
affinely mapped triangle](EXAMPLES.md#an-affinely-mapped-triangle).

See also: [`Gradient`](#gradient).

A concrete subclass of `AbstractExpression` for the pull-back $\psi^* f$ of an
expression $f$ on $Q$ by the parameterisation $\psi$ of a domain $D$, with $P$
the parameter domain and $Q$ the ambient domain of $D$.

| `PullBack(f, D)` |                                                                                                            |
| ---------------- | ---------------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `pull_back(f, D)`                                                                                          |
| Takes            | `f: AbstractExpression`, `D: DomainMap \| SideMap` |
| Properties       | `expression` $f$, and `domain` $D$, a `DomainMap` or a [`SideMap`](#sidemap)                 |
| Requires         | A [valid pair](#valid-pair) $(f, Q)$, and an equal extent for an `Index` free in both $f$ and $D$          |
| Static type      | The value shape of $f$, and the free indices of $f$ followed by those of $D$, with each `Index` taken once |
| Raises           | `TypeError` for $D$, and `ValueError` if another requirement fails                                         |

The second operand is a domain, not an expression, because a domain states
$Q$.

The pull-back is the composition $f \circ \psi$, the expression on $P$ with
the value $f(\psi(p))$ at each point $p$ of $P$. The axes of the value shape
are not transformed. An axis of an expression on $Q$ is an array in the
coordinate frame of $Q$, so the pull-back keeps the frame of $Q$.

A `PullBack` of a downstream terminal or of a
[`Function`](#function) remains, and the library of the terminal evaluates the
terminal at $\psi(p)$. For a function, the library
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

# PullBack(f, e), with the support UnitSimplex(1) and value shape ().
f = X[0] * X[1]
fe = pull_back(f, e)

# The edge of K from (3, 1) to (2, 3), as a domain in RD(2).
E = domain_map(S, pull_back(psi, e), RD(2))
```

- A later pass evaluates `fe` to $(1 - t) t$, and the parameterisation of `E`
  to $(3 - t, 1 + 2 t)$.
- `pull_back(f, domain_map(S, tensor((one - t, t)), RD(2)))`
  raises `ValueError`, because that domain lies in `RD(2)`, not in `T`.

### Jacobian and density

Requires: [`DomainMap`](#domainmap), [`Gradient`](#gradient),
[`Inverse`](#inverse), [Determinant](#determinant), [`Absolute`](#absolute),
[Elementary functions](#elementary-functions-very-incomplete).

See also: [Concept: integral over a domain](#concept-integral-over-a-domain),
[`Integral` over a decomposed domain](#integral-over-a-decomposed-domain).

The free functions `jacobian(D)`, `inverse_jacobian(D)` and `density(D)` are
the Jacobian matrix, the left inverse of the Jacobian matrix on the tangent
space, and the density of a domain $D$ = `DomainMap(P, psi, Q)` of dimension
$k$. $n$ is the dimension of the ambient coordinate domain of $Q$, and $m$ the
dimension of the ambient coordinate domain of $P$, with $m = k$ for a
coordinate domain $P$.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `jacobian(D)` | | `D: DomainMap` | `gradient(D.parameterisation, D.parameter_domain)` |
| `inverse_jacobian(D)` | | `D: DomainMap` | The matrix $K$ below |
| `density(D)` | | `D: DomainMap` | The density $\rho$ below |

The Jacobian matrix $J$ = `jacobian(D)` has value shape `(n, m)`. With
$\Pi$ = `gradient(coordinate(P), P)`, the orthogonal projection onto the
tangent space of $P$ (see [`Gradient`](#gradient)), and $I$ the identity
matrix of value shape `(m, m)`, the metric of $D$ is

$$
G = J^T J + I - \Pi.
$$

$J^T J$ is zero on the normal space of $P$ and $I - \Pi$ is the identity
there, so $G$ is invertible, and on the tangent space $G$ is the metric
$J^T J$ that $\psi$ induces. For $k = m$, $\Pi$ is the identity and
$G = J^T J$. The density $\rho$ is the factor by which $\psi$ scales the
length, the area or the volume of a small region of $P$, and $K$, of value
shape `(m, n)`, is the left inverse of $J$ on the tangent space of $P$:

$$
\rho = \sqrt{\det G}, \qquad K = G^{-1} J^T.
$$

The definitions fix the values. `density(D)` and `inverse_jacobian(D)` read
$k$, $m$ and $n$ to choose the composition, so that a cell gets the closed
form and not the normal equations:

| Case | `density(D)` | `inverse_jacobian(D)` |
| ---- | ------------ | --------------------- |
| $m = 0$ | `Integer(1)` | `inverse(J)`, of value shape `(0, n)` |
| $0 < k = m = n$ | `absolute(determinant(J))` | `inverse(J)` |
| $0 < k = m < n$ | `square_root(determinant(G))` | `inverse(J)`, the pseudo-inverse |
| $k < m$, a decomposed $P$ such as a surface mesh | `square_root(determinant(G))` | `component_tensor(index_sum(inverse(G)[a, c] * J[l, c], (c,)), (a, l))` |

$J^T J$ is `component_tensor(index_sum(J[l, a] * J[l, b], (l,)), (a, b))`, and
the indices in the compositions are derived indices. For $k < m$, $J$ is the
tangential gradient of $\psi$ and has less than full rank, so `inverse(J)` is
undefined, and $G$ is the device that gives the Gram determinant and the
pseudo-inverse of $J$ on the tangent space.

The density for $m > n$ is undefined, and construction does not check. The
core writes $\det J$ only as `determinant(jacobian(D))` and $K$ only as
`inverse_jacobian(D)`, so every use of $\det J$ or $K$ for one $D$ is one
expression, and a pass that computes the geometry of a cell once computes
each once. A library that needs $J$, $\det J$ or $K$, such as for a Piola
map, uses `jacobian(D)`, `determinant(jacobian(D))` and
`inverse_jacobian(D)`.

### Concept: integral over a domain

Requires: [`PullBack`](#pullback), [`Integral`](#integral), [Jacobian and
density](#jacobian-and-density).

The integral of an expression $f$ on $Q$ over a domain $D$ =
`DomainMap(P, psi, Q)` is

$$
\int_D f = \int_P (f \circ \psi) \, \rho \, \mathrm{d}X,
$$

which is `integral(pull_back(f, D) * density(D), P)` for a scalar $f$. The
integral does not depend on the direction in which $\psi$ traverses $D$.

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
cells of a decomposed domain must be oriented consistently, such as for the
flux through a surface in $R^3$, the mesh library provides a sign $+1$ or $-1$
for each cell as an expression (see [The flux through a
surface](EXAMPLES.md#the-flux-through-a-surface)).

### Points and point evaluation

Requires: [`RD`](#rd), [`DomainMap`](#domainmap),
[`PullBack`](#pullback), [`Integral`](#integral), [An affinely mapped
triangle](EXAMPLES.md#an-affinely-mapped-triangle).

See also: [Restriction](#restriction), [Side domains](#side-domains).

A point is the image of a domain map. A point in a domain $Q$ is
`domain_map(RD(0), x, Q)`, with $x$ a constant expression of value shape
`(n,)` for the dimension $n$ of the ambient coordinate domain of $Q$. The core
needs no node class and no free function for a point.

#### Point evaluation

The value of an expression $f$ of any value shape at a point $p$ is
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
point in that cell. These are the side tables of the topology of the points (see [Sides of a topology](#sides-of-a-topology)). The
topology of the points has one cell of a reference cell type of dimension 0
for each point, and the side `("inside", T, j)`. A family of points in a mesh is
therefore a topology with a side.

The mesh library locates the cells and inverts the map of each cell outside
the language, for example by Newton's method, when the mesh library builds the
topology $S$ of the points, as the mesh library computes every side table.
`cells` has the shape `(m,)` and the extent `T.num_cells[j]`. The inclusion
of each point is the constant map onto the position of the point in $K_j$,
and the points have no reference normal.

An expression $f$ on $D$ is evaluated at the points through the
[restriction](#restriction) to the points, with $P$ = `side_domain(D, S)` the
decomposed domain of the points:

```python
inside = side_map(D, P, "inside")
values = pull_back(pull_back(f, inside), cell_domain(P, 0))
total = integral(pull_back(f, inside), P)
```

`values` is the family of the values at the points, with the cell index of
$S$, and `total` is their sum. [`expand_pull_backs`](#expanding-pull-backs)
lowers the restriction as on any other topology with a side, to the value on
the cell on the side at the position of the point, so point evaluation needs
no pass outside the core.

The core does not rewrite `pull_back(u, M)` for a function `u` on $D$,
because the library of the space of `u` defines `u` on the cells only (see
[`AbstractFunctionSpace`](#abstractfunctionspace)), and such a pull-back
remains in the core normal form. A family of points known in advance is
therefore written as the topology of the points, not as $M$. Points known only
after a pass of the compiler, such as points that depend on the quadrature
points, are located by a pass after the core, of the compiler with the mesh
library, which rewrites the pull-back as the restriction to the topology of
the points, with the index $i$ of the family in place of the cell index:

```python
g = gathered_index(cells, (i,))
R = domain_map(RD(0), positions[i], K_j)
pull_back(pull_back(u, gathered_cell_domain(D, j, g)), R)
```

with `positions[i]` the position in $K_j$ of the point $i$, which the mesh
library computes.

The pull-back of `u` by a gathered cell domain is in the contract of the
space. The maps below `u` may be a chain, such as a family of maps from a
decomposed domain, pulled back by a cell of that domain and then by the
quadrature points of the reference cell. The pass reads the chain as one
family with all of its free indices (see [A 1D network in a 3D
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

Worked examples for this chapter in EXAMPLES.md: [an affinely mapped triangle](EXAMPLES.md#an-affinely-mapped-triangle), [circulation along an edge](EXAMPLES.md#circulation-along-an-edge), [integrating over an edge](EXAMPLES.md#integrating-over-an-edge), [integrating over the boundary of a triangle](EXAMPLES.md#integrating-over-the-boundary-of-a-triangle), [quadrature](EXAMPLES.md#quadrature), [integrating over a cut cell](EXAMPLES.md#integrating-over-a-cut-cell), [the centroid of a triangle](EXAMPLES.md#the-centroid-of-a-triangle).

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

```{dropdown} Legacy UFL
In legacy UFL, a `Mesh` holds both the topology and the coordinate element, and the measure `dx` hides the sum over the cells. In UFLx, the topology, the parameterisations of the pieces of a decomposed domain and the `IndexSum` over `cell_index(T, j)` are separate.
```

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

On each side, a cell of type $i$ of $S$ is the image of the cell type $F_i$
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
topology](EXAMPLES.md#a-static-topology)):

- For a matching facet, the position is that of the same vertex.
- For a periodic boundary, the position is that of the vertex one period away.
- For a piece, the position is that of a point of the cell.

#### Side tables

A concrete topology with sides also provides `side_tables(i, side)`, the
triple `(cells, inclusion, reference_normal)` for the side named `side` of
group $i$, with the side `(side, T, j)`, the cell type $F_i$ of group $i$ and
`e = cell_index(S, i)`:

- `cells` is an [`AbstractIndexValue`](#abstractindexvalue) of the shape
  `(S.num_cells[i],)` and the extent `T.num_cells[j]`. Entry $e$ is the cell of
  type $j$ of $T$ on the side of the cell $e$ of $S$.
- `inclusion` is a [`DomainMap`](#domainmap) from $F_i$ into $K_j$ with the
  free index `e`. Member $e$ is the inclusion of the cell $e$ into the cell on
  the side.
- `reference_normal` is an expression with the support $F_i$ or no support,
  the value shape `(k_j,)` and the free index `e`, for a group of cells of
  dimension $k_j - 1$, and `None` for any other group. A normal that varies
  along the cell, such as on a curved piece, has the support $F_i$. Member $e$ is a normal in $K_j$ of
  the image of the cell $e$. The normal points out of $K_j$ for an image on the
  boundary of $K_j$, and away from the centroid of $K_j$ for an image inside
  $K_j$.

The inclusions of the sides of one group parameterise each cell $e$ in the same
way: a point of $F_i$ is the same point of the cell $e$ from every side, so the
difference of two restrictions is the jump (see [Restriction](#restriction)).
[Restriction](#restriction) and [`Normal`](#normal) read the tables. The core
does not check the sides or the tables.

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
[Points and point evaluation](#points-and-point-evaluation)). The inclusion of
a piece of a curved cell may be an approximation, as a quadrature rule is. The
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

```{dropdown} Legacy UFL
The legacy UFL measures `ds` and `dS` integrate over a topology of facets with one and two sides, and legacy UFL names the two sides `'+'` and `'-'`. The legacy UFL measure `dx(1)` integrates over a topology of the cells of a subdomain, and the legacy UFL measure `dP` integrates over a topology of vertices.
```

Worked examples for this chapter in EXAMPLES.md: [the simplex as a cell type](EXAMPLES.md#the-simplex-as-a-cell-type), [a static topology](EXAMPLES.md#a-static-topology).

## Reference cells

### Reference cell

Requires: [`AbstractCoordinateDomain`](#abstractcoordinatedomain),
[`AbstractCellType`](#abstractcelltype).

See also: [`AbstractTopology`](#abstracttopology), [Integrating over
the boundary of a triangle](EXAMPLES.md#integrating-over-the-boundary-of-a-triangle).

An abstract subclass of `AbstractCoordinateDomain` and of `AbstractCellType`
for a reference cell `AbstractReferenceCell`, a cell type that is also a
coordinate domain. The core reads no vertices of a reference cell. A library
that builds inclusions, such as from vertex positions, keeps the vertices of
its reference cells (see [The unit simplex as a reference
cell](EXAMPLES.md#the-unit-simplex-as-a-reference-cell)).

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `AbstractCoordinateDomain`, `AbstractCellType` |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](EXAMPLES.md#the-unit-simplex-as-a-reference-cell)                                                   |

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
chi = component_tensor(psi_j[k], (c, k))[g]
```

with the partial indexing `[g]` (see [Partial indexing](#partial-indexing)),
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

```{dropdown} Legacy UFL
In legacy UFL, a `Mesh` is a decomposed domain, and `SpatialCoordinate(mesh)` is `coordinate(D)`.
```

### Restriction

Requires: [`DecomposedDomain`](#decomposeddomain),
[Sides of a topology](#sides-of-a-topology), [`GatheredIndex`](#gatheredindex),
[`PullBack`](#pullback), [Reference cell](#reference-cell).

See also: [`Integral` over a decomposed domain](#integral-over-a-decomposed-domain).

The restriction of an expression $f$ on a decomposed domain $D$ to the cells
of a topology $S$ with sides, from the side $s$, is the pull-back of $f$ along
the inclusion of each cell of $S$ into the cell on the side $s$. The side
tables of $S$ give the cell on the side, through a
[`GatheredIndex`](#gatheredindex) entry, and the inclusion of each cell of $S$
into that cell. The restriction of $f$ to a decomposed domain $E$ of the cells
of $S$ from the side $s$ is `pull_back(f, side_map(D, E, s))` (see
[`SideMap`](#sidemap)).

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_domain(D, S)` | | `D: DecomposedDomain`, `S: AbstractTopology` | The decomposed domain of the cells of $S$, with the geometry of the first side of each group on the topology of $D$ |

Notation, for the group $i$ of $S$ with the side `(s, T, j)`:

- $T$ is the topology of $D$, and $K_j$ is the reference cell of cell type $j$
  of $T$.
- `(cells, inclusion, reference_normal)` are the side tables of the side $s$
  of group $i$.
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
- $A_i$, the inclusion into the cell on the side, is `inclusion` of the side
  tables, a domain map from $F_i$ into $K_j$. $F_i$ must be a reference cell.

The restriction of $f$ to group $i$ from the side $s$ is
`pull_back(pull_back(f, B_i), A_i)`, an expression on $F_i$ with the free index
`e`.

The inclusions of all sides parameterise the cell $e$ in the same way (see
[Side tables](#side-tables)). The restrictions from two sides are therefore
expressions on the same $F_i$, with the same coordinates, and the difference
of the two restrictions is the jump. A cell of $S$ with more
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
| Requires           | Conditions 1 to 4 below                                                                                 |
| Raises             | `ValueError` if a requirement fails                                                                     |

1. $D$ has a topology.
2. Every cell type of $S$ is an `AbstractReferenceCell`.
3. Every group of $S$ has a side named $s$ with the topology of $D$.
4. The parameter domain of piece $j$ of $D$ is the reference cell $K_j$ of
   the topology of $D$, because the inclusions of the side tables map into
   $K_j$.

`pull_back(f, SideMap(D, E, s))` is an expression with the support $E$, the
restriction of $f$ to all the groups at once. The pull-back of that expression
by `cell_domain(E, i)` is `pull_back(pull_back(f, B_i), A_i)`.

#### Side domains

Piece $i$ of `side_domain(D, S)` is
`domain_map(F_i, pull_back(pull_back(coordinate(D), B_i), A_i), Q)`, with
$B_i$ and $A_i$ of the first side of group $i$ whose topology is the
topology of $D$. Piece $i$ has the one free index `cell_index(S, i)`, the free
index `e` of $B_i$ and $A_i$. $D$ must have a topology, every group of $S$ must have a side
on that topology, and $F_i$ must be an `AbstractReferenceCell`. The topology
of the result is $S$. The core does not
check that the other sides have the same geometry.

#### Traces and normals

The restriction of an expression is the value of the expression on the cell
of $S$, from the side:

- The restriction of `gradient(u, D)` is the whole gradient of `u` on the cell
  on the side, as for a normal derivative. The last axis is in the frame of
  the ambient coordinate domain from every side.
- `gradient(pull_back(u, side_map(D, E, s)), E)` is the tangential gradient of the
  restriction along the cells of $E$.

The outward normal of the cell on a side is [`normal(side_map(D, E, s))`](#normal).

#### Example

The example uses `D` and `square` of the
[`DecomposedDomain`](#decomposeddomain) example, `I = ReferenceSimplex(1)`, a
function `u` on `D`, and the diagonal of `square`:

```python
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
E = side_domain(D, diagonal)
jump = pull_back(u, side_map(D, E, "+")) - pull_back(u, side_map(D, E, "-"))
```

- From `"+"`, the side tables are the cells `(0,)` and the inclusion
  $A_+(t) = (1 - t, t)$, through the local vertices 1 and 2 of `K`. From
  `"-"`, the side tables are the cells `(1,)` and the inclusion
  $A_-(t) = (0, t)$, through the local vertices 0 and 2. The coordinates of
  `E` are $(1 - t, t)$.
- Let `pull_back(u, cell_domain(D, 0))` be
  `a[c, 0] + a[c, 1] * X[0] + a[c, 2] * X[1]`. The two restrictions pulled
  back to the diagonal are $a_{0,0} + a_{0,1}(1 - t) + a_{0,2} t$ and
  $a_{1,0} + a_{1,2} t$. The two are equal for all $t$ if `u` is continuous at
  the vertices 1 and 2.
- With the diagonal listed as `(2, 1)`, the inclusion from `"+"` is
  $A_+(t) = (t, 1 - t)$. Both sides use the vertex numbering of the diagonal, so the restrictions from both sides run the
  other way along the diagonal, and the jump at each point of the diagonal is
  unchanged.
- For the facets `S` of a triangle and quadrilateral mesh, with the three
  groups of [Sides of a topology](#sides-of-a-topology), the same two lines
  with `E = side_domain(D, S)` give one jump for all three groups.

```{dropdown} Legacy UFL
The legacy UFL restriction `f('+')` is `pull_back(f, side_map(D, E, "+"))` on the decomposed domain $E$ of the interior facets, whose sides are named `'+'` and `'-'`.
```

### `Integral` over a decomposed domain

Requires: [`Integral`](#integral),
[Concept: integral over a domain](#concept-integral-over-a-domain),
[`DecomposedDomain`](#decomposeddomain), [Jacobian and
density](#jacobian-and-density).

`Integral(f, D)` for a decomposed domain $D$ integrates $f$ against the
density that the pieces induce. $f$ has the support $D$, the support
`D.ambient_domain` or no support.

Notation, for the piece $j$ of $D$:

- $P_j$ is the parameter domain and $\psi_j$ is the parameterisation.
- `indices_j` are the free indices of piece $j$.

#### Expansion into pieces

The integral is the sum over the pieces $j$ of the integrals over the domains
of piece $j$. For a scalar $f$, `Integral(f, D)` has the value of that sum:

```python
D_j = cell_domain(D, j)
index_sum(integral(pull_back(f, D_j) * density(D_j), P_j), indices_j)
```

- A piece with no free indices has no `IndexSum`.
- An array-valued $f$ integrates component by component, with the density
  written into each component in index notation.
- For a decomposed domain $P_j$, the inner integral expands in turn.
- For a parameter domain $P_j$ of dimension $m = 0$, the density is
  `Integer(1)`, so the integral over a set of points is the sum of the values.

The density needs no orientation, so an integral of a positive $f$ is
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
- `integral(pull_back(u, side_map(D, B, "inside")), B)` integrates `u` over the boundary.

```{dropdown} Legacy UFL
`f * dx`, `f * ds` and `f * dS` are `integral(f, D)` over the decomposed domain of the cells, of the boundary facets and of the interior facets, with the fields restricted to the facets in `f`. The three are one integral over three domains.
```

### `Normal`

Requires: [Restriction](#restriction), [Reference cell](#reference-cell),
[Jacobian and density](#jacobian-and-density).

A concrete subclass of `AbstractExpression` for the unit outward normal of the
cell on a side, on the cells of a topology with sides.

| `Normal(sm)` |                                                                                                       |
| ------------ | ----------------------------------------------------------------------------------------------------- |
| Extends      | `AbstractExpression` |
| Free function | `normal(sm)`                                                                                         |
| Takes        | `sm: SideMap` |
| Properties   | `side_map`, a [`SideMap`](#sidemap) `SideMap(D, E, s)`                                               |
| Requires     | $E$ of the dimension $k - 1$ for the dimension $k$ of $D$                                             |
| Static type  | Value shape `(n,)` for the dimension $n$ of the ambient coordinate domain, the support $E$, real      |
| Raises       | `TypeError` for an `sm` that is not a `SideMap`, and `ValueError` if the requirement fails           |

The normal of `sm` belongs to the same side as the restriction
`pull_back(f, sm)`, so a form writes both from one side map.

Notation, for `sm` = `SideMap(D, E, s)` and the piece $i$ of $E$:

- $B_i$ and $A_i$ are the maps of [Restriction](#restriction) for the side $s$.
- $K_j$ is the reference cell of the cell on the side.
- $k$ is the dimension of $D$, and $n$ is the dimension of the ambient
  coordinate domain.

#### The pull-back of `Normal`

The pull-back of `Normal(SideMap(D, E, s))` by `cell_domain(E, i)` is

$$
\nu = \frac{K^T \hat m}{\lvert K^T \hat m \rvert},
$$

where:

- $K$ is `pull_back(inverse_jacobian(B_i), A_i)`, the inverse Jacobian of
  the cell on the side at the points of $A_i$ (see [Jacobian and
  density](#jacobian-and-density)). The parameter domain of $B_i$ is a
  reference cell (see [`SideMap`](#sidemap)), so $K$ is `inverse(J)`.
- $\hat m$ is `reference_normal` of the side tables of the side $s$ of group
  $i$ (see [Side tables](#side-tables)), a normal of the cell of $E$ in $K_j$.

The composition is

```python
w = transpose(K) @ m_hat
nu = w / square_root(w @ w)
```

The core needs no numbering of the facets and no vertices of $K_j$, because
the library of the topology gives $\hat m$.

#### The direction of the normal

- For a cell of $E$ on the boundary of $K_j$, such as a facet or a piece of a
  facet, $\hat m$ points out of $K_j$ (see [Side tables](#side-tables)).
- For a cell of $E$ inside $K_j$, such as a cut through the cell, $\hat m$
  points away from the centroid of $K_j$, by convention.
- For $k = n$, the normals from two sides of one cell of $E$ are opposite.
- For $k < n$, such as on the edges of a surface in $R^3$, $\nu$ is the unit
  conormal. The conormal is tangent to the cell on the side and orthogonal to
  the cell of $E$. The conormals from two sides are opposite only where the
  cells on the two sides have the same tangent space. On a surface, a jump
  therefore weights each side by the normal of that side, as in
  `u_p * n_p + u_m * n_m`.

`Normal` is the normal of a cell on a side, not the unit normal field of a
domain of dimension $n - 1$, such as of a surface in $R^3$. The
unit normal field needs a consistent orientation of the cells of the domain
(see [Concept: integral over a domain](#concept-integral-over-a-domain)) and
belongs to a downstream library.

#### Examples

The examples use the following objects:

- `B` of the [`Integral` over a decomposed
  domain](#integral-over-a-decomposed-domain) examples.
- The decomposed domain `E` of the interior facets of `D`, with the sides
  `"+"` and `"-"`.
- `grad_u = gradient(u, D)` and `i = index("i")`.

```python
inside, plus, minus = side_map(D, B, "inside"), side_map(D, E, "+"), side_map(D, E, "-")
n_B = normal(inside)

# The flux of grad u out of the square, the legacy UFL dot(grad(u), n) * ds.
flux = integral(index_sum(pull_back(grad_u, inside)[i] * n_B[i], (i,)), B)

# The jump of the normal derivative across the diagonal, from the side "+".
n_E = normal(plus)
jump = pull_back(grad_u, plus)[i] - pull_back(grad_u, minus)[i]
penalty = integral(index_sum(jump * n_E[i], (i,)) ** Integer(2), E)
```

Worked examples for this chapter in EXAMPLES.md: [an interior facet integral](EXAMPLES.md#an-interior-facet-integral), [a Bank-Weiser error estimator](EXAMPLES.md#a-bank-weiser-error-estimator), [a fluid and a solid](EXAMPLES.md#a-fluid-and-a-solid), [a 1D network in a 3D mesh](EXAMPLES.md#a-1d-network-in-a-3d-mesh), [a deformed domain](EXAMPLES.md#a-deformed-domain).

## Function spaces

A function is a member of a function space, such as a coefficient or an
argument of a variational form.

### `AbstractFunctionSpace`

Requires: [`AbstractDomain`](#abstractdomain),
[Scalar type](#scalar-type).

An abstract subclass of the base object `AbstractFunctionSpace` for a space of
functions on a domain.

| `AbstractFunctionSpace` |                                                                                                                                                       |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `Base` |
| Provides                | `domain`, an [`AbstractDomain`](#abstractdomain) $D$ or `None`, `value_shape`, `scalar_type` and `argument_spaces`, `()` for a space of fields |
| Defined by              | A subclass outside the core, such as a subclass for the finite element spaces of Basix, and [`DualSpace`](#dualspace) in the core                           |

A space with the domain `None` has functions that do not depend on a domain.

#### Functions in a space

A function in $V$ is an expression with the support $D$ and with the value
shape and the scalar type of $V$.

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
constants. A function in the space is an expression with no support, a
single array with the value shape and the scalar type of $V$.

The integrals of a form decide the cells on which a constant acts.

A space of constants is needed for a constant unknown, such as a Lagrange
multiplier, because an argument comes from a space. Constant data needs no
space. An expression with no support, such as a literal or an array that a
kernel receives, is already a constant.

#### Representation on the reference cell

The pull-back by a cell fixes what a function in $V$ is on the cell, and
keeps the frame of each axis of the value shape (see [`PullBack`](#pullback)).
The library of $V$ fixes how the function is represented on the reference
cell, as a map of a function $\hat u$ on the reference cell.

In the table, $k$ is the dimension of $D$, $n$ is the dimension of the ambient
coordinate domain, $J$ is the Jacobian matrix of the cell, of value shape
`(n, k)`, $K$ is `inverse_jacobian` of the cell, and $\rho$ is the density of the cell (see
[Jacobian and density](#jacobian-and-density)):

| Space                             | Value shape | Map of $\hat u$                                    | Legacy UFL map      |
| --------------------------------- | ----------- | -------------------------------------------------- | ------------------- |
| Lagrange                          | `()`        | $\hat u$                                           | identity            |
| Vector Lagrange                   | `(n,)`      | $\hat u$                                           | identity            |
| Nédélec, $H(\mathrm{curl})$       | `(n,)`      | $K^T \hat u$                                       | covariant Piola     |
| Raviart–Thomas, $H(\mathrm{div})$ | `(n,)`      | $J \hat u / \det J$ for $k = n$, and $J \hat u / \rho$ for $k < n$ | contravariant Piola |
| Discontinuous densities           | `()`        | $\hat u / \det J$ for $k = n$, and $\hat u / \rho$ for $k < n$       | $L^2$ Piola         |
| Regge                             | `(n, n)`    | $K^T \hat u K$                                     | double covariant    |

Each map is a fixed composition of $\hat u$ with `jacobian`, `inverse`,
`determinant` and `density`, so a library writes each map once, or takes the
maps from a language extension (see [A Piola
extension](#a-piola-extension)).

#### Mixed and Sobolev spaces

A mixed problem needs no mixed space, because an argument of a [form](#form)
may be a `tuple` of functions in several spaces.

The Sobolev space of $V$ is not yet covered. The Sobolev space is expected to
be stated through the [restrictions](#restriction) of the functions of $V$ to
the facets of the topology of $D$, from the two sides of each facet.

#### Examples

The examples use a decomposed domain `D` of `ReferenceSimplex(2)` cells in
`RD(2)`, so $k = n = 2$.

- The Lagrange space of polynomial degree 1 on `D` has the value shape `()`.
- The Nédélec and the Raviart–Thomas spaces on `D` both have the value shape
  `(2,)`, and differ in the map of $\hat u$.

```{dropdown} Legacy UFL
The legacy UFL `Real` space is a space of constants, and the legacy UFL `Constant` is constant data.

In legacy UFL, a `FunctionSpace` holds a mesh and an element, and the element states the value shape, the map to the reference cell and the Sobolev space. In UFLx, the pull-back states the meaning of a function on a cell, the library states the representation, and the space holds a decomposed domain, the legacy mesh.
```

### `Function`

Requires: [`AbstractFunctionSpace`](#abstractfunctionspace),
[`AbstractExpression`](#abstractexpression), [Support](#support).

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

A function is not an argument or a coefficient by itself. Whether a function is
an argument of a form is stated by the form (see [`Form`](#form)).

Examples, with a decomposed domain `D` and a space `V` on `D` of the value
shape `()`:

- `u = function(V, "u")` has the support `D` and the value shape `()`, and
  `function(V, "u") == u` is `True`.

```{dropdown} Legacy UFL
The legacy UFL `Coefficient`, `TestFunction` and `TrialFunction` are all functions. A test function and a trial function are the functions that a form states as the arguments at the positions 0 and 1.
```

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
value shape.

A form is linear in an argument with the space $V$, and antilinear in an
argument with the space `conjugate_space(V)` (see
[`ConjugateSpace`](#conjugatespace)). The form states the linearity when the
form binds the argument (see [`Form`](#form)).

The argument spaces are static, as the value shape is, so a form with a
wrong arity is rejected during construction, and a load that is antilinear in
`v` and a form that is linear in `v` do not add.

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
| Derived             | `domain`, `value_shape`, `scalar_type` and `argument_spaces`, those of $V$ |
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
| Free function  | `dual_space(S)`                                                                                                            |
| Takes          | `S: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]` |
| Properties     | `space` $S$, an `AbstractFunctionSpace`, or a `tuple` of at least two spaces                                               |
| Derived        | `domain`, `value_shape`, `scalar_type` and `argument_spaces`, below |
| Raises         | `TypeError` for an $S$ that is not a space or a `tuple` of spaces, and `ValueError` for a `tuple` of fewer than two spaces |

`dual_space(S)` provides:

- The domain `None` and the value shape `()`.
- The argument spaces `(S',)`, with $S'$ the conjugate of $S$ as in
  [`Conjugate`](#conjugate), member by member for a `tuple` $S$.
- The scalar type of $S$, complex if any space of a `tuple` is complex.

#### Antilinear forms

A function in `dual_space(S)` is an antilinear form on $S$, a linear form on
`conjugate_space(S)` (see [`ConjugateSpace`](#conjugatespace)). A form is
antilinear in the same way in an argument bound as `conjugate(v)`, such as the
argument $v$ of $\int u \bar v$. A residual bound as `(conjugate(v),)` and an
assembled load therefore have the same argument spaces, and the difference of
the two is a form.

The space of linear forms on a complex $S$ is
`dual_space(conjugate_space(S))`, with the argument spaces `(S,)`. The core
therefore fixes no convention on the argument in which a form is antilinear.

The pairing of a function `l` in `dual_space(S)` with a function `w` in $S$ is
`action(l, 0, conjugate(w))`, because the argument of `l` is in
`conjugate_space(S)` and every action is linear in its argument (see
[`Action`](#action)). The library that provides the values of `l` evaluates
the pairing.

#### The dual of a product

For a `tuple` $S$, a function in `dual_space(S)` is one antilinear form on the
product of the spaces, such as the load of a mixed problem. No single function
is in the product.

#### The double dual

`dual_space(dual_space(V))` is the space of antilinear forms on
`dual_space(V)`, with the argument spaces
`(conjugate_space(dual_space(V)),)`. `dual_space(dual_space(V))` is not `V`,
because `V` has the argument spaces `()`.

A function in `dual_space(dual_space(V))` is a leaf, such as data of a library.
A form never binds a function in a dual space (see [`Form`](#form)), so the
canonical map from `V` is not a form of the core.

#### Examples

The examples use the spaces `V` and `Q` and the functions `v` in `V` and `q`
in `Q` of the [`Form`](#form) examples.

- `l = function(dual_space(V), "l")` has the argument spaces
  `(conjugate_space(V),)`. `action(l, 0, conjugate(v))` is the pairing
  $\langle l, v \rangle$.
- `m = function(dual_space((V, Q)), "m")` has the argument spaces
  `((conjugate_space(V), conjugate_space(Q)),)`.
  `action(m, 0, (conjugate(v), conjugate(q)))` is the pairing of `m` with
  `(v, q)`.
- With `l` in `dual_space(V)`, `form(conjugate(action(l, 0, conjugate(v))), (conjugate(l),))`
  raises `ValueError`, because `l` has argument spaces.

```{dropdown} Legacy UFL
The legacy UFL `Cofunction` and `Coargument` are a fixed and a bound function in `dual_space(V)`. In legacy UFL, the dual of the dual space is $V$. In UFLx, the spaces keep distinct argument spaces. A form of the core binds no `Coargument`, and an operator with a `Coargument` is a node class of a language extension.
```

### `Form`

Requires: [Argument spaces](#argument-spaces), [`Function`](#function).

A concrete subclass of `AbstractExpression` for a variational form, an
expression with some of its functions bound as arguments.

| `Form(expression, arguments)` |                                                                                                                                                                                                    |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                       | `AbstractExpression` |
| Free function                 | `form(expression, arguments)`                                                                                                                                                                      |
| Takes                         | `expression: AbstractExpression`, `arguments: tuple[Function \| tuple[Function, ...], ...]` |
| Properties                    | `expression`, an expression, and `arguments`, a `tuple` of entries, each a [function](#function) `v`, the [conjugate](#conjugate) `conjugate(v)` of a function, or a `tuple` of at least two of these |
| Requires                      | Every function in `arguments` distinct, with the argument spaces `()` |
| Static type                   | That of the expression, with the argument spaces of the expression followed by the space $V$ of each entry `v`, `conjugate_space(V)` of each entry `conjugate(v)`, or the `tuple` of the spaces of the members of a `tuple` entry, in order |
| Raises                        | `TypeError` for an entry that is not a `Function`, a `Conjugate` of a `Function` or a `tuple` of these, and `ValueError` for a repeated function, a function with argument spaces or a `tuple` of fewer than two members |

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
of one form and fixed in another. The expression may itself be a form, and the
form then binds more arguments, after the arguments of the expression, as the
[Gateaux derivative](#gateaux-derivative) does. No pass of the core flattens
`Form(Form(e, a), b)` into `Form(e, a + b)`, because a product of forms would
need a reordering of the arguments, so the code that assembles reads the
order of the arguments from the static type. A form never binds a function in a dual space,
so forms are first order, and an action never substitutes a form into a pairing.
An operator with an argument in a dual space, such as the adjoint of an
interpolation, is a node class of a language extension that states the
argument spaces (see [Language extensions](#language-extensions)). The scope of a bound function is lexical, as
the scope of an index is (see [`ComponentTensor` and
`IndexSum`](#componenttensor-and-indexsum)): an occurrence of a function refers
to the innermost enclosing `Form` that binds the function, and is fixed
otherwise.

An argument that is a `tuple` of functions gives the forms of mixed problems
with no mixed space. Each member is a function in an ordinary space, so the
members may have different value shapes, such as a Raviart–Thomas and a
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
- With `r = integral(f * conjugate(v), D)` and `l = function(dual_space(V), "l")`,
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

```{dropdown} Legacy UFL
In legacy UFL, the arguments of a form are the test functions and the trial functions in the integrand, with the argument numbers of the functions, and a form is not an expression. In UFLx, a function has no argument number, the form states the arguments, and a functional is any expression with the argument spaces `()`. The legacy UFL `MixedElement` and `split` are an argument that is a `tuple` of functions.
```

### `Action`

Requires: [`Form`](#form).

See also: [Applying actions](#applying-actions).

A concrete subclass of `AbstractExpression` for the action of a form, an
expression given to one argument of the form.

| `Action(f, position, g)` |                                                                                                                                                        |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Extends                  | `AbstractExpression` |
| Free function            | `action(f, position, g)`                                                                                                                               |
| Takes                    | `f: AbstractExpression`, `position: int`, `g: AbstractExpression \| tuple[AbstractExpression, ...]` |
| Properties               | `form`, an expression with more than `position` argument spaces, `position`, an `int`, and `argument`, the expression $g$, or a `tuple` of expressions |
| Requires                 | The conditions below                                                                                                                                   |
| Static type              | That of $f$, with the free indices of $f$ followed by those of $g$, each `Index` taken once, and the argument spaces of $f$ without the entry at `position` |
| Raises                   | `ValueError` if a requirement fails                                                                                                                    |

Let the argument space at `position` be $V$ or `ConjugateSpace(V)`.

- For a space $V$, $g$ has the value shape and the argument spaces of a
  function in $V$, the scalar type of that function or the real scalar type,
  and the support of that function or no support. The argument spaces are `()`, or `(S',)` for $V$ =
  `dual_space(S)` (see [`DualSpace`](#dualspace)).
- For `ConjugateSpace(V)`, $g$ has the same static type, and $g$ is the value
  of an entry `conjugate(v)`, such as `conjugate(w)` for a function `w` in $V$.
- For a `tuple` of spaces, $g$ is a `tuple` with one such member per space.
- An `Index` free in both $f$ and $g$ has the same extent in both.

$g$ may have free indices, and is then a family of arguments, one for each
value of the free indices, as a domain map with free indices is a family of
domains (see [`DomainMap`](#domainmap)). The action is then the family of
actions. The free indices of $g$ are the free indices of the members of a
`tuple`.

For a `position` below the arity $n$ of $e$, the value of
`Action(Form(e, arguments), position, g)` is
`Form(Action(e, position, g), arguments)`. Otherwise, with the entry at
`position` the entry at `position - n` of `arguments`, the value of
`Action(Form(e, arguments), position, g)` is `Form(e, others)`,
with `others` the entries other than the entry at `position`, and with $g$
substituted for the argument at `position` (see [Passes](#passes)). A function
of `others` that is free in $g$ is therefore renamed. If no entries remain, the
value is $e$ with $g$ substituted. For an entry
`conjugate(v)`, `conjugate(g)` is substituted for `v`, so $g$ is the value of
the entry and the action is linear in $g$ at every position. For a `tuple` entry, the members
of $g$ are substituted simultaneously for the matching members. The action of any other form is
defined by linearity, and the action on a leaf, such as a function in a dual
space, is the pairing, which the library of the leaf states, with one pairing
for each value of the free indices of $g$.

#### Examples

The examples use `a`, `u`, `v` and `w` of the [`Form`](#form) examples, and
`l = function(dual_space(V), "l")`.

- `action(form(a, (v, u)), 1, w)` is a linear form in `v`, with the value
  `form(b, (v,))`, with `b` the expression `a` with `u` replaced by `w`.
- `action(action(form(a, (v, u)), 1, w), 0, w)` is a functional.
- `action(form(a, (v, u)), 1, v)` is a linear form in a renamed function, with
  the fixed `v` in place of `u`, because the form binds `v`.
- `action(l, 0, conjugate(w))` is the pairing of `l` with `w`, with no support.
- `action(form(a, (conjugate(v), u)), 0, conjugate(w))` is the linear form `a`
  with `w` in place of `v`.
- `action(form(a, (v, u)) + 2 * form(a, (v, u)), 1, w)` is, by linearity,
  the sum of the actions on both operands, with `2` times the action on the
  form of the product.
- With fixed functions `w0` and `w1` in `V` and an index `k`,
  `action(form(a, (v, u)), 1, tensor((w0, w1))[k])` is a family of two linear
  forms in `v`, with the free index `k`, the actions on `w0` and on `w1`.
- Let `w` in `V` and `r` in `Q` be fixed functions, with the Stokes form of the
  [`Form`](#form) examples. `action(form(e, ((v, q), (u, p))), 1, (w, r))` is
  a linear form in `(v, q)`, the residual at `(w, r)`.

```{dropdown} Legacy UFL
Legacy UFL operations are actions and rebindings. The action on the last argument of a bilinear form is `action(f, 1, w)`. The adjoint of `form(e, (conjugate(v), u))` is `form(conjugate(e), (conjugate(u), v))`, because `conjugate(e)` is linear in `v` and antilinear in `u`. The Gateaux derivative with respect to a function adds an argument (see [Gateaux derivative](#gateaux-derivative)).
```

Worked examples for this chapter in EXAMPLES.md: [the Laplace-Beltrami operator](EXAMPLES.md#the-laplace-beltrami-operator), [the flux through a surface](EXAMPLES.md#the-flux-through-a-surface), [Helmholtz with a perfectly matched layer](EXAMPLES.md#helmholtz-with-a-perfectly-matched-layer), [mixed Poisson with Raviart-Thomas elements](EXAMPLES.md#mixed-poisson-with-raviart-thomas-elements), [dual-weighted residual error estimation](EXAMPLES.md#dual-weighted-residual-error-estimation), [a classic FEniCS layer for the Poisson problem](EXAMPLES.md#a-classic-fenics-layer-for-the-poisson-problem).

### Gateaux derivative

Requires: [`Form`](#form), [`Function`](#function),
[`DirectionalDerivative`](#directionalderivative),
[`ConjugateSpace`](#conjugatespace), [`DualSpace`](#dualspace).

See also: [Applying Gateaux derivatives](#applying-gateaux-derivatives).

A concrete subclass of `AbstractExpression` for the Gateaux derivative of an
expression $F$ with respect to a variable $x$ at a point $a$, which adds an
argument $w$ after the arguments of $F$.

| `GateauxDerivative(F, x, a, w)` |  |
| ------------------------------- | - |
| Extends     | `AbstractExpression` |
| Takes       | `F: AbstractExpression`, `x: Function \| tuple`, `a: AbstractExpression \| tuple`, `w: Function \| Conjugate \| tuple` |
| Properties  | `expression` $F$, `variable` $x$, `point` $a$ and `argument` $w$, an entry as in [`Form`](#form) |
| Requires    | $x$ a `Function` or a `tuple` of functions with the argument spaces `()`, $a$ and the functions of $w$ with the value shapes of $x$, member by member, the members of a `tuple` $w$ all functions or all conjugates, no function of $w$ free in $F$, the members of a `tuple` $x$ all real or all complex, and each member of $a$ real for a real $x$ |
| Static type | The value shape of $F$, the free indices of $F$ followed by those of $a$, the argument spaces of $F$ followed by the space of the entry $w$, and the scalar type of $F$ for a real $x$, and complex for a complex $x$ |
| Raises      | `TypeError` for an $x$ that is not a function, and `ValueError` if another requirement fails |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `gateaux_derivative(F, x, a, w)` | | `F: AbstractExpression`, `x: Function \| tuple`, `a: AbstractExpression \| tuple`, `w: Function \| Conjugate \| tuple` | `GateauxDerivative(F, x, a, w)` |

The members of a `tuple` $x$ have one scalar type, the scalar type of $x$.
The linear and the antilinear parts below are defined only on a complex
space, and the product of a real and a complex space is a real space, so a
`tuple` that mixes the two is rejected.

```{dropdown} Open: mixed real and complex tuples
A complex state with a real unknown, such as a Lagrange multiplier, needs a
mixed `tuple`. A convention could put the part in the real members into
$J_+$. Until then, take the derivatives separately.
```

The real differential in a direction $v$ is
$\delta F[v]$ = `directional_derivative(F, x, a, v)`, with the members of a
`tuple` variable substituted at once (see
[`DirectionalDerivative`](#directionalderivative)). $\delta F$ is
real-linear in $v$, and the value of the Gateaux derivative binds $w$ by the
scalar type of the space of $x$:

| Space and entry | Value |
| --------------- | ----- |
| A real space | `Form(δF[w], (w,))` |
| A complex space, entry `w` | `Form(J₊, (w,))`, with $J_+[w] = \tfrac12 (\delta F[w] - i \, \delta F[i w])$, the linear part |
| A complex space, entry `conjugate(w)` | `Form(J₋, (conjugate(w),))`, with $J_-[w] = \tfrac12 (\delta F[w] + i \, \delta F[i w])$, the antilinear part |

$\delta F = J_+ + J_-$, and the two parts are unique. For a real functional
$\Pi$ on a complex space, $J_-$ with the entry `conjugate(w)` is the
Wirtinger derivative $\partial \Pi / \partial \bar u$, a function in
`dual_space(V)` (see [Antilinear forms](#antilinear-forms)), so $J_-$ minus a
load is a residual. If no occurrence of $x$ in $F$ is under a `Conjugate`, an
`Absolute`, a `RealPart` or an `ImaginaryPart`, $F$ is holomorphic in $x$,
$J_+$ is $\delta F$ and $J_-$ is zero. The Jacobian of the residual of a
sesquilinear problem, such as the Helmholtz equation, is such a $J_+$.
Otherwise the Jacobian is the pair of forms $J_+$ and $J_-$, which do not add
because the argument spaces differ, and Newton's method solves
$J_+ \, \delta u + J_- \, \overline{\delta u} = -R$ with an algorithm that the
solver chooses.

The value is a form that binds $w$ after the arguments of $F$, so the
derivative of a residual form is the tangent form, and an action on the new
argument is the derivative in a fixed direction. The second derivative of a
real functional is symmetric, or Hermitian, and the graph does not record the
symmetry.

## Evaluation

Requires: [Passes](#passes).

A lazy node class is a node class whose section states the value of the node,
and which a pass may rewrite as an expression of other node classes with that
value. The lazy node classes of the core are `GateauxDerivative`, `Action`,
`DirectionalDerivative`, `Integral` over a decomposed domain, `PullBack` and
`Gradient`.

A pass that rewrites a node keeps the value of the node and the static type of
the node, with the free indices possibly in another order, never simplifies,
and substitutes with the renaming of [Passes](#passes). The one exception to
the static type is `Zero`, which is real and has no support: where a rewrite
produces `Zero`, the scalar type of the result may be real for a complex node,
and the result may have no support for a node with a support. The value is
the same, because real arithmetic embeds in complex arithmetic and a constant
is a field on every domain. Because the pass never
simplifies, a term that is zero stays in the graph, and zero times an
undefined value is undefined. So a rule never multiplies a term that may be
zero by an expression that may be undefined where the rewritten node is
defined, such as the logarithm in the derivative of a power (see
[Derivatives of the node classes](#derivatives-of-the-node-classes)).

`Action` and `DirectionalDerivative` are real-linear, and each commutes with
the following node classes in the listed operands, when the other operands do
not depend on what the evaluated node changes:

| Node class | Operands |
| ---------- | -------- |
| `Sum`, `Tensor` | Every operand |
| `Indexed`, `ComponentTensor`, `IndexSum`, `Form`, `DirectionalDerivative` | The expression |
| `Conjugate`, `RealPart`, `ImaginaryPart` | The operand |
| `Division` | The numerator |
| `Gradient`, `PullBack`, `Integral` | The expression |

The value of a commuting node is real-linear in the listed operands, and the
evaluated node acts on the value of an operand alone: an action substitutes
into a bound argument, which is a field like any other, and a directional
derivative differentiates in a real parameter with the domains fixed. Moving
an evaluated node into a binding node first makes the binding node fresh for
the other operands of the evaluated node (see [Passes](#passes)). The partial
derivative of [`apply_gradients`](#applying-gradients) commutes with the same
node classes except `Gradient`, `PullBack` and `Integral`, which the rows of
that pass match first.

The passes of the core run in this order, each on the output grammar of the
one before, from $N_0$ (see [Grammars](#grammars)):

| Pass | Evaluates | Output grammar |
| ---- | --------- | -------------- |
| [`apply_gateaux_derivatives`](#applying-gateaux-derivatives) | `GateauxDerivative` | $N_1$ |
| [`apply_actions`](#applying-actions) | `Action`, except a pairing | $N_2$ |
| [`apply_derivatives`](#applying-derivatives) | `DirectionalDerivative` | $N_3$ |
| [`expand_integrals`](#expanding-integrals) | `Integral` over a decomposed domain | $N_4$ |
| [`apply_gradients`](#applying-gradients) | `Gradient`, except of a terminal, a `PullBack` or a `Gradient` | $N_5$ |
| [`expand_pull_backs`](#expanding-pull-backs) | `PullBack`, except one that reads a value of a library | $N_6$ |
| [`apply_gradients`](#applying-gradients) | The gradients that `expand_pull_backs` creates, such as the Jacobian of a cell | $N_7$ |
| [`expand_pull_backs`](#expanding-pull-backs) | The pull-backs that the second `apply_gradients` uncovers | $N_8$, [the core normal form](#the-core-normal-form) |

The first `apply_gradients` evaluates each gradient on the domain on which the
form states the gradient, before a pull-back by a cell moves the gradient to
the reference cell. A gradient of an expression of the coordinates of a
decomposed domain $D$ is therefore evaluated on $D$, such as
`gradient(coordinate(D), D)`, which is the identity matrix for $D$ of the
dimension of its ambient coordinate domain, and not on the reference cell as
the product of the Jacobian $J$ and its inverse $K$. The
second `apply_gradients` evaluates the gradients that the pull-backs create on
the reference cells, and the second `expand_pull_backs` moves to the leaves
the pull-backs that enclosed those gradients.

A rewrite of index notation that keeps the value, such as the reduction of an
`Indexed` of a `ComponentTensor`, is a pass of the compiler, not of the core.

### Grammars

Requires: [Evaluation](#evaluation), [Argument spaces](#argument-spaces).

See also: [New questions](#new-questions), [Rule tables](#rule-tables).

A grammar is a set of expressions. A grammar is checked by a pass that asks
whether an expression belongs to the grammar (see [New
questions](#new-questions)), built from rule tables as any other pass (see
[Rule tables](#rule-tables)).

- An expression node of an expression $e$ is an expression reached from $e$
  along the edges of the graph, including $e$ (see [The base object
  `Base`](#the-base-object-base)). The grammar $N_0$ is the set of the
  expressions $e$ in which every expression node is an instance of a concrete
  class of the core or a terminal of a downstream library (see [Leaves and
  subclassing](#leaves-and-subclassing)). Every object is valid, so $N_0$
  needs no other condition. A node class of a language extension is not in
  $N_0$, and an extension that adds node classes states the grammars that
  contain them.
- A pass maps a source grammar into a target grammar if, for every
  expression of the source grammar, the result of the pass belongs to the
  target grammar.
- A pattern is a condition on one expression node that reads only the node
  and the objects within a bounded number of edges of the node. A pattern
  never reads the result of a pass.
- For a grammar $N$ and a set of patterns $P$, the delta $N \setminus P$ is
  the set of the expressions of $N$ with no expression node that matches a
  pattern of $P$. An intermediate grammar is written as a delta of an earlier
  grammar.
- The core normal form is the last grammar of the passes of the core, written
  as a list of the node classes and the patterns that the grammar allows. A
  form compiler receives the integrand layer of an expression in the core
  normal form and continues with passes of its own, and the code that
  assembles receives the form layer.

A grammar may restrict two layers of an expression separately. The form layer
consists of every `Action` and every node with argument spaces other than
`()`. The integrand layer consists of the other nodes.

The section of each pass states the input grammar and the output grammar of
the pass, such as $N_1$ of
[`apply_gateaux_derivatives`](#applying-gateaux-derivatives), and [The core
normal form](#the-core-normal-form) lists what the last of them allows.

### Applying Gateaux derivatives

Requires: [Grammars](#grammars), [Gateaux derivative](#gateaux-derivative).

| `apply_gateaux_derivatives(f)` |  |
| ------------------------------ | - |
| Evaluates                      | `GateauxDerivative` |
| Input grammar                  | $N_0$ |
| Output grammar                 | $N_1 = N_0 \setminus \{$a `GateauxDerivative`$\}$ |

The pass rewrites every `GateauxDerivative(F, x, a, w)`, inner derivatives
first, as a `Form` that binds $w$ after the arguments of $F$, with
$\delta F[v]$ the directional derivative of [Gateaux
derivative](#gateaux-derivative):

| Space of $x$ and entry $w$ | Rewrite |
| -------------------------- | ------- |
| A real space | `Form(δF[w], (w,))` |
| A complex space, $F$ holomorphic in $x$, entry `w` | `Form(δF[w], (w,))` |
| A complex space, $F$ holomorphic in $x$, entry `conjugate(w)` | `Zero` with the static type of the derivative |
| A complex space, entry `w` | `Form(Integer(1) / Integer(2) * (δF[w] - ImaginaryUnit() * δF[ImaginaryUnit() * w]), (w,))` |
| A complex space, entry `conjugate(w)` | `Form(Integer(1) / Integer(2) * (δF[w] + ImaginaryUnit() * δF[ImaginaryUnit() * w]), (conjugate(w),))` |

- $F$ is holomorphic in $x$ if no occurrence of a function of $x$ in $F$ is
  under a `Conjugate`, an `Absolute`, a `RealPart` or an `ImaginaryPart`. The
  rows for a holomorphic $F$ give $J_+$ and $J_-$ with no terms that cancel.
  An action is linear in its argument (see [`Action`](#action)), so an
  action never hides a conjugation, and the four node classes are the only
  ones antilinear in an operand.
- For a `tuple` $w$, `ImaginaryUnit() * w` is the `tuple` of the products of
  `ImaginaryUnit()` with each member.
- The pass runs first, so the forms of $F$ may still be under actions, and
  [`apply_actions`](#applying-actions) substitutes an argument for $w$ like
  any other argument. An action on $w$ is then the derivative in a fixed
  direction.
- The terms of the last two rows are each real-linear in $w$, and are linear
  or antilinear only in sum, so a linearity check that classifies each term
  does not confirm the linearity of these rows.

### Applying actions

Requires: [Applying Gateaux derivatives](#applying-gateaux-derivatives), [`Action`](#action), [`Form`](#form),
[Passes](#passes).

| `apply_actions(f)` |  |
| ------------------ | - |
| Evaluates          | `Action` |
| Input grammar      | $N_1$ |
| Output grammar     | $N_2 = N_1 \setminus \{$an `Action` whose form is not a `Function` in a dual space$\}$ |

The pass rewrites every `Action(f, p, g)`, inner actions first, with the first
row that matches the form $f$, and keeps the value of [`Action`](#action).

| Form $f$ | Rewrite of `Action(f, p, g)` |
| -------- | ---------------------------- |
| `Form(e, arguments)`, $p$ at least the arity $n$ of $e$ | The value of `Action`: `Form(e, others)` with $g$ substituted for the argument at $p - n$, or $e$ with $g$ substituted if no entries remain |
| `Product(a, b)`, with $n$ the arity of $a$ | `Product(Action(a, p, g), b)` for $p < n$, otherwise `Product(a, Action(b, p - n, g))` |
| A commuting node class (see [Evaluation](#evaluation)) | The node with the action moved into the listed operands |
| `Zero` | `Zero` with the free indices of $g$ added and the entry at $p$ removed |
| `Function` in a dual space | Remains, a pairing of $N_2$ |
| `Action` | Cannot occur, because the inner action is rewritten first, and an action that remains has the argument spaces `()` |
| Any other class of the core | Cannot occur, because the argument spaces are `()` |

- `DirectionalDerivative` binds its variable, so moving the action into it
  renames the variable if the variable is free in $g$, and the $x$ in $g$ is
  not differentiated.
- Forms are first order (see [`Form`](#form)), so a substitution never creates
  an `Action` on a `Form`. Each rule removes the action from one node or moves
  the action to operands that are smaller, so the pass terminates.

### Derivatives of the node classes

Requires: [Evaluation](#evaluation), [Elementary
functions](#elementary-functions-very-incomplete), [`Inverse`](#inverse),
[`Action`](#action).

The passes that differentiate share one table, the derivative of each node
class in terms of the derivatives of its operands. In the table, $\delta b$ is
the derivative of the operand $b$ in the pass, and an operand is constant for
the pass if the first row of the pass matches the operand. Each pass matches
its own rows first, for the leaves of the expression, and then the rows below.

| Expression | Derivative |
| ---------- | ---------- |
| `Product(b, c)` | $\delta b \, c + b \, \delta c$ |
| `Division(b, c)` | $\delta b / c - (b \, \delta c) / (c c)$ |
| `Power(b, c)`, $c$ constant | $c \, b^{c - 1} \, \delta b$ |
| `Power(b, c)` | $b^c (\log(b) \, \delta c + c \, \delta b / b)$ |
| `Absolute(b)` | $\operatorname{Re}(\bar b \, \delta b) / \lvert b \rvert$, with `real_part`, so the derivative is real |
| `SquareRoot(b)` | $\delta b / (2 \sqrt b)$ |
| `Logarithm(b)` | $\delta b / b$ |
| `Inverse(A)` | $-A^+ \, \delta A \, A^+ + A^+ (A^+)^T (\delta A^T - \delta A^T A A^+)$, with `@`, `transpose` and `-` |
| `Action(l, 0, w)`, a pairing | `Action(l, 0, δw)` |
| A commuting node class (see [Evaluation](#evaluation)) | The node with the derivative moved into the listed operands |

- `Power` has two rows. The general row multiplies by $\log(b)$, which is
  undefined for a negative base, so the row for a constant exponent avoids
  the logarithm. That row is undefined at a zero base when the value of the
  exponent is 0, which occurs when an integer power is differentiated more
  times than its degree, and is the one known exception to the principle of
  [Evaluation](#evaluation).
- A commuting row applies only to an operand that no row of the pass
  matches first, so each pass decides through its own rows which node classes
  the derivative moves into, as [`apply_gradients`](#applying-gradients) does
  for `Gradient`, `PullBack` and `Integral`.
- A language extension that adds a node class adds the derivative of the
  class to this table, and every pass that differentiates then uses the row.
- The derivative of each row is a composition of node classes of the core,
  so a pass that applies the table terminates.

### Applying derivatives

Requires: [Applying actions](#applying-actions),
[`DirectionalDerivative`](#directionalderivative), [Derivatives of the node
classes](#derivatives-of-the-node-classes).

| `apply_derivatives(f)` |  |
| ---------------------- | - |
| Evaluates              | `DirectionalDerivative` |
| Input grammar          | $N_2$ |
| Output grammar         | $N_3 = N_2 \setminus \{$a `DirectionalDerivative`$\}$ |

The pass rewrites every `DirectionalDerivative(f, x, a, h)`, inner derivatives
first, with the first row that matches $f$. $\mathrm{D} b$ is short for
`DirectionalDerivative(b, x, x, h)`. The rules give the derivative at $x$, and
if $a$ is not $x$, the pass then substitutes $a$ for $x$ in the result, every
member of a `tuple` at once.

| Expression $f$ | Rewrite of the derivative at $x$ |
| -------------- | -------------------------------- |
| No function of $x$ free in $f$ | `Zero`, with the value shape and the argument spaces of $f$, and the free indices of $f$ followed by those of the members of $a$ and $h$ |
| $f$ equal to the member $x_k$ of $x$ | $h_k$, or $h$ for a function $x$ |
| Any other expression | The row of [Derivatives of the node classes](#derivatives-of-the-node-classes), with $\delta b = \mathrm{D} b$ |
| `DirectionalDerivative` | Cannot occur, because the inner derivative is rewritten first |

- The first row comes first, so no derivative is built for an operand that
  does not depend on $x$, and a terminal other than a function of $x$ needs no
  row. An operand is constant for the pass if no function of $x$ is free in
  the operand.
- A `Form` that binds $x$ has no free $x$, so the first row gives `Zero`.
- The variable has the argument spaces `()`, so the variable is never the
  function $l$ of a pairing.

### Expanding integrals

Requires: [Applying derivatives](#applying-derivatives), [`Integral` over a
decomposed domain](#integral-over-a-decomposed-domain), [Jacobian and
density](#jacobian-and-density).

| `expand_integrals(f)` |  |
| --------------------- | - |
| Evaluates             | `Integral` over a decomposed domain |
| Input grammar         | $N_3$ |
| Output grammar        | $N_4 = N_3 \setminus \{$an `Integral` whose domain is a `DecomposedDomain`$\}$ |

The pass rewrites every `Integral(f, D)` with a decomposed domain $D$, inner
integrals first, with the rows below, and applies the rows again to each
`Integral` that a rewrite creates. For piece $j$ of $D$, $P_j$ is the parameter
domain, `indices_j` are the free indices, and `D_j = cell_domain(D, j)`.

| Integral | Rewrite |
| -------- | ------- |
| `Integral(f, D)`, $f$ of rank $r \ge 1$ | `component_tensor(Integral(f[K], D), K)`, with `K` one derived `Index` for each axis |
| `Integral(f, D)`, $f$ scalar | The sum over the pieces $j$ of `index_sum(integral(g_j * density(D_j), P_j), indices_j)`, with no `IndexSum` for a piece with no free indices |

$g_j$ is `pull_back(f, piece_j)` for $f$ with the support `D.ambient_domain`,
and `pull_back(f, D_j)` otherwise. An $f$ with no support is pulled back too,
so that $g_j$ has the free indices of piece $j$ when the density does not,
such as for a family of points, whose density is `Integer(1)`.

- Each application moves to a parameter domain that is nested less deeply, so
  the pass terminates.
- A free index of $f$ shared with piece $j$ is identified by the static type of
  `PullBack`, and the `IndexSum` pairs each value of the index with a member of
  the family, so the result has the static type of the `Integral`.
- The density is the same for `cell_domain(D, j)` and for piece $j$, because
  both have the parameterisation $\psi_j$.
- The integrals that remain are over coordinate domains, which a form compiler
  integrates by quadrature.

### Applying gradients

Requires: [Expanding integrals](#expanding-integrals),
[`Gradient`](#gradient), [Derivatives of the node
classes](#derivatives-of-the-node-classes).

| `apply_gradients(f)` |  |
| -------------------- | - |
| Evaluates            | `Gradient` |
| Input grammar        | $N_4$, and $N_6$ in the second run |
| Output grammar       | $N_5 = N_4 \setminus G$, and $N_7 = N_6 \setminus G$ in the second run, with $G = \{$a `Gradient(f, P)` whose expression $f$ is not a terminal, a `PullBack` or a `Gradient`, or is `Coordinate(P)` with $P$ of the dimension of its ambient coordinate domain$\}$ |

The pass rewrites every `Gradient(f, P)` that $G$ matches, inner
gradients first, as
`component_tensor(∂_j f[a_1, …, a_r], (a_1, …, a_r, j))`, with $r$ the rank
of $f$ and derived indices $a_1, \ldots, a_r$ and $j$. $\partial_j b$, entry
$j$ of the last axis of the gradient of $b$ on $P$, has the value shape of
$b$ and the free index $j$ of the extent $n$, the dimension of the ambient
coordinate domain of $P$. $\partial_j$ applies the first row that matches:

| Expression $b$ | $\partial_j b$ |
| -------------- | -------------- |
| $b$ with no support | `Zero`, with the value shape and the argument spaces of $b$, and the free indices of $b$ followed by $j$ |
| `Coordinate(P)`, $P$ of the dimension of its ambient coordinate domain | `component_tensor(I[i, j], (i,))`, with $I$ the identity matrix of value shape `(n, n)` and a derived index $i$ |
| A terminal, a `PullBack` or a `Gradient` | `component_tensor(gradient(b, P)[a_1, …, a_s, j], (a_1, …, a_s))`, with $s$ the rank of $b$, so the gradient remains |
| Any other expression | The row of [Derivatives of the node classes](#derivatives-of-the-node-classes), with $\delta b = \partial_j b$ |

- An operand is constant for the pass if the operand has no support.
- The gradient of a terminal remains for the library of the terminal, so the
  Hessian of a function is a gradient of a gradient.
- On a decomposed domain $P$ of lower dimension than its ambient coordinate
  domain, `Gradient` is the tangential gradient (see [`Gradient`](#gradient)),
  which obeys the same rows. The `Coordinate` row does not match, so
  `gradient(coordinate(P), P)`, the projection $\Pi$, remains, and a
  pull-back by a cell, such as through a restriction to points, evaluates
  it.
- In the first run, a gradient that $G$ does not match remains, such as the
  gradient of a function on a decomposed domain, and
  [`expand_pull_backs`](#expanding-pull-backs) rewrites the gradient on the
  reference cell. In the second run, the gradients that $G$ matches are the
  gradients that `expand_pull_backs` creates on the reference cells, such as
  the Jacobian of a cell in the inverse Jacobian and, in a Hessian, the
  gradient of the pulled-back gradient.
- A gradient on a cell under a pull-back to the points of a facet, or of
  another cell through the inclusion of a side, remains in $N_6$, because
  the gradient does not commute with that pull-back (see [Expanding
  pull-backs](#expanding-pull-backs)). Once the second run rewrites the
  gradient, the pull-back is pointwise again, as in the inverse Jacobian of
  the cell on a side and in the inverse Jacobian of [`Normal`](#normal), and
  the second `expand_pull_backs` moves the pull-back to the leaves.

### Expanding pull-backs

Requires: [Applying gradients](#applying-gradients), [`PullBack`](#pullback),
[Restriction](#restriction), [`Normal`](#normal), [Jacobian and
density](#jacobian-and-density), [`Inverse`](#inverse).

| `expand_pull_backs(f)` |  |
| ---------------------- | - |
| Evaluates              | `PullBack` |
| Input grammar          | $N_5$, and $N_7$ in the second run |
| Output grammar         | $N_6 = N_4 \setminus R$, and $N_8 = N_7 \setminus R$ in the second run, with $R = \{$a `PullBack` by a domain map whose expression is not a terminal, a `PullBack` or a `Gradient`, a `PullBack` of a `Gradient` by a `cell_domain` or a `gathered_cell_domain`, and a `PullBack` by a domain map of a `PullBack` by a domain map that is not a `cell_domain` or a `gathered_cell_domain`$\}$ |

The pass rewrites every `PullBack(f, ψ)`, inner pull-backs first, with the
first row that matches, and applies the rows again to each `PullBack` that a
rewrite creates. $P$ and $Q$ are the parameter domain and the ambient domain
of $\psi$, $J$ = `jacobian(ψ)`, and for $\psi$ = `cell_domain(E, i)` and a side
$s$, $B_i$ is the cell on the side and $A_i$ the inclusion of the cell of $E$
into it (see [Restriction](#restriction)).

| `PullBack(f, ψ)` | Rewrite |
| ---------------- | ------- |
| $\psi$ a `SideMap`, $f$ with a support | Remains until a pull-back by a cell of $E$ |
| $f$ with no support | $f$, or `f + Zero(…)` with the free indices of $\psi$ that $f$ lacks |
| `Coordinate(Q)` | $\psi$ |
| `PullBack(g, SideMap(D, E, s))`, $\psi$ = `cell_domain(E, i)` | `PullBack(PullBack(g, B_i), A_i)` |
| `PullBack(g, φ)`, φ a `SideMap`, a `cell_domain` or a `gathered_cell_domain` | Remains, so a pull-back nests only inside a read on a cell |
| `PullBack(g, φ)` | `PullBack(g, domain_map(P, pull_back(φ.parameterisation, ψ), φ.ambient_domain))`, the pull-back by the image of $\psi$ under φ (see [`PullBack`](#pullback)) |
| `Gradient(g, Q)`, $\psi$ = `cell_domain(Q, j)` or a `gathered_cell_domain` of $Q$ | $K^T \nabla_P(\psi^* g)$, with $K$ = `inverse_jacobian(ψ)` |
| `Gradient(g, Q)` | Remains, the gradient on $Q$ read at the points of $P$ |
| `Normal(SideMap(D, E, s))`, $\psi$ = `cell_domain(E, i)` | $K^T \hat m / \lvert K^T \hat m \rvert$ of [`Normal`](#normal) |
| A terminal | Remains, for the library of the terminal |
| Any other node class | The node of the pull-backs of the operands that are expressions |

- A pull-back is pointwise, so the last row moves the pull-back into every
  other node class, with a binding node made fresh for $\psi$ first.
- `Zero` is the one node class that carries free indices without operands, so
  `f + Zero(…)` is $f$ for each member of the family $\psi$.
- A restriction is a pull-back by a `SideMap` (see
  [Restriction](#restriction)), so this pass also lowers restrictions. A
  pull-back by a `SideMap` that no pull-back by a cell of $E$ encloses, such
  as on an expression on $E$ that is not integrated, remains.
- `inverse_jacobian` chooses the composition of $K$ for the dimensions of
  $\psi$, including a parameter domain of lower dimension than its ambient
  coordinate domain (see [Jacobian and density](#jacobian-and-density)).
- Only a pull-back by a cell moves a gradient on a decomposed domain to the
  reference cell, so the pass inverts the Jacobian of a cell only. A gradient
  under a pull-back by any other map remains on its domain. The inclusion
  $A_i$ of a cell of a submesh into a cell of the same dimension therefore
  leaves the reference gradient at the points of the cell of $E$, as the
  inclusion of a facet does.
- A pull-back by any other map composes with the pull-back that it encloses,
  so a value that a library reads on a cell is under at most one further
  domain map, such as the inclusion $A_i$, the points of a quadrature rule,
  or the two composed.
- Every pull-back that remains reads a value of a library on a cell, possibly
  at the points of a cell of $E$ through $A_i$, such as the reference
  gradient of a function on the cell at the points of a facet. A pull-back of
  a function or of its gradient by a map other than a cell domain, such as by
  a point of a mesh given by its coordinates, remains, because the library of
  the function defines the function on the cells. A point of a mesh known in
  advance is evaluated through a restriction instead, and other points are
  located by a pass after the core (see [Points and point
  evaluation](#points-and-point-evaluation)).
- The rows create gradients on the reference cells, such as the Jacobian in
  $K$, so $N_6$ is a delta of $N_4$ and not of $N_5$. In the second run,
  every gradient is of a terminal, a `PullBack` or a `Gradient`, and every
  pull-back by a cell already reads a value of a library, so the pass creates
  no gradient that $G$ matches and $N_8$ is contained in $N_7$.
- The pass leaves an expression of $N_6$ unchanged, because each pull-back of
  $N_6$ matches a row that remains, so the second run rewrites only the
  pull-backs that the second `apply_gradients` uncovers.

### The core normal form

Requires: [Grammars](#grammars), [Expanding pull-backs](#expanding-pull-backs).

The core normal form is $N_8$, the output grammar of the second run of
[`expand_pull_backs`](#expanding-pull-backs), the last pass of the core. Every
expression node of an expression in the core normal form is an instance of a
concrete class of the core or a terminal of a downstream library. The table
gives every node class that the core normal form allows, with the condition
on each, and no other node class occurs in the core normal form:

| Chapter | Node class | Allowed in the core normal form |
| ------- | ---------- | ------------------------------- |
| [Algebra](#algebra) | `Integer`, `Float`, `ImaginaryUnit`, `Zero` | Always |
| [Algebra](#algebra) | `Sum`, `Product`, `Division`, `Power`, `Absolute`, `Conjugate`, `RealPart`, `ImaginaryPart`, `SquareRoot`, `Logarithm` | Always |
| [Algebra](#algebra) | `Tensor`, `Indexed`, `ComponentTensor`, `IndexSum` | Always, including index notation that a rewrite could reduce, such as an `Indexed` of a `ComponentTensor`, which the compiler rewrites (see [Evaluation](#evaluation)) |
| [Algebra](#algebra) | `Inverse` | Always, for the compiler to evaluate |
| [Calculus on a domain](#calculus-on-a-domain) | `Coordinate` | Always |
| [Calculus on a domain](#calculus-on-a-domain) | `Gradient` | Of a terminal, a `PullBack` or a `Gradient`, other than `Coordinate(P)` on its own domain $P$ of the dimension of its ambient coordinate domain |
| [Calculus on a domain](#calculus-on-a-domain) | `Integral` | Over a domain that is not a `DecomposedDomain`, such as a reference cell or `RD(0)` |
| [Maps between domains](#maps-between-domains) | `PullBack` | By a `SideMap`, or by a domain map of a terminal, of a `PullBack` by a `SideMap`, a `cell_domain` or a `gathered_cell_domain`, or of a `Gradient` by a domain map that is not a `cell_domain` or a `gathered_cell_domain` |
| [Decomposed domains](#decomposed-domains) | `Normal` | Always |
| [Function spaces](#function-spaces) | `Function` | Always |
| [Variational forms](#variational-forms) | `Form` | Always, including a form whose expression is a form (see [`Form`](#form)) |
| [Variational forms](#variational-forms) | `Action` | A pairing, an `Action` whose form is a `Function` in a dual space |
| | A terminal of a downstream library | Always |

The objects of an expression that are not expressions, such as an `Index`, an
`EllipsisIndex`, a `GatheredIndex`, a domain, a `DomainMap`, a `SideMap` and a
space, occur as in any expression, because a grammar restricts only the
expression nodes.

The code that assembles receives the form layer (see [Grammars](#grammars)):
sums and products of forms, `Zero`, forms whose expression may be a form,
with the order of the arguments from the static type (see [`Form`](#form)),
and pairings, which the library of the function in the dual space evaluates.

A form compiler receives the integrand layer: integrals over reference cells
and over points, whose integrands are arithmetic and index notation of
constants, of the coordinates of the reference cell, and of the values that
the pull-backs read from libraries on the cells. Such a value is a function
pulled back by a cell, the reference gradient of one, or either at the points
of a cell of a topology with sides through the inclusion of a side (see
[Expanding pull-backs](#expanding-pull-backs)).

In particular, the core normal form may still contain:

- A pull-back of a function by points known only after a pass of the
  compiler, which a pass after the core locates (see [Points and point
  evaluation](#points-and-point-evaluation)).
- A pull-back by a `SideMap` that no pull-back by a cell encloses, such as a
  restriction that is not integrated.
- A `Gradient` of the coordinates of a decomposed domain of lower dimension
  than its ambient coordinate domain, the projection $\Pi$, and a `Gradient`
  of a function on a decomposed domain, both under a pull-back by points or
  by another map that is not a cell, such as the coordinates of another
  mesh.

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

```{dropdown} Legacy UFL
In legacy UFL, properties are stored on the node classes, and the set of properties is fixed. Every new property touches every node class, so features such as complex numbers and dual spaces were patched in across the node classes.
```

### Language extensions

Requires: [Leaves and subclassing](#leaves-and-subclassing), [Passes](#passes),
[Evaluation](#evaluation), [Transpose](EXAMPLES.md#transpose).

A language extension adds node classes that are not leaves, in a package
outside the core. A node class of a language extension subclasses
`AbstractExpression`, never a concrete class of the core (see [Leaves and
subclassing](#leaves-and-subclassing)). Each node class of a language
extension provides, as part of the class:

- the value of the node, as the section of a lazy node class of the core
  states the value (see [Evaluation](#evaluation)), and
- a lowering to an equal expression of node classes of the core, if such an
  expression exists.

A node class of a language extension may state argument spaces, as a function
in a dual space does, so that the node is an operator with arguments. The
argument spaces may include a dual space, which a `Form` of the core never
binds. The table of the extension for a pass that evaluates actions gives the
action on each argument.

#### Rule tables

A pass is built from an ordered `tuple` of rule tables. A rule table is a list
of rules, each a pattern and a rewrite. A pattern may match the class of a
node, the classes of the operands, a condition on the static type, or the
result of the pass on an operand. The pass applies to a node the first rule,
in the order of the tables and of the rules in each table, whose pattern
matches the node. Each rule receives the pass, so a rule applies the whole
pass, with every table, to the operands.

- A table of a language extension may match only some nodes. The table of
  the core comes last, and has a rule for every node class of the core.
- A pass applied to a node that no rule matches raises `TypeError`.
- The lowering table of a node class with a lowering matches the nodes of the
  class, lowers the node and applies the pass to the lowering.
- No global registry holds the tables. The code that builds a pass, such as a
  solver library, chooses the tables and the order, including the tables for
  the language extensions that the code loads. The pass is a function of the
  ordered tables, so the result does not depend on the order in which the
  extensions are loaded.

A language extension gives a rule table for each pass that the extension
supports. For a pass of the core or of another library, such as `units`, the
extension gives a table of rules or the lowering table. A rule of an extension
may match a pattern that a rule of the core also matches, and the earlier
table wins. Every rule keeps the value of the node, and a rule must agree with
the lowering followed by the pass, up to mathematical equality. An identity
that holds for a composition of the core, such as a cancellation in the
divergence of a Piola map, needs node classes of an extension for the
composition (see [A Piola extension](#a-piola-extension)).

#### Examples

- A language extension for local external operators adds a node class for an
  operation known only by name, such as `exp`, a material law or a neural
  network. The node class has no lowering. The table of the extension for a
  derivative gives the derivative of the operation, and the table for a
  pull-back gives the node of the pulled-back operands, because the operation
  is pointwise. A pass with no table of the extension raises `TypeError` for
  the node.
- A language extension for interpolation adds a node class for the
  interpolation of functions of a space `U` into a space `V`, with the
  argument spaces `(conjugate_space(dual_space(V)), U)`, the operator from
  `U` to `V` as a form. The action at position 0 on a function in
  `dual_space(V)` is the adjoint interpolation of that function, and the
  action at position 1 on a function in `U` is the interpolant. The table of
  the extension for a pass that evaluates actions gives both.
- A language extension for shape derivatives adds a node class for the
  derivative of an expression when the coordinates of the vertices of a mesh
  move in a direction. After the integrals and the pull-backs are expanded,
  the coordinates occur explicitly in the parameterisations of the cells, and
  the table of the extension lowers the node to a `DirectionalDerivative` with
  respect to the coordinates, which holds the domains fixed. A function keeps
  its values on the reference cell, so the function moves with the mesh.
- A named composition that must survive a pass, or print as written, is a node
  class of a language extension whose lowering is the composition. A
  `Transpose` node class lowers to `component_tensor(A[k, l], (l, k))`. The
  table of the extension for a pull-back gives the `Transpose` of the
  pull-back, because the transpose is pointwise, so `Transpose` survives the
  pull-back. A pass such as a derivative takes the lowering table of
  `Transpose`.
- For a pass that reduces an `Indexed` of a `ComponentTensor`, the table of
  the `Transpose` extension matches an `Indexed` of a `Transpose`, and
  rewrites `Transpose(A)[i, j]` as `A[j, i]`. The table comes before the table
  of the core, which has a rule for every `Indexed`.

### New questions

Requires: [Passes](#passes), [Language extensions](#language-extensions),
[`Coordinate`](#coordinate), [`Gradient`](#gradient).

A downstream library asks a new question of an expression with a pass. The pass
has a rule table for the node classes of the core and the leaf classes of the
library. A pass that meets a node class of a language extension also has a
table for that class before the table of the core, such as the lowering table
(see [Rule tables](#rule-tables)).

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
   `gradient(X[0], RD(2))` is in metres divided
   by metres, so the gradient is dimensionless.
   `X[0] + X[0] * X[1]` is an error of the units library.
3. Calls `units` on an expression whenever the library needs the check.

### New operations

Requires: [Node class or free function](#node-class-or-free-function), [`Tensor`
algebra with `EllipsisIndex`](EXAMPLES.md#tensor-algebra-with-ellipsisindex), [Language
extensions](#language-extensions).

A downstream library builds a new operation from existing node classes, as a
free function. As in the core, the free function may read the static type of
the operands, such as a value shape, to choose the composition (see
[Node class or free function](#node-class-or-free-function)). For example,
`cross(a, b)` reads the value shape `(3,)` of `a` and `b` and returns
`tensor((a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]))`,
so `cross(a, b)` requires the value shape `(3,)`. `trace(A)` reads no extent,
because `index_sum(A[k, k], (k,))`, with a derived index `k`, is the trace for
every extent, and `divergence(u, P)` is
`index_sum(gradient(u, P)[k, k], (k,))` for every extent. An operation that an `EllipsisIndex` expresses, such as
`outer(u, v)`, needs no value shape (see
[`EllipsisIndex`](EXAMPLES.md#tensor-algebra-with-ellipsisindex)).

### New subclasses

Requires: [Leaves and subclassing](#leaves-and-subclassing), [New
operations](#new-operations), [Language extensions](#language-extensions),
[`Coordinate`](#coordinate), [An affinely mapped
triangle](EXAMPLES.md#an-affinely-mapped-triangle).

Every abstract class is a contract. The contract states what a subclass
provides, so that every pass handles the subclass through the abstract class
alone. A downstream library may subclass an abstract class by
fulfilling the contract, without changing the core. A contract is part of the
class, never a global registry. For example, a downstream subclass of
`AbstractCoordinateDomain`, such as
[`UnitSimplex`](EXAMPLES.md#the-unit-simplex), provides the dimension, and every
pass handles the coordinate function of the subclass through `Coordinate`.

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

### A Piola extension

Requires: [Language extensions](#language-extensions), [Jacobian and
density](#jacobian-and-density), [`Gradient`](#gradient), [Representation on
the reference cell](#representation-on-the-reference-cell).

> A sketch of an extension, not a rule of the core.

On a cell that is not affine, the gradient of the contravariant Piola map
$J \hat u / \det J$ contains the gradients of $J$ and of $\det J$, which are
second derivatives of the parameterisation. In the divergence these terms
cancel, because the contravariant Piola map preserves the divergence:

$$
\operatorname{div} \left( \frac{J \hat u}{\det J} \right) = \frac{\widehat{\operatorname{div}} \, \hat u}{\det J}.
$$

The core sees only the composition, so a pass that compares graphs
structurally does not find the cancellation. A language extension makes the
cancellation structural with node classes for the maps and the operators:

| Node class                     | Lowering                                                               |
| ------------------------------ | ---------------------------------------------------------------------- |
| `ContravariantPiola(u_hat, D)` | $J \hat u / \det J$, with `jacobian(D)` and `determinant(jacobian(D))` |
| `CovariantPiola(u_hat, D)`     | $K^T \hat u$, with $K$ = `inverse_jacobian(D)`                         |
| `Divergence(u, P)`             | `index_sum(gradient(u, P)[k, k], (k,))`                                |
| `Curl(u, P)`                   | The curl in $R^3$, a fixed composition of `gradient(u, P)`             |

The library of a space with a Piola map gives the pull-back of a function by
a cell as a `ContravariantPiola` or a `CovariantPiola` of the function
$\hat u$ on the reference cell. The extension states the pull-backs, for a
domain $D$ = `DomainMap(P, psi, Q)` with $k = n$:

- If the pull-back of $u$ by $D$ is `ContravariantPiola(u_hat, D)`, the
  pull-back of `Divergence(u, Q)` is
  `Divergence(u_hat, P) / determinant(jacobian(D))`.
- If the pull-back of $u$ by $D$ is `CovariantPiola(u_hat, D)` and $n = 3$, the
  pull-back of `Curl(u, Q)` is `ContravariantPiola(Curl(u_hat, P), D)`.
- Otherwise the pull-back lowers `Divergence` or `Curl`.

The table of the extension for the pull-back holds these rules (see [Rule
tables](#rule-tables)). Each rule applies the pass to `pull_back(u, D)`, which
the library of the space evaluates, and reads the class of the result, a node
class of the extension.

## Appendix: Class hierarchy

Every object in the graph extends `Base`. Each tree gives the classes of a
chapter under the classes that they extend.

<!-- generated: class-hierarchy, run design/appendices.py to update -->
```
# Domains
Base
└── AbstractDomain
    └── AbstractCoordinateDomain
        └── RD

# Expressions and static type
Base
├── AbstractExpression
├── AbstractExtent
└── Index

# Algebra
AbstractExpression
├── Integer
├── Float
├── ImaginaryUnit
├── Zero
├── Sum
├── Product
├── Power
├── Division
├── Absolute
├── Conjugate
├── RealPart
├── ImaginaryPart
├── SquareRoot
├── Logarithm
├── Tensor
├── Indexed
├── ComponentTensor
├── IndexSum
└── Inverse
Base
├── EllipsisIndex
├── AbstractIndexValue
└── GatheredIndex

# Calculus on a domain
AbstractExpression
├── Coordinate
├── Gradient
├── DirectionalDerivative
└── Integral

# Maps between domains
Base
└── DomainMap
AbstractExpression
└── PullBack

# Topology
Base
├── AbstractCellType
└── AbstractTopology

# Reference cells
AbstractCoordinateDomain, AbstractCellType
└── AbstractReferenceCell

# Decomposed domains
AbstractDomain
└── DecomposedDomain
Base
└── SideMap
AbstractExpression
└── Normal

# Function spaces
Base
└── AbstractFunctionSpace
AbstractExpression
└── Function

# Variational forms
AbstractFunctionSpace
├── ConjugateSpace
└── DualSpace
AbstractExpression
├── Form
├── Action
└── GateauxDerivative
```
<!-- end generated: class-hierarchy -->

## Appendix: Free functions

<!-- generated: free-functions, run design/appendices.py to update -->
| Group                       | Free function                                                            | `AbstractExpression` overload | Takes                                                                                                                  | Returns                                                                                                                                                                         |
| --------------------------- | ------------------------------------------------------------------------ | ----------------------------- | ---------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Expressions and static type | [`index(label)`](#index)                                                 |                               | `label: str \| tuple`                                                                                                  | `Index`                                                                                                                                                                         |
| Algebra                     | [`power(first, second)`](#power)                                         | `**`                          | `first, second: AbstractExpression`                                                                                    | `Power`                                                                                                                                                                         |
|                             | [`negate(f)`](#negation-and-subtraction)                                 | `-f`                          | `f: AbstractExpression`                                                                                                | `multiply(Integer(-1), f)`                                                                                                                                                      |
|                             | [`subtract(f, g)`](#negation-and-subtraction)                            | `f - g`                       | `f, g: AbstractExpression`                                                                                             | `add(f, negate(g))`                                                                                                                                                             |
|                             | [`absolute(f)`](#absolute)                                               | `abs(f)`                      | `f: AbstractExpression`                                                                                                | `Absolute`                                                                                                                                                                      |
|                             | [`conjugate(f)`](#conjugate)                                             |                               | `f: AbstractExpression`                                                                                                | `Conjugate`                                                                                                                                                                     |
|                             | [`square_root(f)`](#elementary-functions-very-incomplete)                |                               | `f: AbstractExpression`                                                                                                | `SquareRoot`                                                                                                                                                                    |
|                             | [`logarithm(f)`](#elementary-functions-very-incomplete)                  |                               | `f: AbstractExpression`                                                                                                | `Logarithm`                                                                                                                                                                     |
|                             | [`tensor(entries)`](#tensor)                                             |                               | `entries: tuple[AbstractExpression, ...]`                                                                              | `Tensor`                                                                                                                                                                        |
|                             | [`indexed(expression, indices)`](#indexed)                               |                               | `expression: AbstractExpression`, `indices: tuple`                                                                     | `Indexed`                                                                                                                                                                       |
|                             | [`component(f, entries)`](#partial-indexing)                             | `[]`                          | `f: AbstractExpression`, `entries: tuple`                                                                              | `component_tensor(indexed(f, entries + K), K)` if the entries are fewer than the rank of $f$ and contain no `EllipsisIndex`, otherwise `indexed(f, entries)`                    |
|                             | [`component_tensor(expression, indices)`](#componenttensor-and-indexsum) |                               | `expression: AbstractExpression`, `indices: tuple[Index, ...]`                                                         | `ComponentTensor`                                                                                                                                                               |
|                             | [`index_sum(expression, indices)`](#componenttensor-and-indexsum)        |                               | `expression: AbstractExpression`, `indices: tuple[Index, ...]`                                                         | `IndexSum`                                                                                                                                                                      |
|                             | [`ellipsis_index(label)`](#ellipsisindex)                                |                               | `label: str \| tuple`                                                                                                  | `EllipsisIndex`                                                                                                                                                                 |
|                             | [`gathered_index(table, indices)`](#gatheredindex)                       |                               | `table: AbstractIndexValue`, `indices: tuple`                                                                          | `GatheredIndex`                                                                                                                                                                 |
|                             | [`add(first, second)`](#tensor-arithmetic)                               | `+`                           | `first, second: AbstractExpression`                                                                                    | `component_tensor(first[K] + second[K], K)` if both operands have rank 1 or more, otherwise `Sum(first, second)`                                                                |
|                             | [`multiply(first, second)`](#tensor-arithmetic)                          | `*`                           | `first, second: AbstractExpression`                                                                                    | `component_tensor(first * second[K], K)` if only `first` is a scalar, `component_tensor(first[K] * second, K)` if only `second` is a scalar, otherwise `Product(first, second)` |
|                             | [`divide(first, second)`](#tensor-arithmetic)                            | `/`                           | `first, second: AbstractExpression`                                                                                    | `component_tensor(first[K] / second, K)` if `first` has rank 1 or more and `second` is a scalar, otherwise `Division(first, second)`                                            |
|                             | [`real_part(f)`](#tensor-arithmetic)                                     |                               | `f: AbstractExpression`                                                                                                | `component_tensor(RealPart(f[K]), K)` if $f$ has rank 1 or more, otherwise `RealPart(f)`                                                                                        |
|                             | [`imaginary_part(f)`](#tensor-arithmetic)                                |                               | `f: AbstractExpression`                                                                                                | `component_tensor(ImaginaryPart(f[K]), K)` if $f$ has rank 1 or more, otherwise `ImaginaryPart(f)`                                                                              |
|                             | [`matrix_product(first, second)`](#tensor-arithmetic)                    | `@`                           | `first, second: AbstractExpression`                                                                                    | The contraction of the last axis of `first` with the first axis of `second`, below                                                                                              |
|                             | [`transpose(A)`](#tensor-arithmetic)                                     |                               | `A: AbstractExpression`                                                                                                | `component_tensor(A[k, l], (l, k))`, with `(k, l)` = `K`, for $A$ of rank 2                                                                                                     |
|                             | [`inverse(A)`](#inverse)                                                 |                               | `A: AbstractExpression`                                                                                                | `Inverse`                                                                                                                                                                       |
|                             | [`determinant(A)`](#determinant)                                         |                               | `A: AbstractExpression`                                                                                                | The Leibniz formula for $\det A$                                                                                                                                                |
| Calculus on a domain        | [`coordinate(P)`](#coordinate)                                           |                               | `P: AbstractDomain`                                                                                                    | `Coordinate`                                                                                                                                                                    |
|                             | [`gradient(f, P)`](#gradient)                                            |                               | `f: AbstractExpression`, `P: AbstractDomain`                                                                           | `Gradient`                                                                                                                                                                      |
|                             | [`directional_derivative(f, x, a, h)`](#directionalderivative)           |                               | `f: AbstractExpression`, `x: Function \| tuple`, `a, h: AbstractExpression \| tuple`                                   | `DirectionalDerivative(f, x, a, h)`                                                                                                                                             |
|                             | [`integral(f, P)`](#integral)                                            |                               | `f: AbstractExpression`, `P: AbstractDomain`                                                                           | `Integral`                                                                                                                                                                      |
| Maps between domains        | [`domain_map(P, psi, Q)`](#domainmap)                                    |                               | `P: AbstractDomain`, `psi: AbstractExpression`, `Q: AbstractDomain`                                                    | `DomainMap`                                                                                                                                                                     |
|                             | [`pull_back(f, D)`](#pullback)                                           |                               | `f: AbstractExpression`, `D: DomainMap \| SideMap`                                                                     | `PullBack`                                                                                                                                                                      |
|                             | [`jacobian(D)`](#jacobian-and-density)                                   |                               | `D: DomainMap`                                                                                                         | `gradient(D.parameterisation, D.parameter_domain)`                                                                                                                              |
|                             | [`inverse_jacobian(D)`](#jacobian-and-density)                           |                               | `D: DomainMap`                                                                                                         | The matrix $K$ below                                                                                                                                                            |
|                             | [`density(D)`](#jacobian-and-density)                                    |                               | `D: DomainMap`                                                                                                         | The density $\rho$ below                                                                                                                                                        |
| Topology                    | [`cell_index(T, j)`](#abstracttopology)                                  |                               | `T: AbstractTopology`, `j: int`                                                                                        | `Index`                                                                                                                                                                         |
| Decomposed domains          | [`decomposed_domain(pieces)`](#decomposeddomain)                         |                               | `pieces: tuple[DomainMap, ...]`                                                                                        | `DecomposedDomain`                                                                                                                                                              |
|                             | [`cell_domain(D, j)`](#cell-domains)                                     |                               | `D: DecomposedDomain`, `j: int`                                                                                        | `DomainMap` of the cells of type `j`                                                                                                                                            |
|                             | [`gathered_cell_domain(D, j, g)`](#gathered-cell-domains)                |                               | `D: DecomposedDomain`, `j: int`, `g: GatheredIndex`                                                                    | `DomainMap` of the cells of type `j` selected by `g`                                                                                                                            |
|                             | [`side_domain(D, S)`](#restriction)                                      |                               | `D: DecomposedDomain`, `S: AbstractTopology`                                                                           | The decomposed domain of the cells of $S$, with the geometry of the first side of each group on the topology of $D$                                                             |
|                             | [`side_map(D, E, s)`](#sidemap)                                          |                               | `D, E: DecomposedDomain`, `s: str`                                                                                     | `SideMap`                                                                                                                                                                       |
|                             | [`normal(sm)`](#normal)                                                  |                               | `sm: SideMap`                                                                                                          | `Normal`                                                                                                                                                                        |
| Function spaces             | [`function(V, label)`](#function)                                        |                               | `V: AbstractFunctionSpace`, `label: str \| tuple`                                                                      | `Function`                                                                                                                                                                      |
| Variational forms           | [`conjugate_space(V)`](#conjugatespace)                                  |                               | `V: AbstractFunctionSpace`                                                                                             | `V` for a real `V`, and `ConjugateSpace(V)` for a complex `V`                                                                                                                   |
|                             | [`dual_space(S)`](#dualspace)                                            |                               | `S: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]`                                                        | `DualSpace`                                                                                                                                                                     |
|                             | [`form(expression, arguments)`](#form)                                   |                               | `expression: AbstractExpression`, `arguments: tuple[Function \| tuple[Function, ...], ...]`                            | `Form`                                                                                                                                                                          |
|                             | [`action(f, position, g)`](#action)                                      |                               | `f: AbstractExpression`, `position: int`, `g: AbstractExpression \| tuple[AbstractExpression, ...]`                    | `Action`                                                                                                                                                                        |
|                             | [`gateaux_derivative(F, x, a, w)`](#gateaux-derivative)                  |                               | `F: AbstractExpression`, `x: Function \| tuple`, `a: AbstractExpression \| tuple`, `w: Function \| Conjugate \| tuple` | `GateauxDerivative(F, x, a, w)`                                                                                                                                                 |
<!-- end generated: free-functions -->

A free function and an operator overload accept only the
types in the "Takes" column.

A free function whose "Returns" column names only a node class `C` returns
`C(...)` of the arguments. Any other free function returns the composition in
the "Returns" column.

## Appendix: Static type table

The argument spaces of each node class of the core. An operand not named in
the argument spaces column must have the argument spaces `()`.

| Class                                    | Argument spaces of the node                                    |
| ---------------------------------------- | -------------------------------------------------------------- |
| `Integer`, `Float`, `ImaginaryUnit`      | `()`                                                           |
| `Zero`                                   | The stored argument spaces                                     |
| `Sum`                                    | Those of the operands, which must be equal                     |
| `Product`                                | Those of the operands, concatenated                            |
| `Power`, `Absolute`                      | `()`                                                           |
| `RealPart`, `ImaginaryPart`              | `()`                                                           |
| `Division`                               | Those of the first operand                                     |
| `Conjugate`                              | Those of the operand, each conjugated                          |
| `SquareRoot`, `Logarithm`, `Inverse`     | `()`                                                           |
| `Tensor`                                 | Those of the entries, which must be equal                      |
| `Indexed`, `ComponentTensor`, `IndexSum` | Those of the expression                                        |
| `Coordinate`                             | `()`                                                           |
| `Normal`                                 | `()`                                                           |
| `Function`                               | Those of the space                                             |
| `Form`                                   | Those of the expression, followed by the space of each entry, conjugated for a `conjugate(v)` entry |
| `Action`                                 | Those of the form without the entry at `position`              |
| `Gradient`                               | Those of the operand                                           |
| `DirectionalDerivative`                  | Those of the expression                                        |
| `GateauxDerivative`                      | Those of the expression, followed by the space of the entry $w$ |
| `Integral`                               | Those of the integrand                                         |
| `PullBack`                               | Those of the expression                                        |

The value shape and the free indices of each node class are in the section
of the class. The support and the scalar type of the node classes that differ
from the default rules of [Support](#support) and [Scalar type](#scalar-type):

| Class                    | Support of the node                | Scalar type of the node           |
| ------------------------ | ---------------------------------- | --------------------------------- |
| `ImaginaryUnit`          | None                               | Complex                           |
| `Zero`                   | None                               | Real                              |
| `Absolute`               | The default                        | Real                              |
| `RealPart`, `ImaginaryPart` | The default                     | Real                              |
| `Coordinate`, `Gradient` | The domain $P$                     | The default                       |
| `DirectionalDerivative`  | The default                        | The scalar type of $f$            |
| `GateauxDerivative`      | The default                        | The scalar type of $F$ for a real $x$, and complex otherwise |
| `Integral`               | None                               | The default                       |
| `PullBack`               | The parameter domain of the domain | The default                       |
| `Normal`                 | $E$                                | Real                              |
| `Function`               | The domain of the space, or none   | The scalar type of the space      |
| `Form`                   | The support of the expression      | The scalar type of the expression |
| `Action`                 | The support of the form            | The scalar type of the form       |

## Appendix: Summary

Counts of the core, from the sections.

<!-- generated: summary, run design/appendices.py to update -->
| Quantity                        | Count |
| ------------------------------- | ----- |
| Classes                         | 48    |
| Abstract classes                | 10    |
| Concrete classes                | 38    |
| Expression node classes         | 29    |
| Free functions                  | 46    |
| Free functions with an overload | 9     |
<!-- end generated: summary -->

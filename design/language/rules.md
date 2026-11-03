# Language rules

## Conventions

### Normative wording

| Wording | Meaning |
| ------- | ------- |
| **Must** | An implementation or library has to satisfy the rule. |
| **Should** | An implementation or library satisfies the rule unless the note states a reason not to. |
| **May** | A choice left open to the implementation or library. |
| **Present tense**, such as “the value shape is a `tuple`” | A rule of the language. |
| **Not yet covered** | A decided gap that the authors intend to specify later. |
| **Open** | A question the authors have not decided. |
| **Legacy UFL** dropdown | A comparison with legacy UFL, not a rule. |

### Mathematical admissibility

Construction checks structural validity. Values and rewrites assume the required
mathematical hypotheses; the core neither verifies them nor extends formulas
through singularities.

### Glossary

| Term | Meaning |
| ---- | ------- |
| [Leaf](#leaves-and-subclassing) | A primitive object rather than an operation; it may have graph dependencies. |
| [Terminal](#leaves-and-subclassing) | A leaf that is an expression, such as a literal, `Coordinate`, `Function` or `Normal`. |
| [Construction free function](#node-class-or-free-function) | A pure free function that may validate inputs and builds fixed expression compositions or canonical space descriptions. |
| [Pass](#passes) | A pure free function that walks a graph after construction to compute a result or rewrite the graph. |
| [Evaluation traversal](evaluation.md#evaluation-traversal) | A walk through expression operands, excluding dependencies stored in geometry, spaces and labels. |
| [Grammar pattern](evaluation.md#grammars) | A structural condition on a node and objects within a bounded number of edges, independent of pass results. |
| [Rule pattern](extensions.md#rule-tables) | A condition selecting a rewrite rule; it may inspect a pass result on an operand. |
| [Array rank](#abstractexpression), or **rank** | The number of axes in an expression's value shape. |
| [Arity](forms.md#argument-spaces) | The number of argument spaces of a form, counting a product-space argument once. |
| [Structurally equal](#the-base-object-base) | Objects with the same class and equal construction arguments, as tested by `==`. |
| [Mathematically equal](#the-base-object-base) | Expressions denoting the same value, whether or not their graphs are structurally equal. |

## Scope

### What UFLx is and isn't

UFLx is a typed symbolic language for writing variational forms.

The UFLx core language, henceforth the core, is independent of any particular
discretisation, solver or finite-element library. This draft requires vector calculus on
domains in $R^n$ and on manifolds embedded in $R^n$, in the coordinates of
$R^n$. One gradient serves every dimension, so the statement of each concept
does not depend on the dimension.

UFLx is *not* a symbolic language for reasoning about problems in functional
analysis or differential geometry. The core should therefore contain only the
concepts and relations needed to write variational forms as code for later
algorithmic transformation.

A UFLx program is ordinary Python code that builds a graph of symbolic objects,
usually a bilinear form, linear form or functional.

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
- Use dimension-independent definitions where the mathematics permits.
  Operations such as `Curl` state their dimension-specific cases explicitly.
- The core contains no algorithms solely for a form compiler's convenience.
  Each compiler states the forms it accepts. A pass outside the core must
  lower to those forms, using only core node classes and
  [language extensions](extensions.md#language-extensions) accepted by that compiler.
  For example, a basic form compiler might only accept integrals over a
  reference cell, such as a unit simplex, as in each cell term of a legacy UFL
  `dx` integral.

```{dropdown} Legacy UFL
TSFC lowers a UFL form to GEM, which repeats much of the algebra of UFL, partly because UFL simplifies during construction. The UFLx core does not simplify during construction, so a form compiler may be able to use a lowered UFLx form directly, adding only leaves for tables and kernel inputs. This has not been confirmed with the TSFC or FFCx developers.
```

## Rules of the language

### The base object `Base`

Every object in the directed acyclic graph (DAG) must extend `Base`. A graph
edge means "is defined from"; geometric edges are explicitly called edges of
a triangle or cell. Everything defining the variational problem must be
recoverable from the graph.

| `Base`     |                                                                    |
| ---------- | ------------------------------------------------------------------ |
| Properties | `construction_arguments`, the construction arguments, as in rule 4 |

Objects must satisfy the following rules:

1. A construction argument must be one of a small curated set of types: a base
   object, or exactly an `int`, a finite `float`, a `str`, or a `tuple` whose
   items are of these types. A subclass of these Python types is not accepted,
   so a `bool` is not an `int` and a named tuple is not a `tuple`. Construction
   arguments are deeply immutable. Their hashes must be deterministic across
   processes on the same machine, and they are compared by type and value.
   `1`, `1.0` and `True` are distinct, `0.0` and `-0.0` are distinct, and every
   value is equal to itself. Construction raises `TypeError` for an argument of
   any other type, and `ValueError` for a `float` that is not finite.
2. Construction must depend only on the construction arguments, not on global
   state.
3. The state of an object is immutable after construction. Every result that
   the object produces after construction depends only on the class and the
   construction arguments.
4. An object must be reconstructible from its class and construction arguments,
   with the result equal to the original. The base object records these arguments
   and returns them in order as the `tuple` `construction_arguments`. A subclass
   constructor validates them and may store data derived from them.
5. The edges of an object must be exactly the base objects among its
   construction arguments, including the base objects nested in `tuple`s, in
   order.
6. Any order that the core or a pass relies on must be deterministic across
   processes. A pass, and any other walk of the graph such as hashing,
   serialisation or printing, follows the edges it visits in order (rule 5).
   A pass states which edges it visits. Where unrelated objects must be sorted,
   the sort uses a total order defined by the class and the construction
   arguments alone.
7. Equality (`==`) must be structural and return a `bool`. Whether `is` is
   structural is implementation dependent. An implementation in which `is` is
   structural must intern objects, so that equal objects are the same object.
   Structural equality does not decide mathematical equality: $a + b$ and
   $b + a$ are not structurally equal. Mathematical equalities are applied
   in later graph algorithms.
8. Constructors validate their nodes. Construction free functions may validate
   their documented inputs, but never convert them. Structural requirements
   must be checked by the free function or the constructors it calls. Unless
   stated otherwise, an unacceptable type or class raises `TypeError`, and an
   unacceptable value of an accepted type or class raises `ValueError`.
   Class tables list only additional details or exceptions to these defaults.
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
construction. Every object that exists is therefore structurally valid.

### Leaves and subclassing

Only leaves are subclassed. A leaf is a primitive object rather than an
operation, such as the coordinate domain `RD(n)` or a function whose values an
external library supplies. Leaves may have graph dependencies: a `Function`
depends on its space, and a `Coordinate` on its domain. A terminal is a leaf
that is an expression, such as a literal, `Coordinate`, `Function` or `Normal`
of the core, or a terminal of a downstream library.

An operation is represented by a single core class or a construction free
function returning core objects, never a subclass. For example, a construction
free function builds the map of an affinely mapped triangle as a
`DomainMap`, rather than a subclass of `DomainMap`. A language extension is the
one exception (see [Language extensions](extensions.md#language-extensions)).

Concrete core classes are final; libraries subclass only abstract classes.
A parent class never constructs a subclass, though its properties may return
existing subclass instances. A subclass may narrow a property's return type.

### Construction policy

| Kind | Role and permitted inspection | Validation and exceptions |
| ---- | ----------------------------- | ------------------------- |
| Property | Returns stored or derived data; creates no objects | Derived objects may be built and stored during construction |
| Constructor | Creates a node and computes derived data; does not simplify expressions | Validates construction arguments |
| Construction free function | Builds a fixed composition from expression static types and documented interfaces of explicit non-expression arguments; does not inspect expression structure to choose a composition | Validates inputs itself or through constructors; may inspect structure for validation only. Space factories apply canonical forms; `indexed`, `component`, `component_tensor` and `index_sum` apply their empty-index identities |
| Pass | After construction, inspects graph structure and returns information or a rewritten graph | Enforces its own preconditions; constructors validate nodes it builds. Space substitution preserves canonical forms |

Free functions and passes are pure. A class method immediately calls a free
function. Constructor arguments are positional only; named arguments raise
`TypeError`. Free functions may take named arguments, such as
`domain_map(P=T, psi=psi, Q=RD(2))`.

A stored property returns a construction argument; a derived property depends
on the class and those arguments. Concrete classes inherit the derived
properties of their abstract classes.

### No simplification during construction

The construction policy preserves the user's written operations. For example,
in `(a + b) + c`, the user may need `(a + b)` evaluated before `c` is added
in floating-point arithmetic. Building a fixed composition defines an
operation; it does not simplify its operands. Later passes, chosen by the
user, may rewrite it.

### Passes

A pass is a pure free function that walks a graph after construction and
returns a value, such as a unit, or a new graph. To build a new graph, it
constructs replacement nodes and reconstructs dependent objects up to
the root along the edges it visits (see [The base object](#the-base-object-base)).
The dependency graph includes every construction argument; the
[evaluation traversal](evaluation.md#evaluation-traversal) follows only expression operands.
A substitution follows dependencies, so replacing an object also
reconstructs an [`Index`](#index) whose label contains the object, other than
a renamed name (see [Derived indices](expressions.md#derived-indices)). Substitution
commutes with a fixed expression composition when its static-type choice is
unchanged: replacing `A` by `B` in `transpose(A)` gives `transpose(B)`, with
the derived indices of `B`. A pass
may cache results for each object, because objects are immutable.

Substitution is simultaneous and respects lexical scope, renaming bound names to
avoid capture (see [Substitution and binding](evaluation.md#substitution-and-binding)). The
core pass `replace(e, mapping)` performs this substitution.

A pass does one thing, such as differentiating or lowering. Passes compose,
and the user of the core chooses which passes run and in what order.

A pass fails with `TypeError` for a node class it has no rule for, or with an
error of the pass, such as inconsistent units.

### Operator overloads

Only `AbstractExpression` defines arithmetic and indexing overloads, listed in
the [table of free functions](reference.md#appendix-free-functions).

An overload passes operands unchanged to its free function, except for the
indexing syntax below. There are no reflected overloads, such as `__radd__`.

An overload applies one free function to the operands, so the graph follows
the grouping and the precedence of Python, and construction never regroups
(see [No simplification during
construction](#no-simplification-during-construction)). For example,
`a + b + c` is `Sum(Sum(a, b), c)` and `a + (b + c)` is `Sum(a, Sum(b, c))`,
and the two are not equal. `a ** b ** c` is `power(a, power(b, c))`, and
`-a ** b` is `negate(power(a, b))`.

Python passes `f[i]` to `__getitem__` as `i`, but `f[i, j]` as `(i, j)`.
`__getitem__` wraps a single entry in a `tuple` and calls `component` (see
[Partial indexing](expressions.md#partial-indexing)).

`AbstractExpression` defines `__iter__` to raise `TypeError`, because
expressions are not iterable and Python otherwise iterates through
`__getitem__`.

### Numbers in expressions

Every number in an expression is written with literals, such as `Integer(2)`,
`Float(0.5)`, `Integer(1) / Integer(3)` or `Complex(Integer(0), Integer(1))`. A Python `int`, `float`,
`complex`, `Fraction` or NumPy scalar is never converted and is rejected with
`TypeError`.

With `X = coordinate(RD(2))`, `X[0]**Integer(2)` is accepted, but `X[0]**2`
raises `TypeError`. The same rule applies inside containers:
`tensor((Integer(1), Integer(2)))` is accepted, but `tensor((1, 2))` is not.

## Domains

Function spaces and operators use a domain's dimension, boundary and cell
decomposition. The only metric is the ambient Euclidean metric, which defines
the tangential [gradient](expressions.md#gradient), integration [density](geometry.md#jacobian-and-density)
and unit [normal](geometry.md#normal).

### `AbstractDomain`

`AbstractDomain` represents a domain, the set on which an expression is defined.
A domain is a coordinate domain supplied by a library (see
[`AbstractCoordinateDomain`](#abstractcoordinatedomain)), or a mesh with a
topology and coordinate expression (see [`Mesh`](geometry.md#mesh)). Every node class
taking a domain accepts either kind; its section defines its value on a mesh
through its cells.

| `AbstractDomain` |                                                                             |
| ---------------- | --------------------------------------------------------------------------- |
| Extends          | `Base` |
| Provides         | `dimension`, an `int` $k \ge 0$, and `ambient_domain`, an `AbstractDomain` |
| Defined by       | `AbstractCoordinateDomain` and `Mesh`                           |

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
(see [`Coordinate`](expressions.md#coordinate)), and the integral is the Lebesgue integral in
the coordinates (see [`Integral`](expressions.md#integral)). Every other domain gets these by
pull-back.

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
| Raises     | `TypeError` if $n$ is not an `int`, and `ValueError` if $n < 0$ |

Worked examples for this chapter: [the unit simplex](../examples/library.md#the-unit-simplex).

## Expressions and static type

### `AbstractExpression`

| `AbstractExpression` |                                                                                      |
| -------------------- | ------------------------------------------------------------------------------------ |
| Extends              | `Base`                                                                               |
| Provides             | `value_shape`, `free_indices`, `scalar_type`, `domain`, `argument_spaces` |

These five derived properties form the static type (see
[Argument spaces](forms.md#argument-spaces)). The value shape is a `tuple` of
[extents](#extents), one per axis; `()` denotes a scalar. The array rank,
written simply as rank, is the number of axes.

### Extents

An extent is the length of an axis of a value shape, or the extent of a free
index. An extent is an `int` $n \ge 0$ or an `AbstractExtent`.

`AbstractExtent` represents an extent whose value is not known during
construction. Its supplied value must be a nonnegative integer.
An `AbstractExtent` is not an expression. A concrete
`AbstractExtent` is a subclass outside the core, such as `RuntimeExtent` in the
[cut cell example](../examples/integration.md#integrating-over-a-cut-cell).

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
read from a table is a [`GatheredIndex`](expressions.md#gatheredindex).
```

### `Index`

`Index` is a free index, not an expression. An `Index` names an axis of a
tensor, as `i` and `j` in `A[i, j]` (see [Index notation](expressions.md#index-notation)).

| `Index(label)` |                                                                                                 |
| -------------- | ----------------------------------------------------------------------------------------------- |
| Extends        | `Base` |
| Free function  | `index(label)`                                                                                  |
| Takes          | `label: str \| tuple` |
| Properties     | `label`, a `str` for a user index, or a `tuple` that contains graph objects for a derived index |
| Raises         | `TypeError` if the label is neither a `str` nor a `tuple`                                       |

An `Index` takes its extent from the indexed axis; it has no extent of its own.

A user index has a `str` label, which the user chooses. Two free occurrences
of one `Index` in one scope are the same index, so an index that the user
means to be new needs a new name. The scope of an index is given in
[`ComponentTensor` and `IndexSum`](expressions.md#componenttensor-and-indexsum).

A derived index has a `tuple` label, which a free function or a pass builds
from the operands, so that the index is fresh without a global counter (see
[Derived indices](expressions.md#derived-indices)).

### Free indices

Every `AbstractExpression` derives the free indices, a `tuple` of pairs. Each
pair is an `Index` and an extent. The order is the order of first occurrence, never the
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

The scalar type is real or complex. Unless the [static type
table](reference.md#appendix-static-type-table) states otherwise, an expression is complex
if and only if at least one edge is a complex expression.

The scalar type is the arithmetic in which the expression is evaluated, not a
statement about the value. A real `SquareRoot`, `Logarithm` or `Power` whose
value is not a real number, such as `square_root(Integer(-1))`, is undefined.

### Domain

An expression's `domain` is where it is defined: an
[`AbstractDomain`](#abstractdomain), or an [`AbstractTopology`](geometry.md#abstracttopology)
whose cell types are reference cells. It is `None` for an expression independent
of a domain, such as a literal.

Unless the [static type table](reference.md#appendix-static-type-table) states otherwise,
a node inherits the domain of its expression operands with a domain. All such
operands must have the same domain, or construction raises `ValueError`. A node
with no such operand has domain `None`.

The property always describes the resulting expression. An `Integral` stores
its region as `integration_domain` but has domain `None`; a `PullBack` stores
its `map` and has that map's `parameter_domain`. A `PushForward` stores its
`mesh`, which is also its domain. For `Coordinate`, `Gradient`, `Divergence`
and `Curl`, the stored domain is also the result's domain.

### Expressions as functions

An expression denotes an array-valued function on the domain of the
expression. An expression on a topology denotes a function on the reference
cell of each cell of the topology (see [`CellMap`](geometry.md#cellmap)). An expression
with no domain denotes a single array. An expression with argument spaces other than `()` denotes a multilinear map from
the argument spaces to such functions (see [Variational
forms](forms.md#variational-forms)).

### Node class or free function

An operation uses a construction free function when it has a fixed
composition under the [construction policy](#construction-policy). A node
class is needed when choosing the result requires expression structure, when
the compiler must choose a numerical algorithm, or when a pass must retain
the operation as written. In the last case, the node's section states its
composition; `Integral` over a mesh is an example.

For example, the transpose is always `component_tensor(A[k, l], (l, k))`,
so it is a free function. `Inverse` and `Determinant` keep the algorithm
choice for the compiler.

A construction free function is part of the core only if a definition or a
rule of the core uses it, such as `density` in [`Integral` over a
mesh](geometry.md#integral-over-a-mesh) and `cell_map` in
[Restriction](geometry.md#restriction), or if a Python operator calls the free function
as an overload, such as `matrix_product` for `@` (see [Operator
overloads](#operator-overloads)). `transpose` is in the core because the rule
for `Inverse` of [Applying derivatives](evaluation.md#applying-derivatives) uses it. Every
other construction free function belongs to a library, such as the operations
of [Tensor algebra](../examples/tensors.md#tensor-algebra).

A node class is part of the core only if a definition or a rule of the core
needs the node class, or if the node class is an operation of the scalar
field, such as conjugation. A rule for the node class itself does not count.
Otherwise the node class belongs to a
[language extension](extensions.md#language-extensions). For example, the core writes a
rational number as the `Division` of two `Integer`s, so `Division` is part of
the core.

The [node-class rationale](reference.md#appendix-node-class-rationale) explains the
individual choices.

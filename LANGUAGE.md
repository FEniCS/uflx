# Core UFLx language definition

This note describes the smallest core language behind UFLx. It intentionally
avoids implementation details, compiler algorithms, and broad redesign ideas.
It assumes knowledge of differential geometry, and functional analysis. Some
advanced differential geometry concepts are explained in more detail.

This document is aimed at users and developers with extensive knowledge of
domain specific languages for variational forms (e.g. UFL, ngsolve).

## Scope

### What UFLx is and isn't

In one sentence:

> UFLx is a typed symbolic language for writing variational forms.

The UFLx core language (henceforth the core) should not use finite element or
numerical analysis terminology - differential geometry, and functional
analysis, and their intersection, without finite element (numerical analysis)
concepts, should be all that is required. This draft of the scheme requires
knowledge of calculus on intrinsic manifolds, i.e. a coordinate-free
representation. This provides a very compact and dimension-free representation
of many concepts already in UFL.

Finite element and numerical analysis concepts, when necessary, will be added
on top of the core, or provided by external libraries (Basix being the
'in-house' example). Basix and DOLFINx and their assumptions should not be
visible in the core.

UFLx is *not* a symbolic language for reasoning about problems in functional
analysis or differential geometry. Hence, the core should only contain the
minimal concepts and relations necessary for writing a wide range of
variational forms as code for subsequent algorithmic transformation.

A UFLx program is ordinary Python code that constructs symbolic objects and
connects them through a graph. The primary results are usually variational
forms, such as bilinear forms, linear forms and functionals. Python code will
be used to motivate examples but a full concrete implementation is out of
scope.

UFLx does not handle concrete data or numerical algorithms; e.g. actual meshes,
boundary markers, numerical values of coefficients, or assembly procedures.
Those belong to the surrounding problem-solving environment, such as DOLFINx,
Firedrake, FIAT and Basix. UFLx instead specifies a set of leaf nodes, and the
external library provides the concrete data behind the leaf nodes.

### Basic principles

- All core object definitions should be closely related to their strict
  mathematical definition.
- Extra definitions, properties or terminology outside of the established
  mathematical fields noted must be avoided.
- Strictness and expressivity, not brevity, is key in the core - brevity can be
  added later with higher-level APIs. The cost of writing code against a
  precise lower-level API has massively decreased with LLMs, and so the end
  library could potentially more concise than it was in the past.
- The mathematical fields provide a rich foundation for writing variational
  problems in a uniform way across geometric dimensions, i.e. 'special cases'
  for particular dimensions do not exist, and so should not exist in the core.
- The core does not contain algorithms for the convenience of a form
  compiler. A form compiler states the forms that the form compiler accepts.
  A pass outside the core lowers a form to those forms, and the pass produces
  only node classes of the core.

## Rules of the language

### Basic implementation principles

- Explicit is better than implicit.

- Boilerplate and error-prone patterns should be minimised.

- UFLx should be deterministic on the same machine using the same byte-code.
  Serialisation/pickling should be considered early on.

- Performance is *important* and so any data structures and algorithms in the
  core should be carefully chosen and tested for scalability.

- The minimal language should only need a relatively small number of object
  kinds. When designing the language, consider whether objects can be composed,
  rather than extended. Consider what could be layered on top of the core.

### The base object `Base`

The base object in UFLx must be extended by all objects that are included in
the graph. The graph is a directed acyclic graph (DAG). An edge means "is
defined from".

Everything needed to define the variational problem must be recoverable from
the graph. No information is held outside of the graph.

Objects must satisfy the following rules:

1. A construction argument must be one of a small curated set of types: a base
   object, exactly `int`, a finite `float`, `str`, or a `tuple` whose items are
   of these types. A `bool` is not an `int` for this rule. Construction
   arguments are therefore deeply immutable, have a hash that is deterministic
   across processes on the same machine, and are compared by type and value.
   `1`, `1.0` and `True` are distinct, `0.0` and `-0.0` are distinct, and every
   value is equal to itself.
2. Construction must depend only on the construction arguments, not on global
   state.
3. The state of an object is immutable after construction. Every result that
   the object produces after construction depends only on the class and the
   construction arguments.
4. An object must be reconstructible from its class and construction arguments,
   and the reconstructed object must be equal to the original.
5. The edges of an object must be exactly the base objects among its
   construction arguments, including the base objects nested in `tuple`s, in
   order.
6. Any order that the core or a pass relies on must be deterministic across
   processes. A pass, and any other walk of the graph such as hashing,
   serialisation or printing, follows the edges of each object in order (rule
   5). Where unrelated objects must be sorted, the sort uses a total order
   defined by the class and the construction arguments alone.
7. Equality (`==` or `is`) must be structural and return a `bool`. Equality
   does not decide mathematical equivalence, so $a + b$ and $b + a$ are not
   equal. Mathematical equivalencies are applied in later graph algorithms.
8. All validation of an object happens in the constructor of the class of the
   object. This rule is a hard rule. A free function never validates and never
   converts arguments. A free function only reads derived properties of
   arguments and constructs objects. An invalid argument is rejected by the
   constructor of the first object that the argument makes invalid.

Hence two objects are equal if and only if the two objects have the same class
and equal construction arguments.

Because objects are immutable, an object need only be validated once, at
construction. Every object that exists is therefore valid.

An object cannot be changed. To obtain a modified object, construct a new
object and reconstruct every object that depends on the original object, up to
the root of the graph.

The above design implies that construction arguments should be small and cheap
to hash and compare, and the graph implementation must be performant.

### Leaves and subclassing

Requires: [The base object `Base`](#the-base-object-base).

See also: [RD](#rd), [ParameterisedDomain](#parameteriseddomain), [Language
extensions](#language-extensions).

Only leaves are subclassed. A leaf is an object that the core cannot build
from other objects, such as the coordinate domain `RD(n)` or a function whose
values an external library supplies. An object built from other objects is a
single class of the core or a free function that returns objects of the core,
never a subclass. For example, every domain is a `ParameterisedDomain`, and an
affinely mapped triangle is a free function that returns a
`ParameterisedDomain`. A language extension is the one exception (see
[Language extensions](#language-extensions)).

A parent class never constructs a subclass. A property of a parent class may
return an object whose class is a subclass, because the property only
returns an object that already exists. A subclass may narrow the return
type of a property.

### Properties and free functions

Requires: [The base object `Base`](#the-base-object-base).

See also: [`AbstractExpression`](#abstractexpression),
[ParameterisedDomain](#parameteriseddomain), [Coordinate](#coordinate),
[Transpose and matrix product](#transpose-and-matrix-product).

A property reads and a free function builds. Data that an object stores, or
derives from the class and the construction arguments of the object, is a
property of the object, such as `value_shape` of an expression or
`parameter_domain` of a domain. A property never constructs an object. Any
operation that constructs objects is a pure free function, such as
`coordinate(P)` or `transpose(A)`. A class method, where
needed, immediately calls a free function.

Free functions are not objects in the graph. A free function builds or returns
objects in the graph. Operator overloads call the free functions.

A stored property returns a construction argument. A derived property returns
data derived from the class and the construction arguments. A concrete class
has the derived properties of the abstract classes of the concrete class.

### No simplification during construction

Requires: [The base object `Base`](#the-base-object-base).

Simplifications and rewriting should be completely avoided during user graph
construction, maintaining the user's intent. Only later selectable algorithmic
passes should re-write. A good example is `(a + b) + c` in floating point
arithmetic where `a`, `b` and `c` are unknown terminals (i.e. to be provided
later). The user may have very good reason to have `(a + b)` performed before
adding the answer to `c`.

### Passes

Requires: [The base object `Base`](#the-base-object-base).

A pass is a pure free function that walks a graph after construction of the
graph. The result of a pass depends only on the arguments of the pass. A pass
returns a value, such as a unit, or a new graph. A pass builds a new graph by
constructing an object for each node that the pass replaces, and by
reconstructing every object that depends on a replaced node, up to the root
(see [The base object](#the-base-object-base)). A pass may cache results for
each object, because objects are immutable.

A pass does one thing, such as differentiating or lowering. Passes compose,
and the user of the core chooses which passes run and in what order.

### Kinds of code that produce objects

Requires: [The base object `Base`](#the-base-object-base), [Properties and free
functions](#properties-and-free-functions), [Passes](#passes).

See also: [Transpose and matrix product](#transpose-and-matrix-product),
[Indexed](#indexed), [New questions](#new-questions).

Three kinds of code produce objects:

| Kind                                               | When the code runs              | Inspects arguments to choose what to build                                | Validates                     |
| -------------------------------------------------- | ------------------------------- | ------------------------------------------------------------------------- | ----------------------------- |
| Constructor                                        | Whenever an object is made      | No, a constructor inspects arguments only to validate                     | Yes, and only constructors do |
| Construction free function, such as `transpose`    | While the user builds the graph | Only the static type, to build one fixed composition for each static type | No                            |
| Pass, such as differentiation, lowering or `units` | After construction              | Yes, and a pass may rewrite the graph                                     | No                            |

A construction free function fails when a constructor rejects the
arguments, for example `transpose` of a vector raises `IndexError` from the
constructor of `Indexed`. A pass fails with `TypeError` for a node class that
the pass has no rule for, or with an error of the pass, such as inconsistent
units.

### Enforcement

Requires: [The base object `Base`](#the-base-object-base).

Enforcement must be implemented to help library developers - look at tools such
as `dataclass`, deterministic hashing, hash-consing, `WeakValueDictionary`,
`ruff` and `semgrep`.

https://github.com/semgrep/semgrep

This rules should be put in place at the start of the project.

It is not clear at the moment how best to do this - it may be that custom
static analysis e.g. (semgrep) is more robust than hacking Python to be 'well
behaved' in some places, for example checking the purity of functions is not
possible in Python.

## Coordinate domains

We begin with differential geometry because the function spaces and operators
of UFLx are defined on a domain. The function spaces and operators require the
topological structure of the domain, which is the dimension, boundary,
orientation and decomposition into cells of the domain. The core needs no
metric.

### AbstractCoordinateDomain

Requires: [The base object `Base`](#the-base-object-base).

See also: [Support](#support), [Coordinate](#coordinate),
[ParameterisedDomain](#parameteriseddomain).

An abstract base class for coordinate domains `AbstractCoordinateDomain`.

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$ whose
points are named by the coordinates of the points themselves. A coordinate
domain is a parameter space. A coordinate domain is the parameter domain of a
domain, the ambient coordinate domain of a domain, and the support of an
expression. A coordinate domain is not a
`ParameterisedDomain`. The domain whose points are the points of a coordinate domain
$P$ is `parameterised_domain(P, coordinate(P), P)` (see
[ParameterisedDomain](#parameteriseddomain)).

A new coordinate domain, such as a reference cell, is a subclass of
`AbstractCoordinateDomain`, because the core cannot build the set from smaller
objects.

Every concrete `AbstractCoordinateDomain` $P$ stores:

- The dimension $n$, provided by the property `dimension`.

Examples:

- The unit simplex in $R^n$ and $R^n$ itself are coordinate domains.
- A face of codimension 1 of the unit simplex in $R^n$ is *not* a coordinate
  domain, because the face has empty interior in $R^n$. The face is the image
  of the unit simplex in $R^{n-1}$ under a map, so the face is a parameterised
  domain (defined below).

### RD

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain).

See also: [Example: a paraboloid](#example-a-paraboloid).

A concrete subclass of `AbstractCoordinateDomain` for the whole of $R^n$.

`RD(n)` stores the dimension $n$, an `int` with $n \ge 0$.

Example:

- `RD(2)` is the parameter domain of the paraboloid example below.

## Expressions and static type

### `AbstractExpression`

Requires: [The base object `Base`](#the-base-object-base).

Every `AbstractExpression` derives the value shape, the free indices, the
scalar type, the support and the degree of the expression.

### Index

Requires: [The base object `Base`](#the-base-object-base), [Properties and free
functions](#properties-and-free-functions).

See also: [Indexed](#indexed), [Transpose and matrix
product](#transpose-and-matrix-product).

> LLM generated content.

Index notation writes tensor algebra as fixed compositions, so that operators
such as `transpose` and `matrix_product` are free functions rather than node
classes.

A concrete subclass of the base object for a free index. An `Index` is not an
expression.

`Index(label)` stores the label. The label is either a `str`, for an index
that the user names, or a `tuple` that contains graph objects, for an index
that a free function derives from the free function's operands. An `Index`
has no extent of its own. An `Index` takes the extent of the axis that the
`Index` indexes.

The user names user indices. The user is responsible for not reusing a name
that is already free in an operand.

A free function that needs fresh indices derives the indices from the
operands. The label of a derived index is the identifier of the free function,
followed by the arguments of the free function and a number. The identifier is
the public name of the free function, qualified by the package name, such as
`"uflx.transpose"` for the core or `"library.trace"` for a library. For example,
`transpose(A)` uses `Index(("uflx.transpose", A, 0))` and
`Index(("uflx.transpose", A, 1))`. An index whose label contains `A` has an
edge to `A`. Every node that uses an `Index` has an edge to the `Index`. So an index
whose label contains `A` cannot occur inside `A`, because `A` would then
contain `A`, and the graph is acyclic. Hence
a derived index is fresh with respect to the operands, and the free function
does not inspect the operands to find a fresh index. A global counter is not
allowed, because a counter makes construction depend on global state.

### EllipsisIndex

Requires: [Index](#index).

See also: [Indexed](#indexed), [Tensor algebra with
EllipsisIndex](#tensor-algebra-with-ellipsisindex).

> LLM generated content.

> Legacy UFL: This numpy-like concept does not appear to be in legacy UFL.

A concrete subclass of the base object for an index that covers any number of
axes. A `EllipsisIndex` is not an expression. `EllipsisIndex(label)` stores the label,
with the same rules as the label of an `Index`.

A `EllipsisIndex` covers every axis of an `Indexed` that the other entries do not
cover, so a `EllipsisIndex` may cover no axis. The extent of a `EllipsisIndex` is the
`tuple` of the extents of the covered axes. In every other rule of index
notation, a `EllipsisIndex` behaves as an `Index`.

### Free indices

Requires: [`AbstractExpression`](#abstractexpression), [Index](#index),
[EllipsisIndex](#ellipsisindex).

See also: [Sum](#sum), [Tensor](#tensor), [Indexed](#indexed), [ComponentTensor
and IndexSum](#componenttensor-and-indexsum).

Every `AbstractExpression` derives the free indices, a `tuple` of pairs of an
`Index` and an `int` extent, or of a `EllipsisIndex` and a `tuple` extent. The
order is the order of first occurrence, never the order of a `set`.

Unless a node class states otherwise, the free indices of a node are the free
indices of the edges of the node that are expressions, in edge order, with each
`Index` taken once at the first occurrence. An `Index` that occurs in two edges
must have the same extent in both edges, or construction raises `ValueError`.
Terminals and literals have no free indices.

The following node classes state otherwise:

- `Sum` and `Tensor`: all operands, or all entries, must have the same free
  indices, with the same extents and in the same order, or construction raises
  `ValueError`.
- `Indexed`: see [Indexed](#indexed).
- `ComponentTensor` and `IndexSum`: see
  [ComponentTensor and IndexSum](#componenttensor-and-indexsum).

Summation is never implicit. `a[i] * b[i]` has the free index `i` and denotes
the component-wise product. The sum over `i` is `index_sum(a[i] * b[i], (i,))`.

### Scalar type

Requires: [`AbstractExpression`](#abstractexpression).

See also: [Literals](#literals), [Absolute](#absolute).

The scalar type is real or complex. Unless a node class states otherwise, an
expression is complex if and only if at least one edge is a complex
expression. `ImaginaryUnit` is complex.

### Support

Requires: [`AbstractExpression`](#abstractexpression),
[AbstractCoordinateDomain](#abstractcoordinatedomain).

The support is an `AbstractCoordinateDomain`. The support is absent for an expression
with no support, such as a literal. Unless a node class states otherwise, the support
of a node is the support of the edges of the node that are expressions with a support.
All such edges must have the same support, or construction raises `ValueError`. A
node with no such edge has no support.

### Degree

Requires: [`AbstractExpression`](#abstractexpression).

See also: [Appendix: Degree table](#appendix-degree-table), [New
subclasses](#new-subclasses), [ExteriorDerivative](#exteriorderivative).

The degree is a `tuple` with the form degree of each slot of the expression,
and every entry is at least 1. An expression with the degree `()` has no
slots and is a 0-form, such as a scalar or an array. The axes of the value
shape are not slots. The degree has no default. Every concrete class states
the degree of the class, and a leaf class of a downstream library states the
degree as part of the contract of the class (see
[New subclasses](#new-subclasses)).

The degree of every node class of the core is in the
[degree table](#appendix-degree-table). If the
degree of an edge breaks the rule of the node, construction raises
`ValueError`. A rule that requires the degree `()` of an edge is a deliberate
statement, such as for `SquareRoot`, because the square root of a differential
form is not defined.

Examples of the degree, on $R^3$ with the coordinate function $X$, the
coordinates $x, y, z$ and a scalar function $f$:

| Object                                                                     | Value shape | Degree   |
| -------------------------------------------------------------------------- | ----------- | -------- |
| The scalar $f$, or the literal `Integer(1)`                                | `()`        | `()`     |
| The coordinate function `X`                                                | `(3,)`      | `()`     |
| The 1-form $\mathrm{d}f$                                                   | `()`        | `(1,)`   |
| The 1-form $\mathrm{d}x$                                                   | `()`        | `(1,)`   |
| The 2-form $\mathrm{d}x \wedge \mathrm{d}y$                                | `()`        | `(2,)`   |
| The volume form $\mathrm{d}x \wedge \mathrm{d}y \wedge \mathrm{d}z$        | `()`        | `(3,)`   |
| The vector-valued 1-form $\mathrm{d}X$, with one 1-form for each component | `(3,)`      | `(1,)`   |
| The double form $\mathrm{d}x \otimes \mathrm{d}y$, a covariant 2-tensor    | `()`        | `(1, 1)` |

### Expressions as differential forms

Requires: [Support](#support), [Degree](#degree).

See also: [Zero](#zero).

An expression denotes an array-valued differential form, a multi-form, on the
support of the expression, with the degree of the expression. A 0-form is an
array-valued function. An expression with no support and the degree `()` denotes a
single array. A `Zero` with a non-empty degree has no support and denotes the zero
form of that degree on every support.

### Node class or free function

Requires: [Properties and free functions](#properties-and-free-functions),
[`AbstractExpression`](#abstractexpression), [Free indices](#free-indices),
[Scalar type](#scalar-type), [Support](#support), [Degree](#degree),
[AbstractCoordinateDomain](#abstractcoordinatedomain).

See also: [Transpose and matrix product](#transpose-and-matrix-product), [Volume
form](#volume-form), [ExteriorDerivative](#exteriorderivative),
[Division](#division), [Elementary functions (very
incomplete)](#elementary-functions-very-incomplete),
[Determinant](#determinant), [Inverse](#inverse), [Tensor algebra with
EllipsisIndex](#tensor-algebra-with-ellipsisindex).

A concept must be a free function, not a node class, if the free function can
build the concept as one fixed composition of existing node classes for each
static type of the arguments. The static type of an expression is the value
shape, the free indices, the scalar type, the support and the degree, and the
static type of a coordinate domain is the dimension. A free function may read
the static type to choose the composition, and never reads the structure of
an argument, such as the class, the edges or the values. Otherwise the concept
gets its own node class. A concept that is evaluated by a numerical algorithm
that the user supplies also gets its own node class, even if a fixed
composition exists, because the composition would fix the algorithm. For
example, the transpose is always `component_tensor(A[k, l], (l, k))`, and
`volume_form(P)` is the product of the $n$ coordinate 1-forms of a coordinate
domain of dimension $n$, so both are free functions. No composition of
existing nodes differentiates, so `ExteriorDerivative` is a node class.
`Division`, `SquareRoot`, `Determinant` and `Inverse` are evaluated by a
numerical algorithm that the user supplies, so each is a node class. A free
function that builds a fixed composition is a definition, not a
simplification. The free function does not read the structure of the
operands, so the free function never removes structure that the user built.
A `EllipsisIndex` makes most tensor algebra one fixed composition for every
value shape (see [EllipsisIndex](#tensor-algebra-with-ellipsisindex)), so a free function rarely
needs to read a value shape.

## Algebra

### Literals

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices), [Scalar type](#scalar-type), [Support](#support),
[Degree](#degree).

See also: [Division](#division).

A literal is a number written down in an expression, e.g. `Integer(2) * u`. A
literal's value is part of the expression itself. Every literal has value shape
`()`, no free indices, no support and the degree `()`.

The literals are concrete subclasses of `AbstractExpression`:

```
Integer(2)
Float(3.5)
ImaginaryUnit()
```

Each literal constructor accepts only exact types and raises `TypeError` or
`ValueError` otherwise:

- `Integer(v)` takes an `int` $v$. A `bool` is rejected, although `bool` is a
  subclass of `int`.
- `Float(v)` takes a finite `float` $v$. NaN and infinities are rejected,
  because a NaN is not equal to itself. `Float(-0.0) != Float(0.0)`.
- `ImaginaryUnit()` takes no arguments and denotes $i$.

A rational number is the `Division` of two `Integer`s, such as
`Integer(1) / Integer(3)`. In floating-point arithmetic, the division of two
integers rounds once, the same as the conversion of a rational number, so the
core has no rational literal. A complex number is a composition such as
`a + ImaginaryUnit() * b`, so the core has no complex literal.

### Zero

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices), [Support](#support), [Degree](#degree).

See also: [Indexed](#indexed), [Applying derivatives](#applying-derivatives).

`Zero(value_shape, free_indices, degree)` is the zero with the given value
shape, free indices and degree, and has no support. `apply_derivatives` produces
`Zero` (see [Applying derivatives](#applying-derivatives)). No fixed
composition gives a zero with free indices or with a non-empty degree, because
`Integer(0) * a[i]` has edges to `a`, and $0 \, \mathrm{d}x$ has edges to
$\mathrm{d}x$. So `Zero` is a node class. For example,
`Zero((), (), (1,))` is the zero 1-form, and `Zero((3,), (), ())` is the zero
array of value shape `(3,)`.

### Arithmetic

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices).

See also: [ComponentTensor and IndexSum](#componenttensor-and-indexsum), [Tensor
algebra with EllipsisIndex](#tensor-algebra-with-ellipsisindex).

`Sum`, `Product`, `Division` and `Power` are concrete subclasses of
`AbstractExpression` that each store two expressions. Arithmetic is scalar.
Both operands must have value shape `()`, or construction raises
`ValueError`, and the result has value shape `()`. Free indices pass through
arithmetic (see [Free indices](#free-indices)), so tensor arithmetic is
written in index notation. For example, a scalar `a` times a tensor `T` of any
rank is `component_tensor(a * T[r], (r,))`, with `r` a
[EllipsisIndex](#ellipsisindex).

| Class      | Free function | Overload |
| ---------- | ------------- | -------- |
| `Sum`      | `add`         | `+`      |
| `Product`  | `multiply`    | `*`      |
| `Division` | `divide`      | `/`      |
| `Power`    | `power`       | `**`     |

### Sum

Requires: [Arithmetic](#arithmetic), [Degree](#degree).

See also: [ExteriorDerivative](#exteriorderivative).

`Sum(first, second)` adds two differential forms of the same degree. Both
operands must have the same degree, or construction raises `ValueError`, and
the sum has that degree. For example, $\mathrm{d}x + \mathrm{d}y$ has the
degree `(1,)`, and $f + \mathrm{d}x$ raises `ValueError`.

### Product

Requires: [Arithmetic](#arithmetic), [Free indices](#free-indices),
[Support](#support), [Degree](#degree),
[AbstractCoordinateDomain](#abstractcoordinatedomain).

See also: [TensorProduct](#tensorproduct), [ComponentTensor and
IndexSum](#componenttensor-and-indexsum),
[ExteriorDerivative](#exteriorderivative).

`Product(first, second)` is the wedge product $\alpha \wedge \beta$ of the
two operands. For two 0-forms, the wedge product is the ordinary product.
`Product` adds the degrees of the operands slot by slot, and an operand with
the degree `()` adds 0 to every slot, so $f \wedge \omega = f \omega$. If both
operands have a non-empty degree, the two degrees must have the same number of
slots, or construction raises `ValueError`. A slot of the result may have a
degree greater than the dimension of the support, and then the result is the zero
form of that degree. The degree rules of the other node classes are in
the [degree table](#appendix-degree-table).

The wedge product is graded commutative,
$\alpha \wedge \beta = (-1)^{kl} \beta \wedge \alpha$ for a $k$-form $\alpha$
and an $l$-form $\beta$. Equality is structural, so `Product` never reorders
the operands, and a later pass may apply the identity. `Product` never adds a
slot, so the product of two 1-forms is a 2-form with the degree `(2,)`, not a
double form with the degree `(1, 1)`.

No fixed composition of the other node classes multiplies, so `Product` is a
node class. `Product` and [`TensorProduct`](#tensorproduct) are the node
classes that combine the free indices of two operands that differ, so every
contraction is an `IndexSum` of a `Product` or a `TensorProduct`.

Examples, with a scalar $f$, the 1-forms $\mathrm{d}x$ and $\mathrm{d}y$ and a
double form $S$ with the degree `(1, 1)`:

| Product                                         | Degree              |
| ----------------------------------------------- | ------------------- |
| $f \cdot f$                                     | `()`                |
| $f \, \mathrm{d}x$                              | `(1,)`              |
| $\mathrm{d}x \wedge \mathrm{d}y$                | `(2,)`              |
| $\mathrm{d}x \wedge \mathrm{d}x$, the zero form | `(2,)`              |
| $f \, S$                                        | `(1, 1)`            |
| $\mathrm{d}x \wedge S$                          | raises `ValueError` |

### Power

Requires: [Arithmetic](#arithmetic), [Support](#support).

See also: [Applying derivatives](#applying-derivatives).

The exponent of `Power` has no support, or construction raises `ValueError`. The
derivative of a power whose exponent has a support needs a logarithm, which the
core does not have.

### Negation and subtraction

Requires: [Sum](#sum), [Product](#product), [Literals](#literals).

Negation and subtraction are free functions, not node classes:

- `negate(f)`, with the overload `-f`, is `multiply(Integer(-1), f)`.
- `subtract(f, g)`, with the overload `f - g`, is `add(f, negate(g))`.

Both compositions are exact in floating-point arithmetic, so neither
composition changes the arithmetic that the user wrote.

### Division

Requires: [Arithmetic](#arithmetic), [Product](#product), [Power](#power),
[Literals](#literals).

Division is a node class, because `multiply(f, power(g, Integer(-1)))`
rounds twice and differs from `f / g` in floating-point arithmetic.

### Absolute

Requires: [`AbstractExpression`](#abstractexpression), [Scalar
type](#scalar-type).

`Absolute` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `absolute(f)`, with the overload `abs(f)`, takes
a scalar and returns a scalar. For a complex argument, `Absolute` is the
modulus. `Absolute` is real for every argument.

### Conjugate

Requires: [`AbstractExpression`](#abstractexpression), [Scalar
type](#scalar-type).

See also: [Variational forms](#variational-forms).

`Conjugate` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `conjugate(f)` takes a scalar and returns the
complex conjugate of the scalar. No fixed composition of the other node classes
conjugates, so `Conjugate` is a node class. A sesquilinear variational form,
such as the $L^2$ inner product $\int u \bar{v}$, needs `Conjugate`.

### TensorProduct

Requires: [Arithmetic](#arithmetic), [Product](#product), [Free
indices](#free-indices), [Degree](#degree).

See also: [ExteriorDerivative](#exteriorderivative),
[ParameterisedDomain](#parameteriseddomain).

`TensorProduct(first, second)` is a concrete subclass of `AbstractExpression`
for the tensor product $\alpha \otimes \beta$ of two differential forms. The
free function is `tensor_product(f, g)`, with no overload.

`TensorProduct` concatenates the degrees of the operands, so the degree of the
result is `first.degree + second.degree`. An operand with the degree `()` adds
no slot, so $f \otimes \omega = f \omega$. Both operands must have value shape
`()`, as in [Arithmetic](#arithmetic), or construction raises `ValueError`.
The free indices combine as in `Product`.

The tensor product is not commutative, and `TensorProduct` never reorders the
operands. No fixed composition of the other node classes adds a slot, so
`TensorProduct` is a node class. `Product` adds degrees within a slot, and
`TensorProduct` adds slots.

Examples, with the setup of [ExteriorDerivative](#exteriorderivative) on
`RD(3)` and `i = index("i")`:

| Expression                                        | Degree              | Meaning                                                                                                 |
| ------------------------------------------------- | ------------------- | ------------------------------------------------------------------------------------------------------- |
| `tensor_product(dx, dy)`                          | `(1, 1)`            | $\mathrm{d}x \otimes \mathrm{d}y$                                                                       |
| `index_sum(tensor_product(dX[i], dX[i]), (i,))`   | `(1, 1)`            | The Euclidean metric                                                                                    |
| `tensor_product(dx * dy, dz)`                     | `(2, 1)`            | $(\mathrm{d}x \wedge \mathrm{d}y) \otimes \mathrm{d}z$                                                  |
| `tensor_product(dx, dy) * tensor_product(dy, dz)` | `(2, 2)`            | $(\mathrm{d}x \wedge \mathrm{d}y) \otimes (\mathrm{d}y \wedge \mathrm{d}z)$, the slotwise wedge product |
| `tensor_product(dX, dX)`                          | raises `ValueError` | `dX` has value shape `(3,)`                                                                             |

For a parameterised domain $D$ with the parameter domain $P$ and
`dpsi = exterior_derivative(D.parameterisation, P)`,
`index_sum(tensor_product(dpsi[i], dpsi[i]), (i,))` is the metric that $\psi$
induces on $P$, the pull-back of the Euclidean metric.

### Elementary functions (very incomplete)

Requires: [`AbstractExpression`](#abstractexpression), [Power](#power),
[Division](#division), [Literals](#literals).

`SquareRoot` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `square_root(f)` takes a scalar and returns a
scalar.

`SquareRoot` is a node class, not the composition `power(f, Integer(1) /
Integer(2))`. IEEE 754 requires a correctly rounded square root but does not
require a correctly rounded power, so the two differ in floating-point
arithmetic (i.e. the compiler/machine provides a specific implementation).

### Tensor

Requires: [`AbstractExpression`](#abstractexpression), [Degree](#degree).

See also: [Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative).

A concrete subclass of `AbstractExpression` that stacks expressions along a
new first axis.

`Tensor(entries)` stores the entries, a non-empty `tuple` of
`AbstractExpression`s that all have the same value shape `s` and the same
degree, or construction raises `ValueError`. `Tensor(entries)` derives the
value shape `(n,) + s`, with $n$ the number of entries, and the degree of the
entries. A matrix is a `Tensor` of `Tensor`s. A `Tensor` of differential forms
is an array-valued differential form, so the new axis is an axis of the value
shape and never a slot.

A value shape may contain the extent 0, for example the value shape `(0,)` of
the coordinate function of `RD(0)`. `Tensor` cannot build an expression with
the extent 0, because an empty `tuple` has no entry to give `s`.

Examples, with `zero, one, two = Integer(0), Integer(1), Integer(2)`:

- `tensor((one, two))` has value shape `(2,)` and no support.
- With the 1-forms $\mathrm{d}x$ and $\mathrm{d}y$ on `RD(2)`,
  `tensor((dx, dy))` has value shape `(2,)` and the degree `(1,)`, and is the
  vector-valued 1-form $\mathrm{d}X$.
- `tensor((tensor((one, zero)), tensor((zero, one))))` has value shape
  `(2, 2)` and no support.
- With `r0 = Float(0.0)` and `r1 = Float(1.0)`,
  `tensor((tensor((r1, r0)), tensor((r0, r1))))` has `Float` entries, so the
  matrix is not equal to `tensor((tensor((one, zero)), tensor((zero, one))))`.
- `tensor((tensor((one, two)),))` has value shape `(1, 2)`.
- `tensor((one, tensor((one, two))))` raises `ValueError`, because the value
  shapes of the entries differ.
- `tensor((one, dx))` raises `ValueError`, because the degrees of the entries
  differ.

### Indexed

Requires: [`AbstractExpression`](#abstractexpression), [Index](#index),
[EllipsisIndex](#ellipsisindex), [Free indices](#free-indices),
[Support](#support), [Degree](#degree).

See also: [Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative),
[Example: a paraboloid](#example-a-paraboloid).

> LLM generated content.

A concrete subclass of `AbstractExpression` that selects one component of a
tensor-valued expression.

`Indexed(expression, indices)` stores:

- The expression, an `AbstractExpression` of rank $r$.
- The indices, a `tuple` of entries. Each entry is an `int`, an `Index` or a
  `EllipsisIndex` (see [Index](#index) and [EllipsisIndex](#ellipsisindex)). At most one entry is a
  `EllipsisIndex`.

Indexing is full indexing only. Each `int` and each `Index` entry covers one
axis, and a `EllipsisIndex` entry covers every axis that the other entries do not
cover. Without a `EllipsisIndex`, the number of entries must equal $r$. With a
`EllipsisIndex`, the number of other entries must be at most $r$. Each `int` entry
$i_j$ must satisfy $0 \le i_j < n_j$, where `(n_0, ..., n_{r-1})` is the value
shape of the expression. Each `Index` entry $i_j$ takes the extent $n_j$. A
`EllipsisIndex` takes the `tuple` of the extents of the axes that the `EllipsisIndex`
covers.

`Indexed(expression, indices)` derives the value shape `()`, the support of the
expression, the degree of the expression, and the free indices of the
expression followed by the `Index` and `EllipsisIndex` entries. The rank $r$
counts the axes of the value shape only. The indices select axes of the value
shape and never a slot, so the slots of a differential form cannot be indexed.

`X[0]` and `A[0, 1]` are `indexed(X, (0,))` and `indexed(A, (0, 1))`.
Construction of `Indexed` raises:

| Error        | Cause                                                                                                     |
| ------------ | --------------------------------------------------------------------------------------------------------- |
| `TypeError`  | An entry that is neither exactly `int`, `Index` nor `EllipsisIndex`.                                      |
| `IndexError` | A number of entries that does not cover exactly the rank. An expression with value shape `()` has rank 0. |
| `ValueError` | More than one `EllipsisIndex` entry.                                                                      |
| `IndexError` | A negative `int` entry, or an `int` entry $i_j \ge n_j$.                                                  |
| `ValueError` | An `Index` that occurs twice, or that is already free.                                                    |

A negative index is rejected rather than converted, so the stored indices are
always the indices that the user wrote.

Examples:

- For the coordinate function $X$ of `RD(3)`, `X[2]` has value shape `()`.
- For a matrix `A` of value shape `(2, 2)`, `A[1, 0]` has value shape `()`,
  and `A[1]` raises `IndexError`.
- `X[-1]` and `X[True]` raise errors.
- `tensor((Integer(1) - X[0], X[0]))` has value shape `(2,)` and the support
  `RD(2)`.
- `tensor((tensor((X[0], X[1])),))` has value shape `(1, 2)`.
- For the vector-valued 1-form `dX = tensor((dx, dy))` on `RD(2)`, `dX[0]` is
  the 1-form $\mathrm{d}x$, with value shape `()` and the degree `(1,)`.
- In the paraboloid example, the height `X[0]**Integer(2) + X[1]**Integer(2)`
  is an expression on `RD(2)` with value shape `()`.

### ComponentTensor and IndexSum

Requires: [Index](#index), [EllipsisIndex](#ellipsisindex), [Free
indices](#free-indices), [Degree](#degree), [Indexed](#indexed).

See also: [ExteriorDerivative](#exteriorderivative).

`ComponentTensor` and `IndexSum` are concrete subclasses of
`AbstractExpression` that bind free indices.
`ComponentTensor(expression, indices)` turns the indices into tensor axes.
`IndexSum(expression, indices)` sums over the indices.

Both store a scalar expression and a non-empty `tuple` of distinct `Index` and
`EllipsisIndex` entries that are free in the expression. Construction raises
`ValueError` if the expression is not scalar, if the indices are empty, if an
index is repeated, or if an index is not free in the expression. Both derive the
free indices of the expression without the bound indices, and the degree of the
expression. The indices bind axes of the value shape, so a `ComponentTensor` of
a differential form is an array-valued differential form, and an `IndexSum` of
differential forms is a sum of differential forms of the same degree.

`ComponentTensor` derives the value shape from the extents of the indices, in
the order of the indices. A `EllipsisIndex` contributes every extent in the
`tuple` extent of the `EllipsisIndex`. `IndexSum` derives the value shape `()`.

One `IndexSum` holds a whole contraction, such as the sum over both indices of
`A[i, j] * B[i, j]`. A pass that chooses the order of summation then finds the
contraction in one node, not in a chain of nested nodes.

Examples, with `i, j = index("i"), index("j")`, two vectors `a` and `b` of value
shape `(3,)`, and a matrix `A` of value shape `(2, 2)`:

- `A[i, j]` has value shape `()` and the free indices `((i, 2), (j, 2))`.
- `a[i] * b[i]` has value shape `()` and the free indices `((i, 3),)`.
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
- `a[i] + b[j]` raises `ValueError`, because the free indices differ.
- `A[i, i]` raises `ValueError`, because `i` occurs twice.

> Note: `opt_einsum` can directly plan a contract path for a given contraction,
> but assumes dense operators.

### Tensor algebra with EllipsisIndex

Requires: [EllipsisIndex](#ellipsisindex), [Indexed](#indexed), [ComponentTensor
and IndexSum](#componenttensor-and-indexsum), [Sum](#sum), [Product](#product),
[Conjugate](#conjugate).

A `EllipsisIndex` writes an operation on tensors of every rank as one fixed
composition. With derived ellipsis indices `r` and `s` and an ordinary index
`k`:

| Operation          | Fixed composition                                              |
| ------------------ | -------------------------------------------------------------- |
| `scale(a, T)`      | `component_tensor(a * T[r], (r,))`                             |
| tensor sum `A + B` | `component_tensor(A[r] + B[r], (r,))`                          |
| `outer(u, v)`      | `component_tensor(u[r] * v[s], (r, s))`                        |
| `dot(u, v)`        | `component_tensor(index_sum(u[r, k] * v[k, s], (k,)), (r, s))` |
| `inner(A, B)`      | `index_sum(A[r] * conjugate(B[r]), (r,))`                      |

`A[r] + B[r]` raises `ValueError` if the value shapes of `A` and `B` differ,
because the extents of `r` differ.

### Determinant

Requires: [`AbstractExpression`](#abstractexpression), [Degree](#degree).

`Determinant` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `determinant(f)` takes an expression of value
shape `(n, n)` and returns a scalar. The determinant of an expression of
value shape `(n, n)` is a polynomial of degree $n$ in the entries. For each
$n$, the cofactor expansion is a fixed composition, but the determinant is
evaluated by a numerical algorithm that the user supplies, such as a closed-form
formula or an LU factorisation. So `Determinant` is a node class.

Construction raises `ValueError` if the expression does not have rank 2, or if
the two extents of the value shape differ.

The determinant of an expression of value shape `(0, 0)` is 1, the empty
product.

### Inverse

Requires: [Determinant](#determinant), [Division](#division).

See also: [Components](#components).

`Inverse` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `inverse(f)` takes an expression of value shape
`(n, n)` and returns the matrix inverse, with value shape `(n, n)`. The entries
of the inverse are the cofactors divided by the determinant. For each $n$,
that is a fixed composition, but the inverse is evaluated by a numerical
algorithm that the user supplies, such as a closed form or an LU
factorisation. So `Inverse` is a node class.

Construction raises `ValueError` if the expression does not have rank 2, or if
the two extents of the value shape differ. Both `Determinant` and `Inverse`
require the degree `()`.

The inverse is defined where the determinant is not zero. Construction does
not check invertibility, as construction of `Division` does not check for a
zero divisor. The inverse of an expression of value shape `(0, 0)` has value
shape `(0, 0)`.

Examples, with a matrix `A` of value shape `(2, 2)`:

- `inverse(A) @ A` is the identity matrix where `A` is invertible, after a pass.
- For a metric `g` on $P$, `inverse(components(components(g, P), P))` is the
  matrix of the inverse metric, which raises the index of a 1-form.
- In hyperelasticity, `inverse(F)` is the inverse of the deformation gradient.

### Transpose and matrix product

Requires: [Index](#index), [Indexed](#indexed), [ComponentTensor and
IndexSum](#componenttensor-and-indexsum), [Product](#product).

Transpose and matrix product are free functions built with derived indices:

- `transpose(A)` is `component_tensor(A[k, l], (l, k))`, with
  $k$ = `Index(("uflx.transpose", A, 0))` and
  $l$ = `Index(("uflx.transpose", A, 1))`.
- `matrix_product(A, B)`, with the overload `A @ B`, is
  `component_tensor(index_sum(A[i, k] * B[k, j], (k,)), (i, j))`, with $i$, $k$
  and $j$ the derived indices `Index(("uflx.matrix_product", A, B, m))` for
  $m = 0, 1, 2$.

An operand that does not have rank 2, including a scalar operand, is rejected
with `IndexError` by the constructor of `Indexed`. Unequal inner extents in
`matrix_product` are rejected with `ValueError` by the constructor of `Product`,
because the index $k$ has two extents. So the overload `@` accepts only rank 2
operands. Contractions of other ranks, such as a matrix-vector product, are
written in index notation.

Examples, with a matrix `A` of value shape `(2, 2)`:

- `transpose(A)` is not equal to `component_tensor(A[i, j], (j, i))` for user
  indices `i` and `j`, because the index labels differ.

### Operator overloads

Requires: [`AbstractExpression`](#abstractexpression), [Sum](#sum),
[Product](#product), [Division](#division), [Power](#power), [Negation and
subtraction](#negation-and-subtraction), [Absolute](#absolute),
[Indexed](#indexed), [Transpose and matrix
product](#transpose-and-matrix-product).

See also: [Numbers in expressions](#numbers-in-expressions), [Appendix: Free
functions](#appendix-free-functions).

Every overload in the [table of free functions](#appendix-free-functions) is a method of `AbstractExpression` (see
[`AbstractExpression`](#abstractexpression)). No other class defines
arithmetic, matrix or indexing overloads. Overloading should be considered
carefully outside of standard mathematical meanings.

An operator overload calls the free function with the operands unchanged, and
the constructor validates (see [Numbers in expressions](#numbers-in-expressions)).
`AbstractExpression` defines no reflected overloads, such as `__radd__`.

`__getitem__` is the one exception. Python passes `f[i]` to `__getitem__` as the
single entry `i`, and passes `f[i, j]` as the `tuple` `(i, j)`. `__getitem__`
converts a single entry into a `tuple` of one entry and calls `indexed`. The
conversion exists because of Python indexing syntax.

This base class is sketched in pseudo-code according to the principles outlined
above. Later classes will not be sketched in detail.

```
class AbstractExpression(Base, ABC):
    """Abstract base class for expressions."""
    @property
    @abstractmethod
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""

    @property
    @abstractmethod
    def free_indices(
        self,
    ) -> tuple[tuple[Index | EllipsisIndex, int | tuple[int, ...]], ...]:
        """The free indices and extents, in order of first occurrence."""

    @property
    @abstractmethod
    def scalar_type(self) -> Literal["real", "complex"]:
        """The scalar type of the expression."""

    @property
    @abstractmethod
    def support(self) -> AbstractCoordinateDomain | None:
        """The support of the expression, or None for an expression with no support."""

    @property
    @abstractmethod
    def degree(self) -> tuple[int, ...]:
        """The form degree of each slot, with every entry at least 1."""

    # Operator overloading must be in the class.
    def __add__(self, other: AbstractExpression) -> AbstractExpression:
        """Add."""
        return add(self, other)

    def __sub__(self, other: AbstractExpression) -> AbstractExpression:
        """Subtract."""
        # Calls the free function subtract.
        ...

    def __mul__(self, other: AbstractExpression) -> AbstractExpression:
        """Multiply."""
        # Calls the free function multiply.
        ...

    def __matmul__(self, other: AbstractExpression) -> AbstractExpression:
        """Matrix multiply."""
        # Calls the free function matrix_product.
        ...

    def __truediv__(self, other: AbstractExpression) -> AbstractExpression:
        """Divide."""
        # Calls the free function divide.
        ...

    def __neg__(self) -> AbstractExpression:
        """Negate."""
        # Calls the free function negate.
        ...

    def __abs__(self) -> AbstractExpression:
        """Absolute value."""
        # Calls the free function absolute.
        ...

    def __pow__(self, exponent: AbstractExpression) -> AbstractExpression:
        """Raise to a power."""
        # Calls the free function power.
        ...

    def __getitem__(self, indices: object) -> AbstractExpression:
        """Index."""
        # Converts a single entry into a tuple of one entry, then calls indexed.
        ...

    def __iter__(self) -> NoReturn:
        """Raise TypeError, because expressions are not iterable."""
        # Without __iter__, Python iterates through __getitem__.
        ...

# A free function
def add(first: AbstractExpression, second: AbstractExpression) -> Sum:
    """Add two scalar differential forms of the same degree.

    Args:
        first: The first summand.
        second: The second summand.

    Returns:
        The sum of the two expressions, with the summands in the given order.

    Raises:
        TypeError: From the constructor of Sum, if a summand is not an
            expression.
        ValueError: From the constructor of Sum, if a summand does not have
            value shape (), if the free indices of the summands differ, if
            the summands have different supports, or if the degrees differ.
    """
    # The constructor of Sum validates.
    return Sum(first, second)
```

### Numbers in expressions

Requires: [Literals](#literals), [Operator overloads](#operator-overloads),
[Tensor](#tensor), [Indexed](#indexed), [RD](#rd).

See also: [Coordinate](#coordinate).

Every number in an expression is written with literals, such as `Integer(2)`,
`Float(0.5)`, `Integer(1) / Integer(3)` or `ImaginaryUnit()`. A Python `int`, `float`,
`complex`, `Fraction` or NumPy scalar is never converted and is rejected with
`TypeError`. The type of every literal is a choice that the user makes
explicitly.

Examples, with `X = coordinate(RD(2))`:

- `X[0]**Integer(2)`, `Integer(2) * X[0]` and `Integer(1) - X[0]` are
  accepted.
- `X[0]**2`, `2 * X[0]`, `1 - X[0]` and `0.5 * X[0]` raise `TypeError`.
- `tensor((1, 2))` raises `TypeError`. `tensor((Integer(1), Integer(2)))` is
  accepted.
- `2 + X[0]` raises `TypeError`, because `AbstractExpression` defines no
  reflected overloads.

## Calculus on a coordinate domain

### Coordinate

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain),
[`AbstractExpression`](#abstractexpression), [Support](#support),
[Degree](#degree).

See also: [ExteriorDerivative](#exteriorderivative),
[ParameterisedDomain](#parameteriseddomain).

A concrete subclass of `AbstractExpression` for the coordinate function of a
coordinate domain.

Every concrete coordinate domain $P$ provides a coordinate function
`Coordinate(P)`. The coordinate function $X$ of $P$ is the terminal expression
that maps each point of $P$ to the $n$-tuple of that point:

$X : P \to R^n$, $X(p) = p$.

`Coordinate(P)` stores the coordinate domain $P$ and derives the value shape
`(n,)`, the support $P$ and the degree `()`. The coordinate function is a 0-form.
The 1-forms $\mathrm{d}X_i$ come from the exterior derivative, not from
`Coordinate`.

Only coordinate domains provide a coordinate function. The coordinates of the
points of a domain are given by the parameterisation of the domain (see
[ParameterisedDomain](#parameteriseddomain)), which is built from the coordinate function
of the parameter domain.

Examples:

- `coordinate(RD(3))` is `Coordinate(RD(3))`, with value shape `(3,)` and the
  support `RD(3)`.
- `coordinate(RD(3)) == coordinate(RD(3))` is `True`, and
  `coordinate(RD(2)) == coordinate(RD(3))` is `False`.
- In the paraboloid example, `X = coordinate(RD(2))` has value shape `(2,)`.

### Valid pair

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain),
[Support](#support).

An expression and a coordinate domain $(f, P)$ are a valid pair if $P$ is an
`AbstractCoordinateDomain` and $f$ has the support $P$ or no support. Otherwise
construction raises `TypeError` for $P$ or `ValueError` for $f$. An $f$ with
no support is a constant on $P$.

### ExteriorDerivative

Requires: [Valid pair](#valid-pair), [Coordinate](#coordinate),
[Degree](#degree), [Product](#product), [Indexed](#indexed).

See also: [Applying derivatives](#applying-derivatives).

A concrete subclass of `AbstractExpression` for the exterior derivative
$\mathrm{d}\omega$ of a differential form $\omega$ on a coordinate domain $P$.

`ExteriorDerivative(f, P)` stores $f$ and $P$ and derives the value shape of
$f$ and the support $P$. The free indices and the scalar type follow the general
rule. `ExteriorDerivative(f, P)` requires a valid pair $(f, P)$ (see
[Valid pair](#valid-pair)). An $f$ with no support is a constant on
$P$.

The exterior derivative acts on differential forms with at most one slot. If
$f$ has the degree `()`, the exterior derivative has the degree `(1,)`. If $f$
has the degree `(k,)`, the exterior derivative has the degree `(k + 1,)`. If
$f$ has two or more slots, construction raises `ValueError`, because the
exterior derivative of a double form needs a connection. A degree greater than
the dimension of $P$ gives the zero form, as for [`Product`](#product).

The exterior derivative acts on each component of an array-valued differential
form, so the axes of the value shape pass through unchanged.

No fixed composition of the other node classes differentiates, so
`ExteriorDerivative` is a node class. `ExteriorDerivative` is lazy.
Construction does not evaluate the derivative. The core transformation
`apply_derivatives` evaluates the derivative with linearity, the graded
Leibniz rule
$\mathrm{d}(\alpha \wedge \beta) = \mathrm{d}\alpha \wedge \beta + (-1)^k \alpha \wedge \mathrm{d}\beta$
for a $k$-form $\alpha$, and $\mathrm{d}\mathrm{d}\omega = 0$ (see
[Applying derivatives](#applying-derivatives)). The exterior derivative of
`Coordinate(P)` is the vector-valued 1-form $\mathrm{d}X$ and remains.

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

A pass that applies the rules of `apply_derivatives` gives the familiar
results of vector calculus:

| Expression                                   | Degree | After a pass                                                               |
| -------------------------------------------- | ------ | -------------------------------------------------------------------------- |
| `d(x * y * z)`                               | `(1,)` | $yz \, \mathrm{d}x + xz \, \mathrm{d}y + xy \, \mathrm{d}z$, the gradient  |
| `d(-y * dx + x * dy)`                        | `(2,)` | $2 \, \mathrm{d}x \wedge \mathrm{d}y$, the curl $(0, 0, 2)$                |
| `d(x * dy * dz + y * dz * dx + z * dx * dy)` | `(3,)` | $3 \, \mathrm{d}x \wedge \mathrm{d}y \wedge \mathrm{d}z$, the divergence 3 |
| `d(d(x * y * z))`                            | `(2,)` | The zero 2-form, because the curl of a gradient is zero                    |
| `d(d(-y * dx + x * dy))`                     | `(3,)` | The zero 3-form, because the divergence of a curl is zero                  |
| `d(Integer(1))`                              | `(1,)` | The zero 1-form, because a constant has no gradient                        |

The second row is the reason that half the integral of
$-y \, \mathrm{d}x + x \, \mathrm{d}y$ around a closed curve in the plane is
the enclosed area, by Stokes' theorem.

- `d(dx * dy * dz)` has the degree `(4,)` and is the zero form, because
  `RD(3)` has dimension 3.
- With a double form $S$ with the degree `(1, 1)`, such as a metric, `d(S)`
  raises `ValueError`.
- `exterior_derivative(X, RD(2))` raises `ValueError`, because `X` has the
  support `RD(3)`.

### Components

Requires: [Valid pair](#valid-pair), [ExteriorDerivative](#exteriorderivative),
[Degree](#degree), [ComponentTensor and
IndexSum](#componenttensor-and-indexsum).

See also: [PullBack](#pullback), [Expanding components](#expanding-components),
[ParameterisedDomain](#parameteriseddomain).

A concrete subclass of `AbstractExpression` that turns the first slot of a
differential form into an axis of the value shape, with the coordinate frame of
a coordinate domain $P$ of dimension $n$.

`Components(f, P)` stores $f$ and $P$ and derives the support $P$. The free
indices and the scalar type follow the general rule. `Components(f, P)`
requires a valid pair $(f, P)$ (see [Valid pair](#valid-pair)),
and $f$ must have at least one slot, or construction raises `ValueError`.

With $\partial_j$ the coordinate vector fields of $P$, dual to the coordinate
1-forms $\mathrm{d}X_j$, the entry $j$ of the new axis is the interior product
of $\partial_j$ with the first slot of $f$,

$$
\mathrm{Components}(f, P)_j(v_2, \dots, v_k) = f(\partial_j, v_2, \dots, v_k).
$$

`Components(f, P)` derives the value shape `f.value_shape + (n,)`. If the
first slot of $f$ has the degree $k > 1$, the first slot of the result has the
degree $k - 1$. If the first slot of $f$ has the degree 1, the result does not
have that slot. For example, the degree `(1,)` gives `()`, the degree `(2,)`
gives `(1,)`, and the degree `(1, 1)` gives `(1,)`.

No fixed composition of the other node classes turns a slot into an axis, so
`Components` is a node class. The reverse direction is a fixed composition,
`index_sum(v[i] * dX[i], (i,))` for a vector `v` of value shape `(n,)`.

The components depend on the coordinates of $P$. So a pull-back does not move
inside `Components`. `PullBack(Components(f, Q), D)` is the array of the
components of $f$ on $Q$, evaluated at $\psi(p)$, and remains until
`expand_components` evaluates the `Components` (see
[Expanding components](#expanding-components)).

With `Components`, the arrays of standard differential calculus are fixed
compositions. With the examples of [ExteriorDerivative](#exteriorderivative)
on `RD(3)`, and `f = x * y * z`:

| Expression                              | Value shape | Degree | Meaning                                                      |
| --------------------------------------- | ----------- | ------ | ------------------------------------------------------------ |
| `components(d(f), P)`                   | `(3,)`      | `()`   | The gradient $(yz, xz, xy)$, after evaluation                |
| `components(d(components(d(f), P)), P)` | `(3, 3)`    | `()`   | The Hessian of $f$                                           |
| `components(dX, P)`                     | `(3, 3)`    | `()`   | The identity matrix                                          |
| `components(dx * dy, P)`                | `(3,)`      | `(1,)` | $(\mathrm{d}y, -\mathrm{d}x, 0)$                             |
| `components(components(dx * dy, P), P)` | `(3, 3)`    | `()`   | The antisymmetric matrix of $\mathrm{d}x \wedge \mathrm{d}y$ |

- `components(x, P)` raises `ValueError`, because `x` has no slot.
- For a parameterised domain $D$ with the parameter domain $P$,
  `components(exterior_derivative(D.parameterisation, P), P)` is the Jacobian
  matrix of value shape `(m, n)`, with $m$ the dimension of the ambient
  coordinate domain of $D$.

### Volume form

Requires: [Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative),
[Product](#product), [Literals](#literals).

See also: [Integral](#integral), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

The free function `volume_form(P)` is the coordinate volume form of a
coordinate domain $P$ of dimension $n$,

$$
\mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}.
$$

With `dX = exterior_derivative(coordinate(P), P)`, `volume_form(P)` is
`dX[0] * ... * dX[n - 1]`, with the degree `(n,)`. For $n = 0$,
`volume_form(P)` is `Integer(1)`. The free function reads the dimension of
$P$, so the free function is one fixed composition for each dimension. The
volume form fixes the orientation of $P$ (see [Integral](#integral)). The
volume form needs no metric. A Riemannian volume form is a density times
`volume_form(P)`.

Examples:

- `volume_form(RD(3))` has the degree `(3,)` and is
  $\mathrm{d}x \wedge \mathrm{d}y \wedge \mathrm{d}z$.
- With `T = UnitSimplex(2)`, `integral(volume_form(T), T)` is the area $1/2$
  of `T`.
- `volume_form(RD(0))` is `Integer(1)`.

### Integral

Requires: [Valid pair](#valid-pair), [Volume form](#volume-form),
[Support](#support), [Degree](#degree).

See also: [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

At this point we have everything we need to integrate an expression over a
single domain. Note we do not need to explicitly introduce any
topological notions of the domain to build a basic integral.

`Integral(f, P)` is the integral of a top form $f$ over a coordinate domain
$P$ of dimension $n$, with the orientation
$\mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}$ of the coordinates $X$
of $P$. With $f = f_{0 \dots n-1} \, \mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}$,

$$
\int_P f = \int_P f_{0 \dots n-1} \,\mathrm{d}X,
$$

with the Lebesgue measure $\mathrm{d}X$ on $P$. For $n = 0$, the top form is a
0-form, and the integral is the value of $f$ at the single point of $R^0$.

`Integral(f, P)` requires a valid pair, and $f$
must have the degree `(n,)`, or the degree `()` for $n = 0$, or construction
raises `ValueError`.

`Integral(f, P)` stores $f$ and $P$ and derives the value shape, the free
indices and the scalar type of $f$, and the degree `()`. `Integral(f, P)` has
no support. An array-valued differential form integrates component by component.

`Integral` is lazy. Construction does not evaluate the integral. A later pass
might evaluate the integral, for example by quadrature.

Note that an `Integral` is an expression, so integrals can be added,
multiplied and used wherever an expression is accepted.

#### Examples

With `T = UnitSimplex(2)` and `X = coordinate(T)`:

- `integral(X[0] * volume_form(T), T)` is $\int_T X_0 \,\mathrm{d}X$. A later
  pass evaluates the integral to $1/6$.
- `integral(X[0], T)` raises `ValueError`, because `X[0]` is a 0-form, not a
  2-form.
- With `four, zero = Integer(4), Integer(0)` and the metric
  `G = tensor((tensor((four, zero)), tensor((zero, Integer(1)))))` on `T`,
  `integral(square_root(determinant(G)) * volume_form(T), T)` is the area of
  `T` with the metric `G`, which is $\sqrt{4} \cdot 1/2 = 1$.
- `integral(X[0] * volume_form(T), T) * integral(X[1] * volume_form(T), T)` is a
  product of two integrals, with value shape `()` and no support.

## Maps between domains

### Domains

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain), [Node class or
free function](#node-class-or-free-function).

See also: [ParameterisedDomain](#parameteriseddomain).

A domain $D$ is the image $\psi(P)$ of a coordinate domain $P$ under a
parameterisation $\psi$ into a coordinate domain $Q$. The core has one class of
domain, [`ParameterisedDomain`](#parameteriseddomain), which stores $P$, $\psi$
and $Q$.

A particular kind of domain, such as an affinely mapped triangle or a point, is
a free function that returns a `ParameterisedDomain`, not a subclass. Every
kind of domain is a fixed composition of $P$, $\psi$ and $Q$, so every kind of
domain is a free function (see
[Node class or free function](#node-class-or-free-function)).

### ParameterisedDomain

Requires: [Domains](#domains),
[AbstractCoordinateDomain](#abstractcoordinatedomain),
[Coordinate](#coordinate), [Support](#support), [Free indices](#free-indices),
[Scalar type](#scalar-type).

See also: [PullBack](#pullback), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

A concrete subclass of `Base`.

The domain $D$ = `ParameterisedDomain(P, psi, Q)` stores:

- The parameter domain $P$, an `AbstractCoordinateDomain`, returned by the
  property `parameter_domain`. The dimension $k$ of $P$ is the topological
  dimension of $D$.
- The parameterisation $\psi$, an `AbstractExpression` on $P$ with value shape
  `(n,)`, returned by the property `parameterisation`. The parameterisation
  $\psi$ is built from the coordinate function $X$ of $P$, or is a constant
  with no support. The parameterisation $\psi$ maps each point of $P$ to a point of
  $Q$, and the image $\psi(P)$ is $D$.
- The ambient coordinate domain $Q$, an `AbstractCoordinateDomain`, returned by
  the property `ambient_domain`. The coordinates of $Q$ name the points of $D$.
  The dimension $n$ of $Q$ is the geometric dimension of $D$.

The 'geometric dimension' is `D.ambient_domain.dimension`, and the 'topological
dimension' is `D.parameter_domain.dimension`.

Construction raises `TypeError` if $P$ or $Q$ is not an
`AbstractCoordinateDomain`. Construction raises `ValueError` if $\psi$ does not
have value shape `(n,)` with $n$ the dimension of $Q$, if $\psi$ has a support
other than $P$, if $\psi$ has free indices, if $\psi$ is complex, or if $k >
n$. Construction does not check that $\psi(P)$ lies in $Q$.

$Q$ is usually `RD(n)`. A domain in a smaller $Q$, such as an edge of the unit
triangle in the unit triangle, pulls back an expression on $Q$ (see
[PullBack](#pullback)).

Examples:

- A sphere in $R^3$ has geometric dimension 3, topological dimension 2 and the
  ambient coordinate domain `RD(3)`.
- With `T = UnitSimplex(2)`, `parameterised_domain(T, coordinate(T), T)` is the
  unit triangle as a domain, with $k = n = 2$.
- The hypotenuse of the unit triangle is a parameterised domain with $k = 1$
  and $n = 2$. $P$ is the unit interval in $R^1$, and
  $\psi(X) = (1 - X_0, X_0)$. $Q$ is `RD(2)` for the hypotenuse in the plane,
  and `UnitSimplex(2)` for the hypotenuse in the unit triangle.
- With `X = coordinate(RD(2))` and `h = X[0]**Integer(2) + X[1]**Integer(2)`,
  `parameterised_domain(RD(2), tensor((X[0], X[1], h)), RD(3))` is the
  paraboloid, with $k = 2$ and $n = 3$.
- A physical triangle in $R^2$ is a parameterised domain with $k = n = 2$. $P$
  is the unit triangle, $\psi$ is the affine map from $P$ onto the physical
  triangle, and $Q$ is `RD(2)`. See the affinely mapped triangle example.
- A point $p$ in $R^2$ is a parameterised domain with $k = 0$ and $n = 2$. $P$
  is `RD(0)`, $\psi$ is the constant `tensor((p0, p1))`, with no support, and $Q$
  is `RD(2)`.

### Example: a paraboloid

Requires: [RD](#rd), [Coordinate](#coordinate), [Tensor](#tensor),
[Indexed](#indexed), [ParameterisedDomain](#parameteriseddomain).

The paraboloid over `RD(2)` collects the examples of the sections above.

```python
# AbstractCoordinateDomain: RD(2), no edges.
P = RD(2)

# AbstractExpression: Coordinate(RD(2)), value shape (2,).
X = coordinate(P)

# AbstractExpression: Tensor((Indexed(X, (0,)), Indexed(X, (1,)), Sum(...))),
# value shape (3,).
psi = tensor((X[0], X[1], X[0]**Integer(2) + X[1]**Integer(2)))

# ParameterisedDomain(RD(2), psi, RD(3)), k = 2, n = 3.
D = parameterised_domain(P, psi, RD(3))
```

### Example: an affinely mapped triangle

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain),
[Coordinate](#coordinate), [Tensor](#tensor), [Components](#components),
[ParameterisedDomain](#parameteriseddomain).

The physical triangle with vertices $v_0 = (1, 0)$, $v_1 = (3, 1)$ and
$v_2 = (2, 3)$ is the image of the unit triangle under an affine map.

Important point is that there is no topology here, just geometry.

Define the unit simplex in $R^n$,

$\{p \in R^n : p_i \ge 0 \text{ for every } i, \ \sum_i p_i \le 1\}$,

as a concrete `AbstractCoordinateDomain`.

```python
class UnitSimplex(AbstractCoordinateDomain):
    def __init__(self, n: int):
        if n < 0:
            raise ValueError("Dimension must be non-negative.")
        self._n = n

    @property
    def dimension(self) -> int:
        return self._n
```

`UnitSimplex(2)` is the unit triangle with vertices $(0, 0)$, $(1, 0)$ and
$(0, 1)$. The affine map from `UnitSimplex(2)` onto the physical triangle is

$$\psi(X) = v_0 + X_0 (v_1 - v_0) + X_1 (v_2 - v_0)
= (1 + 2 X_0 + X_1, X_0 + 3 X_1).$$

```python
# AbstractCoordinateDomain: UnitSimplex(2), no edges.
T = UnitSimplex(2)

# AbstractExpression: Coordinate(UnitSimplex(2)), value shape (2,).
X = coordinate(T)

# AbstractExpression: Tensor(..., ...), value shape (2,).
one, two, three = Integer(1), Integer(2), Integer(3)
psi = tensor((one + two * X[0] + X[1], X[0] + three * X[1]))

# ParameterisedDomain(UnitSimplex(2), psi, RD(2)), k = n = 2.
K = parameterised_domain(T, psi, RD(2))
```

The Jacobian matrix of $\psi$ is the constant matrix

$$J = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix},$$

so $\det J = 5$, and the area of the physical triangle is $5$ times the area
$1/2$ of the unit triangle.

The physical triangle is a parameterised domain with $k = n$, not a coordinate
domain, although the physical triangle has non-empty interior in $R^2$. The
class records how the domain is built.

### PullBack

Requires: [ParameterisedDomain](#parameteriseddomain), [Support](#support),
[Degree](#degree), [ExteriorDerivative](#exteriorderivative), [Example: an
affinely mapped triangle](#example-an-affinely-mapped-triangle).

See also: [Expanding pull-backs](#expanding-pull-backs).

A concrete subclass of `AbstractExpression` for the pull-back $\psi^* f$ of a
differential form $f$ on $Q$ by the parameterisation $\psi$ of a domain $D$,
with $P$ the parameter domain and $Q$ the ambient coordinate domain of $D$.

`PullBack(f, D)` stores $f$ and $D$ and derives the value shape of $f$, the
degree of $f$ and the support $P$. The free indices and the scalar type follow the
general rule. The second operand is a domain, not an expression, because a
domain states $Q$. Construction raises `ValueError` unless $f$ has no support or
the support of $f$ is $Q$. Construction raises `TypeError` if $D$ is not a
`ParameterisedDomain`.

For a 0-form $f$, the pull-back is the composition $f \circ \psi$, the
expression on $P$ with the value $f(\psi(p))$ at each point $p$ of $P$. For a
differential form with slots, the pull-back acts on each slot,

$$
(\psi^* f)_p(v_1, \dots, v_k) = f_{\psi(p)}(\mathrm{d}\psi \, v_1, \dots, \mathrm{d}\psi \, v_k),
$$

so `PullBack` accepts any degree, including a double form such as a metric.
The axes of the value shape are not transformed. A slot with a degree greater
than the dimension of $P$ gives the zero form, so a 2-form pulled back to a
curve, and a 1-form pulled back to a point, are zero.

`PullBack` is lazy. The core transformation `expand_pull_backs` moves a
`PullBack` towards the leaves, because the pull-back commutes with pointwise
operations, with the wedge product and with the exterior derivative, and
replaces `PullBack(Coordinate(Q), D)` by $\psi$ (see
[Expanding pull-backs](#expanding-pull-backs)). A `PullBack` of a downstream
terminal remains, and the external library that supplies the terminal
evaluates the terminal at $\psi(p)$ and transforms the slots.

A domain $E$ in the parameter domain of $D$ has the image $\psi(E)$, a subset
of $D$. The image is the parameterised domain with the parameterisation
`pull_back(D.parameterisation, E)` on the parameter domain of $E$.

Example, with `T`, `X`, `psi`, `one` and `K` of the affinely mapped triangle:

```python
# The unit interval, the parameter domain of the edge, with coordinate t.
S = UnitSimplex(1)
t = coordinate(S)[0]

# The edge of T from (1, 0) to (0, 1), as a domain in T.
e = parameterised_domain(S, tensor((one - t, t)), T)

# PullBack(f, e), with the support UnitSimplex(1), value shape () and the degree ().
f = X[0] * X[1]
fe = pull_back(f, e)

# PullBack(dx, e), with the support UnitSimplex(1) and the degree (1,).
dx = exterior_derivative(X[0], T)
dxe = pull_back(dx, e)

# The edge of K from (3, 1) to (2, 3), as a domain in RD(2).
E = parameterised_domain(S, pull_back(psi, e), RD(2))
```

- A later pass evaluates `fe` to $(1 - t) t$ and the parameterisation of `E`
  to $(3 - t, 1 + 2 t)$.
- A later pass evaluates `dxe` to $\mathrm{d}(1 - t) = -\mathrm{d}t$, because
  $X_0 = 1 - t$ on the edge.
- `pull_back(dx * exterior_derivative(X[1], T), e)` has the degree `(2,)` and
  is the zero form, because the edge has dimension 1.
- The Jacobian matrix of the parameterisation of `E` is
  $\begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix} \begin{pmatrix} -1 \\ 1
  \end{pmatrix} = \begin{pmatrix} -1 \\ 2 \end{pmatrix}$, so the edge of `K`
  has the length $\sqrt{5}$.
- `pull_back(f, parameterised_domain(S, tensor((one - t, t)), RD(2)))`
  raises `ValueError`, because that domain lies in `RD(2)`, not in `T`.

### Points and point evaluation

Requires: [RD](#rd), [ParameterisedDomain](#parameteriseddomain),
[PullBack](#pullback), [Integral](#integral), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

A point is a domain! A point in a coordinate domain $Q$ of dimension $n$ is
`parameterised_domain(RD(0), x, Q)`, with $x$ a constant expression of value
shape `(n,)`. The core needs no node class and no free function for a point.

The value of a 0-form $f$ of any value shape at a point $p$ is
`integral(pull_back(f, p), RD(0))`, because the integral of a 0-form over
`RD(0)` is the value at the single point (see [Integral](#integral)).
`pull_back(f, p)` is an expression on `RD(0)`, and the integral has no support, so
the integral combines with any expression.

A point in `T` is given in the coordinates of `T`. For a point $x$ of `K`, the
surrounding environment supplies $\psi^{-1}(x)$. For a non-affine $\psi$, the
surrounding environment supplies a numerical approximation of $\psi^{-1}(x)$,
for example by Newton's method.

Examples, with `T`, `X`, `psi` and `K` of the affinely mapped triangle, and
`f = X[0] * X[1]`:

```python
third = Integer(1) / Integer(3)

# ParameterisedDomain: the centroid of T, as a point in T.
c = parameterised_domain(RD(0), tensor((third, third)), T)

# Integral: the value of f at the centroid, with no support and value shape ().
# A later pass evaluates the integral to 1/9.
fc = integral(pull_back(f, c), RD(0))

# ParameterisedDomain: the centroid of K, as a point in RD(2). A later pass
# evaluates pull_back(psi, c) to (2, 4/3).
C = parameterised_domain(RD(0), pull_back(psi, c), RD(2))
```

- With a second point
  `d = parameterised_domain(RD(0), tensor((third, Integer(2) / Integer(3))), T)`,
  `integral(pull_back(f, c), RD(0)) + integral(pull_back(f, d), RD(0))` is the action of
  a sum of two Dirac measures on $f$.
- `parameterised_domain(RD(0), tensor((third, third)), RD(2))` is the point
  $(1/3, 1/3)$ in the plane. With that point as `q`, `pull_back(f, q)` raises
  `ValueError`, because `q` lies in `RD(2)`, not in `T`.
- `parameterised_domain(RD(0), X, T)` raises `ValueError`, because `X` has the
  support `T`, not `RD(0)`.

## Topology of domains

### Cells

Requires: [Domains](#domains).

A domain can be divided into finitely many cells. A cell is a piece of the
domain with the same dimension as the domain. Two cells meet only on the
boundaries of the two cells.

### Facets

Requires: [Cells](#cells).

A facet is a piece of the boundary of a cell, one dimension lower than the
cell. A facet in the boundary of exactly one cell is an exterior facet. A facet shared
by two cells is an interior facet.

Cells and facets are topological. Defining cells and facets does not require a
metric.

## Function spaces

## Functions

## Variational forms

Requires: [Integral](#integral), [Degree](#degree), [Functions](#functions).

A variational form is an `Integral`, or a combination of `Integral`s. The arity
of a variational form is the number of test and trial functions in the
integrand, and the integrand depends linearly on each of them. A functional has
the arity 0, a linear form the arity 1 and a bilinear form the arity 2. The
integrand of every variational form is a top form, so the degree of a
differential form and the arity of a variational form are independent.

## Core transformations

Requires: [Passes](#passes), [PullBack](#pullback),
[ExteriorDerivative](#exteriorderivative), [Components](#components).

The core has three transformations, each a pass (see
[Passes](#passes)) that
evaluates one lazy node class:

| Core transformation                             | Evaluates            |
| ----------------------------------------------- | -------------------- |
| [`expand_pull_backs(f)`](#expanding-pull-backs) | `PullBack`           |
| [`apply_derivatives(f)`](#applying-derivatives) | `ExteriorDerivative` |
| [`expand_components(f)`](#expanding-components) | `Components`         |

A rule of a core transformation rewrites a node of the evaluated class to a
mathematically equal expression of node classes of the core, and never
simplifies. The rule dispatches on the class of the operand and moves the
evaluated node towards the leaves. The pass rewrites inner nodes first, until
no rule applies. A node without a rule is rebuilt on the rewritten edges.

Each core transformation rewrites only the evaluated class, so the three run
separately or together. Evaluation runs the three together until no rule
applies. A `PullBack` moves inside an `ExteriorDerivative`, and never the
reverse, so evaluation terminates. A downstream lazy node class may join
evaluation if the rules of the class rewrite only that class, move the node
towards the leaves, and never reverse a core rule, as moving an
`ExteriorDerivative` inside a `PullBack` would.

A terminal is a `Coordinate`, or an expression with a support whose class is a leaf
class of a downstream library. After evaluation, every `PullBack`,
`ExteriorDerivative` and `Components` applies to an expression that contains a
downstream terminal, except the coordinate frame $\mathrm{d}X$. An expression
without a downstream terminal evaluates to a sum of products and tensor
products of 0-forms and the coordinate 1-forms $\mathrm{d}X_i$.

### Expanding pull-backs

Requires: [Core transformations](#core-transformations), [PullBack](#pullback),
[Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative).

For a domain $D$ with the parameterisation $\psi$, the parameter domain $P$
and the ambient coordinate domain $Q$, $\psi^*$ commutes with every node class
of the core except `Components` and `PullBack`. Two rules are specific:

- $\psi^* X = \psi$ for `Coordinate(Q)`, and an expression with no support is
  unchanged.
- $\psi^* \mathrm{d}f = \mathrm{d}(\psi^* f)$, with the exterior derivative on
  $P$.

### Applying derivatives

Requires: [Core transformations](#core-transformations),
[ExteriorDerivative](#exteriorderivative), [Product](#product),
[TensorProduct](#tensorproduct), [Zero](#zero), [Conjugate](#conjugate).

The exterior derivative is linear for the linear node classes and follows the
chain rule for the node classes of the degree `()`. Two rules are specific:

- `Product` follows the graded Leibniz rule,
  $\mathrm{d}(\alpha \wedge \beta) = \mathrm{d}\alpha \wedge \beta + (-1)^k \alpha \wedge \mathrm{d}\beta$
  for a $k$-form $\alpha$. A `TensorProduct` under an exterior derivative has at
  most one slot, so the same rule applies.
- $\mathrm{d}\mathrm{d}\omega = 0$, and the exterior derivative of an
  expression with no support is `Zero`, with the degree raised by one.

> Open: the rule for `Conjugate` is the conjugate of the exterior derivative,
> which needs `Conjugate` of a differential form.

### Expanding components

Requires: [Core transformations](#core-transformations),
[Components](#components), [Product](#product), [TensorProduct](#tensorproduct),
[Zero](#zero).

The interior product $\iota_j$ with $\partial_j$ is linear for the linear node
classes, and the new axis stays the last axis of the value shape. Two rules are
specific:

- `Product` follows the graded rule,
  $\iota_j(\alpha \wedge \beta) = \iota_j\alpha \wedge \beta + (-1)^k \alpha \wedge \iota_j\beta$
  for an $\alpha$ whose first slot has the degree $k$, and
  $\iota_j(f \wedge \beta) = f \, \iota_j\beta$ for a 0-form $f$.
  `TensorProduct` acts on the first slot only, so
  $\iota_j(\alpha \otimes \beta) = (\iota_j\alpha) \otimes \beta$ if
  $\alpha$ has a slot, and $\iota_j(f \otimes \beta) = f \, \iota_j\beta$
  for a 0-form $f$.
- $\iota_j \mathrm{d}X_i = \delta_{ij}$, and the components of a `Zero` are a
  `Zero` with the new value shape and degree.

## Extending the core

Requires: [Passes](#passes), [Scalar type](#scalar-type), [Absolute](#absolute).

Extensibility is a requirement of the core. Downstream libraries must be able
to extend UFLx without changing the core.

> Legacy note: In legacy UFL, properties are stored on the node classes, and
> the set of properties is fixed. Every new property touches every node class,
> so many language features were hacked in, for example complex numbers and
> dual spaces.

The core has a small, fixed set of node classes and an open set of passes. A
downstream library extends UFLx mainly by adding passes, not by adding
properties to the node classes of the core. A new pass needs one rule for each
node class of the core and does not change the core. Only a language extension
adds node classes.

The core supports five kinds of extension:

> LLM generated table, could be made more concise.

| Kind of extension                    | Examples                                                  | Section                                     |
| ------------------------------------ | --------------------------------------------------------- | ------------------------------------------- |
| New questions asked of an expression | Physical units, polynomial degree, linearity checks       | [New questions](#new-questions)             |
| New subclasses                       | Elements and functions from Basix, new coordinate domains | [New subclasses](#new-subclasses)           |
| New operations                       | Trace, divergence                                         | [New operations](#new-operations)           |
| New transformations                  | Gateaux derivative                                        | [New transformations](#new-transformations) |
| Language extensions                  | External operators                                        | [Language extensions](#language-extensions) |

A property that the rules of the core depend on is part of the core, not an
extension. For example, the scalar type is part of the core, because
`Absolute` is real for a complex operand. Adding such a property changes the
core.

### Language extensions

Requires: [Leaves and subclassing](#leaves-and-subclassing), [Passes](#passes),
[Core transformations](#core-transformations), [Transpose and matrix
product](#transpose-and-matrix-product).

A language extension adds node classes that are not leaves, in a package
outside the core. Each node class of a language extension provides, as part
of the class:

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

UFLx core is likely to include its own language extensions, such as UFL
compatibility layers and a simpler layer based on classical notions of
differential geometry compatible with most engineering and physics syllabus.

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
[Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative),
[Components](#components).

A downstream library asks a new question of an expression with a pass. The pass
has one rule for each node class of the core and one rule for each leaf class
of the library. The pass needs no rule for a node class of a language extension
that has a lowering (see [Language extensions](#language-extensions)).

For example, a library that checks physical units, like the work of Habera et
al.:

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

The core never knows about units. Objects of the core are immutable and are
compared by content, so the leaves are reliable keys of the table.

An error reported by `units` does not make an object invalid. An object is
valid if the constructors of the core accept the object.

The pass never guesses a rule. A `TypeError` from `units` means that the units
library does not know a node class. An error of the units library means that the
expression has inconsistent units.

### New operations

Requires: [Node class or free function](#node-class-or-free-function), [Tensor
algebra with EllipsisIndex](#tensor-algebra-with-ellipsisindex), [Language
extensions](#language-extensions).

A downstream library builds a new operation from existing node classes, as a
free function. As in the core, the free function may read the static type of
the operands, such as a value shape, to choose the composition (see
[Node class or free function](#node-class-or-free-function)). For
example, `trace(A)` reads the extent $n$ of `A` and returns the sum of the
diagonal entries `A[0, 0] + ... + A[n - 1, n - 1]`. An operation that a
`EllipsisIndex` expresses, such as `scale(a, T)`, needs no value shape (see
[EllipsisIndex](#tensor-algebra-with-ellipsisindex)).

A library node class that a pass replaces by node classes of the core gives
nothing that the free function does not give, unless a core transformation
must keep the structure of the operation. Such an operation, and an operation
that is not a composition of node classes of the core, needs a language
extension.

### New subclasses

Requires: [Leaves and subclassing](#leaves-and-subclassing), [New
operations](#new-operations), [Language extensions](#language-extensions),
[Coordinate](#coordinate), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

Every abstract class is a contract. The contract states what a subclass
provides, so that every core transformation handles the subclass through the
abstract class alone. A downstream library may subclass an
abstract class by fulfilling the contract, without changing the core. Only
leaves are subclassed (see
[Leaves and subclassing](#leaves-and-subclassing)). An
operation built from other expressions is not a subclass but a free function
(see [New operations](#new-operations)) or a node class of a language extension
(see [Language extensions](#language-extensions)). A contract is
part of the class, never a global registry. For example, a
downstream subclass of `AbstractCoordinateDomain`, such as `UnitSimplex` in the
[affinely mapped triangle](#example-an-affinely-mapped-triangle), provides the
dimension, and every core transformation handles the coordinate function of
the subclass through `Coordinate`.

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
a Gateaux derivative rewrites the functions that the derivative is taken with
respect to, so the Gateaux derivative needs a rule for the leaf classes of
those functions. A new transformation must produce only node classes of the
core. A transformation that needs a new node class needs a language extension.

There should be helpers to check whether a transformation covers all classes
of the core.

## Appendix: Class hierarchy

Every object in the graph extends `Base`. The classes defined in this document
form the following hierarchy.

> LLM generated graphs.

```
# Coordinate domains
Base
└── AbstractCoordinateDomain              parameter space, stores n
    └── RD                                R^n

# Expressions and static type
Base
├── AbstractExpression
├── Index                                 named by a label, not an expression
└── EllipsisIndex                         covers any number of axes

# Algebra
AbstractExpression
├── Integer                               literal
├── Float                                 literal
├── ImaginaryUnit                         literal
├── Zero                                  zero of any shape and degree, no support
├── Sum                                   +
├── Product                               *, the wedge product
├── Power                                 **
├── Division                              /
├── Absolute                              abs
├── Conjugate
├── TensorProduct                         ⊗, concatenates slots
├── SquareRoot
├── Tensor                                stacks entries along a new axis
├── Indexed                               full indexing, int or Index entries
├── ComponentTensor                       binds indices to tensor axes
├── IndexSum                              sums over bound indices
├── Determinant
└── Inverse

# Calculus on a coordinate domain
AbstractExpression
├── Coordinate                            terminal, X on P
├── ExteriorDerivative                    d with respect to P
├── Components                            first slot to an axis, with the frame of P
└── Integral                              integral of a top form over P, no support

# Maps between domains
Base
└── ParameterisedDomain                   stores P, ψ and Q
AbstractExpression
└── PullBack                              pull-back of a differential form by a domain, support P of the domain
```

## Appendix: Free functions

> LLM generated table.

| Group                           | Free function                                                   | `AbstractExpression` overload | Takes                                                       | Returns                             |
| ------------------------------- | --------------------------------------------------------------- | ----------------------------- | ----------------------------------------------------------- | ----------------------------------- |
| Expressions and static type     | [`index(label)`](#index)                                        |                               | `label: str \| tuple`                                       | `Index`                             |
|                                 | [`ellipsis_index(label)`](#ellipsisindex)                       |                               | `label: str \| tuple`                                       | `EllipsisIndex`                     |
| Algebra                         | [`add(f, g)`](#sum)                                             | `+`                           | `f, g: AbstractExpression`                                  | `Sum`                               |
|                                 | [`multiply(f, g)`](#product)                                    | `*`                           | `f, g: AbstractExpression`                                  | `Product`                           |
|                                 | [`power(f, g)`](#power)                                         | `**`                          | `f, g: AbstractExpression`                                  | `Power`                             |
|                                 | [`negate(f)`](#negation-and-subtraction)                        | unary `-`                     | `f: AbstractExpression`                                     | `Product` of `Integer(-1)` and `f`  |
|                                 | [`subtract(f, g)`](#negation-and-subtraction)                   | `-`                           | `f, g: AbstractExpression`                                  | `Sum` of `f` and `negate(g)`        |
|                                 | [`divide(f, g)`](#division)                                     | `/`                           | `f, g: AbstractExpression`                                  | `Division`                          |
|                                 | [`absolute(f)`](#absolute)                                      | `abs`                         | `f: AbstractExpression`                                     | `Absolute`                          |
|                                 | [`conjugate(f)`](#conjugate)                                    |                               | `f: AbstractExpression`                                     | `Conjugate`                         |
|                                 | [`tensor_product(f, g)`](#tensorproduct)                        |                               | `f, g: AbstractExpression`                                  | `TensorProduct`                     |
|                                 | [`square_root(f)`](#elementary-functions-very-incomplete)       |                               | `f: AbstractExpression`                                     | `SquareRoot`                        |
|                                 | [`tensor(entries)`](#tensor)                                    |                               | `entries: tuple[AbstractExpression, ...]`                   | `Tensor`                            |
|                                 | [`indexed(f, indices)`](#indexed)                               | `f[i]`                        | `f: AbstractExpression`, `indices: tuple`                   | `Indexed`                           |
|                                 | [`component_tensor(f, indices)`](#componenttensor-and-indexsum) |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`       | `ComponentTensor`                   |
|                                 | [`index_sum(f, indices)`](#componenttensor-and-indexsum)        |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`       | `IndexSum`                          |
|                                 | [`determinant(f)`](#determinant)                                |                               | `f: AbstractExpression`                                     | `Determinant`                       |
|                                 | [`inverse(f)`](#inverse)                                        |                               | `f: AbstractExpression`                                     | `Inverse`                           |
|                                 | [`transpose(f)`](#transpose-and-matrix-product)                 |                               | `f: AbstractExpression`                                     | `ComponentTensor` of an `Indexed`   |
|                                 | [`matrix_product(f, g)`](#transpose-and-matrix-product)         | `@`                           | `f, g: AbstractExpression`                                  | `ComponentTensor` of an `IndexSum`  |
| Calculus on a coordinate domain | [`coordinate(P)`](#coordinate)                                  |                               | `P: AbstractCoordinateDomain`                               | `Coordinate`                        |
|                                 | [`exterior_derivative(f, P)`](#exteriorderivative)              |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain`      | `ExteriorDerivative`                |
|                                 | [`components(f, P)`](#components)                               |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain`      | `Components`                        |
|                                 | [`volume_form(P)`](#volume-form)                                |                               | `P: AbstractCoordinateDomain`                               | `Product` of the coordinate 1-forms |
|                                 | [`integral(f, P)`](#integral)                                   |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain`      | `Integral`                          |
| Maps between domains            | [`parameterised_domain(P, psi, Q)`](#parameteriseddomain)       |                               | `P, Q: AbstractCoordinateDomain`, `psi: AbstractExpression` | `ParameterisedDomain`               |
|                                 | [`pull_back(f, D)`](#pullback)                                  |                               | `f: AbstractExpression`, `D: ParameterisedDomain`           | `PullBack`                          |

The value-shape and other constraints on the arguments are stated in the section
of each free function. A free function and an operator overload accept only the
types in the "Takes" column.

A free function whose "Returns" column names only a node class `C` returns
`C(...)` of the arguments. Any other free function returns the composition in
the "Returns" column. The sections of the node classes do not
repeat this rule.

## Appendix: Properties

A derived property listed for an abstract class is abstract, unless the table
gives the value. Each concrete class implements an abstract property by the
rules in the section of the concrete class.

> LLM generated table.

| Class                                          | Stored properties                                        | Derived properties                                                |
| ---------------------------------------------- | -------------------------------------------------------- | ----------------------------------------------------------------- |
| `AbstractCoordinateDomain`                     | `dimension`, in every concrete subclass such as `RD`     |                                                                   |
| `AbstractExpression`                           |                                                          | `value_shape`, `free_indices`, `scalar_type`, `support`, `degree` |
| `Index`, `EllipsisIndex`                       | `label`                                                  |                                                                   |
| `Integer`, `Float`                             | `value`                                                  |                                                                   |
| `Zero`                                         | `value_shape`, `free_indices`, `degree`                  |                                                                   |
| `Sum`, `Product`, `Power`, `Division`          | `first`, `second`                                        |                                                                   |
| `Absolute`, `Conjugate`                        | `operand`                                                |                                                                   |
| `TensorProduct`                                | `first`, `second`                                        |                                                                   |
| `SquareRoot`                                   | `operand`                                                |                                                                   |
| `Tensor`                                       | `entries`                                                |                                                                   |
| `Indexed`, `ComponentTensor`, `IndexSum`       | `expression`, `indices`                                  |                                                                   |
| `Determinant`, `Inverse`                       | `operand`                                                |                                                                   |
| `Coordinate`                                   | `domain`                                                 |                                                                   |
| `ExteriorDerivative`, `Components`, `Integral` | `expression`, `domain`                                   |                                                                   |
| `ParameterisedDomain`                          | `parameter_domain`, `parameterisation`, `ambient_domain` |                                                                   |
| `PullBack`                                     | `expression`, `domain`                                   |                                                                   |

## Appendix: Degree table

| Class                                        | Edges that are expressions must have                      | Degree of the node                            |
| -------------------------------------------- | --------------------------------------------------------- | --------------------------------------------- |
| `Integer`, `Float`, `ImaginaryUnit`          | No such edges                                             | `()`                                          |
| `Zero`                                       | No such edges                                             | The stored degree                             |
| `Sum`                                        | Equal degrees                                             | The degree of the operands                    |
| `Product`                                    | The degree `()`, or equal numbers of slots                | The slotwise sum (see [Product](#product))    |
| `Power`, `Division`, `Absolute`, `Conjugate` | The degree `()`                                           | `()`                                          |
| `TensorProduct`                              | Any degree                                                | The concatenated degrees                      |
| `SquareRoot`                                 | The degree `()`                                           | `()`                                          |
| `Tensor`                                     | Equal degrees                                             | The degree of the entries                     |
| `Indexed`, `ComponentTensor`, `IndexSum`     | Any degree                                                | The degree of the expression                  |
| `Determinant`, `Inverse`                     | The degree `()`                                           | `()`                                          |
| `Coordinate`                                 | No such edges                                             | `()`                                          |
| `ExteriorDerivative`                         | The degree `()` or `(k,)`                                 | `(1,)` for `()`, and `(k + 1,)` for `(k,)`    |
| `Components`                                 | At least one slot                                         | The first slot lowered by 1, and removed at 0 |
| `Integral`                                   | `(n,)` with $n$ the dimension of $P$, or `()` for $n = 0$ | `()`                                          |
| `PullBack`                                   | Any degree                                                | The degree of the expression                  |

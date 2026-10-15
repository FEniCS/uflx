# Core UFLx language definition

This note describes the smallest core language behind UFLx. It intentionally
avoids implementation details, compiler algorithms, and broad redesign ideas.
It assumes knowledge of differential geometry, and functional analysis. Some
advanced differential geometry concepts are explained in more detail.

This document is aimed at users and developers with extensive knowledge of
domain specific languages for variational forms (e.g. UFL, ngsolve).

## What UFLx is and isn't

In one sentence:

> UFLx is a typed symbolic language for writing variational forms.

The UFLx core language (henceforth the core) should not use finite element or
numerical analysis terminology - differential geometry, and functional analysis,
and their intersection, without finite element (numerical analysis) concepts,
should be all that is required.

Finite element and numerical analysis concepts, when necessary, will be added
on top of the core, or provided by external libraries (Basix being the
'in-house' example). Basix and DOLFINx and their assumptions should not be
visible in the core.

UFLx is *not* a symbolic language for reasoning about problems in functional
analysis or differential geometry. Hence, the core should only contain the
minimal concepts and relations necessary for writing a wide range of
variational forms as code for subsequent algorithmic transformation.

A UFLx program is ordinary Python code that constructs symbolic objects and
connects them through a graph. The primary result(s) are usually Form(s), such
as a bilinear form, linear form, or functional. Python code will be used to
motivate examples but a full concrete implementation is out of scope.

UFLx does not handle concrete data or numerical algorithms; e.g. actual meshes,
boundary markers, numerical values of coefficients, or assembly procedures.
Those belong to the surrounding problem-solving environment, such as DOLFINx,
Firedrake, FIAT and Basix. UFLx instead specifies a set of leaf nodes, and the
external library provides the concrete data behind the leaf nodes.

## Basic principles

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
- The core should not contain any algorithm that is largely focused on
  producing convenient output for form compilers.

## Basic implementation principles

- Explicit is better than implicit.

- Boilerplate and error-prone patterns should be minimised.

- UFLx should be deterministic on the same machine using the same byte-code.
  Serialisation/pickling should be considered early on.

- Performance is *important* and so any data structures and algorithms in the
  core should be carefully chosen and tested for scalability.

- The minimal language should only need a relatively small number of object
  kinds. When designing the language, consider whether objects can be composed,
  rather than extended. Consider what could be layered on top of the core.

- Only leaves are subclassed. A leaf is an object that the core cannot build
  from other objects, such as the coordinate domain `RD(n)` or a function whose
  values an external library supplies. An object built from other objects is a
  single class of the core or a free function that returns objects of the core,
  never a subclass. For example, every domain is a `ParameterisedDomain`, and an
  affinely mapped triangle is a free function that returns a
  `ParameterisedDomain`. A language extension is the one exception (see
  [Language extensions](#language-extensions)).

- A property reads and a free function builds. Data that an object stores, or
  derives from the class and the construction arguments of the object, is a
  property of the object, such as `value_shape` of an expression or
  `parameter_domain` of a domain. A property never constructs an object. Any
  operation that constructs objects is a pure free function, such as
  `coordinate(P)`, `jacobian(D)` or `metric(D)`. A class method, where
  needed, immediately calls a free function.

- A parent class never constructs a subclass. A property of a parent class may
  return an object whose class is a subclass, because the property only
  returns an object that already exists. A subclass may narrow the return
  type of a property.

- A concept must be a free function, not a node class, if the free function can
  build the concept as one fixed composition of existing node classes applied to
  the arguments, *without inspecting the arguments to choose the composition*.
  Otherwise the concept gets its own node class. For example, the metric is
  always `transpose(J) @ J`, so the metric is a free function. No composition of
  existing nodes differentiates, so `Gradient` is currently a node class. A free
  function that builds a fixed composition is a definition, not a
  simplification. The free function does not inspect the operands, so the free
  function never removes structure that the user built. This rule binds the
  free functions of the core. A layer above the core may read value shapes to
  choose what to build, if the layer builds only objects of the core and
  annotates the result. A `EllipsisIndex` makes most tensor algebra a fixed
  composition (see [EllipsisIndex](#ellipsisindex)), so the layer rarely needs to
  read a value shape.

- Simplifications and rewriting should be completely avoided during user graph
  construction, maintaining the user's intent. Only later selectable algorithmic
  passes should re-write. A good example is `(a + b) + c` in floating point
  arithmetic where `a`, `b` and `c` are unknown terminals (i.e. to be provided
  later). The user may have very good reason to have `(a + b)` performed before
  adding the answer to `c`.

- A pass is a pure free function that walks a graph after construction of the
  graph. The result of a pass depends only on the arguments of the pass. A pass
  returns a value, such as a unit, or a new graph. A pass builds a new graph by
  constructing an object for each node that the pass replaces, and by
  reconstructing every object that depends on a replaced node, up to the root
  (see [The base object](#the-base-object-base)). A pass may cache results for
  each object, because objects are immutable.

  A pass does one thing, such as differentiating or lowering. Passes compose,
  and the user of the core chooses which passes run and in what order.

Three kinds of code produce objects:

| Kind                                               | When the code runs              | Inspects arguments to choose what to build            | Validates                     |
| -------------------------------------------------- | ------------------------------- | ----------------------------------------------------- | ----------------------------- |
| Constructor                                        | Whenever an object is made      | No, a constructor inspects arguments only to validate | Yes, and only constructors do |
| Construction free function, such as `transpose`    | While the user builds the graph | No, the free function builds one fixed composition    | No                            |
| Pass, such as differentiation, lowering or `units` | After construction              | Yes, and a pass may rewrite the graph                 | No                            |

A construction free function fails when a constructor rejects the
arguments, for example `transpose` of a vector raises `IndexError` from the
constructor of `Indexed`. A pass fails with `TypeError` for a node class that
the pass has no rule for, or with an error of the pass, such as inconsistent
units.

## The minimal object model

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

### Enforcement

Enforcement must be implemented to help library developers - look at tools such
as `dataclass`, deterministic hashing, hash-consing, `WeakValueDictionary`,
`ruff` and `semgrep`.

https://github.com/semgrep/semgrep

This rules should be put in place at the start of the project.

It is not clear at the moment how best to do this - it may be that custom
static analysis e.g. (semgrep) is more robust than hacking Python to be 'well
behaved' in some places, for example checking the purity of functions is not
possible in Python.

### Class hierarchy

Every object in the graph extends `Base`. The classes defined in this document
form the following hierarchy.

> LLM generated graphs.

```
Base
├── AbstractExpression
│   ├── Integer                           literal
│   ├── Float                             literal
│   ├── ImaginaryUnit                     literal
│   ├── Zero                              zero of any shape, no set
│   ├── Coordinate                        terminal, X on P
│   ├── Tensor                            stacks entries along a new axis
│   ├── Indexed                           full indexing, int or Index entries
│   ├── ComponentTensor                   binds indices to tensor axes
│   ├── IndexSum                          sums over bound indices
│   ├── AnnotatedExpression               identity
│   ├── Gradient                          derivative with respect to X of P
│   ├── Composition                       f composed with a domain, set P of the domain
│   ├── LebesgueIntegral                  Lebesgue integral over P, no set
│   ├── Sum                               +
│   ├── Product                           *
│   ├── Division                          /
│   ├── Power                             **
│   ├── Absolute                          abs
│   ├── Conjugate
│   ├── SquareRoot
│   └── Determinant
│
├── Index                                 named by a label, not an expression
├── EllipsisIndex                         covers any number of axes
│
├── AbstractCoordinateDomain              parameter space, stores n
│   └── RD                                R^n
│
└── ParameterisedDomain                   stores P, ψ and Q
```

### Free functions

Free functions are not objects in the graph. A free function builds or returns
objects in the graph. Operator overloads call the free functions.

> LLM generated table.

| Group                 | Free function                                                   | `AbstractExpression` overload | Takes                                                       | Returns                                                             |
| --------------------- | --------------------------------------------------------------- | ----------------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------- |
| Arithmetic            | [`add(f, g)`](#arithmetic)                                      | `+`                           | `f, g: AbstractExpression`                                  | `Sum`                                                               |
|                       | [`subtract(f, g)`](#arithmetic)                                 | `-`                           | `f, g: AbstractExpression`                                  | `AnnotatedExpression`: `Sum` of `f` and `negate(g)`                 |
|                       | [`multiply(f, g)`](#arithmetic)                                 | `*`                           | `f, g: AbstractExpression`                                  | `Product`                                                           |
|                       | [`divide(f, g)`](#arithmetic)                                   | `/`                           | `f, g: AbstractExpression`                                  | `Division`                                                          |
|                       | [`negate(f)`](#arithmetic)                                      | unary `-`                     | `f: AbstractExpression`                                     | `AnnotatedExpression`: `Product` of `Integer(-1)` and `f`           |
|                       | [`absolute(f)`](#arithmetic)                                    | `abs`                         | `f: AbstractExpression`                                     | `Absolute`                                                          |
|                       | [`conjugate(f)`](#arithmetic)                                   |                               | `f: AbstractExpression`                                     | `Conjugate`                                                         |
|                       | [`power(f, g)`](#arithmetic)                                    | `**`                          | `f, g: AbstractExpression`                                  | `Power`                                                             |
| Elementary functions  | [`square_root(f)`](#elementary-functions-very-incomplete)       |                               | `f: AbstractExpression`                                     | `SquareRoot`                                                        |
| Coordinate domains    | [`coordinate(P)`](#coordinate)                                  |                               | `P: AbstractCoordinateDomain`                               | `Coordinate`                                                        |
| Tensors and indexing  | [`tensor(entries)`](#tensor)                                    |                               | `entries: tuple[AbstractExpression, ...]`                   | `Tensor`                                                            |
|                       | [`indexed(f, indices)`](#indexed)                               | `f[i]`                        | `f: AbstractExpression`, `indices: tuple`                   | `Indexed`                                                           |
| Index notation        | [`index(label)`](#index)                                        |                               | `label: str \| tuple`                                       | `Index`                                                             |
|                       | [`ellipsis_index(label)`](#ellipsisindex)                       |                               | `label: str \| tuple`                                       | `EllipsisIndex`                                                     |
|                       | [`component_tensor(f, indices)`](#componenttensor-and-indexsum) |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`       | `ComponentTensor`                                                   |
|                       | [`index_sum(f, indices)`](#componenttensor-and-indexsum)        |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`       | `IndexSum`                                                          |
| Annotation            | [`annotate(f, annotation)`](#annotation)                        |                               | `f: AbstractExpression`, `annotation: tuple`                | `AnnotatedExpression`                                               |
| Tensor algebra        | [`determinant(f)`](#determinant)                                |                               | `f: AbstractExpression`                                     | `Determinant`                                                       |
|                       | [`transpose(f)`](#transpose-and-matrix-product)                 |                               | `f: AbstractExpression`                                     | `AnnotatedExpression`: `ComponentTensor`                            |
|                       | [`matrix_product(f, g)`](#transpose-and-matrix-product)         | `@`                           | `f, g: AbstractExpression`                                  | `AnnotatedExpression`: `ComponentTensor` of an `IndexSum`           |
| Differentiation       | [`gradient(f, P)`](#gradient)                                   |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain`      | `Gradient`                                                          |
|                       | [`identity(P)`](#identity)                                      |                               | `P: AbstractCoordinateDomain`                               | `AnnotatedExpression`: `Gradient` of `coordinate(P)`                |
| Parameterised domains | [`parameterised_domain(P, psi, Q)`](#parameteriseddomain)       |                               | `P, Q: AbstractCoordinateDomain`, `psi: AbstractExpression` | `ParameterisedDomain`                                               |
| Geometry              | [`jacobian(D)`](#geometry)                                      |                               | `D: ParameterisedDomain`                                    | `AnnotatedExpression`: `Gradient` of `D.parameterisation`           |
|                       | [`metric(D)`](#geometry)                                        |                               | `D: ParameterisedDomain`                                    | `AnnotatedExpression`: `transpose(J) @ J`                           |
|                       | [`measure_density(D)`](#geometry)                               |                               | `D: ParameterisedDomain`                                    | `AnnotatedExpression`: `SquareRoot` of the determinant of $g$       |
| Composition           | [`compose(f, D)`](#composition)                                 |                               | `f: AbstractExpression`, `D: ParameterisedDomain`           | `Composition`                                                       |
| Integration           | [`lebesgue_integral(f, P)`](#lebesgueintegral)                  |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain`      | `LebesgueIntegral`                                                  |
|                       | [`integrate(f, D)`](#integrate)                                 |                               | `f: AbstractExpression`, `D: ParameterisedDomain`           | `AnnotatedExpression`: `LebesgueIntegral` over `D.parameter_domain` |

Every overload in the table is a method of `AbstractExpression` (see
[`AbstractExpression`](#abstractexpression)). No other class defines
arithmetic, matrix or indexing overloads. Overloading should be considered
carefully outside of standard mathematical meanings.

The value-shape and other constraints on the arguments are stated in the section
of each free function. A free function and an operator overload accept only the
types in the "Takes" column.

A free function whose "Returns" column names only a node class `C` returns
`C(...)` of the arguments, unchanged. The sections of the node classes do not
repeat this rule.

### Properties

A stored property returns a construction argument. A derived property returns
data derived from the class and the construction arguments. A concrete class
has the derived properties of the abstract classes of the concrete class.

A derived property listed for an abstract class is abstract, unless the table
gives the value. Each concrete class implements an abstract property by the
rules in the section of the concrete class.

> LLM generated table.

| Class                                                | Stored properties                                        | Derived properties                                  |
| ---------------------------------------------------- | -------------------------------------------------------- | --------------------------------------------------- |
| `AbstractExpression`                                 |                                                          | `value_shape`, `free_indices`, `scalar_type`, `set` |
| `Integer`, `Float`                                   | `value`                                                  |                                                     |
| `Zero`                                               | `value_shape`, `free_indices`                            |                                                     |
| `Coordinate`                                         | `domain`                                                 |                                                     |
| `Tensor`                                             | `entries`                                                |                                                     |
| `Indexed`, `ComponentTensor`, `IndexSum`             | `expression`, `indices`                                  |                                                     |
| `AnnotatedExpression`                                | `expression`, `annotation`                               |                                                     |
| `Gradient`, `LebesgueIntegral`, `Composition`        | `expression`, `domain`                                   |                                                     |
| `Absolute`, `Conjugate`, `SquareRoot`, `Determinant` | `operand`                                                |                                                     |
| `Sum`, `Product`, `Division`, `Power`                | `first`, `second`                                        |                                                     |
| `Index`, `EllipsisIndex`                             | `label`                                                  |                                                     |
| `AbstractCoordinateDomain`                           | `dimension`, in every concrete subclass such as `RD`     |                                                     |
| `ParameterisedDomain`                                | `parameter_domain`, `parameterisation`, `ambient_domain` |                                                     |

## Extending the core

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

| Kind of extension             | Examples                                                  | Section                                     |
| ----------------------------- | --------------------------------------------------------- | ------------------------------------------- |
| New questions asked of a form | Physical units, polynomial degree, linearity checks       | [New questions](#new-questions)             |
| New subclasses                | Elements and functions from Basix, new coordinate domains | [New subclasses](#new-subclasses)           |
| New operations                | Trace, divergence                                         | [New operations](#new-operations)           |
| New transformations           | Gateaux derivative                                        | [New transformations](#new-transformations) |
| Language extensions           | External operators                                        | [Language extensions](#language-extensions) |

A property that the rules of the core depend on is part of the core, not an
extension. For example, the scalar type is part of the core, because
`Absolute` is real for a complex operand. Adding such a property changes the
core.

### New questions

A downstream library asks a new question of a form with a pass. The pass has
one rule for each node class of the core and one rule for each leaf class of
the library. The pass needs no rule for a node class of a language extension
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
   `X[0] * X[1]` is in square metres. `gradient(X[0], RD(2))` is in metres
   divided by metres, so `gradient(X[0], RD(2))` is dimensionless.
   `X[0] + X[0] * X[1]` is an error of the units library.
3. Calls `units` on a form whenever the library needs the check.

The core never knows about units. Objects of the core are immutable and are
compared by content, so the leaves are reliable keys of the table.

An error reported by `units` does not make an object invalid. An object is
valid if the constructors of the core accept the object.

The pass never guesses a rule. A `TypeError` from `units` means that the units
library does not know a node class. An error of the units library means that the
expression has inconsistent units.

### New subclasses

Every abstract class is a contract. The contract states what a subclass
provides, so that every core transformation handles the subclass through the
abstract class alone. A downstream library may subclass an
abstract class by fulfilling the contract, without changing the core. Only
leaves are subclassed (see
[Basic implementation principles](#basic-implementation-principles)). An
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
  changes the meaning of a form, so the space is part of the leaf.
- Information that only checks or analyses the variational problem may stay in
  a table of the library, keyed by leaves. For example, units check a form but
  do not change the meaning of the form.

Information never goes in an annotation, because any pass may drop an
annotation.

### New operations

A downstream library builds a new operation from existing node classes, as a
free function. A free function above the core may read value shapes to choose
the composition (see
[Basic implementation principles](#basic-implementation-principles)), and
annotates the result. For example, `trace(A)` reads the extent $n$ of `A` and
returns the sum of the diagonal entries `A[0, 0] + ... + A[n - 1, n - 1]`,
annotated with `("library.trace", A)`. An operation that a `EllipsisIndex`
expresses, such as `scale(a, T)`, needs no value shape (see [EllipsisIndex](#ellipsisindex)).

A library node class that a pass replaces by node classes of the core gives
nothing that the free function does not give, unless a core transformation
must keep the structure of the operation. Such an operation, and an operation
that is not a composition of node classes of the core, needs a language
extension.

### New transformations

A new transformation is a pass with one rule for each node class of the coreV
and one rule for each leaf class that the transformation rewrites. For example,
a Gateaux derivative rewrites the functions that the derivative is taken with
respect to, so the Gateaux derivative needs a rule for the leaf classes of
those functions. A new transformation must produce only node classes of the
core. A transformation that needs a new node class needs a language extension.

There should be helpers to check whether a transformation covers all classes
of the core.

### Language extensions

A language extension adds node classes that are not leaves, in a package
outside the core. Each node class of a language extension provides, as part
of the class:

- a rule for each [core transformation](#core-transformations), which
  produces only node classes of the core or of the language extension, and
- a lowering to an equal expression of node classes of the core, if such an
  expression exists.

A pass without a rule for a node class of a language extension lowers the node
and continues, so a pass such as `units` needs no rule for the node class. For
a node class without a lowering, every pass without a rule raises
`TypeError`. A rule must agree with the lowering followed by the rules for the
lowered expression, up to mathematical equality.

For example, a language extension for external operators adds a node class
for an operation known only by name, such as `exp`, a material law or a neural
network. The rule for `apply_derivatives` differentiates the node class. The
rule for `expand_compositions` composes the operands, because the operation is
pointwise. The node class has no lowering.

## The lowest objects

### Expressions

An expression denotes an array-valued function defined on the set of the
expression. An expression with no set denotes a single array.

#### `AbstractExpression`

Every `AbstractExpression` derives the value shape, the free indices, the
scalar type and the set of the expression.

- The set is an `AbstractCoordinateDomain`. The set is absent for an expression
  with no set, such as a literal. Unless a node class states otherwise, the set
  of a node is the set of the edges of the node that are expressions with a set.
  All such edges must have the same set, or construction raises `ValueError`. A
  node with no such edge has no set.
- The scalar type is real or complex. Unless a node class states otherwise, an
  expression is complex if and only if at least one edge is a complex
  expression. `ImaginaryUnit` is complex.

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
# expressions.py
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
    def set(self) -> AbstractCoordinateDomain | None:
        """The set of the expression, or None for an expression with no set."""

    # Operator overloading must be in the class.
    def __add__(self, other: AbstractExpression) -> AbstractExpression:
        """Add."""
        return add(self, other)

    def __sub__(self, other: AbstractExpression) -> AbstractExpression:
        """Subtract."""
        # Calls the free function arithmetic.subtract.
        ...

    def __mul__(self, other: AbstractExpression) -> AbstractExpression:
        """Multiply."""
        # Calls the free function arithmetic.multiply.
        ...

    def __matmul__(self, other: AbstractExpression) -> AbstractExpression:
        """Matrix multiply."""
        # Calls the free function tensor_algebra.matrix_product.
        ...

    def __truediv__(self, other: AbstractExpression) -> AbstractExpression:
        """Divide."""
        # Calls the free function arithmetic.divide.
        ...

    def __neg__(self) -> AbstractExpression:
        """Negate."""
        # Calls the free function arithmetic.negate.
        ...

    def __abs__(self) -> AbstractExpression:
        """Absolute value."""
        # Calls the free function arithmetic.absolute.
        ...

    def __pow__(self, exponent: AbstractExpression) -> AbstractExpression:
        """Raise to a power."""
        # Calls the free function arithmetic.power.
        ...

    def __getitem__(self, indices: object) -> AbstractExpression:
        """Index."""
        # Converts a single entry into a tuple of one entry, then calls indexed.
        ...

    def __iter__(self) -> NoReturn:
        """Raise TypeError, because expressions are not iterable."""
        # Without __iter__, Python iterates through __getitem__.
        ...

# arithmetic.py

# A free function
def add(
    first: AbstractExpression, second: AbstractExpression
) -> AbstractExpression:
    """Add two expressions.

    Args:
        first: The first summand.
        second: The second summand.

    Returns:
        The sum of the two expressions, with the summands in the given order.

    Raises:
        TypeError: From the constructor of Sum, if a summand is not an
            expression.
        ValueError: From the constructor of Sum, if the value shapes or the free
            indices of the summands differ.
    """
    # The constructor of Sum validates.
    return Sum(first, second)
```

### Literals

A literal is a number written down in an expression, e.g. `Integer(2) * u`. A
literal's value is part of the expression itself. Every literal has value shape
`()`, no free indices and no set.

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

`Zero(value_shape, free_indices)` is the zero with the given value shape and
free indices, and has no set. `apply_derivatives` produces `Zero` (see
[Applying derivatives](#applying-derivatives)). No fixed
composition gives a zero with free indices, because `Integer(0) * a[i]` has
edges to `a`, so `Zero` is a node class.

#### Numbers in expressions

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

### Arithmetic

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

The exponent of `Power` has no set, or construction raises `ValueError`. The
derivative of a power whose exponent has a set needs a logarithm, which the
core does not have.

`Absolute` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `absolute(f)`, with the overload `abs(f)`, takes
a scalar and returns a scalar. For a complex argument, `Absolute` is the
modulus. `Absolute` is real for every argument.

`Conjugate` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `conjugate(f)` takes a scalar and returns the
complex conjugate of the scalar. No fixed composition of the other node classes
conjugates, so `Conjugate` is a node class. A sesquilinear form, such as the
$L^2$ inner product $\int u \bar{v}$, needs `Conjugate`.

Negation and subtraction are free functions, not node classes:

- `negate(f)`, with the overload `-f`, is `multiply(Integer(-1), f)`,
  annotated with `("uflx.negate", f)`.
- `subtract(f, g)`, with the overload `f - g`, is `add(f, negate(g))`,
  annotated with `("uflx.subtract", f, g)`.

Both compositions are exact in floating-point arithmetic, so neither
composition changes the arithmetic that the user wrote. Division is a node
class, because `multiply(f, power(g, Integer(-1)))` rounds twice and differs
from `f / g` in floating-point arithmetic.

### Elementary functions (very incomplete)

`SquareRoot` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `square_root(f)` takes a scalar and returns a
scalar.

`SquareRoot` is a node class, not the composition
`power(f, Integer(1) / Integer(2))`.
IEEE 754 requires a correctly rounded square root but does not require a
correctly rounded power, so the two differ in floating-point arithmetic.

## Differential geometry

We begin with differential geometry because the function spaces and operators
of UFLx are defined on a domain. The function spaces and operators require two
kinds of structure on the domain:

- topological structure, which is the dimension, boundary, orientation and
  decomposition into cells of the domain,
- geometric structure, which is the metric of the domain and the measure
  induced by the metric.

### Domains

A domain $D$ is the image $\psi(P)$ of a coordinate domain $P$ under a
parameterisation $\psi$ into a coordinate domain $Q$. The core has one class of
domain, [`ParameterisedDomain`](#parameteriseddomain), which stores $P$, $\psi$
and $Q$.

A particular kind of domain, such as an affinely mapped triangle or a point, is
a free function that returns a `ParameterisedDomain`, not a subclass. Every
kind of domain is a fixed composition of $P$, $\psi$ and $Q$, so every kind of
domain is a free function (see
[Basic implementation principles](#basic-implementation-principles)).

#### AbstractCoordinateDomain

An abstract base class for coordinate domains `AbstractCoordinateDomain`.

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$ whose
points are named by the coordinates of the points themselves. A coordinate
domain is a parameter space. A coordinate domain is the parameter domain of a
domain, the ambient coordinate domain of a domain, and the set of an
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

#### RD

A concrete subclass of `AbstractCoordinateDomain` for the whole of $R^n$.

`RD(n)` stores the dimension $n$, an `int` with $n \ge 0$. Everything else
comes from `AbstractCoordinateDomain`.

Example:

- `RD(2)` is the parameter domain of the paraboloid example below.

#### Coordinate

A concrete subclass of `AbstractExpression` for the coordinate function of a
coordinate domain.

Every concrete coordinate domain $P$ provides a coordinate function
`Coordinate(P)`. The coordinate function $X$ of $P$ is the terminal expression that maps each point
of $P$ to the $n$-tuple of that point:

$X : P \to R^n$, $X(p) = p$.

`Coordinate(P)` stores the coordinate domain $P$ and derives the value shape `(n,)` and the set $P$.

Only coordinate domains provide a coordinate function. The coordinates of the
points of a domain are given by the parameterisation of the domain (see
[ParameterisedDomain](#parameteriseddomain)), which is built from the coordinate function
of the parameter domain.

Examples:

- `coordinate(RD(3))` is `Coordinate(RD(3))`, with value shape `(3,)` and the
  set `RD(3)`.
- `coordinate(RD(3)) == coordinate(RD(3))` is `True`, and
  `coordinate(RD(2)) == coordinate(RD(3))` is `False`.
- In the paraboloid example, `X = coordinate(RD(2))` has value shape `(2,)`.

## Tensors and indexing

### Tensor

A concrete subclass of `AbstractExpression` that stacks expressions along a
new first axis.

`Tensor(entries)` stores the entries, a non-empty `tuple` of
`AbstractExpression`s that all have the same value shape `s`, or construction
raises `ValueError`. `Tensor(entries)` derives the value shape `(n,) + s`, with
$n$ the number of entries. A matrix is a `Tensor` of `Tensor`s.

A value shape may contain the extent 0, for example the value shape `(0,)` of
the coordinate function of `RD(0)`. `Tensor` cannot build an expression with
the extent 0, because an empty `tuple` has no entry to give `s`.

Examples, with `zero, one, two = Integer(0), Integer(1), Integer(2)`:

- `tensor((one, two))` has value shape `(2,)` and no set.
- `tensor((tensor((one, zero)), tensor((zero, one))))` has value shape
  `(2, 2)` and no set.
- With `r0 = Float(0.0)` and `r1 = Float(1.0)`,
  `tensor((tensor((r1, r0)), tensor((r0, r1))))` has `Float` entries, so the
  matrix is not equal to `tensor((tensor((one, zero)), tensor((zero, one))))`.
- `tensor((tensor((one, two)),))` has value shape `(1, 2)`.
- `tensor((one, tensor((one, two))))` raises `ValueError`, because the value
  shapes of the entries differ.

### Indexed

> LLM generated content.

A concrete subclass of `AbstractExpression` that selects one component of a
tensor-valued expression.

`Indexed(expression, indices)` stores:

- The expression, an `AbstractExpression` of rank $r$.
- The indices, a `tuple` of entries. Each entry is an `int`, an `Index` or a
  `EllipsisIndex` (see [Index notation](#index-notation)). At most one entry is a
  `EllipsisIndex`.

Indexing is full indexing only. Each `int` and each `Index` entry covers one
axis, and a `EllipsisIndex` entry covers every axis that the other entries do not
cover. Without a `EllipsisIndex`, the number of entries must equal $r$. With a
`EllipsisIndex`, the number of other entries must be at most $r$. Each `int` entry
$i_j$ must satisfy $0 \le i_j < n_j$, where `(n_0, ..., n_{r-1})` is the value
shape of the expression. Each `Index` entry $i_j$ takes the extent $n_j$. A
`EllipsisIndex` takes the `tuple` of the extents of the axes that the `EllipsisIndex`
covers.

`Indexed(expression, indices)` derives the value shape `()`, the set of the
expression, and the free indices of the expression followed by the `Index` and
`EllipsisIndex` entries.

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
- `tensor((Integer(1) - X[0], X[0]))` has value shape `(2,)` and the set
  `RD(2)`.
- `tensor((tensor((X[0], X[1])),))` has value shape `(1, 2)`.
- In the paraboloid example, the height `X[0]**Integer(2) + X[1]**Integer(2)`
  is an expression on `RD(2)` with value shape `()`.

## Index notation

> LLM generated content.

Index notation writes tensor algebra as fixed compositions, so that operators
such as `transpose` and `matrix_product` are free functions rather than node
classes.

### Index

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
operands. The label of a derived index is the annotation of the free function
followed by a number (see [Annotation](#annotation)). For example,
`transpose(A)` uses `Index(("uflx.transpose", A, 0))` and
`Index(("uflx.transpose", A, 1))`. An index whose label contains `A` has an
edge to `A`. Every node that uses an `Index` has an edge to the `Index`. So an index
whose label contains `A` cannot occur inside `A`, because `A` would then
contain `A`, and the graph is acyclic. Hence
a derived index is fresh with respect to the operands, and the free function
does not inspect the operands to find a fresh index. A global counter is not
allowed, because a counter makes construction depend on global state.

### EllipsisIndex

> LLM generated content.

A concrete subclass of the base object for an index that covers any number of
axes. A `EllipsisIndex` is not an expression. `EllipsisIndex(label)` stores the label,
with the same rules as the label of an `Index`.

A `EllipsisIndex` covers every axis of an `Indexed` that the other entries do not
cover, so a `EllipsisIndex` may cover no axis. The extent of a `EllipsisIndex` is the
`tuple` of the extents of the covered axes. In every other rule of index
notation, a `EllipsisIndex` behaves as an `Index`.

A `EllipsisIndex` writes an operation on tensors of every rank as one fixed
composition. With derived ellipsis indices `r` and `s`, an ordinary index `k`, and
annotations omitted:

| Operation          | Fixed composition                                              |
| ------------------ | -------------------------------------------------------------- |
| `scale(a, T)`      | `component_tensor(a * T[r], (r,))`                             |
| tensor sum `A + B` | `component_tensor(A[r] + B[r], (r,))`                          |
| `outer(u, v)`      | `component_tensor(u[r] * v[s], (r, s))`                        |
| `dot(u, v)`        | `component_tensor(index_sum(u[r, k] * v[k, s], (k,)), (r, s))` |
| `inner(A, B)`      | `index_sum(A[r] * conjugate(B[r]), (r,))`                      |

`A[r] + B[r]` raises `ValueError` if the value shapes of `A` and `B` differ,
because the extents of `r` differ.

*Downsides:* expressions that differ only by bound-index labels are not equal
and are not shared, every substitution must avoid capturing indices, and
derived indices enlarge the graph. There are well-established algorithms
for simplifying these expressions

### Free indices

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

### ComponentTensor and IndexSum

`ComponentTensor` and `IndexSum` are concrete subclasses of
`AbstractExpression` that bind free indices.
`ComponentTensor(expression, indices)` turns the indices into tensor axes.
`IndexSum(expression, indices)` sums over the indices.

Both store a scalar expression and a non-empty `tuple` of distinct `Index` and
`EllipsisIndex` entries that are free in the expression. Construction raises
`ValueError` if the expression is not scalar, if the indices are empty, if an
index is repeated, or if an index is not free in the expression. Both derive the
free indices of the expression without the bound indices.

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
- `a[i] + b[j]` raises `ValueError`, because the free indices differ.
- `A[i, i]` raises `ValueError`, because `i` occurs twice.

> Note: `opt_einsum` can directly plan a contract path for a given contraction,
> but assumes dense operators.

## Annotation

> LLM generated content.

An annotation records how an expression was produced or incidental metadata.
An annotation is a `tuple` of a `str` name followed by the arguments that
produced the expression, for example `("uflx.metric", D)`. Every argument must be an
allowed construction argument (see [The base object](#the-base-object-base)).

`AnnotatedExpression(expression, annotation)` is a concrete subclass of
`AbstractExpression` that stores the wrapped expression and the annotation.

- Identity. `AnnotatedExpression` derives the value shape, the free indices,
  the scalar type and the set of the wrapped expression, so an annotated
  expression stands in wherever an expression is accepted.
- Equality. Equality is structural, so an annotated expression is not equal to
  the wrapped expression. A later algorithmic pass may compare expressions
  while ignoring annotations.
- Passes. An annotation adds structure and never removes structure, so
  annotation is not a simplification. An algorithmic pass may drop any
  annotation. A pass that rewrites inside an annotated expression drops the
  annotation, unless the pass shows that the annotation still holds.

Every free function that is defined as a composition annotates the result with
the identifier of the free function, followed by the arguments of the free
function in the order of the signature. An argument passed by keyword or
filled by a default value gives the same annotation as an argument passed by
position. For example,
`metric(D)` is `AnnotatedExpression(transpose(J) @ J, ("uflx.metric", D))`.

The identifier of a free function is the public name of the free function,
qualified by the package name, such as `"uflx.metric"` for the core or
`"library.trace"` for a library with the package name `library`. The package
name makes the identifier unique without a registry. The identifier names the
public interface, not the module that defines the free function, so moving the
free function to another module does not change the identifier. Importing the
package and looking up the public name gives back the free function, so a pass
can check an annotation by calling the free function on the arguments and
comparing the result with the wrapped expression.

A user may annotate an expression with any `str` name. A name with a dot is
reserved for the identifier of a free function.

Examples:

- `annotate(X[0]**Integer(2) + X[1]**Integer(2), ("height", X))` is an
  `AnnotatedExpression` with value shape `()`, and prints as `height(X)`.
- `metric(D)` prints as `metric(D)`, because the printer shows the last
  component of the identifier.

## Tensor algebra

### Determinant

`Determinant` is a concrete subclass of `AbstractExpression` that stores one
expression. The free function `determinant(f)` takes an expression of value
shape `(n, n)` and returns a scalar. The determinant of an expression of
value shape `(n, n)` is a polynomial of degree $n$ in the entries, so no fixed
composition expresses the determinant, and `Determinant` is a node class.

Construction raises `ValueError` if the expression does not have rank 2, or if
the two extents of the value shape differ.

The determinant of an expression of value shape `(0, 0)` is 1, the empty
product.

### Transpose and matrix product

Transpose and matrix product are free functions built with derived indices:

- `transpose(A)` is `component_tensor(A[k, l], (l, k))`, with
  $k$ = `Index(("uflx.transpose", A, 0))` and
  $l$ = `Index(("uflx.transpose", A, 1))`, annotated with `("uflx.transpose", A)`.
- `matrix_product(A, B)`, with the overload `A @ B`, is
  `component_tensor(index_sum(A[i, k] * B[k, j], (k,)), (i, j))`, with $i$, $k$
  and $j$ the derived indices `Index(("uflx.matrix_product", A, B, m))` for
  $m = 0, 1, 2$, annotated with `("uflx.matrix_product", A, B)`.

An operand that does not have rank 2, including a scalar operand, is rejected
with `IndexError` by the constructor of `Indexed`. Unequal inner extents in
`matrix_product` are rejected with `ValueError` by the constructor of `Product`,
because the index $k$ has two extents. So the overload `@` accepts only rank 2
operands. Contractions of other ranks, such as a matrix-vector product, are
written in index notation.

Examples, with a matrix `A` of value shape `(2, 2)`:

- `transpose(A)` is not equal to `component_tensor(A[i, j], (j, i))` for user
  indices `i` and `j`, because the index labels differ and because
  `transpose(A)` is annotated.

## Differentiation

### Gradient

> LLM generated content.

A concrete subclass of `AbstractExpression` for the gradient of an expression
$f$ with respect to the coordinate function $X$ of a coordinate domain $P$ of
dimension $k$.

`Gradient(f, P)` stores $f$ and $P$ and derives the value shape
`f.value_shape + (k,)` and the set $P$.

An expression and a coordinate domain $(f, P)$ are a valid pair if $P$ is an
`AbstractCoordinateDomain` and $f$ has the set $P$ or no set. Otherwise
construction raises `TypeError` for $P$ or `ValueError` for $f$. An $f$ with
no set is a constant on $P$. `Gradient(f, P)` requires a valid pair.

`Gradient` is lazy. Construction does not evaluate the derivative. The core
transformation `apply_derivatives` evaluates the derivative (see
[Applying derivatives](#applying-derivatives)). `gradient(X, P)` is already the
identity on $P$, so `apply_derivatives` leaves `gradient(X, P)` unchanged (see
[Identity](#identity)).

Differentiation on a parameterised domain is not a primitive. Differentiation
on a parameterised domain will be a composite of `gradient` on the parameter
domain and the geometry of the parameterised domain (see
[Geometry](#geometry)).

Examples, with `X = coordinate(RD(2))`:

- `gradient(X, RD(2))` has value shape `(2, 2)` and is the identity on
  `RD(2)`.
- `gradient(X[0]**Integer(2) + X[1]**Integer(2), RD(2))` has value shape
  `(2,)`. A later pass evaluates the gradient to $(2 X_0, 2 X_1)$.
- `gradient(Integer(1), RD(2))` has value shape `(2,)`, and a later pass
  evaluates the gradient to zero.
- `gradient(X, RD(3))` raises `ValueError`, because `X` has the set `RD(2)`.

### Identity

The free function `identity(P)` is the identity on the tangent space of a
coordinate domain $P$ of dimension $k$. `identity(P)` is
`gradient(coordinate(P), P)`, annotated with `("uflx.identity", P)`. With $X$
the coordinate function of $P$, $\partial X_i / \partial X_j = \delta_{ij}$, so
the identity is a fixed composition, and the core has no identity node class.

`identity(P)` has value shape `(k, k)` and the set $P$. The identity depends on
$P$, not on $k$ alone, because the identity acts on the tangent space of $P$.
`identity(RD(0))` has value shape `(0, 0)` and needs no special case.

The Kronecker delta $\delta_{ij}$ is `identity(P)[i, j]`.

Examples, with `X = coordinate(RD(2))`:

- `identity(RD(2))` has value shape `(2, 2)` and the set `RD(2)`.
- `identity(RD(2))` is not equal to `gradient(X, RD(2))`, because
  `identity(RD(2))` is annotated.
- With `i, j = index("i"), index("j")`, `identity(RD(2))[i, j] * X[j]` has the
  free indices `((i, 2), (j, 2))`, and
  `index_sum(identity(RD(2))[i, j] * X[j], (j,))` is $X_i$.

## Geometry of domains

### ParameterisedDomain

A concrete subclass of `Base`.

The domain $D$ = `ParameterisedDomain(P, psi, Q)` stores:

- The parameter domain $P$, an `AbstractCoordinateDomain`, returned by the
  property `parameter_domain`. The dimension $k$ of $P$ is the topological
  dimension of $D$.
- The parameterisation $\psi$, an `AbstractExpression` on $P$ with value shape
  `(n,)`, returned by the property `parameterisation`. The parameterisation
  $\psi$ is built from the coordinate function $X$ of $P$, or is a constant
  with no set. The parameterisation $\psi$ maps each point of $P$ to a point of
  $Q$, and the image $\psi(P)$ is $D$.
- The ambient coordinate domain $Q$, an `AbstractCoordinateDomain`, returned by
  the property `ambient_domain`. The coordinates of $Q$ name the points of $D$.
  The dimension $n$ of $Q$ is the geometric dimension of $D$.

The 'geometric dimension' is `D.ambient_domain.dimension`, and the 'topological
dimension' is `D.parameter_domain.dimension`.

Construction raises `TypeError` if $P$ or $Q$
is not an `AbstractCoordinateDomain`. Construction raises `ValueError` if
$\psi$ does not have value shape `(n,)` with $n$ the dimension of $Q$, if
$\psi$ has a set other than $P$, if $\psi$ has free indices, if $\psi$ is
complex, or if $k > n$. Construction does not check that $\psi(P)$ lies in
$Q$.

$Q$ is usually `RD(n)`. A domain in a smaller $Q$, such as an edge of the unit
triangle in the unit triangle, can be composed with an expression on $Q$ (see
[Composition](#composition)). The geometry of $D$ is derived from the
parameterisation (see [Geometry](#geometry)).

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
  is `RD(0)`, $\psi$ is the constant `tensor((p0, p1))`, with no set, and $Q$
  is `RD(2)`.

### Geometry

The geometry of every domain $D$ comes from the parameterisation $\psi$ =
`D.parameterisation` on the parameter domain $P$ = `D.parameter_domain` of
dimension $k$:

| Free function        | Definition                                               | Value shape |
| -------------------- | -------------------------------------------------------- | ----------- |
| `jacobian(D)`        | $J$ = `Gradient(D.parameterisation, D.parameter_domain)` | `(n, k)`    |
| `metric(D)`          | $g$ = `transpose(J) @ J`                                 | `(k, k)`    |
| `measure_density(D)` | `square_root(determinant(g))`                            | `()`        |

For `parameterised_domain(P, coordinate(P), P)`, the parameterisation is the
coordinate function $X$ of $P$. The general rules give
$J$ = `gradient(X, P)`, which is the identity on $P$, $g = I$ and a measure
density of 1, so the domain of the points of $P$ needs no special case.

For a point, $k = 0$. $J$ has value shape `(n, 0)`, $g$ has value shape
`(0, 0)`, and the measure density is $\sqrt{1} = 1$. So the measure of a point
is the measure at the point, and a point needs no special case either.

### Example: a paraboloid

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

# AbstractExpression: AnnotatedExpression(Gradient(psi, P), ("uflx.jacobian", D)),
# value shape (3, 2).
J = jacobian(D)

# AbstractExpression: AnnotatedExpression(transpose(J) @ J, ("uflx.metric", D)),
# value shape (2, 2).
g = metric(D)

# AbstractExpression: AnnotatedExpression(SquareRoot(Determinant(g)),
# ("uflx.measure_density", D)), value shape ().
rho = measure_density(D)
```

Note that at this stage we do not need the concept of a `Point` as we have not
evaluated anything at a `Point` so far.

### Example: an affinely mapped triangle

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

# AbstractExpression: AnnotatedExpression(Gradient(psi, T), ("uflx.jacobian", K)),
# value shape (2, 2).
J = jacobian(K)

# AbstractExpression: AnnotatedExpression(SquareRoot(Determinant(metric(K))),
# ("uflx.measure_density", K)), value shape ().
rho = measure_density(K)
```

A later pass evaluates the Jacobian to the constant matrix

$$J = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix},$$

so $\det J = 5$ and $\sqrt{\det g} = |\det J| = 5$. The area of the physical
triangle is $5$ times the area $1/2$ of the unit triangle.

The physical triangle is a parameterised domain with $k = n$, not a coordinate
domain, although the physical triangle has non-empty interior in $R^2$. The
class records how the domain is built.

### Composition

A concrete subclass of `AbstractExpression` for $f \circ \psi$, the composition
of an expression $f$ with the parameterisation $\psi$ of a domain $D$.

`Composition(f, D)` stores $f$ and $D$. With $P$ the parameter domain and $Q$
the ambient coordinate domain of $D$, `Composition(f, D)` is the expression on
$P$ with the value $f(\psi(p))$ at each point $p$ of $P$. `Composition(f, D)`
derives the value shape of $f$ and the set $P$. The free indices and the scalar
type follow the general rule.

The second operand is a domain, not an expression, because a domain states $Q$.
Construction raises `ValueError` unless $f$ has no set or the set of $f$ is
$Q$. Construction raises `TypeError` if $D$ is not a `ParameterisedDomain`.

`Composition` is lazy. The core transformation `expand_compositions` moves a
`Composition` towards the leaves and replaces `Composition(Coordinate(Q), D)`
by $\psi$ (see [Expanding compositions](#expanding-compositions)). A
`Composition` of a downstream terminal remains, and the external library that
supplies the terminal evaluates the terminal at $\psi(p)$. The core
transformation `apply_derivatives` applies the chain rule.
`Gradient(Composition(f, D), P)` is `Composition(Gradient(f, Q), D)`
contracted with `jacobian(D)`.

A domain $E$ in the parameter domain of $D$ has the image $\psi(E)$, a subset
of $D$. The image is the parameterised domain with the parameterisation
`compose(D.parameterisation, E)` on the parameter domain of $E$.

Example, with `T`, `X`, `psi`, `one` and `K` of the affinely mapped triangle:

```python
I = UnitSimplex(1)
U = coordinate(I)

# The edge of T from (1, 0) to (0, 1), as a domain in T.
e = parameterised_domain(I, tensor((one - U[0], U[0])), T)

# Composition(f, e), with the set UnitSimplex(1) and value shape ().
f = X[0] * X[1]
fe = compose(f, e)

# The edge of K from (3, 1) to (2, 3), as a domain in RD(2).
E = parameterised_domain(I, compose(psi, e), RD(2))
```

- A later pass evaluates `fe` to $(1 - U_0) U_0$ and the parameterisation of
  `E` to $(3 - U_0, 1 + 2 U_0)$.
- The chain rule gives `jacobian(E)` $= \begin{pmatrix} 2 & 1 \\ 1 & 3
  \end{pmatrix} \begin{pmatrix} -1 \\ 1 \end{pmatrix} = \begin{pmatrix} -1 \\ 2
  \end{pmatrix}$, so the measure density of `E` is $\sqrt{5}$, the length of
  the edge of `K`.
- `compose(f, parameterised_domain(I, tensor((one - U[0], U[0])), RD(2)))`
  raises `ValueError`, because that domain lies in `RD(2)`, not in `T`.

## Integrating over the domain

At this point we have everything we need to integrate an expression over a
single domain. Note we do not need to explicitly introduce any
topological notions of the domain to build a basic integral.

### LebesgueIntegral

`LebesgueIntegral(f, P)` is the integral of $f$ over a coordinate domain $P$
with the Lebesgue measure $\mathrm{d}X$ on $P$, in the coordinates $X$ of $P$:

$$
\int_P f \,\mathrm{d}X.
$$

`LebesgueIntegral(f, P)` requires a valid pair $(f, P)$ (see
[Gradient](#gradient)). `LebesgueIntegral(f, P)` stores $f$ and $P$ and derives
the value shape, the free indices and the scalar type of $f$.
`LebesgueIntegral(f, P)` has no set.

`LebesgueIntegral` is lazy. Construction does not evaluate the integral. A later
pass might evaluate the integral, for example by quadrature.

Note that a `LebesgueIntegral` is an expression, so integrals can be added,
multiplied and used wherever an expression is accepted.

The naming is slightly pretentious for good reason here - the vast majority of
users will want to use the `integrate` primitive below which takes a domain.

### Measures

A measure $\mu$ on a domain $D$, with parameter domain $P$, is given by a
density $\rho$ with respect to the Lebesgue measure $\mathrm{d}X$ on $P$:

$$
\int f \,\mathrm{d}\mu = \int_P f \rho \,\mathrm{d}X,
$$

where $f$ is a scalar expression on $P$ or a constant. $\rho \,\mathrm{d}X$ is the
volume element of $\mu$. The core has no measure object. The integral with any
such measure is `lebesgue_integral(f * rho, P)`, with $\rho$ an ordinary
expression.

The induced measure $\mu_D$ of $D$ is the measure that the parameterisation
$\psi$ induces from the Euclidean metric of $R^n$, with density
$\sqrt{\det g}$, $g = J^T J$ and $J = \nabla_P \psi$. The density is
`measure_density(D)`. The induced measure needs no orientation, because
$\sqrt{\det g} \ge 0$.

### integrate

The free function `integrate(f, D)` integrates $f$ over the domain $D$ with the
induced measure $\mu_D$. `integrate(f, D)` is
`lebesgue_integral(f * measure_density(D), D.parameter_domain)`, annotated with
`("uflx.integrate", f, D)`. The free function is the definition above, written as a
fixed composition. Arithmetic is scalar, so $f$ is a scalar, possibly with free
indices. A tensor-valued integral is written in index notation, such as
`component_tensor(integrate(f[r], D), (r,))`, with `r` a
[EllipsisIndex](#ellipsisindex).

Most users will use this method.

#### Examples

With `T = UnitSimplex(2)`, `X = coordinate(T)`, and the physical triangle `K`
of the affinely mapped triangle example:

- `lebesgue_integral(X[0], T)` is $\int_T X_0 \,\mathrm{d}X$. A later pass
  evaluates the integral to $1/6$.
- `integrate(Integer(1), K)` is the area of the physical triangle. The measure
  density of `K` is $\lvert \det J \rvert = 5$, so a later pass evaluates the
  integral to $5/2$.
- With `h = X[0]**Integer(2) + X[1]**Integer(2)` and
  `S = parameterised_domain(T, tensor((X[0], X[1], h)), RD(3))`,
  `integrate(Integer(1), S)` is the area of the paraboloid over the unit
  triangle. The measure density of `S` is $\sqrt{1 + 4 X_0^2 + 4 X_1^2}$.
- With `four, zero, one = Integer(4), Integer(0), Integer(1)` and the metric
  `G = tensor((tensor((four, zero)), tensor((zero, one))))` on `T`,
  `lebesgue_integral(one * square_root(determinant(G)), T)` is the integral of
  1 with the measure of `G`, which is $\sqrt{4} \cdot 1/2 = 1$.
- With `c = Float(0.5)` and
  `p = parameterised_domain(RD(0), tensor((Float(1.0), Float(2.0))), RD(2))`,
  a later pass evaluates `integrate(c, p)` to `c`, because the Lebesgue
  measure on $R^0$ is the counting measure at the single point of $R^0$.
- `integrate(X[0], K) * integrate(X[1], K)` is a product of two integrals, with
  value shape `()` and no set.
- With `fe` and `E` of the [Composition](#composition) example,
  `integrate(fe, E)` is the integral of $X_0 X_1$ over the edge of `K`, which
  is $\sqrt{5}/6$.

### But now don't we need a point?

A point is a domain! A point in a coordinate domain $Q$ of dimension $n$ is
`parameterised_domain(RD(0), x, Q)`, with $x$ a constant expression of value
shape `(n,)`. The core needs no node class and no free function for a point.

The value of a scalar expression $f$ at a point $p$ is
`integrate(compose(f, p), p)`, the integral of $f$ with the counting measure of
$p$. `compose(f, p)` is an expression on `RD(0)`, and the integral has no set,
so the integral combines with any expression.

The value of a tensor-valued expression at a point is written in index
notation, because `integrate` takes a scalar. For $f$ of any value shape, the
value at $p$ is `component_tensor(integrate(compose(f, p)[r], p), (r,))`, with
`r` a [EllipsisIndex](#ellipsisindex).

> Alternative: `lebesgue_integral(compose(f, p), RD(0))` gives the
> tensor-valued value directly, because the measure density of a point is 1.
> The composition is fixed, so a core free function `evaluate(f, p)` could
> build the composition. If $p$ is not a point, the constructor of
> `LebesgueIntegral` raises `ValueError`. The cost is a second way to write
> point evaluation beside `integrate`.

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

# AnnotatedExpression: the value of f at the centroid, with no set and value
# shape (). A later pass evaluates the integral to 1/9.
fc = integrate(compose(f, c), c)

# ParameterisedDomain: the centroid of K, as a point in RD(2). A later pass
# evaluates compose(psi, c) to (2, 4/3).
C = parameterised_domain(RD(0), compose(psi, c), RD(2))
```

- With a second point
  `d = parameterised_domain(RD(0), tensor((third, Integer(2) / Integer(3))), T)`,
  `integrate(compose(f, c), c) + integrate(compose(f, d), d)` is the action of
  a sum of two Dirac measures on $f$.
- `parameterised_domain(RD(0), tensor((third, third)), RD(2))` is the point
  $(1/3, 1/3)$ in the plane. With that point as `q`, `compose(f, q)` raises
  `ValueError`, because `q` lies in `RD(2)`, not in `T`.
- `parameterised_domain(RD(0), X, T)` raises `ValueError`, because `X` has the
  set `T`, not `RD(0)`.

## Core transformations

The core has two transformations, `apply_derivatives` and
`expand_compositions`. Each core transformation is a pass with one rule for
each node class of the core, and produces only node classes of the core. A rule
replaces a node by a mathematically equal expression and never simplifies.

A terminal is a `Coordinate`, or an expression with a set whose class is a leaf
class of a downstream library. Neither core transformation rewrites a terminal.
Every other pass, such as simplification or `units`, is not a core
transformation.

### Applying derivatives

`apply_derivatives(f)` evaluates every `Gradient` in $f$ with the standard
rules of differentiation, including the chain rule for a `Composition`. The
gradient of an expression with no set is `Zero`. `Gradient(Coordinate(P), P)`
is the identity on $P$ and remains (see [Identity](#identity)). Afterwards,
every `Gradient` is applied, one or more times, to a terminal.

### Expanding compositions

`expand_compositions(f)` moves every `Composition` in $f$ towards the leaves,
because composition commutes with pointwise operations, and replaces
`Composition(Coordinate(Q), D)` by the parameterisation of $D$. A
`Composition` of a downstream terminal or of a `Gradient` remains.
`expand_compositions` is the pull-back to the parameter domain.

## Topology of domains

### Cells

A domain can be divided into finitely many cells. A cell is a piece of the
domain with the same dimension as the domain. Two cells meet only on the
boundaries of the two cells.

A facet is a piece of the boundary of a cell, one dimension lower than the
cell. A facet in the boundary of exactly one cell is an exterior facet. A facet shared
by two cells is an interior facet.

Cells and facets are topological. Defining cells and facets does not require a
metric.

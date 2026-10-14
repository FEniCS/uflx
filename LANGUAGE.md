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
- Pure free functions should be strong preferred over classes. Class members,
  where needed, should immediately call free functions wherever possible. Pure
  free functions maximise the chance of code re-use and are easier to reason
  about.
- A parent class should not return a subclass. A method of a parent class that
  needs a subclass calls a free function, and the free function constructs the
  subclass.
- A concept must be a free function, not a node class, if the free function can
  build the concept as one fixed composition of existing node classes applied to
  the arguments, *without inspecting the arguments to choose the composition*.
  Otherwise the concept gets its own node class. For example, the metric is
  always `transpose(J) @ J`, so the metric is a free function. No composition of
  existing nodes differentiates, so `Gradient` is currently a node class. A free
  function that builds a fixed composition is a definition, not a
  simplification. The free function does not inspect the operands, so the free
  function never removes structure that the user built.
- Under the previous rule, a node class may become a free function as UFLx core
  grows, once the node class can be written as a fixed composition of those
  newer nodes *without inspecting the arguments to choose the composition*. For
  example, if subtraction were a node class `Difference`, `subtract(f, g)`
  should become the free function `add(f, negate(g))` once `negate` exists. The
  composition is the same for every `f` and `g`, so the composition is allowed.
  In contrast, `transpose(A)`, assuming no index notation, cannot be implemented
  as:

  ```python
  m, n = A.value_shape
  tensor(tuple(tuple(A[j, i] for j in range(m)) for i in range(n)))
  ```

  because the composition reads the value shape of `A` to choose the entries.
  With free indices, `transpose(A)` is the fixed composition
  `component_tensor(A[i, j], (j, i))`, so `transpose` is a free function.
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

  A pass differs from a free function used during construction in two ways. A
  pass may read any property of an object, such as the value shape, to choose
  what the pass builds. A pass may rewrite the graph, for example by
  simplifying.

  A pass does one thing, such as differentiating or lowering. Passes compose,
  and the user of the core chooses which passes run and in what order.

Three kinds of code produce objects:

| Kind                                                  | When the code runs                | Inspects arguments to choose what to build                | Validates                      |
| ----------------------------------------------------- | --------------------------------- | --------------------------------------------------------- | ------------------------------ |
| Constructor                                           | Whenever an object is made        | No, a constructor inspects arguments only to validate     | Yes, and only constructors do  |
| Construction free function, such as `transpose`       | While the user builds the graph   | No, the free function builds one fixed composition        | No                             |
| Pass, such as differentiation, lowering or `units`    | After construction                | Yes, and a pass may rewrite the graph                     | No                             |

Every object that a free function or a pass returns comes from a constructor, so
every object that a free function or a pass returns is valid. A call can still
fail. A construction free function fails when a constructor rejects the
arguments, for example `transpose` of a vector raises `IndexError` from the
constructor of `Indexed`. A pass fails with `TypeError` for a node class that
the pass has no rule for, or with an error of the pass, such as inconsistent
units.

## The minimal object model

### The base object `Base`

The base object in UFLx must be extended by all objects that are included in
the graph. The graph is a directed acyclic graph (DAG).

Everything needed to define the variational problem must be recoverable from
the graph. No information is held outside of the graph.

Objects must satisfy the following rules:

1. All construction arguments must be deeply immutable values. For example, a
   `tuple` of `list`s is not allowed.
2. All construction arguments must be compared by type and value, so that `1`,
   `1.0` and `True` are distinct, `0.0` and `-0.0` are also distinct, and every
   value is equal (`==` or `is`, depending on implementation) to itself.
3. Every construction argument must have an identity (hash) that is
   deterministic across processes on the same machine.
4. Construction must depend only on the construction arguments, not on global
   state.
5. An object must be reconstructible from its class and construction arguments,
   and the reconstructed object must be equal to the original.
6. Any results produced by the class after construction must be a result of
   only the class and the construction arguments.
7. The state of the class is immutable after construction.
8. Any order that the core or a pass relies on must be deterministic across
   processes. A pass, and any other walk of the graph such as hashing,
   serialisation or printing, follows the edges of each object in order (rule
   9). Where unrelated objects must be sorted, the sort uses a total order
   defined by the class and the construction arguments alone. Rules on sorting
   comes in later passes, e.g. deciding whether `a + b` and `b + a` are
   mathematically equivalent.
9. The edges of an object must be exactly the base objects among its
   construction arguments, including the base objects nested in `tuple`s, in
   order.
10. Equality (`==` or `is`) must be structural and return a `bool`. Equality
    does not decide mathematical equivalence, so $a + b$ and $b + a$ are not
    equal. Mathematical equivalencies are applied in later graph algorithms.
11. All validation of an object happens in the constructor of the class of the
    object. This rule is a hard rule. A free function, including a free
    function that builds a composition, never validates. A free function never
    converts arguments. A free function only reads derived properties of
    arguments and constructs objects. An invalid argument is rejected by the
    constructor of the first object that the argument makes invalid.

Hence two objects are equal if and only if the two objects have the same class
and equal construction arguments.

To satisfy the above, a construction argument must be one of a small curated
set of types: a base object, exactly `int`, a finite `float`, `str`, or a
`tuple` whose items are of these types. A `bool` is not an `int` for this rule.

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
│   ├── AbstractScalar                    literals, no set
│   │   ├── Integer
│   │   ├── Rational
│   │   ├── RealScalar
│   │   └── ComplexScalar
│   ├── Coordinate                        terminal, X on P
│   ├── Tensor                            built from scalar entries
│   ├── Indexed                           full indexing, int or Index entries
│   ├── ComponentTensor                   binds indices to tensor axes
│   ├── IndexSum                          sums over bound indices
│   ├── AnnotatedExpression               Annotation mixin, identity
│   ├── Gradient                          derivative with respect to X of P
│   ├── Integral                          Lebesgue integral over P, no set
│   ├── UnaryOperator
│   │   ├── Determinant
│   │   ├── SquareRoot
│   │   └── Absolute
│   └── BinaryOperator
│       ├── Sum                           +
│       ├── Product                       *
│       ├── Division                      /
│       └── Power                         **
│
├── Index                                 named by a label, not an expression
│
└── AbstractDomain                        parameter_domain, parameterisation, ambient_domain
    ├── AnnotatedDomain                   Annotation mixin, identity
    ├── AbstractCoordinateDomain          stores n, parameterisation is X
    │   └── RD                            R^n
    └── ParameterisedDomain               stores P and ψ, parameterisation is ψ
```

### Free functions

Free functions are not objects in the graph. A free function builds or returns
objects in the graph. Operator overloads call the free functions.

> LLM generated table.

| Group                 | Free function                                           | `AbstractExpression` overload | Takes                                                                | Returns                                                       |
| --------------------- | ------------------------------------------------------- | ----------------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------- |
| Arithmetic            | [`add(f, g)`](#arithmetic)                              | `+`                           | `f, g: AbstractExpression`                                           | `Sum`                                                         |
|                       | [`subtract(f, g)`](#arithmetic)                         | `-`                           | `f, g: AbstractExpression`                                           | `AnnotatedExpression`: `Sum` of `f` and `negate(g)`           |
|                       | [`multiply(f, g)`](#arithmetic)                         | `*`                           | `f, g: AbstractExpression`                                           | `Product`                                                     |
|                       | [`divide(f, g)`](#arithmetic)                           | `/`                           | `f, g: AbstractExpression`                                           | `Division`                                                    |
|                       | [`negate(f)`](#arithmetic)                              | unary `-`                     | `f: AbstractExpression`                                              | `AnnotatedExpression`: `Product` of `Integer(-1)` and `f`     |
|                       | [`absolute(f)`](#arithmetic)                            | `abs`                         | `f: AbstractExpression`                                              | `Absolute`                                                    |
|                       | [`power(f, g)`](#arithmetic)                            | `**`                          | `f, g: AbstractExpression`                                           | `Power`                                                       |
| Elementary functions  | [`square_root(f)`](#elementary-functions)               |                               | `f: AbstractExpression`                                              | `SquareRoot`                                                  |
| Domains               | [`ambient_domain(D)`](#abstractdomain)                  |                               | `D: AbstractDomain`                                                  | `RD(n)`                                                       |
|                       | [`parameter_domain(D)`](#abstractdomain)                |                               | `D: AbstractDomain`                                                  | `AbstractCoordinateDomain`                                    |
|                       | [`parameterisation(D)`](#abstractdomain)                |                               | `D: AbstractDomain`                                                  | `AbstractExpression` with shape `(n,)`                        |
|                       | [`coordinate(P)`](#coordinate)                          |                               | `P: AbstractCoordinateDomain`                                        | `Coordinate`                                                  |
| Tensors and indexing  | [`tensor(entries)`](#tensor)                            |                               | `entries: tuple` (nested, of `AbstractExpression`)                   | `Tensor`                                                      |
|                       | [`tensor(entries)`](#tensor)                            |                               | `entries: tuple[AbstractExpression, ...]`                            | `Tensor` of rank 1                                            |
|                       | [`tensor(entries)`](#tensor)                            |                               | `entries: tuple[tuple[AbstractExpression, ...], ...]`                | `Tensor` of rank 2                                            |
|                       | [`indexed(f, indices)`](#indexed)                       | `f[i]`                        | `f: AbstractExpression`, `indices: tuple[int \| Index, ...]`         | `Indexed`                                                     |
| Index notation        | [`index(label)`](#index)                                |                               | `label: str \| tuple`                                                | `Index`                                                       |
|                       | [`indices(*labels)`](#index)                            |                               | `*labels: str \| tuple`                                              | `tuple` of `Index`                                            |
|                       | [`component_tensor(f, indices)`](#componenttensor)      |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`                | `ComponentTensor`                                             |
|                       | [`index_sum(f, indices)`](#indexsum)                    |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`                | `IndexSum`                                                    |
| Annotation            | [`annotate(x, name, *arguments)`](#annotation)          |                               | `x: AbstractExpression \| AbstractDomain`, `name: str`, `*arguments` | `AnnotatedExpression` or `AnnotatedDomain`                    |
| Tensor algebra        | [`determinant(f)`](#determinant)                        |                               | `f: AbstractExpression`                                              | `Determinant`                                                 |
|                       | [`transpose(f)`](#transpose-and-matrix-product)         |                               | `f: AbstractExpression`                                              | `AnnotatedExpression`: `ComponentTensor`                      |
|                       | [`matrix_product(f, g)`](#transpose-and-matrix-product) | `@`                           | `f, g: AbstractExpression`                                           | `AnnotatedExpression`: `ComponentTensor` of an `IndexSum`     |
| Differentiation       | [`grad(f)`](#gradient)                                  |                               | `f: AbstractExpression`                                              | `Gradient`                                                    |
| Parameterised domains | [`parameterised_domain(P, psi)`](#parameteriseddomain)  |                               | `P: AbstractCoordinateDomain`, `psi: AbstractExpression`             | `ParameterisedDomain`                                         |
| Geometry              | [`jacobian(D)`](#geometry)                              |                               | `D: AbstractDomain`                                                  | `AnnotatedExpression`: `Gradient` of `parameterisation(D)`    |
|                       | [`metric(D)`](#geometry)                                |                               | `D: AbstractDomain`                                                  | `AnnotatedExpression`: `transpose(J) @ J`                     |
|                       | [`measure_density(D)`](#geometry)                       |                               | `D: AbstractDomain`                                                  | `AnnotatedExpression`: `SquareRoot` of the determinant of $g$ |
| Integration           | [`integral(f, P)`](#integral)                           |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain`               | `Integral`                                                    |
|                       | [`integrate(f, D)`](#integrate)                         |                               | `f: AbstractExpression`, `D: AbstractDomain`                         | `AnnotatedExpression`: `Integral` over `parameter_domain(D)`  |

Every overload in the table is a method of `AbstractExpression` (see
[`AbstractExpression`](#abstractexpression)). No other class defines
arithmetic, matrix or indexing overloads. Overloading should be considered
carefully outside of standard mathematical meanings.

The value-shape and other constraints on the arguments are stated in the section
of each free function. A free function and an operator overload accept only the
types in the "Takes" column.

## Extending the core

Extensibility is a requirement of the core. Downstream libraries must be able
to extend UFLx without changing the core.

> Legacy note: In legacy UFL, properties are stored on the node classes, and
> the set of properties is fixed. Every new property touches every node class,
> so many language features were hacked in, for example complex numbers and
> dual spaces.

The core has a small enough, fixed set of node classes and an open set of passes. A
downstream library extends UFLx by adding passes, not by adding properties to
the node classes of the core. A new pass needs one rule for each node class of
the core and does not change the core.

The core supports five kinds of extension:

> LLM generated table, could be made more concise.

| Kind of extension                  | Examples                                                | Section                                                   |
| ---------------------------------- | ------------------------------------------------------- | --------------------------------------------------------- |
| New questions asked of a form      | Physical units, polynomial degree, linearity checks     | [New questions](#new-questions)                           |
| New kinds of leaf                  | Elements and functions from Basix, new geometric leaves | [New leaves](#new-leaves)                                 |
| New operations                     | Trace, external operators                               | [New operations](#new-operations)                         |
| New transformations                | Pull-back, Gateaux derivative                           | [New transformations](#new-transformations)               |
| New properties that the core needs | Complex numbers                                         | [New properties of the core](#new-properties-of-the-core) |

### New questions

A downstream library asks a new question of a form with a pass. The pass has
one rule for each node class of the core and one rule for each leaf class of
the library.

For example, a library that checks physical units, like the work of Habera et
al.:

1. Assigns a unit to leaves in a table of the library. For example, with
   `X = coordinate(RD(2))`, the table `{X: "m"}` states that the coordinates
   are in metres. A leaf that is not in the table, such as a literal, is
   dimensionless.
1. Writes passes, such as `units(expression, table)`, which returns the unit of
   any expression. `units` derives the unit of an expression from the units of
   the leaves, in the same way as the core derives the value shape. `X[0] *
   X[1]` is in square metres. `grad(X[0])` is in metres divided by metres, so
   `grad(X[0])` is dimensionless. `X[0] + X[0] * X[1]` is an error of the units
   library.
1. Calls `units` on a form whenever the library needs the check.

The core never knows about units. Objects of the core are immutable and are
compared by content, so the leaves are reliable keys of the table.

An error reported by `units` does not make an object invalid. An object is
valid if the constructors of the core accept the object.

The pass never guesses a rule. A `TypeError` from `units` means that the units
library does not know a node class. An error of the units library means that the
expression has inconsistent units.

### New leaves

Each kind of leaf has an abstract class in the core, such as `AbstractScalar` or
`AbstractCoordinateDomain`. A downstream leaf class is a subclass of the
abstract class. The abstract class is the contract between the downstream leaf
class and the transformations of the core. For example, differentiation gives
zero for every `AbstractScalar`, so differentiation handles a downstream
subclass of `AbstractScalar` without a new rule.

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

A downstream library adds a new operation in one of two ways:

- Builds the operation from existing node classes, as a free function.
- Defines the operation as a node class of the library, together with a pass
  that replaces every node of that class by node classes of the core. The
  library runs the pass before differentiation or any other transformation of
  the core. For example, a node class `Trace` is replaced by the sum of the
  diagonal entries of the matrix. The pass reads the value shape of the matrix
  to choose the entries.

An operation that fits neither way needs a language extension. A language
extension adds node classes and tells every transformation of the core how to
handle the new node classes. Every transformation of the core must still
produce only node classes of the core or of the language extension. For
example, a language extension for external operators adds a node class for an
operation known only by name, such as `exp`, a material law or a neural
network, together with the rule that differentiates the node class.

### New transformations

A new transformation, such as a pull-back, is a pass with one rule for each
node class of the core. A new transformation
must produce only node classes of the core. A transformation that needs a new
node class is a proposal to change the core.

There should be helpers to check whether a transformation covers all classes
of the core.

### New properties of the core

A property that the rules of the core depend on is part of the core, not an
extension. For example, the scalar type is part of the core, because
`Absolute` is real for a complex operand. Adding such a property changes the
core.

## The lowest objects

### Expressions

An expression denotes an array-valued function defined on the set of the
expression. An expression with no set denotes a single array.

#### `AbstractExpression`

Every `AbstractExpression` derives the value shape, the free indices, the
scalar type and the set of the expression.

- The set is an `AbstractDomain`. The set is absent for an expression with no
  set, such as a literal.
- The scalar type is real or complex. Unless a node class states otherwise, an
  expression is complex if and only if at least one edge is a complex
  expression. `ComplexScalar` is complex.

An operator overload calls the free function with the operands unchanged, and
the constructor validates. An operand that is not an `AbstractExpression`, such
as the Python `int` in `X[0] + 2`, is rejected with `TypeError` by the
constructor. `AbstractExpression` defines no reflected overloads, such as
`__radd__`, so `2 + X[0]` also raises `TypeError`.

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
    def free_indices(self) -> tuple[tuple[Index, int], ...]:
        """The free indices and extents, in order of first occurrence."""

    @property
    @abstractmethod
    def scalar_type(self) -> Literal["real", "complex"]:
        """The scalar type of the expression."""

    @property
    @abstractmethod
    def set(self) -> AbstractDomain | None:
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
literal's value is part of the expression itself.

Base:

```
AbstractScalar(AbstractExpression)
```

Concrete:

```
Integer(2)
Rational(2, 4)
RealScalar(3.5)
ComplexScalar(RealScalar(3.0), RealScalar(4.0))
```

Each literal constructor accepts only exact types and raises `TypeError` or
`ValueError` otherwise:

- `Integer(v)` takes an `int` $v$. A `bool` is rejected, although `bool` is a
  subclass of `int`.
- `Rational(p, q)` takes two `int`s $p$ and $q$ with $q \ne 0$. $p$ and $q$
  are stored as given, without normalisation, so
  `Rational(2, 4) != Rational(1, 2)`.
- `RealScalar(v)` takes a finite `float` $v$. NaN and infinities are rejected,
  because a NaN is not equal to itself. `RealScalar(-0.0) != RealScalar(0.0)`.
- `ComplexScalar(re, im)` takes two `RealScalar`s.

Literals should not fold on construction, e.g.

```
Integer(2) + Integer(3) != Integer(5)
Integer(2) + Integer(3) == Sum(Integer(2), Integer(3))
```

#### Numbers in expressions

Every number in an expression is written as a literal, such as `Integer(2)`,
`RealScalar(0.5)`, `Rational(1, 3)` or
`ComplexScalar(RealScalar(0.0), RealScalar(1.0))`. A Python `int`, `float`,
`complex`, `Fraction` or NumPy scalar is never converted and is rejected with
`TypeError`. The type of every literal is a choice that the user makes
explicitly.

Examples, with `X = coordinate(RD(2))`:

- `X[0]**Integer(2)`, `Integer(2) * X[0]` and `Integer(1) - X[0]` are
  accepted.
- `X[0]**2`, `2 * X[0]`, `1 - X[0]` and `0.5 * X[0]` raise `TypeError`.
- `tensor((1, 2))` raises `TypeError`. `tensor((Integer(1), Integer(2)))` is
  accepted.

### Operators

```UnaryOperator(AbstractExpression)```
```BinaryOperator(AbstractExpression)```

Downstream basic ops are relatively obvious.

### Arithmetic

Concrete subclasses of `BinaryOperator`. Each class stores two expressions,
which are the two edges in order. The set is the set of the operands. Operands
with a set must have the same set, or construction raises `ValueError`.
Construction raises `ValueError` for wrong value shapes.

| Class      | Free function | Overload | Operand shapes           | Result shape |
| ---------- | ------------- | -------- | ------------------------ | ------------ |
| `Sum`      | `add`         | `+`      | `s`, `s`                 | `s`          |
| `Product`  | `multiply`    | `*`      | `()`, `s` or `s`, `()`   | `s`          |
| `Division` | `divide`      | `/`      | `s`, `()`                | `s`          |
| `Power`    | `power`       | `**`     | `()`, `()`               | `()`         |

The exponent of `Power` has no set, or construction raises `ValueError`. The
derivative of a power whose exponent has a set needs a logarithm, which the
core does not have.

The free indices of every operator follow the general rule of
[Free indices](#free-indices).

`Absolute` is a concrete subclass of `UnaryOperator`. The free function
`absolute(f)`, with the overload `abs(f)`, takes a scalar and returns a scalar.
For a complex argument, `Absolute` is the modulus. `Absolute` is real for
every argument.

Negation and subtraction are free functions, not node classes:

- `negate(f)`, with the overload `-f`, is `multiply(Integer(-1), f)`,
  annotated with `("negate", f)`.
- `subtract(f, g)`, with the overload `f - g`, is `add(f, negate(g))`,
  annotated with `("subtract", f, g)`.

Both compositions are exact in floating-point arithmetic, so neither
composition changes the arithmetic that the user wrote. Division is a node
class, because `multiply(f, power(g, Integer(-1)))` rounds twice and differs
from `f / g` in floating-point arithmetic.

### Elementary functions (very incomplete)

`SquareRoot` is a concrete subclass of `UnaryOperator`. The free function
`square_root(f)` takes a scalar and returns a scalar. Like every unary
operator, `SquareRoot` stores one expression, which is the only edge, derives
the set of that expression and performs no simplification.

`SquareRoot` is a node class, not the composition `power(f, Rational(1, 2))`.
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

#### AbstractDomain

An abstract base class for domains `AbstractDomain`.

Every `AbstractDomain` $D$ has:

- The geometric dimension $n$, which is the dimension of the Euclidean space
  $R^n$ whose coordinates name the points of $D$.
- The topological dimension $k$, which is the dimension of the parameter domain
  of $D$.

Example:

- A sphere in $R^3$ has geometric dimension 3 and topological dimension 2.

Three free functions accept every `AbstractDomain` $D$ with topological
dimension $k$ and geometric dimension $n$:

- `parameter_domain(D)` returns an `AbstractCoordinateDomain` $P$ of dimension
  $k$.
- `parameterisation(D)` returns an `AbstractExpression` $\psi$ on $P$ with value
  shape `(n,)`. The parameterisation $\psi$ maps each point of $P$ to a point of
  $D$, and the image $\psi(P)$ is $D$.
- `ambient_domain(D)` returns `RD(n)`.

The geometry of $D$ is derived from the parameterisation (see
[Geometry](#geometry)).

#### AbstractCoordinateDomain

An abstract base class for coordinate domains `AbstractCoordinateDomain`, and a
subclass of `AbstractDomain`.

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$ whose
points are named by the coordinates of the points themselves.

`AbstractCoordinateDomain` is distinguished from other subclasses of
`AbstractDomain` by how the domain is built, not by the set of points. A
parameterised domain with $k = n$ can have non-empty interior in $R^n$, but a
parameterised domain is not a coordinate domain (see
[Example: an affinely mapped triangle](#example-an-affinely-mapped-triangle)).

Every concrete `AbstractCoordinateDomain` $P$ stores:

- The dimension $n$, provided by the property `dimension`.

An `AbstractCoordinateDomain` $P$ derives:

- The geometric dimension and the topological dimension, which are both $n$.
- The ambient domain `ambient_domain(P)`, which is `RD(n)`.
- The coordinate function `coordinate(P)` (see [Coordinate](#coordinate)).
- The parameter domain `parameter_domain(P)`, which is $P$.
- The parameterisation `parameterisation(P)`, which is the identity map
  `coordinate(P)`.

Examples:

- The unit simplex in $R^n$ and $R^n$ itself are coordinate domains.
- A face of codimension 1 of the unit simplex in $R^n$ is *not* a coordinate
  domain, because the face has empty interior in $R^n$. The face is the image
  of the unit simplex in $R^{n-1}$ under a map, so the face is a parameterised
  domain (defined below).

#### RD

A concrete subclass of `AbstractCoordinateDomain` for the whole of $R^n$.

`RD(n)` stores the dimension $n$, an `int` with $n \ge 0$. The integer $n$ is
not a base object, so `RD(n)` has no edges in the graph. `RD(m) == RD(n)` if
and only if $m = n$. Everything else comes from `AbstractCoordinateDomain`, so
`ambient_domain(RD(n)) == RD(n)`.

Example:

- `RD(2)` is the parameter domain of the paraboloid example below.

#### Coordinate

A concrete subclass of `AbstractExpression` for the coordinate function of a
coordinate domain.

Every concrete coordinate domain $P$ provides a coordinate function
`Coordinate(P)`, obtained with the free function `coordinate(P)`. The
coordinate function $X$ of $P$ is the terminal expression that maps each point
of $P$ to the $n$-tuple of that point:

$X : P \to R^n$, $X(p) = p$.

`Coordinate(P)` stores the coordinate domain $P$, which is the only edge.
`Coordinate(P)` derives the value shape `(n,)` and the set $P$.

Only coordinate domains provide a coordinate function. The coordinates of the
points of any other domain are given by `parameterisation` (see
[AbstractDomain](#abstractdomain)), which is built from the coordinate function
of the parameter domain.

Examples:

- `coordinate(RD(3))` is `Coordinate(RD(3))`, with value shape `(3,)` and the
  set `RD(3)`.
- `coordinate(RD(3)) == coordinate(RD(3))` is `True`, and
  `coordinate(RD(2)) == coordinate(RD(3))` is `False`.
- In the paraboloid example, `X = coordinate(RD(2))` has value shape `(2,)`.

## Tensors and indexing

### Tensor

A concrete subclass of `AbstractExpression` that builds a tensor-valued
expression from scalar expressions.

`Tensor(entries)` stores:

- The entries, a nested `tuple` of scalar `AbstractExpression`s. The nesting
  depth $r \ge 1$ is the rank of the tensor.

The nested `tuple` must be rectangular. Every `tuple` at the same depth must
have the same length, and every length must be at least 1. Every entry must
have value shape `()`.

`Tensor(entries)` derives:

- The value shape `(n_0, ..., n_{r-1})`, where $n_i$ is the length of every
  `tuple` at depth $i$.
- The set of the tensor. The set of the tensor is the set of the entries. All
  entries with a set must have the same set, or construction raises
  `ValueError`. If no entry has a set, the tensor has no set.

The edges of `Tensor(entries)` are the entries in row-major order, so the
entry `entries[i_0]...[i_{r-1}]` comes before every later entry in
lexicographic order of the indices $(i_0, \ldots, i_{r-1})$.

Construction performs no simplification. For example,
`tensor((Integer(1) + Integer(1),)) != tensor((Integer(2),))`.

The free function `tensor(entries)` returns `Tensor(entries)`.

A value shape may contain the extent 0, for example the value shape `(0,)` of
the coordinate function of `RD(0)`. `Tensor` cannot build an expression with
the extent 0, because a `tuple` of length 0 has no entries to give the nesting
depth.

Examples, with `zero, one, two = Integer(0), Integer(1), Integer(2)`:

- `tensor((one, two))` has value shape `(2,)` and no set.
- `tensor(((one, zero), (zero, one)))` has value shape `(2, 2)` and no set.
- With `r0 = RealScalar(0.0)` and `r1 = RealScalar(1.0)`,
  `tensor(((r1, r0), (r0, r1)))` has `RealScalar` entries, so the matrix is not
  equal to `tensor(((one, zero), (zero, one)))`.
- `tensor(((one, two),))` has value shape `(1, 2)`, and `tensor((one, two))`
  has value shape `(2,)`.
- `tensor((((one,),),))` has value shape `(1, 1, 1)`.

### Indexed

> LLM generated content.

A concrete subclass of `AbstractExpression` that selects one component of a
tensor-valued expression.

`Indexed(expression, indices)` stores:

- The expression, an `AbstractExpression` of rank $r \ge 1$.
- The indices, a `tuple` of $r$ entries. Each entry is an `int` or an `Index`
  (see [Index notation](#index-notation)).

The edges of `Indexed(expression, indices)` are the expression followed by the
`Index` entries, in order.

Indexing is full indexing only. The number of indices must equal the rank $r$
of the expression. Each `int` entry $i_j$ must satisfy $0 \le i_j < n_j$, where
`(n_0, ..., n_{r-1})` is the value shape of the expression. Each `Index` entry
$i_j$ takes the extent $n_j$.

`Indexed(expression, indices)` derives the value shape `()`, the set of the
expression, and the free indices of the expression followed by the `Index`
entries.

The free function `indexed(expression, indices)` takes a `tuple` of entries
and returns `Indexed(expression, indices)`. The overload `expression[indices]`
converts a single entry into a `tuple` of one entry and calls `indexed`, so
`X[0]` and `A[0, 1]` are `indexed(X, (0,))` and `indexed(A, (0, 1))`.
Construction of `Indexed` raises:

| Error        | Cause                                                        |
| ------------ | ------------------------------------------------------------ |
| `TypeError`  | An entry that is neither exactly `int` nor `Index`.          |
| `IndexError` | A number of indices that differs from the rank. An expression with value shape `()` has rank 0. |
| `IndexError` | A negative `int` entry, or an `int` entry $i_j \ge n_j$.     |
| `ValueError` | An `Index` that occurs twice, or that is already free.       |

A negative index is rejected rather than converted, so the stored indices are
always the indices that the user wrote.

Construction performs no simplification. For example, `tensor((a, b))[0]` is
`Indexed(Tensor((a, b)), (0,))` and is not equal to `a`.

Examples:

- For the coordinate function $X$ of `RD(3)`, `X[2]` has value shape `()`.
- For a matrix `A` of value shape `(2, 2)`, `A[1, 0]` has value shape `()`,
  and `A[1]` raises `IndexError`.
- `X[-1]` and `X[True]` raise errors.
- For the coordinate function $X$ of `RD(2)`, `Tensor((X[0], X[1])) != X`,
  because construction performs no simplification.
- `tensor((Integer(1) - X[0], X[0]))` has value shape `(2,)` and the set
  `RD(2)`.
- `tensor(((X[0], X[1]),))` has value shape `(1, 2)`.
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
that a free function derives from the free function's operands. The edges of
`Index(label)` are the graph objects in the label, in order. Two indices are
equal if and only if the two labels are equal.

The free function `index(label)` returns `Index(label)`. The free function
`indices(*labels)` returns a `tuple` of `Index`, one for each label. An `Index`
has no extent of its own. An `Index` takes the extent of the axis that the
`Index` indexes.

The user names user indices. The user is responsible for not reusing a name
that is already free in an operand.

A free function that needs fresh indices derives the indices from the
operands. For example, `transpose(A)` uses `Index(("transpose", A, 0))` and
`Index(("transpose", A, 1))`. An index whose label contains `A` has an edge to
`A`. Every node that uses an `Index` has an edge to the `Index`. So an index
whose label contains `A` cannot occur inside `A`, because `A` would then
contain `A`, and the graph is acyclic. Hence
a derived index is fresh with respect to the operands, and the free function
does not inspect the operands to find a fresh index. A global counter is not
allowed, because a counter makes construction depend on global state.

Note that `transpose(A)` `is` or `==` `transpose(A)` with this design.

*Downsides:* expressions that differ only by bound-index labels are not equal
and are not shared, every substitution must avoid capturing indices, and
derived indices enlarge the graph. There are well-established algorithms
for simplifying these expressions

### Free indices

Every `AbstractExpression` derives the free indices, a `tuple` of pairs of an
`Index` and an extent. The order is the order of first occurrence, never the
order of a `set`.

Unless a node class states otherwise, the free indices of a node are the free
indices of the edges of the node that are expressions, in edge order, with each
`Index` taken once at the first occurrence. An `Index` that occurs in two edges
must have the same extent in both edges, or construction raises `ValueError`.
Terminals and literals have no free indices.

The following node classes state otherwise:

- `Sum`: both operands must have the same free indices with the same extents,
  in any order, or construction raises `ValueError`. The order of the free
  indices of the sum is the order of the first operand.
- `Indexed`: see [Indexed](#indexed).
- `ComponentTensor` and `IndexSum`: see below.

Summation is never implicit. `a[i] * b[i]` has the free index `i` and denotes
the component-wise product. The sum over `i` is `index_sum(a[i] * b[i], (i,))`.

### ComponentTensor

A concrete subclass of `AbstractExpression` that turns free indices into
tensor axes.

`ComponentTensor(expression, indices)` stores a scalar expression and a
non-empty `tuple` of distinct `Index`es that are free in the expression. The
edges are the expression followed by the indices, in order.
`ComponentTensor(expression, indices)` derives the value shape from
the extents of the indices, in the order of the indices, and the free indices
of the expression without the bound indices.

The free function `component_tensor(expression, indices)` returns
`ComponentTensor(expression, indices)`. Construction raises `ValueError` if the
expression is not scalar, if the indices are empty, if an index is repeated, or
if an index is not free in the expression.

### IndexSum

A concrete subclass of `AbstractExpression` for the sum of a scalar expression
over free indices.

`IndexSum(expression, indices)` stores a scalar expression and a non-empty
`tuple` of distinct `Index`es that are free in the expression. The edges are
the expression followed by the indices, in order.
`IndexSum(expression, indices)` derives the value shape `()` and the free
indices of the expression without the summation indices.

One `IndexSum` holds a whole contraction, such as the sum over both indices of
`A[i, j] * B[i, j]`. A pass that chooses the order of summation then finds the
contraction in one node, not in a chain of nested nodes.

The free function `index_sum(expression, indices)` returns
`IndexSum(expression, indices)`. Construction raises `ValueError` if the
expression is not scalar, if the indices are empty, if an index is repeated, or
if an index is not free in the expression.

Examples, with `i, j = indices("i", "j")`, two vectors `a` and `b` of value
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
but assumes dense operators.

## Annotation

> LLM generated content.

An annotation records how an object was produced or incidental metadata. An
annotation is a `tuple` of a `str` name followed by the arguments that produced
the object, for example `("metric", D)`. Every argument must be an allowed
construction argument (see [The base object](#the-base-object-base)).

> Matthew expressed dislike for MixIns, rethink?

The mixin `Annotation` holds what every annotated class shares. Each kind of
object that can be annotated has one concrete annotated class, which combines
`Annotation` with the abstract class of that kind. An annotated object is an
object of the same abstract kind as the wrapped object, so an annotated object
stands in wherever that abstract kind is accepted. An annotated object does not
stand in where a more specific class is required (see
[AnnotatedDomain](#annotateddomain)).

The `Annotation` mixin specifies:

- Storage. An annotated class stores the wrapped object and the annotation. The
  edges are the wrapped object, followed by the graph objects in the
  annotation, in order.
- Identity. An annotated object derives every attribute of the kind from the
  wrapped object.
- Equality. Equality is structural, so an annotated object is not equal to the
  wrapped object. A later algorithmic pass may compare objects while ignoring
  annotations.
- Passes. An annotation adds structure and never removes structure, so
  annotation is not a simplification. An algorithmic pass may drop any
  annotation. A pass that rewrites inside an annotated object drops the
  annotation, unless the pass shows that the annotation still holds.

The free function `annotate(x, name, *arguments)` returns the annotated class
for the abstract kind of `x`, with the annotation `(name, *arguments)`.

Every free function in the core that is defined as a composition annotates the
result with the name and arguments of the free function. For example,
`metric(D)` is `AnnotatedExpression(transpose(J) @ J, ("metric", D))`.

### AnnotatedExpression

`AnnotatedExpression(expression, annotation)` combines `Annotation` with
`AbstractExpression`. `AnnotatedExpression` derives the value shape, the free
indices, the scalar type and the set of the expression. Construction of
`Gradient` keeps an annotated expression as the edge. The later pass that
evaluates the `Gradient` of an annotated expression differentiates the wrapped
expression and drops the annotation.

### AnnotatedDomain

`AnnotatedDomain(domain, annotation)` combines `Annotation` with
`AbstractDomain`. `AnnotatedDomain` derives the dimensions, the parameter
domain, the parameterisation and the ambient domain of the domain.

An annotated domain is an `AbstractDomain` but not an
`AbstractCoordinateDomain`, even when the domain is a coordinate domain. A
coordinate function, a gradient and a parameterised domain need the plain
coordinate domain.

Examples:

- `annotate(X[0]**Integer(2) + X[1]**Integer(2), "height", X)` is an
  `AnnotatedExpression` with value shape `()`, and prints as `height(X)`.
- `annotate(parameterised_domain(RD(2), psi), "paraboloid")` is an
  `AnnotatedDomain` with $k = 2$ and $n = 3$.

## Tensor algebra

### Determinant

`Determinant` is a concrete subclass of `UnaryOperator`. The free function
`determinant(f)` takes an expression of value shape `(n, n)` and returns a
scalar. `Determinant` stores one expression, which is the only edge, and
derives the set of that expression. The determinant of an expression of value
shape `(n, n)` is a polynomial of degree $n$ in the entries, so no fixed
composition expresses the determinant, and `Determinant` is a node class.

Construction raises `ValueError` if the expression does not have rank 2, or if
the two extents of the value shape differ.

The determinant of an expression of value shape `(0, 0)` is 1, the empty
product. The measure density of a point (the Lebesgue couting measure) relies
on this value (see [Geometry](#geometry)).

### Transpose and matrix product

Transpose and matrix product are free functions built with derived indices:

- `transpose(A)` is `component_tensor(A[k, l], (l, k))`, with
  $k$ = `Index(("transpose", A, 0))` and $l$ = `Index(("transpose", A, 1))`,
  annotated with `("transpose", A)`.
- `matrix_product(A, B)`, with the overload `A @ B`, is
  `component_tensor(index_sum(A[i, k] * B[k, j], (k,)), (i, j))`, with $i$, $k$
  and $j$ the derived indices `Index(("matrix_product", A, B, m))` for
  $m = 0, 1, 2$, annotated with `("matrix_product", A, B)`.

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

## Geometry of domains

### Differentiation

#### Gradient

> LLM generated content.

A concrete subclass of `AbstractExpression` for the gradient of an expression
$f$ with respect to the coordinate function $X$ of a coordinate domain $P$ of
dimension $k$.

`Gradient(f, P)` stores $f$ and $P$, which are the two edges, in that order.
`Gradient(f, P)` derives the value shape `f.value_shape + (k,)` and the set
$P$.

Construction raises `TypeError` if $P$ is not an `AbstractCoordinateDomain`,
and raises `ValueError` if $f$ has a set other than $P$. An $f$ with no set is
a constant on $P$, and the gradient of a constant is zero.

The free function `grad(f)` returns `Gradient(f, P)`, with $P$ the set of $f$.
If $f$ has no set, `grad(f)` passes no coordinate domain, and the constructor
raises `TypeError`.

`Gradient` is lazy. Construction does not evaluate the derivative. A later pass
evaluates the derivative with the standard rules, including `grad(X)` $= I$
and a zero derivative for every expression with no set.

Differentiation on a parameterised domain is not a primitive. Differentiation
on a parameterised domain will be a composite of `grad` on the parameter domain
and the geometry of the parameterised domain (see [Geometry](#geometry)).

Examples, with `X = coordinate(RD(2))`:

- `grad(X)` has value shape `(2, 2)`. A later pass evaluates `grad(X)` to the
  identity.
- `grad(X[0]**Integer(2) + X[1]**Integer(2))` has value shape `(2,)`. A later
  pass evaluates the gradient to $(2 X_0, 2 X_1)$.
- `grad(Integer(1))` raises `TypeError`, because the literal `Integer(1)` has
  no set. `Gradient(Integer(1), RD(2))` has value shape `(2,)`, and a later
  pass evaluates the gradient to zero.

### Parameterised domains

#### ParameterisedDomain

A concrete subclass of `AbstractDomain` for parameterised domains.

A parameterised domain $D$ is the image $\psi(P)$ of a coordinate domain
$P \subset R^k$ under a parameterisation $\psi : P \to R^n$.

`ParameterisedDomain(P, psi)` stores:

- The parameter domain $P$, an `AbstractCoordinateDomain` of dimension $k$,
  returned by `parameter_domain(D)`.
- The parameterisation $\psi$, an `AbstractExpression` on $P$ with value shape
  `(n,)`, returned by `parameterisation(D)`. The parameterisation $\psi$ is
  built from the coordinate function $X$ of $P$, or is a constant with no
  set.

$P$ and $\psi$ are the two edges, in that order.

`ParameterisedDomain(P, psi)` derives:

- The topological dimension $k$, which is the dimension of $P$.
- The geometric dimension $n$, which is taken from the value shape `(n,)` of
  $\psi$. The geometric dimension is derived, not stored, so the geometric
  dimension cannot disagree with the value shape of $\psi$.

The free function `parameterised_domain(P, psi)` returns
`ParameterisedDomain(P, psi)`. Construction raises `TypeError` if $P$ is not an
`AbstractCoordinateDomain`. Construction raises `ValueError` if $\psi$ does not
have rank 1, if $\psi$ has a set other than $P$, if $\psi$ has free indices,
or if $k > n$.

Examples:

- The hypotenuse of the unit triangle is a parameterised domain with $k = 1$
  and $n = 2$. $P$ is the unit interval in $R^1$, and
  $\psi(X) = (1 - X_0, X_0)$.
- With `X = coordinate(RD(2))` and `h = X[0]**Integer(2) + X[1]**Integer(2)`,
  `parameterised_domain(RD(2), tensor((X[0], X[1], h)))` is the paraboloid,
  with $k = 2$ and $n = 3$.
  `ambient_domain` of the paraboloid is `RD(3)`.
- A physical triangle in $R^2$ is a parameterised domain with $k = n = 2$. $P$
  is the unit triangle, and $\psi$ is the affine map from $P$ onto the
  physical triangle. See the affinely mapped triangle example.
- A point $p$ in $R^2$ is a parameterised domain with $k = 0$ and $n = 2$. $P$
  is `RD(0)`, and $\psi$ is the constant `tensor((p0, p1))`, with no set.

### Geometry

The geometry of every domain $D$, coordinate or parameterised, comes from the
parameterisation $\psi$ = `parameterisation(D)` on the parameter domain $P$ =
`parameter_domain(D)` of dimension $k$:

| Free function        | Definition                                                 | Value shape |
| -------------------- | ---------------------------------------------------------- | ----------- |
| `jacobian(D)`        | $J$ = `Gradient(parameterisation(D), parameter_domain(D))` | `(n, k)`    |
| `metric(D)`          | $g$ = `transpose(J) @ J`                                   | `(k, k)`    |
| `measure_density(D)` | `square_root(determinant(g))`                              | `()`        |

For a coordinate domain, the parameterisation is the coordinate function $X$.
The general rules give $J$ = `grad(X)` $= I$, $g = I$ and a measure density of
1, so a coordinate domain needs no special case.

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

# AbstractDomain: ParameterisedDomain(RD(2), psi), k = 2, n = 3.
D = parameterised_domain(P, psi)

# AbstractExpression: AnnotatedExpression(Gradient(psi, P), ("jacobian", D)),
# value shape (3, 2).
J = jacobian(D)

# AbstractExpression: AnnotatedExpression(transpose(J) @ J, ("metric", D)),
# value shape (2, 2).
g = metric(D)

# AbstractExpression: AnnotatedExpression(SquareRoot(Determinant(g)),
# ("measure_density", D)), value shape ().
rho = measure_density(D)
```

Construction of every object in the example performs no evaluation and no
simplification. `J`, `g` and `rho` contain the lazy node `Gradient(psi, P)`,
which a later pass evaluates.

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
        super().__init__(n)
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

# AbstractDomain: ParameterisedDomain(UnitSimplex(2), psi), k = n = 2.
K = parameterised_domain(T, psi)

# AbstractExpression: AnnotatedExpression(Gradient(psi, T), ("jacobian", K)),
# value shape (2, 2).
J = jacobian(K)

# AbstractExpression: AnnotatedExpression(SquareRoot(Determinant(metric(K))),
# ("measure_density", K)), value shape ().
rho = measure_density(K)
```

A later pass evaluates the Jacobian to the constant matrix

$$J = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix},$$

so $\det J = 5$ and $\sqrt{\det g} = |\det J| = 5$. The area of the physical
triangle is $5$ times the area $1/2$ of the unit triangle.

The physical triangle is a parameterised domain with $k = n$, not a coordinate
domain, although the physical triangle has non-empty interior in $R^2$. The
class records how the domain is built.

## Integrating over the domain

At this point we have everything we need to integrate a scalar or tensor-valued
expression over a single domain. Note we do not need to introduce any
topological notions of the domain to build a basic integral.

### Integral

A domain $D$ with parameter domain $P \subset R^k$ and parameterisation $\psi$
has the measure $\mu_D$ that $\psi$ induces from the Lebesgue measure
$\lambda^k$ on $P$:

$$
\int_D f \,\mathrm{d}\mu_D = \int_P f \sqrt{\det g} \,\mathrm{d}\lambda^k,
\qquad g = J^T J, \qquad J = \nabla_P \psi,
$$

where $f$ is an expression on $P$ or a constant. $\sqrt{\det g}
\,\mathrm{d}\lambda^k$ is the volume element of $D$, and $\sqrt{\det g}$ is
`measure_density(D)`. The measure $\mu_D$ needs no orientation, because
$\sqrt{\det g} \ge 0$.

A concrete subclass of `AbstractExpression` for the integral of an expression
$f$ over a coordinate domain $P$ of dimension $k$ with the Lebesgue measure
$\lambda^k$:

$$
\mathrm{Integral}(f, P) = \int_P f \,\mathrm{d}\lambda^k.
$$

`Integral(f, P)` stores $f$ and $P$, which are the two edges, in that order.
`Integral(f, P)` derives the value shape, the free indices and the scalar type
of $f$. `Integral(f, P)` has no set.

Construction raises `TypeError` if $P$ is not an `AbstractCoordinateDomain`,
and raises `ValueError` if $f$ has a set other than $P$. An $f$ with no set is
a constant on $P$.

The free function `integral(f, P)` returns `Integral(f, P)`.

`Integral` is lazy. Construction does not evaluate the integral. A later pass
might evaluate the integral, for example by quadrature.

Note that an `Integral` is an expression, so integrals can be added, multiplied
and used wherever an expression is accepted.

### integrate

The free function `integrate(f, D)` integrates $f$ over any domain $D$ with the
measure $\mu_D$. `integrate(f, D)` is
`integral(f * measure_density(D), parameter_domain(D))`, annotated with
`("integrate", f, D)`. The free function is the definition above, written as
a fixed composition.

### Examples

With `T = UnitSimplex(2)`, `X = coordinate(T)`, and the physical triangle `K`
of the affinely mapped triangle example:

- `integral(X[0], T)` is $\int_T X_0 \,\mathrm{d}\lambda^2$. A later pass
  evaluates the integral to $1/6$.
- `integrate(Integer(1), K)` is the area of the physical triangle. The measure
  density of `K` is $\lvert \det J \rvert = 5$, so a later pass evaluates the
  integral to $5/2$.
- With `h = X[0]**Integer(2) + X[1]**Integer(2)` and
  `S = parameterised_domain(T, tensor((X[0], X[1], h)))`,
  `integrate(Integer(1), S)` is the area of the paraboloid over the unit
  triangle. The measure density of `S` is $\sqrt{1 + 4 X_0^2 + 4 X_1^2}$.
- With `c = RealScalar(0.5)` and
  `p = parameterised_domain(RD(0), tensor((RealScalar(1.0), RealScalar(2.0))))`,
  `integrate(c, p)` is `c`. The measure density of the point `p` is 1, and the
  Lebesgue measure on $R^0$ is the counting measure at the single point of
  $R^0$.
- `integrate(X[0], K) * integrate(X[1], K)` is a product of two integrals, with
  value shape `()` and no set.

Integrating over a point works only for an expression with no set. An
expression `f` on another coordinate domain `Q`, such as a function on a
triangle that contains the point, must first be composed with the
parameterisation of the point, giving $f \circ \psi$ on `RD(0)`. At this point
in the document the core has no composition yet.

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

## Target language

> The following is LLM generated. It should be updated and inspected as the
> language progresses and can be used for formal proofs later.

The target language states, without reference to UFLx, which mathematical
expressions the core can write. The core is complete if every well-typed
expression of the target language translates to an object of the core with the
same meaning.

### Target grammar

A target expression $e$ has a shape $s$, a set (a coordinate domain $P$, or
none), a scalar type (real or complex), and free indices.

$$
\begin{aligned}
e \;::=\;\; & c && \text{literal: } n,\ p/q,\ r,\ r + \mathrm{i}\, r' \\
  \mid\;\; & \ell && \text{a scalar leaf from } \Sigma \\
  \mid\;\; & X_P && \text{the coordinate function of } P \\
  \mid\;\; & e_1 + e_2 \;\mid\; e_1 e_2 \;\mid\; e_1 / e_2 \;\mid\; e_1^{e_2} && \text{arithmetic} \\
  \mid\;\; & \sqrt{e} \;\mid\; \lvert e \rvert \;\mid\; \det e && \\
  \mid\;\; & \nabla_P e && \text{gradient with respect to } X_P \\
  \mid\;\; & \textstyle\int_P e \,\mathrm{d}\lambda && \text{Lebesgue integral over } P \\
  \mid\;\; & [e_1, \ldots, e_m] && \text{tensor from entries} \\
  \mid\;\; & e_{a_1 \ldots a_r} && a_j \text{ an integer or an index symbol} \\
  \mid\;\; & [e]_{i_1 \ldots i_m} && \text{the tensor whose } (i_1, \ldots, i_m) \text{ component is } e \\
  \mid\;\; & \textstyle\sum_{i_1 \ldots i_m} e && \text{sum over index symbols} \\[4pt]
D \;::=\;\; & P && \text{a coordinate domain: } \mathbb{R}^n \text{, or a leaf from } \Sigma \\
  \mid\;\; & \psi(P) && \psi \text{ an expression on } P \text{ of shape } (n,)
\end{aligned}
$$

### Core grammar

The node classes of the core form the abstract syntax of UFLx:

```
Expression ::= Integer | Rational | RealScalar | ComplexScalar
             | <leaf subclass of AbstractScalar>
             | Coordinate(CoordinateDomain)
             | Sum(Expression, Expression) | Product(Expression, Expression)
             | Division(Expression, Expression) | Power(Expression, Expression)
             | SquareRoot(Expression) | Absolute(Expression)
             | Determinant(Expression) | Gradient(Expression, CoordinateDomain)
             | Integral(Expression, CoordinateDomain)
             | Tensor(entries) | Indexed(Expression, indices)
             | ComponentTensor(Expression, indices) | IndexSum(Expression, indices)
             | AnnotatedExpression(Expression, annotation)

Domain     ::= CoordinateDomain
             | ParameterisedDomain(CoordinateDomain, Expression)
             | AnnotatedDomain(Domain, annotation)

CoordinateDomain ::= RD(n) | <leaf subclass of AbstractCoordinateDomain>
```

### Translation

Each row gives the condition under which a target expression is well typed, and
the node class that the target expression translates to. The condition of each
row is exactly the condition that the constructor of the node class checks.

| Target expression              | Well typed when                                                                  | Shape                | Node class                          |
| ------------------------------ | -------------------------------------------------------------------------------- | -------------------- | ----------------------------------- |
| $n$, $p/q$, $r$                | $q \ne 0$, $r$ finite                                                            | $()$                 | `Integer`, `Rational`, `RealScalar` |
| $r + \mathrm{i}\, r'$          | $r$, $r'$ finite                                                                 | $()$                 | `ComplexScalar`                     |
| $\ell$                         | always                                                                           | from the leaf        | the leaf class                      |
| $X_P$                          | always                                                                           | $(n,)$               | `Coordinate`                        |
| $e_1 + e_2$                    | equal shapes, equal free indices, compatible sets                                | $s$                  | `Sum`                               |
| $e_1 e_2$                      | one operand scalar, compatible sets                                              | $s$                  | `Product`                           |
| $e_1 / e_2$                    | $e_2$ scalar, compatible sets                                                    | $s$                  | `Division`                          |
| $e_1^{e_2}$                    | both scalar, $e_2$ has no set                                                    | $()$                 | `Power`                             |
| $\sqrt{e}$, $\lvert e \rvert$  | $e$ scalar                                                                       | $()$                 | `SquareRoot`, `Absolute`            |
| $\det e$                       | shape $(m, m)$                                                                   | $()$                 | `Determinant`                       |
| $\nabla_P e$                   | $P$ of dimension $k$, $e$ has no set or the set $P$                              | $s + (k,)$           | `Gradient`                          |
| $\int_P e \,\mathrm{d}\lambda$ | $e$ has no set or the set $P$                                                    | $s$, no set          | `Integral`                          |
| $[e_1, \ldots, e_m]$           | entries scalar, rectangular, $m \ge 1$, compatible sets                          | $(m, \ldots)$        | `Tensor`                            |
| $e_{a_1 \ldots a_r}$           | $r$ equals the rank, $0 \le a_j < n_j$, no index symbol repeated or already free | $()$                 | `Indexed`                           |
| $[e]_{i_1 \ldots i_m}$         | $e$ scalar, $m \ge 1$, indices distinct and free in $e$                          | extents of the $i_j$ | `ComponentTensor`                   |
| $\sum_{i_1 \ldots i_m} e$      | $e$ scalar, $m \ge 1$, indices distinct and free in $e$                          | $()$                 | `IndexSum`                          |
| $\psi(P)$                      | $\psi$ on $P$, shape $(n,)$, no free indices, $k \le n$                          | a domain             | `ParameterisedDomain`               |

"Compatible sets" means that every operand with a set has the same set.

Every other construction of the core is a definition in terms of the rows
above. For example, $-e$ is $(-1)\, e$, $A^T$ is $[A_{kl}]_{lk}$, and $AB$ is
$[\sum_k A_{ik} B_{kj}]_{ij}$. An annotation has no counterpart in the target
language, because an annotation does not change the meaning of an object.

# JSH worked to here.


## 5. What belongs to the language

The minimal language includes:

- symbolic domains/function spaces,
- symbolic terminals,
- tensor-valued expressions,
- index notation,
- differential operators,
- integration measures,
- integrals,
- and forms.

This is enough to describe ordinary weak forms such as mass, stiffness,
elasticity, Stokes, mixed Poisson, and '.

## 6. What does not belong to the minimal language

The following are not primitive language concepts:

- actual mesh connectivity,
- actual coefficient vectors,
- basis tabulation,
- quadrature-loop generation,
- form compilation,
- assembly,
- solver setup,
- boundary-condition application,
- expression simplification strategy,
- code generation.

These are consumers or transformations of the symbolic language, not the core
language itself.

## 8. Examples

### Poisson equation

```python
coordinate_element = FiniteElement("Lagrange", triangle, 1, shape=(2,))
element = FiniteElement("Lagrange", triangle, 1)

mesh = Mesh(coordinate_element)
V = FunctionSpace(mesh, element)

u = TrialFunction(V, name="u")
v = TestFunction(V, name="v")
f = Coefficient(V, name="f")

a = inner(grad(v), grad(u)) * dx
L = v * f * dx
```

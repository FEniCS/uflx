# Core UFLx language definition

This note describes the smallest core language behind UFLx. It intentionally
avoids implementation details, compiler algorithms, and broad redesign ideas.
It assumes knowledge of differentiable geometry and functional analysis. Some
advanced differential geometry concepts are explained in more detail.

This document is aimed at users and developers with extensive knowledge of
domain specific languages for variational forms (e.g. UFL, ngsolve).

## What UFLx is and isn't

In one sentence:

> UFLx is a typed symbolic language writing variational forms.

UFLx is a symbolic language embedded in Python for writing variational forms.

The UFLx core language (henceforth the core) should not use finite element or
numerical analysis terminology - differential geometry and functional analysis
,and their intersection, without finite element (numerical analysis) concepts,
should be all that should is required.

Finite element and numerical analysis concepts, when necessary, will be added
on top of the core, or provided by external libraries (Basix being the
'in-house' example). Basix and DOLFINx and their assumptions should not be
visible in the core.

The core requires an understanding of modern differential geometry and
functional analysis. UFLx will then provide higher-level functionality based on
this core which mimics the 'easier' interface presented by other competing
packages.

UFLx is *not* a symbolic language for reasoning about problems in functional
analysis or differential geometry. Hence, the core should only contain the
minimal concepts and relations necessary for writing a wide range of
variational forms as code for subsequent algorithmic transformation.

A UFLx program is ordinary Python code that constructs symbolic objects and
connects them through a graph. The primary result is usually a Form, such as a
bilinear form, linear form, or functional. Python code will be used to motivate
examples but a full concrete implementation is out of scope.

UFLx does not handle concrete data or numerical algorithms; e.g. actual meshes,
boundary markers, numerical values of coefficients, or assembly procedures.
Those belong to the surrounding problem-solving environment, such as DOLFINx,
Firedrake, FIAT and Basix. UFLx instead provides a specified set of leaf nodes
which the external library must provide.

## Basic principles

- All core object definitions should be closely related to their strict
  mathematical definition.
- Extra definitions, properties or terminology outside of the established
  mathematical fields noted must be avoided, in particular ideas from more
  applied fields such as mechanics and engineering-oriented computational
  sciences.
- Strictness and expressivity, not brevity, is key in the core - brevity can be
  added later with higher-level APIs.
- The mathematical fields provide a rich foundation for writing variational
  problems in a uniform way across geometric dimensions, i.e. 'special cases'
  for particular dimensions do not exist, and so should not exist in the core.
- The core should not contain any algorithm that is designed to produce
  convenient output for form compilers.

## Basic implementation principles

- Explicit is better than implicit.
- Boilerplate and error-prone patterns should be minimised.
- UFLx should deterministic on the same machine using the same byte-code.
- Performance is *important* and so any data structures and algorithms in the
  core should be carefully chosen and tested for scalability.
- The minimal language should only need a relatively small number of object
  kinds. When designing the language, consider whether objects can be composed,
  rather than extended.
- Free functions should be preferred over class members. Class members should
  immediately call free functions wherever possible. Free functions maximise
  the chance of code re-use and are easier to reason about.
- A parent class should avoid returning a subclass.
- A concept is a free function, not a node class, if the free function can
  build the concept as one fixed composition of existing nodes applied to the
  arguments, *without inspecting the arguments to choose the composition*. A
  free function may inspect the arguments to validate the arguments. Otherwise
  the concept gets its own node class. For example, the metric is always
  `transpose(J) @ J`, so the metric is a free function. No composition of
  existing nodes differentiates, so `Gradient` is a node class. A free
  function that builds a fixed composition is a definition, not a
  simplification. The free function does not inspect the operands, so the free
  function never removes structure that the user built. The free function
  annotates the result with the name and arguments of the free function (see
  [Annotation](#annotation)), so the graph stays readable.
- Under the previous rule, a node class is likely to become a free function as
  UFLx grows, once the node class can be written as a fixed composition of
  newer nodes *without inspecting the arguments to choose the composition*. For
  example, if subtraction were a node class `Difference`, `subtract(f, g)`
  should become the free function `add(f, negate(g))` once `negate` exists. The
  composition is the same for every `f` and `g`, so the composition is allowed.
  In contrast, `transpose(A)`, assuming no `Indexing`-related functionality,
  cannot be:

  ```python
  m, n = A.value_shape
  tensor(tuple(tuple(A[j, i] for j in range(m)) for i in range(n)))
  ```

  because the composition reads the value shape of `A` to choose the entries.
  With free indices, `transpose(A)` is the fixed composition
  `component_tensor(A[i, j], (j, i))`, so `transpose` is a free function.
  `Gradient` is a primitive, and `Gradient` stays a node class.
- Simplifications and rewriting should be avoided at all costs during user
  graph construction, maintaining the user's intent. Only later controllable
  algorithmic passes should re-write. A good example is `(a + b) + c` in
  floating point arithmetic where `a`, `b` and `c` are unknown terminals (i.e.
  to be provided later. The user may have very good reason to have `(a + b)`
  performed before adding the answer to `c`.

## Tradeoffs of the above decisions

Good things:

- Much smaller algebraic lowering pass.
- Fewer node types. Every algorithm, including differentiation, will handle a
  smaller and more stable set of classes.
- Repeated construction gives equal, shared graphs, for example
  `transpose(A) == transpose(A)`. Spellings with different bound-index labels
  are equal only up to renaming, which a later pass decides.
- Graphs stay readable. Every composition is wrapped in an `AnnotatedExpression` node
  that records the name and arguments of the free function, so a printer can
  show `metric(D)` instead of the expanded composition (see
  [Annotation](#annotation)).

Bad things:

- An annotated result is not structurally equal to the same composition
  written by hand. Comparison that ignores annotations is a later pass.
- Graph size. *Requires highly performant graph/hashing implementation*.
- Shape-dependent operations stay nodes (until lowering). For example,
  `Determinant` and `Inverse` cannot be lowered without inspecting the value
  shape. 

## The minimal object model

### The base object `Base`

The base object in UFLx must be extended by all objects that are included in
the graph. The graph is a directed acyclic graph (DAG).

Everything needed to define the variational problem must be recoverable from
the graph. No information is held outside of the graph.

Objects must satisfy the following rules:

1. All construction arguments must be deeply immutable values. For example, a
   `tuple` of `list`s is not allowed.
1. All construction arguments must be compared by type and value, so that `1`,
   `1.0` and `True` are distinct, `0.0` and `-0.0` are also distinct, and every
   value is equal (`==` or `is`, depending on implementation) to itself.
1. Every construction argument must have an identity (hash) that is
   deterministic across processes on the same machine.
1. Construction must depend only on the construction arguments, not on global
   state.
1. An object must be reconstructible from its class and construction arguments,
   and the reconstructed object must be equal to the original.
1. Any results produced by the class after construction must be a result of
   only the class and the construction arguments.
1. The state of the class must be immutable after construction.
1. Any ordering of objects must come from a total order that is deterministic
   across processes.
1. The edges of an object must be exactly the base objects among its
   construction arguments, in order.
1. Equality (`==`) must be structural and return a `bool`. Equality does not
   decide mathematical equivalence, so $a + b$ and $b + a$ are not equal.
   Mathematical equivalencies are applied in later graph algorithms.

Hence two objects are equal if and only if the two objects have the same class
and equal construction arguments.

To satisfy the above, the base object must only accept a small curated set of
types.

Because objects are immutable, an object need only be validated once, at
construction. Every object that exists is therefore valid.

An object cannot be changed. To obtain a modified object, construct a new
object and reconstruct every object that depends on the original object, up to
the root of the graph.

The above design implies that construction arguments should be small and cheap
to hash and compare, and the graph implementation must be performant.

The base object also fixes two Python protocols for every object in the graph:

- Truth testing raises `TypeError`. An object has no truth value, so
  `if a == b` cannot be mistaken for a test of mathematical equality, and
  `if x` cannot replace `if x is not None`.
- The base object sets `__array_ufunc__ = None`. NumPy then defers to the
  object rather than broadcasting over the object, so `np.float64(0.5) * X`
  raises `TypeError` and comparison with a NumPy array does not become
  element-wise.

### Enforcement

Enforcement must be implemented statically and dynamically - look at tools such as
`dataclass`, deterministic hashing, hash-consing, `WeakValueDictionary`, `ruff`
and `semgrep`.

https://github.com/semgrep/semgrep

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
│   ├── IndexSum                          sums over one bound index
│   ├── AnnotatedExpression               Annotation mixin, identity
│   ├── UnaryOperator
│   │   ├── Gradient                      grad, on coordinate domains only
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
└── AbstractDomain                        parameter_domain(D), mapping(D)
    ├── AnnotatedDomain                   Annotation mixin, identity
    ├── AbstractCoordinateDomain          stores n, mapping is X
    │   └── RD                            R^n
    └── ParameterisedDomain               stores P and ψ, mapping is ψ
```

The edges of the graph point from an object to the objects that the object
stores, and never back:

```
ParameterisedDomain ──► P (AbstractCoordinateDomain)
        │
        └─────────────► ψ (AbstractExpression) ──► ... ──► Coordinate ──► P

Gradient ──► f        Indexed ──► expression        Tensor ──► entries
```

### Free functions

Free functions are not objects in the graph. A free function builds or returns
objects in the graph. Operator overloads call the free functions.

> LLM generated table.

| Group                 | Free function                                           | `AbstractExpression` overload | Takes                                                                        | Returns                                                       |
| --------------------- | ------------------------------------------------------- | ----------------------------- | ---------------------------------------------------------------------------- | ------------------------------------------------------------- |
| Arithmetic            | [`add(f, g)`](#arithmetic)                              | `+`                           | `f, g: AbstractExpression`                                                   | `Sum`                                                         |
|                       | [`subtract(f, g)`](#arithmetic)                         | `-`                           | `f, g: AbstractExpression`                                                   | `AnnotatedExpression`: `Sum` of `f` and `negate(g)`           |
|                       | [`multiply(f, g)`](#arithmetic)                         | `*`                           | `f, g: AbstractExpression`                                                   | `Product`                                                     |
|                       | [`divide(f, g)`](#arithmetic)                           | `/`                           | `f, g: AbstractExpression`                                                   | `Division`                                                    |
|                       | [`negate(f)`](#arithmetic)                              | unary `-`                     | `f: AbstractExpression`                                                      | `AnnotatedExpression`: `Product` of `-1` and `f`              |
|                       | [`absolute(f)`](#arithmetic)                            | `abs`                         | `f: AbstractExpression`                                                      | `Absolute`                                                    |
|                       | [`power(f, g)`](#arithmetic)                            | `**`                          | `f, g: AbstractExpression`                                                   | `Power`                                                       |
| Elementary functions  | [`square_root(f)`](#elementary-functions)               |                               | `f: AbstractExpression`                                                      | `SquareRoot`                                                  |
| Domains               | [`ambient_domain(D)`](#abstractcoordinatedomain)        |                               | `D: AbstractDomain`                                                          | `RD(n)`                                                       |
|                       | [`parameter_domain(D)`](#abstractdomain)                |                               | `D: AbstractDomain`                                                          | `AbstractCoordinateDomain`                                    |
|                       | [`mapping(D)`](#abstractdomain)                         |                               | `D: AbstractDomain`                                                          | `AbstractExpression` with shape `(n,)`                        |
|                       | [`coordinate(P)`](#coordinate)                          |                               | `P: AbstractCoordinateDomain`                                                | `Coordinate`                                                  |
| Tensors and indexing  | [`tensor(entries)`](#tensor)                            |                               | `entries: tuple` (nested, of `AbstractExpression`)                           | `Tensor`                                                      |
|                       | [`vector(entries)`](#tensor)                            |                               | `entries: tuple[AbstractExpression, ...]`                                    | `Tensor` of rank 1                                            |
|                       | [`matrix(entries)`](#tensor)                            |                               | `entries: tuple[tuple[AbstractExpression, ...], ...]`                        | `Tensor` of rank 2                                            |
|                       | [`indexed(f, indices)`](#indexed)                       | `f[i]`                        | `f: AbstractExpression`, `indices: int \| Index \| tuple[int \| Index, ...]` | `Indexed`                                                     |
| Index notation        | [`index(label)`](#index)                                |                               | `label: str \| tuple`                                                        | `Index`                                                       |
|                       | [`indices(*labels)`](#index)                            |                               | `*labels: str \| tuple`                                                      | `tuple` of `Index`                                            |
|                       | [`component_tensor(f, indices)`](#componenttensor)      |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`                        | `ComponentTensor`                                             |
|                       | [`index_sum(f, i)`](#indexsum)                          |                               | `f: AbstractExpression`, `i: Index`                                          | `IndexSum`                                                    |
| Annotation            | [`annotate(x, name, *arguments)`](#annotation)          |                               | `x: AbstractExpression \| AbstractDomain`, `name: str`, `*arguments`         | `AnnotatedExpression` or `AnnotatedDomain`                    |
| Tensor algebra        | [`determinant(f)`](#determinant)                        |                               | `f: AbstractExpression`                                                      | `Determinant`                                                 |
|                       | [`transpose(f)`](#transpose-and-matrix-product)         |                               | `f: AbstractExpression`                                                      | `AnnotatedExpression`: `ComponentTensor`                      |
|                       | [`matrix_product(f, g)`](#transpose-and-matrix-product) | `@`                           | `f, g: AbstractExpression`                                                   | `AnnotatedExpression`: `ComponentTensor` of an `IndexSum`     |
| Differentiation       | [`grad(f)`](#gradient)                                  |                               | `f: AbstractExpression`                                                      | `Gradient`                                                    |
| Parameterised domains | [`parameterised_domain(P, psi)`](#parameteriseddomain)  |                               | `P: AbstractCoordinateDomain`, `psi: AbstractExpression`                     | `ParameterisedDomain`                                         |
| Geometry              | [`jacobian(D)`](#geometry)                              |                               | `D: AbstractDomain`                                                          | `AnnotatedExpression`: `Gradient` of `mapping(D)`             |
|                       | [`metric(D)`](#geometry)                                |                               | `D: AbstractDomain`                                                          | `AnnotatedExpression`: `transpose(J) @ J`                     |
|                       | [`measure_density(D)`](#geometry)                       |                               | `D: AbstractDomain`                                                          | `AnnotatedExpression`: `SquareRoot` of the determinant of $g$ |

Every overload in the table is a method of `AbstractExpression`, together with
the reflected method where one exists, such as `__add__` and `__radd__` for
`+` (see [`AbstractExpression`](#abstractexpression)). No other class defines
operator overloads.

The value-shape and other constraints on the arguments are stated in the
section of each free function. A Python `int` is accepted wherever an
expression is (see
[Python integers in expressions](#python-integers-in-expressions)).

## The lowest objects

### Expressions

An expression denotes an array-valued function defined on a set of the
expression. An expression with no set denotes a single array. An
expression can be real or complex valued.

#### `AbstractExpression`

Every `AbstractExpression` derives the value shape, the free indices and the
set of the expression. The set is absent for an expression with no set, such
as a literal.

An operator overload returns `NotImplemented` for an operand type that the
overload does not accept, rather than raising an error. Python then tries the
reflected overload of the other operand and raises `TypeError` if neither
overload accepts the operands.

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

    # Operator overloading must be in the class.
    def __add__(self, other: object) -> AbstractExpression:
        """Add."""
        # Calls the free function arithmetic.add.
        ...

    def __radd__(self, other: object) -> AbstractExpression:
        """Add."""
        # Calls add(other, self), keeping the order that the user wrote.
        ...

    def __sub__(self, other: object) -> AbstractExpression:
        """Subtract."""
        # Calls the free function arithmetic.subtract.
        ...

    def __rsub__(self, other: object) -> AbstractExpression:
        """Subtract."""
        # Calls subtract(other, self), keeping the order that the user wrote.
        ...

    def __mul__(self, other: object) -> AbstractExpression:
        """Multiply."""
        # Calls the free function arithmetic.multiply.
        ...

    def __rmul__(self, other: object) -> AbstractExpression:
        """Multiply."""
        # Calls multiply(other, self), keeping the order that the user wrote.
        ...

    def __matmul__(self, other: object) -> AbstractExpression:
        """Matrix multiply."""
        # Calls the free function tensor_algebra.matrix_product.
        ...

    def __truediv__(self, other: object) -> AbstractExpression:
        """Divide."""
        # Calls the free function arithmetic.divide.
        ...

    def __rtruediv__(self, other: object) -> AbstractExpression:
        """Divide."""
        # Calls divide(other, self), keeping the order that the user wrote.
        ...

    def __neg__(self) -> AbstractExpression:
        """Negate."""
        # Calls the free function arithmetic.negate.
        ...

    def __abs__(self) -> AbstractExpression:
        """Absolute value."""
        # Calls the free function arithmetic.absolute.
        ...

    def __pow__(self, exponent: object) -> AbstractExpression:
        """Raise to a power."""
        # Calls the free function arithmetic.power.
        ...

    def __rpow__(self, base: object) -> AbstractExpression:
        """Raise a Python int to a power."""
        # Calls power(base, self).
        ...

    def __getitem__(self, indices: object) -> AbstractExpression:
        """Index."""
        # Calls the free function indexed.
        ...

    def __iter__(self) -> NoReturn:
        """Raise TypeError, because expressions are not iterable."""
        # Without __iter__, Python iterates through __getitem__.
        ...

# arithmetic.py

# A free function
def add(
    first: AbstractExpression | int, second: AbstractExpression | int
) -> AbstractExpression:
    """Add two expressions.

    Args:
        first: The first summand. A Python int becomes an Integer.
        second: The second summand. A Python int becomes an Integer.

    Returns:
        The sum of the two expressions, with the summands in the given order.

    Raises:
        TypeError: If a summand is neither an expression nor an int, or is a bool.
        ValueError: If the value shapes or the free indices of the summands differ.
    """
    # Converts ints, then checks value_shape and free_indices.
    return Sum(first, second)
```

### Literals

A literal is a number written down in an expression, e.g. `Integer(2) * u`. A
literal's value is part of the expression itself.

> Literals are Expressions, but they do not depend on a set (no domain).

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
- `Rational(p, q)` takes two `int`s $p$ and $q$ with $q \ne 0$.
- `RealScalar(v)` takes a finite `float` $v$. NaN and infinities are rejected,
  because a NaN is not equal to itself. `RealScalar(-0.0) != RealScalar(0.0)`.
- `ComplexScalar(re, im)` takes two `RealScalar`s.

Literals should not fold on construction, e.g.

```
Integer(2) + Integer(3) != Integer(5)
Integer(2) + Integer(3) == Sum(Integer(2), Integer(3))
```

#### Python integers in expressions

>>> LLM generated content.

Wherever an operator overload or a free function accepts an expression, the
operator overload or free function also accepts a Python `int`, which becomes
`Integer`. A `bool` is rejected with `TypeError`, although `bool` is a subclass
of `int`. Class constructors never convert.

Every other number must be written as a literal, such as `RealScalar(0.5)`,
`Rational(1, 3)` or `ComplexScalar(RealScalar(0.0), RealScalar(1.0))`. A
`float`, `complex`, `Fraction` or NumPy scalar is rejected with `TypeError`.
Only `int` converts because only the conversion of `int` is unambiguous. The
type of every other literal is a choice that the user makes explicitly.

Examples, with `X = coordinate(RD(2))`:

- `X[0]**2`, `2 * X[0]` and `1 - X[0]` are accepted.
- `0.5 * X[0]` raises `TypeError`. `RealScalar(0.5) * X[0]` is accepted.
- `vector((1, 2))` is `Tensor((Integer(1), Integer(2)))`.

### Operators

```UnaryOperator(AbstractExpression)```
```BinaryOperator(AbstractExpression)```

Downstream basic ops are relatively obvious.

### Arithmetic

Concrete subclasses of `BinaryOperator`. Each class stores two expressions,
which are the two edges in order. The set is the set of the operands. Operands
with a set must have the same set. Each free function accepts Python `int`s as
operands, raises `ValueError` for wrong value shapes and performs no
simplification.

| Class      | Free function | Overload | Operand shapes | Result shape |
| ---------- | ------------- | -------- | -------------- | ------------ |
| `Sum`      | `add`         | `+`      | `s`, `s`       | `s`          |
| `Product`  | `multiply`    | `*`      | `()`, `s`      | `s`          |
|            |               |          | `s`, `()`      | `s`          |
| `Division` | `divide`      | `/`      | `s`, `()`      | `s`          |
| `Power`    | `power`       | `**`     | `()`, `()`     | `()`         |

At least one operand of `Product` is scalar. The free indices of every operator
follow the general rule of [Free indices](#free-indices).

`Absolute` is a concrete subclass of `UnaryOperator`. The free function
`absolute(f)`, with the overload `abs(f)`, takes a scalar and returns a scalar.
For a complex-valued argument, `Absolute` is the modulus and is real-valued.

Negation and subtraction are free functions, not node classes:

- `negate(f)`, with the overload `-f`, is `multiply(-1, f)`, annotated with
  `("negate", f)`.
- `subtract(f, g)`, with the overload `f - g`, is `add(f, negate(g))`,
  annotated with `("subtract", f, g)`.

Both compositions are exact in floating-point arithmetic, so neither
composition changes the arithmetic that the user wrote. Division is a node
class, because `multiply(f, power(g, -1))` rounds twice and differs from
`f / g` in floating-point arithmetic.

### Elementary functions

`SquareRoot` is a concrete subclass of `UnaryOperator`. The free function
`square_root(f)` takes a scalar and returns a scalar. Like every unary
operator, `SquareRoot` stores one expression, which is the only edge, derives
the set of that expression and performs no simplification.

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

Defined by:

- The ambient dimension $n$ of the Euclidean space $R^n$ whose coordinates name the
  points of this domain (geometric dimension). 
- The topological dimension, which is the maximum dimension of the domains
  cells (standard Lebesgue definition). 

Examples:
- A domain consisting of triangle and interval domain has topological dimension
  2.
- A sphere in $R^3$ has a geometric dimension 3 and topological dimension 2.

Every `AbstractDomain` $D$ with topological dimension $k$ and geometric
dimension $n$ provides two free functions:

- `parameter_domain(D)`, an `AbstractCoordinateDomain` $P$ of dimension $k$.
- `mapping(D)`, an `AbstractExpression` $\psi$ on $P$ with value shape `(n,)`.
  The mapping $\psi$ maps each point of $P$ to a point of $D$, and the image
  $\psi(P)$ is $D$.

The mapping of $D$ gives the coordinates in $R^n$ of the points of $D$. The
geometry of $D$ is derived from the mapping (see
[Geometry](#geometry)).

#### AbstractCoordinateDomain

An abstract base class for coordinate domains `AbstractCoordinateDomain`, and a
subclass of `AbstractDomain`.

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$.
Hence the topological dimension of a coordinate domain equals the geometric
dimension $n$. The main purpose of this class over `AbstractDomain` is to
enforce that the topological dimension equals the geometric dimension.

Every concrete `AbstractCoordinateDomain` $P$ stores:

- The dimension $n$, provided by the property `dimension`.

An `AbstractCoordinateDomain` $P$ derives:

- The geometric dimension and the topological dimension, which are both $n$.
- The ambient domain `ambient_domain(P)`, which is `RD(n)`.
- The coordinate function `coordinate(P)` (see [Coordinate](#coordinate)).
- The parameter domain `parameter_domain(P)`, which is $P$.
- The mapping `mapping(P)`, which is the identity map `coordinate(P)`.

Examples:

- A reference cell in $R^n$ and $R^n$ itself are coordinate domains.
- A facet of a reference cell in $R^n$ is *not* a coordinate domain, because the
  facet has empty interior in $R^n$. The facet is the image of a reference cell
  in $R^{n-1}$ under a map, so the facet is a parameterised domain (defined
  below).

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

The coordinate function $X$ is a function on $P$, not a point of $P$.
Evaluating $X$ at a point $p$ of $P$ gives the $n$-tuple of $p$.

`Coordinate(P)` stores the coordinate domain $P$, which is the only edge.
`Coordinate(P)` derives the value shape `(n,)` and the set $P$. Two coordinate
functions are equal if and only if the two coordinate domains are equal.

Only coordinate domains provide a coordinate function. The coordinates of the
points of any other domain are given by `mapping` (see
[AbstractDomain](#abstractdomain)), which is built from the coordinate function
of the parameter domain.

Examples:

- `coordinate(RD(3))` is `Coordinate(RD(3))`, with value shape `(3,)` and the
  set `RD(3)`.
- `coordinate(RD(3)) == coordinate(RD(3))` is `True`, and
  `coordinate(RD(2)) == coordinate(RD(3))` is `False`.
- In the paraboloid example, `X = coordinate(RD(2))` has value shape `(2,)`.

### Tensors and indexing

#### Tensor

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
  entries with a set must have the same set. If no entry has a set, the tensor
  has no set.

The edges of `Tensor(entries)` are the entries in row-major order, so the
entry `entries[i_0]...[i_{r-1}]` comes before every later entry in
lexicographic order of the indices $(i_0, \ldots, i_{r-1})$.

Construction performs no simplification. For example,
`vector((Integer(1) + Integer(1),)) != vector((2,))`.

The free function `tensor(entries)` accepts Python `int`s as entries and
returns `Tensor(entries)`. The free functions `vector` and `matrix` are
`tensor` restricted to nesting depth 1 and nesting depth 2, and raise
`ValueError` for any other nesting depth.

Examples:

- `vector((1, 2))` has value shape `(2,)` and no set.
- `matrix(((1, 0), (0, 1)))` has value shape `(2, 2)` and no set. The entries
  are `Integer`s.
- With `one = RealScalar(1.0)` and `zero = RealScalar(0.0)`,
  `matrix(((one, zero), (zero, one)))` has `RealScalar` entries, so the matrix
  is not equal to `matrix(((1, 0), (0, 1)))`.
- `vector((0.5, 1))` raises `TypeError`.
- `matrix(((1, 2),))` has value shape `(1, 2)`.
- `tensor((((1,),),))` has value shape `(1, 1, 1)`.

#### Indexed

> LLM generated content.

A concrete subclass of `AbstractExpression` that selects one component of a
tensor-valued expression.

`Indexed(expression, indices)` stores:

- The expression, an `AbstractExpression` of rank $r \ge 1$. The expression is
  the only edge of `Indexed(expression, indices)`.
- The indices, a `tuple` of $r$ entries. Each entry is an `int` or an `Index`
  (see [Index notation](#index-notation)).

Indexing is full indexing only. The number of indices must equal the rank $r$
of the expression. Each `int` entry $i_j$ must satisfy $0 \le i_j < n_j$, where
`(n_0, ..., n_{r-1})` is the value shape of the expression. Each `Index` entry
$i_j$ takes the extent $n_j$.

`Indexed(expression, indices)` derives the value shape `()`, the set of the
expression, and the free indices of the expression followed by the `Index`
entries.

The free function `indexed(expression, indices)` accepts a single entry, which
the free function converts into a `tuple` of one entry, or a `tuple` of
entries, and returns `Indexed(expression, indices)`. The overload
`expression[indices]` calls `indexed`, so `X[0]` and `A[0, 1]` are
`indexed(X, 0)` and `indexed(A, (0, 1))`. The free function `indexed` raises:

| Error        | Cause                                                        |
| ------------ | ------------------------------------------------------------ |
| `TypeError`  | An entry that is neither exactly `int` nor `Index`.          |
| `TypeError`  | An expression with value shape `()`.                         |
| `IndexError` | A number of indices that differs from the rank.              |
| `IndexError` | A negative `int` entry, or an `int` entry $i_j \ge n_j$.     |
| `ValueError` | An `Index` that occurs twice, or that is already free.       |

A negative index is rejected rather than converted, so the stored indices are
always the indices that the user wrote. Expressions are not iterable (see
`AbstractExpression.__iter__`).

Construction performs no simplification. For example, `vector((a, b))[0]` is
`Indexed(Tensor((a, b)), (0,))` and is not equal to `a`.

Partial indexing and slicing are not part of `Indexed`. A repeated `Index`
does not imply summation (see [Index notation](#index-notation)).

Examples:

- For the coordinate function $X$ of `RD(3)`, `X[2]` has value shape `()`.
- For `A = matrix(((1, 2), (3, 4)))`, `A[1, 0]` has value shape `()`, and
  `A[1]` raises `IndexError`.
- `X[-1]` and `X[True]` raise errors.
- For the coordinate function $X$ of `RD(2)`, `Tensor((X[0], X[1])) != X`,
  because construction performs no simplification.
- `vector((1 - X[0], X[0]))` has value shape `(2,)` and the set `RD(2)`.
- `matrix(((X[0], X[1]),))` has value shape `(1, 2)`.
- In the paraboloid example, the height `X[0]**2 + X[1]**2` is an expression
  on `RD(2)` with value shape `()`.

### Index notation

>>> LLM generated content.

Index notation writes tensor algebra as fixed compositions, so that operators
such as `transpose`, `inner` and `dot` are free functions rather than node
classes.

#### Index

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
`Index(("transpose", A, 1))`. An index whose label contains `A` cannot occur
inside `A`, because `A` would then contain `A`, and the graph is acyclic. Hence
a derived index is fresh with respect to the operands, and the free function
does not inspect the operands to find a fresh index. A global counter is not
allowed, because a counter makes construction depend on global state.

Note that `transpose(A)` `is` or `==` `transpose(A)` with this design.

*Downsides:* expressions that differ only by bound-index labels are not equal
and are not shared, every substitution must avoid capturing indices, and
derived indices enlarge the graph. Later passes, not construction, could
address these downsides by hashing modulo renaming of bound indices (Maziarz et
al., PLDI 2021) and by removing index-notation scaffolding.

#### Free indices

Every `AbstractExpression` derives the free indices, a `tuple` of pairs of an
`Index` and an extent. The order is the order of first occurrence, never the
order of a `set`.

Unless a node class states otherwise, the free indices of a node are the free
indices of the edges of the node, in edge order, with each `Index` taken once
at the first occurrence. An `Index` that occurs in two edges must have the same
extent in both edges, or construction raises `ValueError`. Terminals and
literals have no free indices.

The following node classes state otherwise:

- `Sum`: both operands must have the same free indices with the same extents,
  in any order, or construction raises `ValueError`. The order of the free
  indices of the sum is the order of the first operand.
- `Indexed`: see [Indexed](#indexed).
- `ComponentTensor` and `IndexSum`: see below.

Summation is never implicit. `a[i] * b[i]` has the free index `i` and denotes
the component-wise product. The sum over `i` is `index_sum(a[i] * b[i], i)`.

Equality is structural. Two expressions that differ only by the names of bound
indices are not equal. A later algorithmic pass may compare expressions up to
renaming of bound indices.

A form integrand must have no free indices.

#### ComponentTensor

A concrete subclass of `AbstractExpression` that turns free indices into
tensor axes.

`ComponentTensor(expression, indices)` stores a scalar expression and a `tuple`
of distinct `Index`es that are free in the expression. The expression is the
only edge. `ComponentTensor(expression, indices)` derives the value shape from
the extents of the indices, in the order of the indices, and the free indices
of the expression without the bound indices.

The free function `component_tensor(expression, indices)` returns
`ComponentTensor(expression, indices)`. The free function raises `ValueError`
if the expression is not scalar, if an index is repeated, or if an index is not
free in the expression.

#### IndexSum

A concrete subclass of `AbstractExpression` for the sum of a scalar expression
over one free index.

`IndexSum(expression, index)` stores a scalar expression and an `Index` that is
free in the expression. The expression is the only edge.
`IndexSum(expression, index)` derives the value shape `()` and the free indices
of the expression without the summation index.

The free function `index_sum(expression, index)` returns
`IndexSum(expression, index)`. The free function raises `ValueError` if the
expression is not scalar, or if the index is not free in the expression.

Examples, with `i, j = indices("i", "j")`, two vectors `a` and `b` of value
shape `(3,)`, and `A = matrix(((1, 2), (3, 4)))`:

- `A[i, j]` has value shape `()` and the free indices `((i, 2), (j, 2))`.
- `a[i] * b[i]` has value shape `()` and the free indices `((i, 3),)`.
- `index_sum(a[i] * b[i], i)` is the inner product of `a` and `b`, with no
  free indices.
- `component_tensor(A[i, j], (j, i))` has value shape `(2, 2)` and is the
  transpose of `A`, written with user indices.
- `a[i] + b[j]` raises `ValueError`, because the free indices differ.
- `A[i, i]` raises `ValueError`, because `i` occurs twice.

### Annotation

An annotation records how an object was produced. An annotation is a `tuple`
of a `str` name followed by the arguments that produced the object, for example
`("metric", D)`.

The mixin `Annotation` holds what every annotated class shares. Each kind of
object that can be annotated has one concrete annotated class, which combines
`Annotation` with the abstract class of that kind. An annotated object is an
object of the same kind as the wrapped object, so an annotated object stands in
wherever the wrapped object does.

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
for the kind of `x`, with the annotation `(name, *arguments)`.

Every free function in the core that is defined as a composition annotates the
result with the name and arguments of the free function. For example,
`metric(D)` is `AnnotatedExpression(transpose(J) @ J, ("metric", D))`.

#### AnnotatedExpression

`AnnotatedExpression(expression, annotation)` combines `Annotation` with
`AbstractExpression`. `AnnotatedExpression` derives the value shape, the set
and the free indices of the expression. `grad` of an annotated expression
differentiates the expression and drops the annotation.

#### AnnotatedDomain

`AnnotatedDomain(domain, annotation)` combines `Annotation` with
`AbstractDomain`. `AnnotatedDomain` derives the dimensions, the parameter
domain, the mapping and the ambient domain of the domain.

An annotated domain is an `AbstractDomain` but not an
`AbstractCoordinateDomain`, even when the domain is a coordinate domain. A
coordinate function, a gradient and a parameterised domain need the plain
coordinate domain.

Examples:

- `annotate(X[0]**2 + X[1]**2, "height", X)` is an `AnnotatedExpression` with
  value shape `()`, and prints as `height(X)`.
- `annotate(parameterised_domain(RD(2), psi), "paraboloid")` is an
  `AnnotatedDomain` with $k = 2$ and $n = 3$.

### Tensor algebra

#### Determinant

`Determinant` is a concrete subclass of `UnaryOperator`. The free function
`determinant(f)` takes an expression of value shape `(n, n)` and returns a
scalar. `Determinant` stores one expression, which is the only edge, and
derives the set of that expression. A determinant of value shape `(n, n)` is a
polynomial of degree $n$ in the entries, so no fixed composition expresses the
determinant, and `Determinant` is a node class.

#### Transpose and matrix product

Transpose and matrix product are free functions built with derived indices:

- `transpose(A)` is `component_tensor(A[k, l], (l, k))`, with
  $k$ = `Index(("transpose", A, 0))` and $l$ = `Index(("transpose", A, 1))`,
  annotated with `("transpose", A)`.
- `matrix_product(A, B)`, with the overload `A @ B`, is
  `component_tensor(index_sum(A[i, k] * B[k, j], k), (i, j))`, with $i$, $k$
  and $j$ the derived indices `Index(("matrix_product", A, B, m))` for
  $m = 0, 1, 2$, annotated with `("matrix_product", A, B)`.

`transpose` raises `ValueError` unless `A` has rank 2. `matrix_product` raises
`ValueError` unless `A` and `B` have rank 2 and the second extent of `A` equals
the first extent of `B`. The overload `@` accepts only rank 2 operands.
Contractions of other ranks, such as a matrix-vector product, are written in
index notation.

Examples, with `A = matrix(((1, 2), (3, 4)))`:

- `transpose(A)` is `component_tensor(A[k, l], (l, k))`, annotated with
  `("transpose", A)`, with the derived indices
  `k = Index(("transpose", A, 0))` and `l = Index(("transpose", A, 1))`.
- `transpose(A)` is not equal to `component_tensor(A[i, j], (j, i))` for user
  indices `i` and `j`, because the index labels differ and because
  `transpose(A)` is annotated.

### Differentiation

#### Gradient

> LLM generated content.

A concrete subclass of `UnaryOperator` for the gradient of an expression $f$
on a coordinate domain $P$ of dimension $n$. The gradient of $f$ is the
derivative of $f$ with respect to the coordinate function $X$ of $P$.

`Gradient(f)` stores $f$, which is the only edge. The set of $f$ must be an
`AbstractCoordinateDomain` $P$. `Gradient(f)` derives the value shape
`f.value_shape + (n,)` and the set $P$.

The free function `grad(f)` returns `Gradient(f)`. The free function raises `ValueError` if $f$ has no set, or if
the set of $f$ is not a coordinate domain.

`Gradient` is lazy. Construction does not evaluate the derivative. A later pass
evaluates the derivative with the standard rules, including `grad(X)` $= I$
and a zero derivative for every literal.

The gradient on a parameterised domain is not a primitive. The gradient on a
parameterised domain is a composite of `grad` on the parameter domain and the
geometry of the parameterised domain (see
[Geometry](#geometry)).

Examples, with `X = coordinate(RD(2))`:

- `grad(X)` has value shape `(2, 2)`. A later pass evaluates `grad(X)` to the
  identity.
- `grad(X[0]**2 + X[1]**2)` has value shape `(2,)`. A later pass evaluates the
  gradient to $(2 X_0, 2 X_1)$.
- `grad(1)` raises `ValueError`, because the literal `Integer(1)` has no set.

### Parameterised domains

#### ParameterisedDomain

A concrete subclass of `AbstractDomain` for parameterised domains.

A parameterised domain $D$ is the image $\psi(P)$ of a coordinate domain
$P \subset R^k$ under a map $\psi : P \to R^n$.

`ParameterisedDomain(P, psi)` stores:

- The parameter domain $P$, an `AbstractCoordinateDomain` of dimension $k$,
  returned by `parameter_domain(D)`.
- The mapping $\psi$, an `AbstractExpression` on $P$ with value shape `(n,)`,
  returned by `mapping(D)`. The mapping $\psi$ is built from the coordinate
  function $X$ of $P$.

$P$ and $\psi$ are the two edges, in that order.

`ParameterisedDomain(P, psi)` derives:

- The topological dimension $k$, which is the dimension of $P$.
- The geometric dimension $n$, which is taken from the value shape `(n,)` of
  $\psi$. The geometric dimension is derived, not stored, so the geometric
  dimension cannot disagree with the value shape of $\psi$.
- The ambient domain `ambient_domain(D)`, which is `RD(n)`.

A parameterised domain does not provide a coordinate function.

The free function `parameterised_domain(P, psi)` returns
`ParameterisedDomain(P, psi)`. The free function raises
`ValueError` if $\psi$ does not have rank 1, if $\psi$ has a set other than
$P$, or if $k > n$.

A coordinate domain is *not* a parameterised domain. `AbstractCoordinateDomain`
and `ParameterisedDomain` are both subclasses of `AbstractDomain`. Both classes
provide `parameter_domain` and `mapping`.

Examples:

- The hypotenuse of the reference triangle is a parameterised domain with
  $k = 1$ and $n = 2$. $P$ is the reference interval in $R^1$, and
  $\psi(X) = (1 - X_0, X_0)$.
- `parameterised_domain(RD(2), vector((X[0], X[1], X[0]**2 + X[1]**2)))`, with
  `X = coordinate(RD(2))`, is the paraboloid, with $k = 2$ and $n = 3$.
  `ambient_domain` of the paraboloid is `RD(3)`.
- A physical triangle in $R^2$ is a parameterised domain with $k = n = 2$. $P$
  is a reference triangle, and $\psi$ is the affine map from $P$ onto the
  physical triangle. See the affinely mapped triangle example.

### Geometry

The geometry of every domain $D$, coordinate or parameterised, comes from the
mapping $\psi$ = `mapping(D)` on the parameter domain $P$ =
`parameter_domain(D)` of dimension $k$:

| Free function        | Definition                       | Value shape |
| -------------------- | -------------------------------- | ----------- |
| `jacobian(D)`        | $J$ = `grad(mapping(D))`         | `(n, k)`    |
| `metric(D)`          | $g$ = `transpose(J) @ J`         | `(k, k)`    |
| `measure_density(D)` | `square_root(determinant(g))`    | `()`        |

Each free function annotates the result with the name of the free function
and $D$, for example `("jacobian", D)`.

The Jacobian, metric and measure density are lazy expressions. Construction
performs no simplification.

For a coordinate domain, the mapping is the coordinate function $X$. The
general rules give $J$ = `grad(X)` $= I$, $g = I$ and a measure density of 1,
so a coordinate domain needs no special case.

### Example: a paraboloid

The paraboloid over `RD(2)` collects the examples of the sections above.

```python
# AbstractCoordinateDomain: RD(2), no edges.
P = RD(2)

# AbstractExpression: Coordinate(RD(2)), value shape (2,).
X = coordinate(P)

# AbstractExpression: Tensor((Indexed(X, (0,)), Indexed(X, (1,)), Sum(...))),
# value shape (3,).
psi = vector((X[0], X[1], X[0]**2 + X[1]**2))

# AbstractDomain: ParameterisedDomain(RD(2), psi), k = 2, n = 3.
D = parameterised_domain(P, psi)

# AbstractExpression: AnnotatedExpression(Gradient(psi), ("jacobian", D)),
# value shape (3, 2).
J = jacobian(D)

# AbstractExpression: AnnotatedExpression(transpose(J) @ J, ("metric", D)),
# value shape (2, 2).
g = metric(D)

# AbstractExpression: AnnotatedExpression(SquareRoot(Determinant(g)),
# ("measure_density", D)), value shape ().
rho = measure_density(D)
```

Every object in the example is lazy.

Note that at this stage we do not need the concept of a `Point` as we have not
evaluated anything at a `Point` so far.

### Example: an affinely mapped triangle

The physical triangle with vertices $v_0 = (1, 0)$, $v_1 = (3, 1)$ and
$v_2 = (2, 3)$ is the image of a reference triangle under an affine map.

Important point is that there is no topology here, just geometry.

$\{p \in R^n : p_i \ge 0 \text{ for every } i, \ \sum_i p_i \le 1\}$

as a concrete `AbstractCoordinateDomain`. A concrete coordinate domain provides
only the dimension $n$. Everything else comes from `AbstractCoordinateDomain`.

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

`UnitSimplex(2)` is the reference triangle with vertices $(0, 0)$, $(1, 0)$ and
$(0, 1)$. The affine map from `UnitSimplex(2)` onto the physical triangle is

$$\psi(X) = v_0 + X_0 (v_1 - v_0) + X_1 (v_2 - v_0)
= (1 + 2 X_0 + X_1, X_0 + 3 X_1).$$

```python
# AbstractCoordinateDomain: UnitSimplex(2), no edges.
T = UnitSimplex(2)

# AbstractExpression: Coordinate(UnitSimplex(2)), value shape (2,).
X = coordinate(T)

# AbstractExpression: Tensor(..., ...), value shape (2,).
psi = vector((1 + 2 * X[0] + X[1], X[0] + 3 * X[1]))

# AbstractDomain: ParameterisedDomain(UnitSimplex(2), psi), k = n = 2.
K = parameterised_domain(T, psi)

# AbstractExpression: AnnotatedExpression(Gradient(psi), ("jacobian", K)),
# value shape (2, 2).
J = jacobian(K)

# AbstractExpression: AnnotatedExpression(SquareRoot(Determinant(metric(K))),
# ("measure_density", K)), value shape ().
rho = measure_density(K)
```

A later pass evaluates the Jacobian to the constant matrix

$$J = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix},$$

so $\det J = 5$ and $\sqrt{\det g} = |\det J| = 5$. The area of the physical
triangle is $5$ times the area $1/2$ of the reference triangle.

The physical triangle is a parameterised domain with $k = n$, not a coordinate
domain, although the physical triangle has non-empty interior in $R^2$. The
class records how the domain is built.

### Cells

A domain can be divided into finitely many cells. A cell is a piece of the
domain with the same dimension as the domain. Two cells meet only on the
boundaries of the two cells.

A facet is a piece of the boundary of a cell, one dimension lower than the
cell. A facet in the boundary of one cell is an exterior facet. A facet shared
by two cells is an interior facet.

Cells and facets are topological. Defining cells and facets does not require a
metric.

# JSH worked to here.

### Finite elements

A finite element describes the local finite-dimensional field used on a cell.

At minimum, an element is defined by:

- the cell type,
- the element family,
- the polynomial degree or equivalent order information,
- the value shape,
- and, for mapped elements, the required pullback or mapping behaviour.

Examples:

```python
P1 = FiniteElement(family="Lagrange", cell=triangle, degree=1, shape=())
vP2 = FiniteElement(family="Lagrange", cell=triangle, degree=2, shape=(2,))
```

The element describes what kind of discrete field is being represented. It does
not store actual coefficient values.

### Domains

A domain is the integration domain of a form. For a minimal language
definition, a domain only needs to provide:

- topological dimension,
- geometric dimension,
- and a coordinate element, if geometry is represented using a finite element.

The actual mesh is still external to UFLx.

### Function spaces

A function space defines a finite-dimensional space of functions. There are
three types of function spaces:

1. standard finite-element spaces, defined by a domain and an element,
1. constant function spaces, defined by a shape and a scalar type,
1. non-finite-element function spaces, defined by domain, shape and scalar
   type, but no element.

Examples:

```python
V = FunctionSpace(domain=domain, element=element)  # standard finite-element space
C = FunctionSpace(shape=(2,), scalar_type="real")  # constant function space
N = FunctionSpace(
    domain=domain, shape=(2,), scalar_type="real"
)  # non-finite-element function space
```

In older UFL examples, `Argument` and `Coefficient` are sometimes created
directly from an element. This is not allowed.

### Terminals

Terminals are the leaves of UFLx expression trees.

The essential terminals are:

- `Argument`, representing an arbitrary basis function of a function space,

  It has a `number` attribute to represent the position of the `Argument` in
  the form. For example, consider a bilinear form `a(u, v)`. The `Argument`
  representing `u` has `number=1`, and the `Argument` representing `v` has
  `number=0`.

  > Experimental (API breaking) note: Existence of an `Argument` and a
  > `Coefficient` is rooted deep into the legacy UFL. However, more useful
  > concept would be to only allow for `Function` terminals, which can
  > represent both known and unknown functions. It would be a responsibility of
  > the user, when creating a form, to define which of the `Function` terminals
  > are to be considered as "arguments" and which not. This would be much
  > closer to a mathematical notation and would solve many user error of having
  > to treat `Argument` and `Coefficient` differently.

- `Coefficient`, representing a known function in a function space,

- literals such as numbers,

- and geometric quantities such as coordinates, normals, Jacobians, cell
  diameter, facet area, etc.

Example:

```python
v = Argument(function_space=V, name="v", number=0)
u = Argument(function_space=V, name="u", number=1)
f = Coefficient(function_space=V, name="f")
```

> Legacy note: There is no `Constant` in UFLx as a means to represent literal
> value which is not known at compile time.

UFLx supports named symbolic terminals. This allows for better error messages
and better postprocessing of expressions. In addition, it improves the
stability of the expression tree.

> Legacy note: In legacy UFL, symbolic objects do not have names. It requires
> expensive and error prone renumbering of terminals in expressions before
> lowering and code generation.

### Expressions

An expression is either a terminal or an operator applied to other expressions.

Examples of operators are addition, multiplication, indexing, `grad`, `div`,
`dot`, `inner`, `jump`, and `avg`.

Expressions are symbolic. They are not evaluated when they are created.

A minimal expression carries these static properties:

- tensor shape,
- free indices,
- associated domain or domains,
- scalar type (real or complex),
- continuity information,
- (possibly) physical dimension.

These properties are enough to reject many meaningless expressions before
lowering or code generation.

### Measures

A measure describes where and how an integrand is integrated.

TODO.

The standard measures are:

- `dx` for cell integrals,
- `ds` for exterior-facet integrals,
- `dS` for interior-facet integrals.

Examples:

```python
dx
ds
dS
```

### Integrals

An integral is made by multiplying a scalar expression by a measure.

Examples:

```python
v * f * dx
inner(grad(v), grad(u)) * dx
jump(v) * dS
```

A valid integral has a scalar-valued integrand with no unresolved free indices.

### Forms

A form is a multi-linear functional
$a(u, v, w, \ldots) : U \times V \times W \times \ldots \to \mathbb{R}$ or
$\mathbb{C}$.

Arguments determine the arity of a form.

> Legacy note: In contrast to the legacy UFL, in UFLx a form does not need to
> be defined as a sum of integrals. For example, the following is a valid UFLx
> form: `a = (u * dx) * (v * dx)`, or
> `a = PointEvaluation(u, (0.5, 0.5)) * PointEvaluation(v, (0.5, 0.5))`.

Examples:

```python
a = inner(grad(v), grad(u)) * dx
L = v * f * dx
F = a - L
```

> Experimental (API breaking) note: Similar to the `Argument` and `Coefficient`
> terminals, the `Form` in the legacy UFL is defined implicitly, by the
> presence of `Argument` terminals. In UFLx, a `Form` could be an explicit
> object that can be created by the user, e.g.

```python
a = Form(inner(grad(v), grad(u)) * dx, arguments=[u, v])
```

## 3. Minimal validity rules

The language is small, but it is typed.

- Addition and subtraction require compatible shapes and compatible free
  indices.
- Multiplication follows scalar, tensor, matrix-vector, matrix-matrix, and
  indexed-expression rules.
- Tensor contractions such as `dot` and `inner` require compatible tensor
  shapes.
- Spatial derivatives require a spatial domain.
- Facet quantities such as `FacetNormal` are only meaningful in facet
  integration contexts.
- Interior-facet expressions may use side restrictions:
  ```python
  u("+")
  u("-")
  ```
- Operators such as `jump` and `avg` are shorthand for combinations of such
  restrictions.
- An integral is valid only when the final integrand is scalar and has no free
  indices.
- A form is valid when all of its integrals are valid and their domains,
  arguments, and coefficients are mutually consistent.

## 4. What transformation procedures does the UFLx core provide

There are two key transformation procedures that the UFLx core provides:

1. Automatic/symbolic differentiation.

   Similar to UFL, UFLx provides a symbolic differentiation via a directional
   (Gateaux) derivative. The derivative of an Expression, Integral or a Form is
   lazy/eagerly evaluated. That means the derivative is stored as an expression
   tree operator, and the actual differentiation is performed when the
   derivative is used in a context where it needs to be evaluated, such as
   lowering or code generation.

   Example:

   ```python
    a = inner(grad(v), grad(u)) * dx
    da = derivative(a, u, du)
   ```

1. Expression and integral transformation between different configurations.

   > Legacy note: In legacy UFL, there is a very central perspective on the
   > configuration of where objects are defined. Usually, all expressions are
   > defined in the physical configuration, and UFL applies pullbacks to the
   > entire expression tree.

   In UFLx, configuration where expressions are defined are a consequence of
   the domain. Since the domain is part of the Expression's static attributes,
   it is propagated through the expression tree based on the language rules.

   TODO.

## 5. What belongs to the language

The minimal language includes:

- symbolic finite-element objects,
- symbolic domains/function spaces,
- symbolic terminals,
- tensor-valued expressions,
- index notation,
- differential operators,
- DG side restrictions,
- integration measures,
- integrals,
- and forms.

This is enough to describe ordinary weak forms such as mass, stiffness,
elasticity, Stokes, mixed Poisson, and DG forms.

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

## 7. Extensibility of the attributes

> Legacy note: In legacy UFL, the set of static attributes is fixed and
> hardcoded. This led to many language features being hacked in (e.g. complex
> numbers support, dual spaces, etc.).

In UFLx, the set of static attributes is extensible. This allows for cleaner
language extensions in the future, such as support for physical dimensions,
dual spaces, error estimates, custom quadrature rules, custom
integrability/continuity rules, etc.

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
